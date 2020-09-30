#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    pyinknife-run.py
#
#    Run the PyInKnife pipeline for Protein Structure Network
#    analysis with PyInteraph.
#
#    Copyright (C) 2020 Valentina Sora 
#                       <sora.valentina1@gmail.com>
#                       Juan Salamanca Viloria 
#                       <juan.salamanca.viloria@gmail.com> 
#                       Elena Papaleo
#                       <elenap@cancer.dk>
#
#    This program is free software: you can redistribute it and/or
#    modify it under the terms of the GNU General Public License as
#    published by the Free Software Foundation, either version 3 of
#    the License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public
#    License along with this program. 
#    If not, see <http://www.gnu.org/licenses/>.


# standard library
import argparse
import logging
import os
import os.path
# third-party packages
import dask
from dask import distributed
from distributed import wait
import MDAnalysis as mda

from . import util



def main():


    ######################### ARGUMENT PARSER #########################


    # description of the script
    description = "\nRun the PyInKnife pipeline for Protein " \
                  "Structure Network Analysis with PyInteraph.\n"
    
    # create the argument parser
    parser = argparse.ArgumentParser(description = description)
    
    # add arguments to the parser
    f_help = "Trajectory."
    parser.add_argument("-f", "--trj", \
                        type = str, \
                        required = True, \
                        help = f_help)

    s_help = "Topology."
    parser.add_argument("-s", "--top", \
                        type = str, \
                        required = True, \
                        help = s_help)

    r_help = "Reference structure."
    parser.add_argument("-r", "--ref", \
                        type = str, \
                        default = None, \
                        help = r_help)
    
    c_help = f"Configuration file used to run the pipeline. Default " \
             f"is {util.DEFCONFIG}."
    parser.add_argument("-c", "--configfile", \
                        type = str, \
                        default = util.DEFCONFIG, \
                        help = c_help)

    d_help = "Directory where to run the pipeline. " \
             "Default is the current working directory."
    parser.add_argument("-d", "--rundir", \
                        type = str, \
                        default = os.getcwd(), \
                        help = d_help)

    n_default = 1
    n_help = f"Number of processes to be started in parallel. " \
             f"Default is {n_default} process(es)."
    parser.add_argument("-n", "--nproc", \
                        type = int, \
                        default = n_default, 
                        help = n_help)

    ncaa_help = "Noncanonical residues present in your system."
    parser.add_argument("-ncaa", "--noncanonical-residues", \
                        type = str, \
                        nargs = "*", \
                        default = None, 
                        help = ncaa_help)

    # parse the arguments
    args = parser.parse_args()
    # get single arguments
    trj = util.get_abspath(args.trj)
    top = util.get_abspath(args.top)
    ref = util.get_abspath(args.ref)
    configfile = util.get_abspath(args.configfile)
    rundir = util.get_abspath(args.rundir)
    nproc = args.nproc
    ncaa = args.noncanonical_residues


    ########################## CONFIGURATION ##########################

    # get the configuration
    config = util.get_pyinknife_configuration(configfile = configfile)
    # get the logging configuration
    logging.basicConfig(level = logging.INFO)


    ######################## PYINKNIFE PIPELINE #######################

  
    # change scheduler if you want to use threads or a single core
    with dask.config.set(scheduler = "processes"):

        # NB: Universe creation must be done internally by any
        # function requiring it since Universes are unpicklable
        # objects, therefore not serializable

        #---------------------------- Dask ---------------------------#

        # create a local cluster with the desired number of workers
        cluster = distributed.LocalCluster(n_workers = nproc, \
                                           silence_logs = logging.INFO, \
                                           processes = True, \
                                           threads_per_worker = 1)
        
        # set a client to submit the jobs to
        client = distributed.Client(cluster)

        # create a list for graph_analysis futures
        futures = []

        #------------------------ Configuration ----------------------#

        # get the configurations for the different runs
        # and for the resampling
        resconfig = config["resampling"]
        pyinconfig = config["pyinteraph"]
        fgconfig = config["filter_graph"]
        gaconfig = config["graph_analysis"]

        #---------------------- Topology/Reference -------------------#

        # create a MDAnalysis Universe object from the topology
        # and the trajectory
        u = mda.Universe(top, trj)
        # get the trajectory
        trajectory = u.trajectory
        # set the selection string for atoms
        selstring = client.submit(util.get_protein_selection_string, \
                                  ncaa = ncaa)
        
        # if the user provided a reference structure
        if ref:
            # the PyInteraph topology will be the topology passed
            # by the user
            pyintop = top
            # the PyInteraph reference will be the reference structure
            pyinref = ref
        # if no reference structure was provided
        else:
            # create the file that will serve as topology 
            # for PyInteraph
            pyintop = client.submit(util.write_trj_subset, \
                                    filename = config["pyintop"], \
                                    trj = trj, \
                                    top = top, 
                                    selstring = selstring,
                                    interval = ((0,1),), \
                                    wd = rundir)
            # create the file that will serve as reference
            # for PyInteraph
            pyinref = client.submit(util.write_trj_subset, \
                                    filename = config["pyinref"], \
                                    trj = trj, \
                                    top = top, \
                                    selstring = selstring, 
                                    interval = ((0,1),), \
                                    wd = rundir)

        #------------------------- Resamplings -----------------------#

        # get the sets of frames composing the full trajectory
        # and possible resampled trajectories
        resamplings = client.submit(util.get_frames_sets, \
                                    trajectory = trajectory, \
                                    config = resconfig).result()

        # for each trajectory (full trajectory and resamplings)
        for i, (trjname, interval) in enumerate(resamplings):

            #------------------- Set the trajectory ------------------#

            # set a subdirectory for the trajectory
            resdir = os.path.join(rundir, trjname)
            # the first trajectory is the full trajectory,
            if i == 0:
                # the sub-trajectory will be the full trajectory
                subtrj = trj
            else:
                # write the sub-trajectory
                subtrj = client.submit(util.write_trj_subset, \
                                       filename = trjname + ".xtc", \
                                       trj = trj, \
                                       top = top, \
                                       selstring = selstring, \
                                       interval = interval, \
                                       wd = resdir).result()

            #-------------------- Set the analysis -------------------#

            # for each analysis
            for analysis in pyinconfig.keys():
                anconfig = pyinconfig[analysis]
                # initialize to None the modes available
                modes = [None]
                # continue if it was not requested
                if not anconfig["run"]:
                    continue
                # get the options
                pyinargs = anconfig["args"]
                # set a subdirectory for the analysis and set the 
                # subdirectory for alternative analysis modes equal
                # to the parent analysis subdirectory
                andir = os.path.join(resdir, analysis)
                modedir = andir
                # check if analyses modes are present/requested
                if "modes" in anconfig.keys():
                    # set the current modes to the ones found
                    modes = anconfig["modes"]
                
                # for each analysis mode
                for mode in modes:
                    # if it is not None (alternative analysis
                    # mode for the analysis)
                    if mode:
                        # set the subdirectory for the mode
                        modedir = os.path.join(andir, mode)

                    #------------------- pyinteraph ------------------#

                    # for each distance cut-off  
                    for dcut in anconfig["dcuts"]:
                        # set the subdirectory for the cut-off 
                        dcutdir = os.path.join(modedir, str(dcut))
                        # set the path to the output file
                        pyinout = os.path.join(dcutdir, \
                                               anconfig["out"])
                        # run the PSN calculation
                        pyinmat = client.submit(util.run_pyinteraph, \
                                                trj = subtrj, \
                                                top = pyintop, \
                                                ref = pyinref, \
                                                analysis = analysis, \
                                                mode = mode, \
                                                co = dcut, \
                                                wd = dcutdir, 
                                                out = pyinout, \
                                                otherargs = pyinargs)
                        
                        #---------------- filter_graph ---------------#
                        
                        # get network filtering options
                        fgargs = fgconfig["args"]

                        # for each persistence cut-off
                        for pcut in fgconfig["pcuts"]:
                            # set the subdirectory for the cut-off
                            pcutdir = os.path.join(dcutdir, str(pcut))
                            # set the output/log file
                            fgout = os.path.join(pcutdir, \
                                                 fgconfig["out"])
                            # run network filtering
                            fgmat = \
                                client.submit(util.run_filter_graph, \
                                              mat = pyinmat, \
                                              out = fgout, \
                                              perco = pcut, \
                                              wd = pcutdir, \
                                              otherargs = fgargs)

                            #------------- graph_analysis ------------#

                            # for each type of network analysis
                            for netanalysis in gaconfig.keys():
                                # get the configuration of the
                                # network analysis
                                netanconfig = gaconfig[netanalysis]
                                # set the subdirectory for the network
                                # analysis
                                netandir = os.path.join(pcutdir, \
                                                        netanalysis)
                                # set the output/log file
                                gaout = os.path.join(netandir, \
                                                     netanconfig["out"])
                                gaargs = netanconfig["args"]
                                # run network analysis     
                                futures.append(\
                                    client.submit(\
                                        util.run_graph_analysis, \
                                            mat = fgmat, \
                                            ref = pyinref, \
                                            out = gaout, \
                                            wd = netandir, \
                                            otherargs = gaargs))

        # wait on the last active futures (graph_analysis futures,
        # since they have no futures depending on them)
        wait(futures)


if __name__ == "__main__":
    main()
