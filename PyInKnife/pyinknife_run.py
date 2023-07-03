#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    pyinknife_run.py
#
#    Run the PyInKnife pipeline for Protein Structure Network
#    analysis with PyInteraph2.
#
#    Copyright (C) 2023 Valentina Sora 
#                       <sora.valentina1@gmail.com>
#                       Juan Salamanca Viloria 
#                       <juan.salamanca.viloria@gmail.com>
#                       Matteo Tiberti 
#                       <matteo.tiberti@gmail.com>
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


# Standard library
import argparse
import logging as log
import os
import os.path
import sys
# Third-party packages
import dask
from dask import distributed
from distributed import wait
import MDAnalysis as mda
# PyInKnife2
from . import util


# Get the module logger
logger = log.getLogger(__name__)


def main():


    #------------------------ Argument parser ------------------------#


    # Description of the script
    description = \
        "\nRun the PyInKnife pipeline for Protein " \
        "Structure Network analysis with PyInteraph.\n"
    
    # Create the argument parser
    parser = argparse.ArgumentParser(description = description)
    
    # Add arguments to the parser
    f_help = "Trajectory."
    parser.add_argument("-f", "--trj",
                        type = str,
                        required = True,
                        help = f_help)

    s_help = "Topology."
    parser.add_argument("-s", "--top",
                        type = str,
                        required = True,
                        help = s_help)

    r_help = "Reference structure."
    parser.add_argument("-r", "--ref",
                        type = str,
                        default = None,
                        help = r_help)

    c_help = \
        f"Configuration file used to run the pipeline. " \
        f"Default is {util.DEF_CONFIG}."
    parser.add_argument("-c", "--configfile",
                        type = str,
                        default = util.DEF_CONFIG,
                        help = c_help)

    d_help = "Directory where to run the pipeline."
    parser.add_argument("-d", "--rundir",
                        type = str,
                        required = True,
                        help = d_help)

    n_default = 1
    n_help = \
        f"Number of processes to be started in parallel. " \
        f"Default is {n_default} process(es)."
    parser.add_argument("-n", "--nproc",
                        type = int,
                        default = n_default,
                        help = n_help)

    ncaa_help = "Noncanonical residues present in your system."
    parser.add_argument("-ncaa", "--noncanonical-residues",
                        type = str,
                        nargs = "*",
                        default = None, 
                        help = ncaa_help)

    # Parse the arguments
    args = parser.parse_args()
    
    # Get the arguments
    trj = util.get_abspath(args.trj)
    top = util.get_abspath(args.top)
    ref = util.get_abspath(args.ref)
    config_file = util.get_abspath(args.configfile)
    run_dir = util.get_abspath(args.rundir)
    n_proc = args.nproc
    ncaa = args.noncanonical_residues


    #------------------------- Configuration -------------------------#


    # Try to load the configuration
    try:

        config = \
            util.get_pyinknife_configuration(\
                config_file = config_file)

    # If something went wrong
    except Exception as e:

        # Log an error and exit
        errstr = \
            f"It was not possible to load the configuration from " \
            f"'{config_file}'. Error: {e}"
        logger.exception(errstr)
        sys.exit(errstr)


    #---------------------------- Logging ----------------------------#


    # Get the module logger
    logger = log.getLogger(__name__)

    # Configure the logger
    log.basicConfig(level = log.INFO)


    #---------------------- PyInKnife pipeline -----------------------#

  
    # Change scheduler if you want to use threads or a single core
    with dask.config.set(scheduler = "processes"):

        # NB: Universe creation must be done internally by any
        # function requiring it since Universes are unpicklable
        # objects, therefore not serializable


        #---------------------------- Dask ---------------------------#


        # Create a local cluster with the desired number of workers
        cluster = \
            distributed.LocalCluster(\
                # Number of workers in the cluster
                n_workers = n_proc,
                # Log messages below this level will be suppressed
                silence_logs = log.INFO,
                # Whether to use processes or nannies
                processes = True,
                # How many threads to use for each worker
                threads_per_worker = 1)
        
        # Create a client to submit the jobs to
        client = distributed.Client(cluster)

        # Create a list to store pending futures
        futures = []


        #------------------------ Configuration ----------------------#


        # Get the configurations for the resampling procedure
        res_config = config["resampling"]

        # Get the configuration to use for 'pyinteraph' runs
        pyin_config = config["pyinteraph"]

        # Get the configuration to use for 'filter_graph' runs
        fg_config = config["filter_graph"]

        # Get the configuration to use for 'graph_analysis' runs
        ga_config = config["graph_analysis"]


        #--------------------- Running directory ---------------------#


        # If the running directory does not exist yet
        if not os.path.exists(run_dir):

            # Create it
            os.makedirs(run_dir, exist_ok = True)

            # Inform the user that the directory was created
            infostr = \
                f"The directory '{run_dir}' was created."
            logger.info(infostr)

        # Otherwise
        else:

            # Inform the user that the directory was found
            infostr = \
                f"The directory '{run_dir}' was found."
            logger.info(infostr)


        #--------------------------- System --------------------------#


        # Try to create a MDAnalysis Universe object from the topology
        # and the trajectory
        try:
            
            u = mda.Universe(top, trj)

        # If something went wrong
        except Exception as e:

            # Log an error and exit
            errstr = \
                f"It was not possible to load the system. Error: {e}"
            logger.error(errstr)
            sys.exit(errstr)
        
        # Get the trajectory
        trajectory = u.trajectory
        
        # Set the selection string for atoms to be considered
        sel_string = util.get_protein_selection_string(ncaa = ncaa)


        #------------------------- Reference -------------------------#
        

        # If the user provided a reference structure
        if ref is not None:
            
            # The PyInteraph topology will be the topology passed
            # by the user
            pyin_top = top
            
            # The PyInteraph reference will be the reference structure
            pyin_ref = ref
        
        # If no reference structure was provided
        else:
            
            # Create the file that will serve as topology
            # for PyInteraph
            pyin_top = \
                util.write_trj_subset(filename = config["pyintop"],
                                      trj = trj,
                                      top = top, 
                                      sel_string = sel_string,
                                      interval = ((0,1),),
                                      wd = run_dir)
            
            # Create the file that will serve as reference
            # for PyInteraph
            pyin_ref = \
                util.write_trj_subset(filename = config["pyinref"],
                                      trj = trj,
                                      top = top,
                                      sel_string = sel_string, 
                                      interval = ((0,1),),
                                      wd = run_dir)

        # Inform the user about which topology will be used
        infostr = \
            f"'{pyin_top}' will be used as topology."
        logger.info(infostr)

        # Inform the user about which reference structure will be used
        infostr = \
            f"'{pyin_ref}' will be used as reference structure."
        logger.info(infostr)


        #------------------------- Resamplings -----------------------#


        # Get the sets of frames composing the full trajectory
        # and possible resampled trajectories
        resamplings = \
            util.get_frames_sets(trajectory = trajectory,
                                 config = res_config)

        # for each trajectory (full trajectory and resamplings)
        for i, (trj_name, interval) in enumerate(resamplings):


            #------------------- Set the trajectory ------------------#


            # Set a sub-directory for the trajectory
            res_dir = os.path.join(run_dir, trj_name)
            
            # If the trajectory is the full trajectory (put first
            # in the list)
            if i == 0:
                
                # The sub-trajectory will be the full trajectory
                sub_trj = trj
            
            # Otherwise
            else:
                
                # Write the sub-trajectory
                sub_trj = \
                    client.submit(util.write_trj_subset,
                                  filename = trj_name + ".xtc",
                                  trj = trj,
                                  top = top,
                                  sel_string = sel_string,
                                  interval = interval,
                                  wd = res_dir).result()


            #-------------------- Set the analysis -------------------#


            # For each analysis
            for analysis in pyin_config.keys():

                # Get the configuration for the current analysis
                an_config = pyin_config[analysis]

                # If the analysis was not requested
                if not an_config["run"]:

                    # Continue
                    continue
                
                # Get the options to run pyinteraph
                pyin_args = an_config["args"]
                
                # Set a sub-directory for the analysis
                an_dir = os.path.join(res_dir, analysis)


                #-------------- Modes/Imins/Corrections --------------#

            
                # If different analysis modes are requested
                if "modes" in an_config.keys():

                    # Set the 'modes' to the ones found
                    modes = an_config["modes"]

                # If different interaction strength cut-offs
                # were requested (acPSN)
                elif "imins" in an_config.keys():

                    # Set the 'modes' to the interaction
                    # strength cut-offs found 
                    modes = an_config["imins"]

                # If different corrections were requested
                # (cmPSN)
                elif "corrections" in an_config.keys():

                    # Set the 'modes' to the corrections
                    # found
                    modes = an_config["corrections"]

                # Otherwise
                else:

                    # Set it to an empty list containing 'None'
                    modes = [None]
                
                # For each analysis mode
                for mode in modes:
                    
                    # If it is not None (alternative analysis
                    # mode for the analysis were requested)
                    if mode is not None:

                        # Set the subdirectory for the mode
                        mode_dir = os.path.join(an_dir, str(mode))

                    # Otherwise
                    else:
                        
                        # Set the sub-directory for alternative
                        # analysis modes equal to the parent analysis
                        # sub-directory
                        mode_dir = an_dir
                    

                    #------------------- pyinteraph ------------------#


                    # For each distance cut-off  
                    for d_cut in an_config["dcuts"]:
                        
                        # Set the subdirectory for the cut-off 
                        d_cut_dir = os.path.join(mode_dir, str(d_cut))
                        
                        # Set the path to the pyinteraph output file
                        pyin_out = \
                            os.path.join(d_cut_dir,
                                         an_config["out"])
                        
                        # Run the PSN calculation with pyinteraph
                        # (it returns the path to the file containing
                        # the matrix representing the PSN)
                        pyin_mat = \
                            client.submit(util.run_pyinteraph,
                                          trj = sub_trj,
                                          top = pyin_top,
                                          ref = pyin_ref,
                                          analysis = analysis,
                                          mode = mode,
                                          co = d_cut,
                                          wd = d_cut_dir, 
                                          out = pyin_out,
                                          otherargs = pyin_args)


                        #---------------- filter_graph ---------------#

                        
                        # Get the network filtering options
                        fg_args = fg_config["args"]

                        # For each occurrence cut-off
                        for p_cut in fg_config["pcuts"]:
                            
                            # Set the sub-directory for the cut-off
                            p_cut_dir = \
                                os.path.join(d_cut_dir, str(p_cut))
                            
                            # Set the path to the output/log file
                            fg_out = \
                                os.path.join(p_cut_dir, 
                                             fg_config["out"])
                            
                            # Run the network filtering (it returns
                            # the path to the file containing the
                            # matrix representing the filtered PSN)
                            fg_mat = \
                                client.submit(util.run_filter_graph,
                                              mat = pyin_mat,
                                              out = fg_out,
                                              perco = p_cut,
                                              wd = p_cut_dir,
                                              otherargs = fg_args)


                            #------------- graph_analysis ------------#


                            # For each type of network analysis
                            for net_analysis in ga_config.keys():
                                
                                # Get the configuration of the
                                # network analysis
                                net_an_config = \
                                    ga_config[net_analysis]
                                
                                # Set the sub-directory for the
                                # network analysis
                                net_an_dir = \
                                    os.path.join(p_cut_dir,
                                                 net_analysis)
                                
                                # Set the path to the output/log file
                                ga_out = \
                                    os.path.join(net_an_dir,
                                                 net_an_config["out"])
                                
                                # Get the options for the current
                                # analysis
                                ga_args = net_an_config["args"]
                                
                                # Run the analysis     
                                futures.append(\
                                    client.submit(\
                                        util.run_graph_analysis,
                                            mat = fg_mat,
                                            ref = pyin_ref,
                                            out = ga_out,
                                            wd = net_an_dir,
                                            otherargs = ga_args))

        # Wait on the last active futures ('graph_analysis' futures,
        # since they have no futures depending on them)
        wait(futures)


if __name__ == "__main__":
    main()
