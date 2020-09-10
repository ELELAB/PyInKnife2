#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    pyinknife.py
#
#    PyInKnife pipeline for Protein Structure Network analysis
#    with PyInteraph.
#
#    Copyright (C) 2020 Juan Salamanca Viloria 
#                       <juan.salamanca.viloria@gmail.com> 
#                       Valentina Sora 
#                       <sora.valentina1@gmail.com>
#                       Elena Papaleo
#                       <elenap@cancer.dk>
#
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
from pathlib import Path
import subprocess
# dask
import dask
from dask import distributed
from distributed import fire_and_forget
from distributed.config import initialize_logging
import distributed.utils
# others
import MDAnalysis as mda
import numpy as np
import yaml


# to address a bug that resets the distributed.worker
# logger to WARNING level when a task is launched on
# the worker, no matter what the configuration was
def reset_worker_logger():
    """Utility function to reset a Dask logger handlers
    and level to desired values.
    """

    # new level
    NEWLEVEL = logging.INFO
    # get the logger
    logger = logging.getLogger("distributed.worker")
    # define the handlers to keep
    htokeep = [h for h in logger.handlers if type(h).__name__ == \
               distributed.utils.DequeHandler.__name__]
    # remove all the handlers
    for h in logger.handlers:
        logger.removeHandler(h)
    # add the handlers to keep
    for h in htokeep:
        # set the new level
        h.setLevel(NEWLEVEL)
        # add the handler to the logger
        logger.addHandler(h)
    # reset the logger level to the new level
    logger.setLevel(NEWLEVEL)
    # return the new logger
    return logger


############################# RUN COMMANDS ############################


def run_pyinteraph(trj, \
                   top, \
                   ref, \
                   analysis, \
                   mode, \
                   co, \
                   out, \
                   wd, \
                   otherargs):
    """Run the `pyinteraph` executable to generate the Protein
    Structure Network.
    """

    # reset the distributed.worker logger
    logger = reset_worker_logger()

    # make sure that the directory exists. If it does not, create it.
    os.makedirs(wd, exist_ok = True)
    # set the trajectory, topology, reference, analysis and
    # distance cut-off arguments
    args = ["pyinteraph", "--trj", trj, "--top", top, \
            "--ref", ref, "--" + analysis + "-co", str(co)]
    # add the analysis mode/class if the analysis is either
    # hydrogen bonds or salt bridges
    if analysis == "hb":
        args.extend(["--hb-class", mode])
    elif analysis == "sb":
        args.extend(["--sb-mode", mode])
    # launch the process
    p = subprocess.Popen(args = args + otherargs, \
                         stdout = open(out, "w"), \
                         stderr = subprocess.STDOUT, \
                         cwd = wd)
    # wait for the process to complete
    p.wait()
    # set the log string
    logstr = \
        f"pyinteraph run in {wd} exited with code {p.returncode}. " \
        f"Command line: {' '.join(p.args)}"   
    # if the command completed successfully
    if p.returncode == 0:
        # inform the user only for debug purposes
        logger.info(logstr)
    # if an error occurred
    else:
        # inform the user in any case 
        logger.error(logstr)
    # return the path to output matrix
    outmatrix = "--" + analysis + "-graph"
    return os.path.join(wd, otherargs[otherargs.index(outmatrix)+1])


def run_filter_graph(mat, perco, out, wd, otherargs):
    """Run the `filter_graph` executable to filter the Protein
    Structure Network.
    """

    # reset the distributed.worker logger
    logger = reset_worker_logger()

    # make sure that the directory exists. If it does not, create it.
    os.makedirs(wd, exist_ok = True)
    # set the matrix and persistence cutoff arguments
    args = ["filter_graph", "--input-dat", mat, \
            "--filter-threshold", str(perco)]
    # launch the process
    p = subprocess.Popen(args = args + otherargs, \
                         stdout = open(out, "w"), \
                         stderr = subprocess.STDOUT, \
                         cwd = wd)
    # wait for the process to complete
    p.wait()
    # set the log string
    logstr = \
        f"filter_graph run in {wd} exited with code {p.returncode}. " \
        f"Command line: {' '.join(p.args)}"   
    # if the command completed successfully
    if p.returncode == 0:
        # inform the user only for debug purposes
        logger.info(logstr)
    # if an error occurred
    else:
        # inform the user in any case 
        logger.error(logstr)
    # return the path to the output matrix
    outmatrix = "--output-dat"
    return os.path.join(wd, otherargs[otherargs.index(outmatrix)+1])


def run_graph_analysis(mat, ref, out, wd, otherargs):
    """Run the `graph_analysis` executable to analyze the
    Protein Structure Network.
    """

    # reset the distributed.worker logger
    logger = reset_worker_logger()

    # make sure that the directory exists. If it does not, create it.
    os.makedirs(wd, exist_ok = True)
    # set the adjacency matrix and the reference structure arguments
    args = ["graph_analysis", "--adj-matrix", mat, "--reference", ref]
    # launch the process
    p = subprocess.Popen(args = args + otherargs, \
                         stdout = open(out, "w"), \
                         stderr = subprocess.STDOUT, \
                         cwd = wd)
    # wait for the process to complete
    p.wait()
    # set the log string
    logstr = \
        f"graph_analysis run in {wd} exited with code {p.returncode}. " \
        f"Command line: {' '.join(p.args)}"   
    # if the command completed successfully
    if p.returncode == 0:
        # inform the user only for debug purposes
        logger.info(logstr)
    # if an error occurred
    else:
        # inform the user in any case 
        logger.error(logstr)
    # return the output file
    return out


######################### GET CONFIGURATIONS ##########################


def get_pyinknife_configuration(configfile):

    # recursively add the result of a function
    # taking as inputs entries from a dictionary as
    # keyword arguments as a new key in the parent
    # dictionary
    def recursively_add(data, \
                        keys, \
                        newkey, \
                        func):
        
        # if data is a dictionary
        if isinstance(data, dict):
            # for each key, value pair
            for k, v in list(data.items()):
                if k in keys:
                    # if value is a dictionary
                    if isinstance(v, dict):
                        # pass it as keyword arguments
                        data[newkey] = func(**v)
                    else:
                        # pass it as argument
                        data[newkey] = func(v)
                else:
                    # recursively check the
                    # sub-dictionaries
                    recursively_add(data = data[k], \
                                    keys = keys, \
                                    newkey = newkey, \
                                    func = func)

    # create a list of arguments made to be passed
    # to subprocess as arguments of a command
    def create_arguments_list(**kwargs):
        # initialize the list of arguments to an empty list
        args = []
        # for each option passed as keyword argument
        for opt, val in kwargs.items():
            # if the value of the argument is True
            if val is True:
                # only append the option name since it is
                # a flag
                args.append(opt)
            # if the value of the argument is False
            elif val is False:
                # ignore the option
                continue
            # if the value of the argument is a list
            elif isinstance(val, list):
                # convert each item into a string and join
                # them into a single string
                args.extend([opt, ",".join([str(v) for v in val])])
            # if the value is a string, int or float
            elif isinstance(val, (str, int, float)):
                # convert it to a string
                args.extend([opt, str(val)])
        # return the list of arguments
        return args

    # load the configuration
    data = yaml.safe_load(open(configfile, "r"))
    # convert all options given as key, value pairs to a command
    # to a list of command line arguments and add an extra entry
    # to the dictionary for each of these lists
    recursively_add(data = data, \
                    keys = {"options"}, \
                    newkey = "args", \
                    func = create_arguments_list)
    # return the new configuration
    return data


########################## READ/PROCESS DATA ##########################


def get_protein_selection_string(ncaa):
    """Select atoms belonging to protein residues in a Universe.
    """

    # reset the distributed.worker logger
    logger = reset_worker_logger()

    # the starting poin of the selection string is the 'protein'
    # keyword
    selstring = "protein"
    # if the system contains noncanonical residues
    if ncaa:
        # add the residue names to the selection string
        selstring += " or resname ".join(ncaa)
    # return the selection string
    return selstring


def get_frames_sets(trajectory, config):
    """Get the frame sets (full trajectory and possible
    subsets of frames) to perform the analysis on.
    """

    # reset the distributed.worker logger
    logger = reset_worker_logger()

    # utility function that gets intervals corresponding
    # to subsets of frames obtained by jackknife resampling
    def jackknife_resampling(nsamplings, ff, lf, step):
        intervals = []
        for ns in range(nsamplings):
            # get the frame intervals
            interval = ((ff, (step*ns)), (step*(ns+1), lf))
            # if the first chunck of data is excluded
            if ns == 0:
                # exclude the first half of the interval
                # since it is empty
                interval = (interval[1],)
            # if the last chunk of data is excluded
            elif ns == nsamplings-1:
                # exclude the second half of the interval
                # since it is either empty or contains extra
                # frames (in case of uneven splitting the last
                # chunk is the longest)
                interval = (interval[0],)
            # save the interval
            intervals.append(interval)
        # return the intervals
        return intervals 

    # get the directory names and the number of samplings
    # from the configuration
    dirnames = config["dirnames"]
    nsamplings = config["nsamplings"]

    # first frame and last frame are indexed from 0 (frame
    # attribute of MDAnalysis Timestep objects)
    ff, lf = trajectory[0].frame, trajectory[-1].frame
    step = int(np.floor(len(trajectory) / nsamplings))
    #print(step)
    # initialize the intervals to the interval corresponding
    # to the full trajectory
    intervals = [((ff, lf),)]
    # initialize the trajectory names to the name of the
    # full trajectory (name of the directory where analyses of
    # the full trajectory will be stored)
    trjnames = [dirnames["trj"]]
    # if no resampling requested
    if not config["run"]:
        # simply return the name and interval of the full trajectory
        return zip(trjnames, intervals)
    # if jackknife resampling has been requested
    if config["method"] == "jackknife":
        # get the intervals corresponding to the subsets of frames
        subintervals = jackknife_resampling(nsamplings, ff, lf, step)
        # get the names corresponding to the subsets of frames
        subtrjnames = []
        for ns in range(nsamplings):
            subtrjnames.append(\
                dirnames["subtrj"].replace(r"{nsampling}", str(ns)))
        # add intervals and names to the corresponding lists
        intervals.extend(subintervals)
        trjnames.extend(subtrjnames)
    # return names and intervals as a zipped list (cannot return
    # a generator/iterator when passing the function to a Client)
    return list(zip(trjnames, intervals))


############################# WRITE FILES #############################


def write_trj_subset(filename, trj, top, selstring, interval, wd):
    """Write a subset of a full trajectory containing only
    a subset of frames given as an interval. Can be used to
    write a trajectory or single frames.
    """

    # reset the distributed.worker logger
    logger = reset_worker_logger()
    
    # make sure that the directory exists. If it does not, create it.
    os.makedirs(wd, exist_ok = True) 
    # create a new Universe from the topology and the trajectory
    u = mda.Universe(top, trj)
    # select the atoms
    atoms = u.select_atoms(selstring)
    # get the trajectory
    trajectory = u.trajectory
    # set the file path
    filepath = os.path.join(wd, filename)
    # create the writer
    with mda.Writer(filepath, atoms.n_atoms) as Wset:
        # for each portion of the interval
        for start, end in interval:
            # for each frame in the trajectory
            for ts in trajectory[start:end]:
                Wset.write(atoms)
    # return the file path
    return filepath



if __name__ == "__main__":


    ######################### ARGUMENT PARSER #########################


    # description of the script
    description = "\nPyInKnife pipeline for Protein Structure " \
                  "Network Analysis with PyInteraph.\n"
    
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

    c_help = "Configuration file."
    parser.add_argument("-c", "--configfile", \
                        type = str, \
                        required = True, \
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
    trj = os.path.abspath(args.trj)
    top = os.path.abspath(args.top)
    ref = os.path.abspath(args.ref) if args.ref else None
    configfile = os.path.abspath(args.configfile)
    rundir = os.path.abspath(args.rundir)
    nproc = args.nproc
    ncaa = args.noncanonical_residues


    ########################## CONFIGURATION ##########################

    # get the configuration
    config = get_pyinknife_configuration(configfile = configfile)
    # get the logging configuration
    logging.basicConfig(level = logging.INFO)


    ######################### RUN THE PIPELINE ########################

  
    # change scheduler if you want to use threads or a single core
    with dask.config.set(scheduler = "processes"):

        # NB: Universe creation must be done internally by any
        # function requiring it since Universes are unpicklable
        # objects, therefore not serializable

        # create a local cluster with the desired number of workers
        cluster = distributed.LocalCluster(n_workers = nproc, \
                                           silence_logs = logging.INFO, \
                                           processes = True, \
                                           threads_per_worker = 1)
        
        # set a client to submit the jobs to
        client = distributed.Client(cluster)
        
        # get the configurations for the different runs
        # and for the resampling
        resconfig = config["resampling"]
        pyinconfig = config["pyinteraph"]
        fgconfig = config["filter_graph"]
        gaconfig = config["graph_analysis"]

        # create a MDAnalysis Universe object from the topology
        # and the trajectory
        u = mda.Universe(top, trj)
        # get the trajectory
        trajectory = u.trajectory
        # set the selection string for atoms
        selstring = client.submit(get_protein_selection_string, \
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
            pyintop = client.submit(write_trj_subset, \
                                    filename = config["pyintop"], \
                                    trj = trj, \
                                    top = top, 
                                    selstring = selstring,
                                    interval = ((0,1),), \
                                    wd = rundir)
            # create the file that will serve as reference
            # for PyInteraph
            pyinref = client.submit(write_trj_subset, \
                                    filename = config["pyinref"], \
                                    trj = trj, \
                                    top = top, \
                                    selstring = selstring, 
                                    interval = ((0,1),), \
                                    wd = rundir)

        # get the sets of frames composing the full trajectory
        # and possible resampled trajectories
        resamplings = client.submit(get_frames_sets, \
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
            elif i != 0:
                # write the sub-trajectory
                subtrj = client.submit(write_trj_subset, \
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
                        pyinmat = client.submit(run_pyinteraph, \
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
                            fgmat = client.submit(run_filter_graph, \
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
                                fire_and_forget(\
                                    client.submit(\
                                        run_graph_analysis, \
                                            mat = fgmat, \
                                            ref = pyinref, \
                                            out = gaout, \
                                            wd = netandir, \
                                            otherargs = gaargs))


