#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    util.py
#
#    Utility functions for the pyinknife scripts.
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
import collections
import logging
import os
import os.path
from pathlib import Path
import subprocess
# third-party packages
import dask
from dask import distributed
from distributed import fire_and_forget
from distributed.config import initialize_logging
from pkg_resources import resource_filename, Requirement
#import distributed.utils
import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
import pandas as pd
import seaborn as sns
import yaml

from .dask_patches import reset_worker_logger


# directory where the default configuration files are stored
CONFIGDIR = resource_filename(Requirement("PyInKnife"), \
                              "PyInKnife/config")

# default configuration file for pyinknife_run
DEFCONFIG = os.path.join(CONFIGDIR, "run.yaml")

# default configuration file for pyinknife_aggregate
DEFCONFIGAGGR = os.path.join(CONFIGDIR, "aggregate.yaml")

# default configuration files for pyinknife_plot
DEFCONFIGCCS = os.path.join(CONFIGDIR, "plot_ccs_barplot.yaml")
DEFCONFIGHUBS = os.path.join(CONFIGDIR, "plot_hubs_barplot.yaml")

# pyinteraph/filter_graph/graph_analysis options
# that should not be specified in the configuration
NOOPTSDICT = \
    {"pyinteraph": \
        {"-h", "--help", "-s", "--top", "-t", "--trj", \
         "-r", "--ref", "--sb-mode", "--hb-class", \
         "-hc-co", "--hc-cutoff", "-sb-co", "--sb-cutoff", \
         "-hb-co", "--hb-cutoff"}, \
    "filter_graph" : \
        {"-d", "--input-dat", "-t", "--filter-threshold"}, \
    "graph_analysis" : \
        {"-r", "--reference", "-a", "--adj-matrix"}}


########################### UTILITY FUNCTIONS #########################


def get_abspath(path):
    """Return the absolute path of a path if the path is not None,
    else return None.
    """
    
    return os.path.abspath(path) if path is not None else path



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


    # reset the distributed.worker logger
    logger = reset_worker_logger()

    # recursively perform an action on items of
    # a dictionary
    def recursively_traverse(data, \
                             action, \
                             keys, \
                             newkey = None, \
                             func = None):
        
        # if data is a dictionary
        if isinstance(data, dict):
            # for each key, value pair
            for k, v in list(data.items()):
                if k in keys:
                    if action == "pop":
                        data.pop(k)
                    # if value is a dictionary
                    if isinstance(v, dict):
                        if action == "add":
                            # pass it as keyword arguments
                            data[newkey] = func(**v)
                    else:
                        if action == "add":
                            # pass it as argument
                            data[newkey] = func(v)
                else:
                    # recursively check the
                    # sub-dictionaries
                    recursively_traverse(data = data[k], \
                                         action = action, \
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

    # for each set of undesired options
    for command, noopts in NOOPTSDICT.items():
        # pop them from the list of options of the corresponding
        # command
        recursively_traverse(data = data[command], \
                             action = "pop", \
                             keys = noopts)

    # convert all dictionaries storing options into lists of
    # options and add this list as an extra key
    recursively_traverse(data = data, \
                         action = "add", \
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
    # inform the user about the selection string
    logger.debug(f"Atom selection: {selstring}")
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


def parse_cc_out(outfile, firstccs):
    """Parse the output file containing the list of
    connected components.
    """
    
    # line that indicates the beginning of the list
    # of connected components
    START = "Connected component"
    # how data are separated in a line
    DATASEP = ")"
    # how the different nodes are separated in the same
    # connected component
    NODESEP = ","
    # index of the output Series
    INDEX = "cc_num"
    
    with open(outfile, "r") as f:
        # create a list to store the raw data
        rawdata = []
        # initialize the flag to start parsing to False
        parse = False
        # initialize the counter for connected components found to 0
        cccounter = 0
        # for each line 
        for l in f:
            # if you found the beginning of the list
            if l.startswith(START):
                # start parsing
                parse = True
                # go to the next line
                continue
            # if parsing is enabled the the counter is still below
            # the maximum number of connected components to be reported
            if parse and cccounter <= firstccs:
                # update the counter
                cccounter += 1
                # get the nodes composing the connected component
                cc = [i for i in l.split(DATASEP)[1].split(NODESEP) if i]
                # save the number of nodes in the current connected
                # component
                rawdata.append(len(cc))
        # sort the connected components by decreasing size and only
        # take the first 'firstccs' components
        ccs = sorted(rawdata, reverse = True)[:firstccs]
        # create a Series with the number of each connected
        # component associated to the number of nodes in it
        series = pd.Series(data = ccs, \
                           index = range(1, len(ccs)+1), \
                           dtype = np.int64)
        # rename the index of the Series and return it
        return series.rename_axis(index = INDEX)


def parse_hubs_out(outfile):
    """Parse the output file containing the list of hubs."""

    # beginning of the line that indicates that no hubs
    # have been found
    NOHUBS = "WARNING:root:No hubs"
    # line that indicates the beginning of the list of hubs    
    START = "Hubs:"
    # how data are separated in a line
    DATASEP = "\t"
    # index of the output Series
    INDEX = "k"
    
    with open(outfile, "r") as f:
        # create a default dictionary to store the hubs
        hubs = collections.defaultdict(int)
        # initialize the flag to start parsing to False       
        parse = False
        # for each line
        for l in f:
            # if no hubs were found
            if l.startswith(NOHUBS):
                # return an empty Series
                return pd.Series()
            # if you found the beginning of the list
            if l.startswith(START):
                # turn on the flag to start parsing
                parse = True
                # skip the next line
                next(f)
                # go to two lines after
                continue
            # if parsing is enabled
            if parse:
                # get the degree of the hub
                k = [i for i in l.rstrip("\n").split(DATASEP) if i][1]
                # update the counter for the hubs having the
                # same degree
                hubs[int(k)] += 1
        # create a Series with the degree of a hub associated
        # to the number of hubs having that degree (degree
        # distribution)
        series = pd.Series(data = hubs, dtype = np.int64)
        # rename the index of the Series and return it
        return series.rename_axis(index = INDEX)


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
    # inform the user about the file written
    logger.debug(f"Trajectory {filepath} written. " \
                 f"Frames: {start}-{end}.")
    # return the file path
    return filepath


############################ GENERATE PLOTS ###########################


def get_interval(values, config):
    """Get the numeric interval of a plot axis."""


    # reset the distributed.worker logger
    logger = reset_worker_logger()
    
    # default to an empty configuration if no configuration was
    # provided for the interval
    inconfig = config.get("interval", {})
    # default to a discrete interval if the type of interval
    # was not provided
    inttype = inconfig.get("type", "discrete")

    # if it is a discrete interval
    if inttype == "discrete":
        # default to rounding to the nearest integer
        defrtn = 1
        # default top value is the maximum of the values provided
        deftop = int(max(values))
        # default bottom value is the minimum of the values
        # provided
        defbottom = int(min(values))
        # default number of steps is the number of values provided
        defsteps = len(values)
        # default spacing is the one between two steps
        defspacing = int(np.linspace(defbottom, \
                                     deftop, \
                                     defsteps, \
                                     retstep = True)[1])
        # the interval is not centered in zero by default
        defciz = False

    # if it is a continuous interval
    elif inttype == "continuous":
        # default to rounding to the nearest 0.5
        defrtn = 0.5
        # default top value is the rounded up maximum of the values
        deftop = np.ceil(max(values)*(1/defrtn))/(1/defrtn)
        # default bottom value is the rounded down minimum of
        # the values
        defbottom = np.floor(min(values)*(1/defrtn))/(1/defrtn)
        # default is 5 steps 
        defsteps = 5
        # default spacing is the one between two steps
        defspacing = np.linspace(defbottom, \
                                 deftop, \
                                 defsteps, \
                                 retstep = True)[1]
        # the interval is not centered in zero by default
        defciz = False

    # get the configurations
    rtn = inconfig.get("round_to_nearest", defrtn)
    top = inconfig.get("top", deftop)
    bottom = inconfig.get("bottom", defbottom)
    steps = inconfig.get("steps", defsteps)
    spacing = inconfig.get("spacing", defspacing)
    ciz = inconfig.get("center_in_zero", defciz)
    
    # if the interval should be centered in zero
    if ciz:
        # get the absolute values of top and bottom
        abstop, absbottom = abs(top), abs(bottom)
        # get the highest absolute value
        absval = abstop if abstop > absbottom else absbottom
        # top and bottom will be opposite numbers with
        # absolute value equal to absval
        top, bottom = absval, -absval

    # return the interval
    return np.arange(bottom, top + spacing, spacing)


def set_barplot(ax, \
                df, \
                y, \
                yerr, \
                config):
    """Set up a bar plot."""
    
    # number of bars to be drawn
    numdata = len(df)
    # get as many colors from the color palette as the
    # number of bars
    color = sns.color_palette(config["palette"])[:numdata]
    # plot the bars
    bars = ax.bar(x = range(numdata), \
                  height = df[y], \
                  color = color, \
                  **config["bars"])
    # plot error bars
    ax.errorbar(x = range(numdata), \
                y = df[y], \
                yerr = df[yerr], \
                **config["yerr"])
    # return the axis
    return ax


def set_axis(ax, \
             axis, \
             config, \
             ticks = None, \
             ticklabels = None):
    """Set up the x- or y-axis."""
    
    # get the tick locations already on the axis 
    ticklocs = plt.xticks()[0]
    # default to the tick locations already present if
    # no custom ones were provided
    ticks = ticks if ticks is not None else ticklocs
    # default to the string representations of the tick
    # locations as labels if no custom labels were provided
    ticklabels = ticklabels if ticklabels is not None \
                 else [str(t) for t in ticks]
    
    # if it is the x-axis
    if axis == "x":
        # set the axis label
        ax.set_xlabel(**config["label"])
        # set the ticks
        ax.set_xticks(ticks = ticks)
        # set the tick labels
        ax.set_xticklabels(labels = ticklabels, \
                           **config["ticklabels"])
        # set the axis boundaries
        ax.spines["bottom"].set_bounds(ticks[0], ticks[-1])
    
    # if it is the y-axis
    elif axis == "y":
        # set the axis label
        ax.set_ylabel(**config["label"])
        # set the ticks
        ax.set_yticks(ticks = ticks)
        # set the tick labels
        ax.set_yticklabels(labels = ticklabels, \
                           **config["ticklabels"])
        # set the axis boundaries
        ax.spines["left"].set_bounds(ticks[0], ticks[-1])

    # if a configuration for the tick parameters was provided
    if config.get("tick_params"):
        # apply the configuration to the ticks
        ax.tick_params(axis = axis, \
                       **config["tick_params"])

    # return the axis
    return ax