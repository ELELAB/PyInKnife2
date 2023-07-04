#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    util.py
#
#    Utility functions for the PyInKnife2 executables.
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
import collections
import logging
import os
import os.path
from pathlib import Path
import subprocess
# Third-party packages
import dask
from dask import distributed
from distributed import fire_and_forget
from distributed.config import initialize_logging
from pkg_resources import resource_filename, Requirement
import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
import pandas as pd
import seaborn as sns
import yaml
# PyInKnife
from .dask_patches import reset_worker_logger


#-------------------------- Default values ---------------------------#


# Directory where the default configuration files are stored
CONFIG_DIR = \
    resource_filename(Requirement("PyInKnife"), \
                      "PyInKnife/config")

# Default configuration file for pyinknife_run
DEF_CONFIG = os.path.join(CONFIG_DIR, "run.yaml")

# Default configuration file for pyinknife_aggregate
DEF_CONFIG_AGGR = os.path.join(CONFIG_DIR, "aggregate.yaml")

# Default configuration files for pyinknife_plot
DEF_CONFIG_CCS = os.path.join(CONFIG_DIR, "plot_ccs_barplot.yaml")
DEF_CONFIG_HUBS = os.path.join(CONFIG_DIR, "plot_hubs_barplot.yaml")

# pyinteraph/filter_graph/graph_analysis options
# that should not be specified in the configuration
NO_OPTS_DICT = \
    {"pyinteraph": \
        {"-h", "--help",
         "-s", "--top", "-t", "--trj", "-r", "--ref",
         "--cmpsn-correction",
         "--acpsn-imin",
         "--sb-mode",
         "--hb-class",
         "--cmpsn-co", "--cmpsn-cutoff",
         "--acpsn-co", "--acpsn-cutoff",
         "--hc-co", "--hc-cutoff",
         "--sb-co", "--sb-cutoff",
         "--hb-co", "--hb-cutoff"},
    "filter_graph" : \
        {"-d", "--input-dat", "-t", "--filter-threshold"},
    "graph_analysis" : \
        {"-r", "--reference", "-a", "--adj-matrix"}}


#---------------------------- Miscellanea ----------------------------#


def get_abspath(path):
    """Return the absolute path of a path if the path is not None,
    else return None.
    """

    # If the path is None
    if path is None:

        # Return None
        return path

    # If the path is a symbolic link
    if os.path.islink(path):

        # The path is the one the symbolic link points to
        path = os.readlink(path)

    # Return the absolute path
    return os.path.abspath(path)


#--------------------------- Run utilities ---------------------------#


def run_pyinteraph(trj,
                   top,
                   ref,
                   analysis,
                   mode,
                   co,
                   out,
                   wd,
                   otherargs):
    """Run the ``pyinteraph`` executable to generate the PSN.
    """


    #---------------------------- Logging ----------------------------#

    # Reset the distributed.worker logger
    logger = reset_worker_logger()

    #----------------------- Working directory -----------------------#

    # Make sure that the working directory exists. If it does not,
    # create it.
    os.makedirs(wd,
                exist_ok = True)

    #--------------------------- Arguments ---------------------------#

    # Set the trajectory, topology, reference, analysis and
    # distance cut-off arguments
    args = \
        ["pyinteraph", "--trj", trj, "--top", top,
         "--ref", ref, f"--{analysis}-co", str(co)]

    # Add the analysis mode/class if the analysis requires it
    if analysis == "cmpsn":
        args.extend(["--cmpsn-correction", mode])
    elif analysis == "acpsn":
        args.extend(["--acpsn-imin", str(mode)])
    elif analysis == "hb":
        args.extend(["--hb-class", mode])
    elif analysis == "sb":
        args.extend(["--sb-mode", mode])

    #---------------------------- Process ----------------------------#

    # Launch the process
    p = subprocess.Popen(args = args + otherargs,
                         stdout = open(out, "w"),
                         stderr = subprocess.STDOUT,
                         cwd = wd)
    
    # Wait for the process to complete
    p.wait()
    
    # Set the log string
    logstr = \
        f"pyinteraph run in {wd} exited with code {p.returncode}. " \
        f"Command line: {' '.join(p.args)}" 
    
    # If the command completed successfully
    if p.returncode == 0:
        
        # Log as an information
        logger.info(logstr)
    
    # If an error occurred
    else:
        
        # Log as an error
        logger.error(logstr)

    #------------------------- Output matrix -------------------------#
    
    # Get the name of the option defining the output matrix
    out_matrix = "--" + analysis + "-graph"
    
    # Get the matrix full name
    out_matrix_fullname = otherargs[otherargs.index(out_matrix)+1]
    
    # Divide the matrix name from the file extension
    out_matrix_name, out_matrix_ext = out_matrix_fullname.split(".")
    
    # Return the path to output matrix
    return os.path.join(wd,
                        out_matrix_name + "_all." + out_matrix_ext)


def run_filter_graph(mat,
                     perco,
                     out,
                     wd,
                     otherargs):
    """Run the ``filter_graph`` executable to filter the PSN.
    """

    #--------------------------- Logging -----------------------------#

    # Reset the distributed.worker logger
    logger = reset_worker_logger()

    #----------------------- Working directory -----------------------#

    # Make sure that the directory exists. If it does not, create it.
    os.makedirs(wd,
                exist_ok = True)

    #--------------------------- Arguments ---------------------------#

    # set the matrix and persistence cutoff arguments
    args = \
        ["filter_graph", "--input-dat", mat,
         "--filter-threshold", str(perco)]

    #---------------------------- Process ----------------------------#
    
    # Launch the process
    p = subprocess.Popen(args = args + otherargs,
                         stdout = open(out, "w"),
                         stderr = subprocess.STDOUT,
                         cwd = wd)

    # Wait for the process to complete
    p.wait()
    
    # Set the log string
    logstr = \
        f"filter_graph run in {wd} exited with code {p.returncode}. " \
        f"Command line: {' '.join(p.args)}"   
    
    # If the command completed successfully
    if p.returncode == 0:
        
        # Log as an information
        logger.info(logstr)
    
    # If an error occurred
    else:
        
        # Log as an error
        logger.error(logstr)

    #------------------------- Output matrix -------------------------#
    
    # Set the option containing the output matrix
    out_matrix_opt = "--output-dat"
    
    # Return the path to the output matrix
    return os.path.join(wd,
                        otherargs[otherargs.index(out_matrix_opt)+1])


def run_graph_analysis(mat,
                       ref,
                       out,
                       wd,
                       otherargs):
    """Run the ``graph_analysis`` executable to analyze the PSN.
    """

    #---------------------------- Logging ----------------------------#

    # Reset the distributed.worker logger
    logger = reset_worker_logger()

    #----------------------- Working directory -----------------------#

    # Make sure that the directory exists. If it does not, create it.
    os.makedirs(wd,
                exist_ok = True)

    #--------------------------- Arguments ---------------------------#

    # Set the adjacency matrix and the reference structure arguments
    args = \
        ["graph_analysis", "--adj-matrix", mat, "--reference", ref]

    #---------------------------- Process ----------------------------#

    # Launch the process
    p = subprocess.Popen(args = args + otherargs,
                         stdout = open(out, "w"),
                         stderr = subprocess.STDOUT,
                         cwd = wd)
    
    # Wait for the process to complete
    p.wait()
    
    # Set the log string
    logstr = \
        f"graph_analysis run in {wd} exited with code {p.returncode}. " \
        f"Command line: {' '.join(p.args)}"   
    
    # If the command completed successfully
    if p.returncode == 0:
        
        # Log as an information
        logger.info(logstr)
    
    # If an error occurred
    else:
        
        # Log as an error
        logger.error(logstr)

    #-------------------------- Output file --------------------------#
    
    # Return the output file
    return out


#--------------------------- Configuration ---------------------------#


def get_pyinknife_configuration(config_file):
    """Get the configuration to run the PyInKnife pipeline.
    """

    #---------------------------- Logging ----------------------------#

    # Reset the distributed.worker logger
    logger = reset_worker_logger()

    #----------------------- Helper functions ------------------------#

    # Recursively perform an action on items of a dictionary
    def recursive_traverse(data,
                           action,
                           keys,
                           new_key = None,
                           func = None):
        
        # If 'data' is a dictionary
        if isinstance(data, dict):
            
            # For each key, value pair
            for k, v in list(data.items()):

                # If the current key is in the list of key of
                # interest
                if k in keys:

                    # If the 'remove' action was selected
                    if action == "remove":

                        # Remove the key and associated value
                        # from the dictionary 
                        data.pop(k)

                    # If 'v' is a dictionary
                    if isinstance(v, dict):
                        
                        # If the 'add' action was requested
                        if action == "add":

                            # Add a new key 'new_key' whose
                            # associated value is the return
                            # value of the 'func' function,
                            # which takes the key, value pairs
                            # in 'v' as keyword arguments                            
                            data[new_key] = func(**v)

                    # If 'v' is not a dictionary
                    else:

                        # If the 'add' action was requested
                        if action == "add":

                            # Add a new key 'new_key' whose
                            # associated value is the return
                            # value of the 'func' function,
                            # which takes 'v' as a positional
                            # argument
                            data[new_key] = func(v)

                # If the current key is not in the list of
                # keys of interest 
                else:
                    
                    # Recursively check the sub-dictionaries
                    recursive_traverse(data = data[k],
                                       action = action,
                                       keys = keys,
                                       new_key = new_key,
                                       func = func)

    # Create a list of arguments made to be passed
    # to subprocess.Popen as arguments of a command
    def create_arguments_list(**kwargs):
        
        # Initialize the list of arguments to an empty list
        args = []
        
        # For each option passed as a keyword argument
        for opt, val in kwargs.items():

            # If the value of the argument is True (we use the
            # 'is' operator instead of '==' because we want to
            # check if 'val' is the boolean 'True', not just
            # whether it evaluates to 'True')
            if val is True:
                
                # Only append the option name since it is a
                # boolean flag
                args.append(opt)

            # If the value of the argument is False
            elif val is False:

                # Ignore the option
                continue

            # If the value of the argument is a list
            elif isinstance(val, list):
                
                # Convert each item into a string, join them
                # into a single string, and add the argument
                # and its value to the list of arguments
                args.extend(\
                    [opt, ",".join([str(v) for v in val])])
            
            # If the value is a string, int or float
            elif isinstance(val, (str, int, float)):
                
                # Convert it to a string and add the argument
                # and its value the list of arguments
                args.extend([opt, str(val)])
        
        # Return the list of arguments
        return args

    #--------------------- Configuration loading ---------------------#

    # Load the configuration
    config = yaml.safe_load(open(config_file, "r"))

    #----------------------- Undesired options -----------------------#

    # For each set of undesired options
    for command, no_opts in NO_OPTS_DICT.items():
        
        # Remove them from the list of options of the corresponding
        # command
        recursive_traverse(data = config[command],
                           action = "remove",
                           keys = no_opts)

    #----------------------- List of arguments -----------------------#

    # Convert all dictionaries storing options into lists of
    # arguments and add this list as an extra key
    recursive_traverse(data = config,
                       action = "add",
                       keys = {"options"},
                       new_key = "args",
                       func = create_arguments_list)

    # Return the updated configuration
    return config


#------------------------- Read/process data -------------------------#


def get_protein_selection_string(ncaa):
    """Select atoms belonging to protein residues in a Universe.
    """

    #---------------------------- Logging ----------------------------#

    # Reset the distributed.worker logger
    logger = reset_worker_logger()

    #----------------------- Selection string ------------------------#

    # The starting point of the selection string is the 'protein'
    # keyword
    sel_string = "protein"
    
    # If the system contains noncanonical residues
    if ncaa:
        
        # Add the residue names to the selection string
        sel_string += " or resname ".join(ncaa)
    
    # Inform the user about the selection string
    logger.info(f"Atom selection: '{sel_string}'.")
    
    # Return the selection string
    return sel_string


def get_frames_sets(trajectory, config):
    """Get the frame sets (full trajectory and possible
    subsets of frames) to perform the analysis on.
    """

    #---------------------------- Logging ----------------------------#
    
    # Reset the distributed.worker logger
    logger = reset_worker_logger()

    #----------------------- Helper functions ------------------------#

    # Utility function that gets intervals corresponding
    # to subsets of frames obtained by jackknife resampling
    def jackknife_resampling(n_samplings,
                             ff,
                             lf,
                             step):
        
        # Create an empty list to store the intervals        
        intervals = []

        # For each sampling  
        for ns in range(n_samplings):

            # Get the frame interval
            interval = ((ff, (step*ns)), (step*(ns+1), lf))
            
            # If the first chunck of data is excluded
            if ns == 0:

                # Exclude the first half of the interval
                # since it is empty
                interval = (interval[1],)

            # If the last chunk of data is excluded
            elif ns == n_samplings-1:

                # Exclude the second half of the interval
                # since it is either empty or contains extra
                # frames (in case of uneven splitting the last
                # chunk is the longest)
                interval = (interval[0],)
            
            # Append the interval to the list
            intervals.append(interval)
        
        # Return the intervals
        return intervals

    #------------------------- Configuration -------------------------#

    # Get the directory names and the number of samplings
    # from the configuration
    dir_names = config["dirnames"]
    n_samplings = config["nsamplings"]

    # Available resampling methods
    AVALILABLE_RESAMPLING_METHODS = ["jackknife"]

    #--------------------------- Intervals ---------------------------#

    # The first and last frame are indexed from 0 (frame
    # attribute of MDAnalysis Timestep objects)
    ff, lf = trajectory[0].frame, trajectory[-1].frame

    # Get the step size
    step = int(np.floor(len(trajectory) / n_samplings))

    # Initialize the intervals as a list containing as the
    # first and only element to the interval corresponding
    # to the full trajectory
    intervals = [((ff, lf),)]

    # Initialize the trajectory name as a list containing
    # as the first and only element the name of the full
    # trajectory
    trj_names = [dir_names["trj"]]

    # If no resampling is requested
    if not config["run"]:

        # Simply return the name and interval of the full trajectory
        return zip(trj_names, intervals)

    # If the resampling method is not among those available
    if config["method"] not in AVALILABLE_RESAMPLING_METHODS:

        # Raise an error
        available_methods_str = \
            ", ".join(f"'{m}'" for m in AVALILABLE_RESAMPLING_METHODS)
        errstr = \
            f"The '{config['method']}' resampling method is " \
            f"not supported. So far, PyInKnife supports the " \
            f"following resampling methods: {available_methods_str}."
        raise ValueError(errstr)

    # If jackknife resampling was requested
    if config["method"] == "jackknife":

        # Get the intervals corresponding to the subsets of frames
        sub_intervals = \
            jackknife_resampling(n_samplings = n_samplings,
                                 ff = ff,
                                 lf = lf,
                                 step = step)
        
        # Get the names corresponding to the subsets of frames
        sub_trj_names = []
        for ns in range(n_samplings):
            sub_trj_names.append(\
                dir_names["subtrj"].replace(r"{nsampling}", str(ns)))
        
        # Add intervals and names to the corresponding lists
        intervals.extend(sub_intervals)
        trj_names.extend(sub_trj_names)

    # Return names and intervals as a zipped list (cannot return
    # a generator/iterator when passing the function to a Client)
    return list(zip(trj_names, intervals))


def parse_cc_out(out_file,
                 first_ccs):
    """Parse the output file containing the list of
    connected components.
    """

    #--------------------------- Constants ---------------------------#
    
    # Line that indicates the beginning of the list of connected
    # components
    START = "Connected component"
    
    # Start of the line where the nodes belonging to the connected
    # component are stored
    DATA_START = " ("
    
    # How data are separated in a line
    DATA_SEP = ")"
    
    # How the different nodes are separated in the same connected
    # component
    NODE_SEP = ","
    
    # Index of the output Series
    INDEX = "cc_num"

    #------------------------- File parsing --------------------------#
    
    # Open the output file
    with open(out_file, "r") as f:
        
        # Initialize an empty list to store the raw data
        raw_data = []

        # Initialize the flag to start parsing to False
        parse = False
        
        # For each line 
        for l in f:

            # If we found the beginning of the list
            if l.startswith(START):

                # Turn on the flag to start parsing
                parse = True

                # Go to the next line
                continue

            # If parsing is enabled and the line is the
            # first line containing data of interest
            if parse and l.startswith(DATA_START):
                
                # Get the nodes in the connected component
                cc = \
                    [i for i in l.split(DATA_SEP)[1].split(NODE_SEP) \
                     if i]

                # Save the number of nodes in the current connected
                # component
                raw_data.append(len(cc))
        
        # Sort the connected components by decreasing size and only
        # take the first 'first_ccs' components
        ccs = sorted(raw_data,
                     reverse = True)[:first_ccs]

        # Create a Series with the number of each connected
        # component associated to the number of nodes in it
        series = pd.Series(data = ccs,
                           index = range(1, len(ccs)+1),
                           dtype = np.int64)

        # Rename the index of the Series
        series = series.rename_axis(index = INDEX)

        # Return the Series
        return series


def parse_hubs_out(out_file):
    """Parse the output file containing the list of hubs.
    """

    #--------------------------- Constants ---------------------------#

    # Beginning of the line that indicates that no hubs
    # have been found
    NO_HUBS = "WARNING:root:No hubs"
    
    # Line that indicates the beginning of the list of hubs    
    START = "Hubs:"
    
    # How data are separated in a line
    DATA_SEP = "\t"
    
    # index of the output Series
    INDEX = "k"

    #------------------------- File parsing --------------------------#
    
    with open(out_file, "r") as f:
        
        # Create a default dictionary to store the hubs
        hubs = collections.defaultdict(int)
        
        # Set the flag indicating when to start parsing to False     
        parse = False

        # For each line
        for l in f:

            # If no hubs were found
            if l.startswith(NO_HUBS):
                
                # Return an empty Series
                return pd.Series()
            
            # If you found the beginning of the list
            if l.startswith(START):
                
                # Turn on the flag to start parsing
                parse = True
                
                # Skip the next line
                next(f)
                
                # Go to two lines after
                continue

            # If parsing is enabled
            if parse:
                
                # Get the degree of the hub
                k = [i for i in l.rstrip("\n").split(DATA_SEP) if i][1]
                
                # Update the counter for the hubs having the
                # same degree
                hubs[int(k)] += 1
        
        # Create a Series with the hub degreeassociated
        # to the number of hubs having that degree (degree
        # distribution)
        series = pd.Series(data = hubs,
                           dtype = np.int64)
        
        # Rename the index of the Series
        series = series.rename_axis(index = INDEX)

        # Return the Series
        return series


#---------------------------- Write files ----------------------------#


def write_trj_subset(filename,
                     trj,
                     top,
                     sel_string,
                     interval,
                     wd):
    """Write a subset of a full trajectory containing only
    a subset of frames given as an interval. Can be used to
    write a trajectory or single frames.
    """

    #---------------------------- Logging ----------------------------#

    # reset the distributed.worker logger
    logger = reset_worker_logger()
    
    #----------------------- Working directory -----------------------#

    # Make sure that the directory exists. If it does not, create it.
    os.makedirs(wd,
                exist_ok = True)

    #---------------------------- System -----------------------------#

    # Create a new Universe from the topology and the trajectory
    u = mda.Universe(top, trj)
    
    # Select the atoms
    atoms = u.select_atoms(sel_string)
    
    # Get the trajectory
    trajectory = u.trajectory
    
    # Set the path to the file that will contain the trajectory
    filepath = os.path.join(wd, filename)

    # Create the writer
    with mda.Writer(filepath, atoms.n_atoms) as Wset:
        
        # For each portion of the interval
        for start, end in interval:
            
            # For each frame in the interval
            for ts in trajectory[start:end]:

                # Write out the coordinates
                Wset.write(atoms)
    
    # Inform the user about the file written
    infostr = \
        f"Trajectory '{filepath}' written. Frames: {start}-{end}."
    logger.info(infostr)

    # Return the file path
    return filepath


#-------------------------- Generate plots ---------------------------#


def get_axis_interval(values,
                      config):
    """Get the numeric interval of a plot axis.
    """

    #---------------------------- Values -----------------------------#

    # Convert the values to a NumPy array
    values = np.array(values)

    # Filter out NaN values
    values = values[np.logical_not(np.isnan(values))]

    #------------------------- Configuration -------------------------#

    # Get the top configuration
    config = config.get("interval")

    # If there is no configuration for the interval
    if config is None:

        # Raise an error
        errstr = \
            f"No 'interval' section was found in the " \
            f"configuration for the {item}."
        raise KeyError(errstr)

    # Get the configurations
    int_type = config.get("type")
    rtn = config.get("round_to_nearest")
    top = config.get("top")
    bottom = config.get("bottom")
    steps = config.get("steps")
    spacing = config.get("spacing")
    ciz = config.get("center_in_zero")

    #--------------------------- Rounding ---------------------------#

    # If no rounding was specified
    if rtn is None:

        # If the interval is discrete
        if int_type == "discrete":

            # Default to rounding to the nearest 1
            rtn = 1

        # If the interval is continuous
        elif int_type == "continuous":
        
            # Default to rounding to the nearest 0.5
            rtn = 0.5

    #----------------------- Top/bottom values -----------------------#

    # If the maximum of the ticks interval was not specified
    if top is None:

        # If the interval is discrete
        if int_type == "discrete":
            
            # The default top value will be the
            # maximum of the values provided
            top = int(max(values))
        
        # If the interval is continuous
        elif int_type == "continuous":
            
            # The default top value will be the
            # rounded-up maximum of the values
            top = \
                np.ceil(max(values)*(1/rtn)) / (1/rtn)

    # If the minimum of the ticks interval was not specified
    if bottom is None:
        
        # If the interval is discrete
        if int_type == "discrete":
            
            # The default bottom value is the
            # minimim of the values provided
            bottom = int(min(values))

        # If the interval is continuous
        elif int_type == "continuous":
            
            # The default bottom value is the rounded
            # down minimum of the values
            bottom = \
                np.floor(min(values)*(1/rtn)) / (1/rtn)

    # If the two extremes of the interval coincide
    if top == bottom:
        
        # Return only one value
        return np.array([bottom])

    #----------------------------- Steps -----------------------------#

    # If the number of steps the interval should have
    # was not specified
    if steps is None:

        # The default number of steps is 10
        steps = 10

    #---------------------------- Spacing ----------------------------#

    # If the interval spacing was not specified
    if spacing is None:
        
        # If the interval is discrete
        if int_type == "discrete":

            # The default spacing is the one between two steps,
            # rounded up
            spacing = \
                int(np.ceil(np.linspace(bottom,
                                        top,
                                        steps,
                                        retstep = True)[1]))

        # If the interval is continuous
        elif int_type == "continuous":
            
            # The default spacing is the one between two steps
            spacing = round(np.linspace(bottom,
                                        top,
                                        steps,
                                        retstep = True)[1], 2)

    #------------------------ Center in zero -------------------------#

    # If the interval should be centered in zero
    if ciz:
        
        # Get the highest absolute value
        absval = \
            np.ceil(top) if top > bottom else np.floor(bottom)
        
        # Top and bottom will be opposite numbers with
        # absolute value equal to absval
        top, bottom = absval, -absval

        # Get an evenly-spaced interval between the bottom
        # and top value
        interval = np.linspace(bottom, top, steps)
        
        # Return the interval
        return interval

    # Get the interval
    interval = np.arange(bottom, top + spacing, spacing)

    # Return the interval
    return interval


def set_barplot(ax,
                df,
                y,
                yerr,
                config):
    """Set up a bar plot.
    """
    
    # Get the number of bars to be drawn
    num_bars = len(df)
    
    # Get as many colors from the color palette as the
    # number of bars
    color = sns.color_palette(config["palette"])[:num_bars]
    
    # Plot the bars
    bars = ax.bar(x = range(num_bars),
                  height = df[y],
                  color = color,
                  **config["bars"])
    
    # Plot error bars
    ax.errorbar(x = range(num_bars),
                y = df[y],
                yerr = df[yerr],
                **config["yerr"])
    
    # Return the axis
    return ax


def set_axis(ax,
             axis,
             config,
             ticks = None,
             ticklabels = None):
    """Set up the x- or y-axis.
    """
    
    # Default to the tick locations already present if
    # no custom ones were provided
    ticks = \
        ticks if ticks is not None else plt.xticks()[0]
    
    # Default to the string representations of the tick
    # locations as labels if no custom labels were provided
    ticklabels = \
        ticklabels if ticklabels is not None \
        else [str(t) for t in ticks]
    
    # If it is the x-axis
    if axis == "x":
        
        # Set the axis label
        ax.set_xlabel(**config["label"])
        
        # Set the ticks
        ax.set_xticks(ticks = ticks)
        
        # Set the tick labels
        ax.set_xticklabels(labels = ticklabels,
                           **config["ticklabels"])
        
        # Set the axis boundaries
        ax.spines["bottom"].set_bounds(ticks[0], ticks[-1])
    
    # If it is the y-axis
    elif axis == "y":
        
        # Set the axis label
        ax.set_ylabel(**config["label"])
        
        # Set the ticks
        ax.set_yticks(ticks = ticks)
        
        # Set the tick labels
        ax.set_yticklabels(labels = ticklabels, \
                           **config["ticklabels"])
        
        # Set the axis boundaries
        ax.spines["left"].set_bounds(ticks[0], ticks[-1])

    # If a configuration for the tick parameters was provided
    if config.get("tick_params"):
        
        # Apply the configuration to the ticks
        ax.tick_params(axis = axis,
                       **config["tick_params"])

    # Return the axis
    return ax