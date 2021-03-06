#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    pyinknife-aggregate.py
#
#    Aggregate data generated by pyinknife-run.
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
import os
import os.path
# third-party packages
import pandas as pd
import yaml

from . import util


def main():


    ######################### ARGUMENT PARSER #########################


    # description of the script
    description = "\nAggregate data generated " \
                  "by the PyinKnife pipeline.\n"
    
    # create the argument parser
    parser = argparse.ArgumentParser(description = description)
    
    # add arguments
    c_help = f"Configuration file used to run the pipeline. " \
             f"Default is {util.DEFCONFIG}."
    parser.add_argument("-c", "--configfile", \
                        type = str, \
                        default = util.DEFCONFIG, \
                        help = c_help)

    ca_help = f"Configuration file for the data aggregation. " \
              f"Default is {util.DEFCONFIGAGGR}."
    parser.add_argument("-ca", "--configfile-aggregate", \
                        type = str, \
                        default = util.DEFCONFIGAGGR, \
                        help = ca_help)

    d_help = "Directory where the pipeline was run."
    parser.add_argument("-d", "--rundir", \
                        type = str, \
                        required = True, \
                        help = d_help)

    od_help = "Directory where to save the output files."
    parser.add_argument("-od", "--outdir", \
                        type = str, \
                        required = True, \
                        help = od_help)
    
    firstccs_def = 5 
    firstccs_help = f"First # most populated connected components " \
                    f"to be considered. Default is {firstccs_def}."
    parser.add_argument("--firstccs", \
                        type = int, \
                        default = firstccs_def, \
                        help = firstccs_help)
    
    # parse the arguments
    args = parser.parse_args()
    # get single arguments
    configfile = util.get_abspath(args.configfile)
    configfile_aggregate = util.get_abspath(args.configfile_aggregate)
    rundir = util.get_abspath(args.rundir)
    outdir = util.get_abspath(args.outdir)
    firstccs = args.firstccs


    ########################## CONFIGURATION ##########################  


    # load the configuration
    config = yaml.safe_load(open(configfile, "r"))
    aggrconfig = yaml.safe_load(open(configfile_aggregate, "r"))
    
    # output file name for the hubs calculation
    outhubs = config["graph_analysis"]["hubs"]["out"]
    # output file name for the connected components calculation
    outcc = config["graph_analysis"]["ccs"]["out"]
    # directory name of the full trajectory results
    trjdirname = config["resampling"]["dirnames"]["trj"]

    # format specifications for the output dataframes
    extension = aggrconfig["output"]["extension"]
    sep = aggrconfig["output"]["sep"]

    # column names for the statistics to be computed
    meancol = aggrconfig["columns"]["mean"]
    mediancol = aggrconfig["columns"]["median"]
    mincol = aggrconfig["columns"]["min"]
    maxcol = aggrconfig["columns"]["max"]
    stdcol = aggrconfig["columns"]["std"]
    semcol = aggrconfig["columns"]["sem"]


    ######################### DATA AGGREGATION ########################


    # create a defaultdict defaulting to a pandas DataFrame
    # to store the results
    results = collections.defaultdict(pd.DataFrame)
    
    # get the running directory base name so that the paths reported
    # by os.walk will be rooted there
    rundirbasename = os.path.basename(rundir)
    
    # traverse the directory tree rooted in the running directory
    for path, dirs, files in os.walk(rundirbasename):
        
        # split the path into its components
        splitpath = path.split("/")
        
        # check if the path correspond to a leaf
        if dirs:
            # if not, keep traversing (results are stored in the leafs)
            continue
        
        # if the path is 5-directories deep
        if len(splitpath) == 6:
            # the analysis had only one possible mode
            trj, analysis, dcut, pcut, graphanalysis = \
                *splitpath[1:4], *splitpath[-2:]
            # build the key that will identify the results for
            # this analysis
            key = f"{analysis}_{dcut}_{pcut}_{graphanalysis}"
        
        # if the path is 6-directories deep
        elif len(splitpath) == 7:
            # the analysis had multiple modes
            trj, analysis, mode, dcut, pcut, graphanalysis = \
                *splitpath[1:5], *splitpath[-2:]
            # build the key that will identify the results for
            # this analysis
            key = f"{analysis}_{mode}_{dcut}_{pcut}_{graphanalysis}"

        else:
            # ignore the path since it does not contain the
            # data we are interested in
            continue
        
        # if this leaf corresponds to the results for hubs 
        if graphanalysis == "hubs":
            # get the output file path
            outfile = os.path.join(path, outhubs)
            # parse the output file
            result = util.parse_hubs_out(outfile = outfile)
        
        # if this leaf corresponds to the results for connected
        # components
        elif graphanalysis == "ccs":
            # get the output file path
            outfile = os.path.join(path, outcc)
            # parse the output file
            result = util.parse_cc_out(outfile = outfile, \
                                       firstccs = firstccs)
        
        # the name of the Series will be the name of the trajectory
        result.name = trj
        # add the result for the current trajectory as a new column
        # in the dataframe storing data for this analysis
        results[key][trj] = result


    # for each dataframe
    for outname, df in results.items():

        # set the output file path
        outaggrfile = os.path.join(outdir, outname + extension)
        # create the output directory if it does not exist yet
        os.makedirs(outdir, exist_ok = True)
        # select only those columns containing results from
        # resamplings to compute statistics
        dfresamplings = df.loc[:, df.columns != trjdirname]
        # put the column identifying the full trajectory first
        dffulltrj = df[trjdirname]
        df.drop(labels = [trjdirname], axis = 1, inplace = True)
        df.insert(0, trjdirname, dffulltrj)
        # compute statistics and add them as columns
        df[meancol] = dfresamplings.mean(axis = 1)
        df[mediancol] = dfresamplings.median(axis = 1)
        df[mincol] = dfresamplings.min(axis = 1)
        df[maxcol] = dfresamplings.max(axis = 1)
        df[stdcol] = dfresamplings.std(axis = 1)
        df[semcol] = dfresamplings.sem(axis = 1)
        # sort the index so that more populated connected components
        # and hubs with lower degrees go first
        df = df.sort_index()
        # write the output file
        df.to_csv(outaggrfile, sep = sep)


if __name__ == "__main__":
    main()