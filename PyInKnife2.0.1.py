#!/usr/bin/env python
# -*- coding: utf-8 -*-

#    Pipeline for graph analysis with PyInteraph [1]_, a software suite 
#    to analyze interactions and interaction networks in
#    structural ensembles.
#
#    Copyright (C) 2016 Juan Salamanca Viloria 
#                        <juan.salamanca.viloria@gmail.com>, 
#                        Elena Papaleo 
#                        <elenap@cancer.dk>
#    Copyright (C) 2018 Juan Salamanca Viloria
#                        <juan.salamanca.viloria@gmail.com>, 
#                        Elena Papaleo 
#                        <elenap@cancer.dk>
#                        Valentina Sora
#                        <vaso@cancer.dk>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#    
#
#    .. [1] Tiberti, Matteo, et al. "PyInteraph: a framework for the
#       analysis of interaction networks in structural ensembles of 
#       proteins." Journal of chemical information and modeling 
#       54.5 (2014): 1537-1551.

import argparse
import sys
import os
import os.path
import subprocess
import multiprocessing
import logging
import functools
import re

import numpy as np
import MDAnalysis as mda

__author__ = ["Valentina Sora", "Juan Salamanca Viloria", "Elena Papaleo"]
__credits__ = ["Valentina Sora", "Juan Salamanca Viloria", "Elena Papaleo"]
__license__ = "GPL"
__version__ = "2.0.1"
__date__ = "12/6/2019"
__maintainer__ = "Valentina Sora"
__email__ = "vaso@cancer.dk"
__status__ = "Development"

class PyInteraphCaller:    
    """Class implementing a wrapper for calling `PyInteraph`
    tools.
    So far, the following PyInteraph commands are supported:
    -   `pyinteraph`
    -   `filter_graph`
    -   `graph_analysis`

    Attributes
    ----------
    ARGS : `dict`
        Dictionary containing all arguments available
        for each command, divided according to their
        being flags, double-dash flags, single-dash 
        or double-dash arguments.
    """

    ARGS = \
        {
         "pyinteraph" : \
            {
                "FLAGS" : \
                ["h", "v", "b", "f", "y", "p"], \

                "DDFLAGS" : \
                ["help", "verbose", "salt_bridges", \
                 "hydrophobic", "hydrogen_bonds", "potential"], \

                "SDARGS" : \
                    ["s", "t", "r"], \

                "DDARGS" : \
                    ["top", "trj", "ref", \
                     "sb_co", "sb_cutoff", "sb_perco", \
                     "sb_persistence_cutoff", "sb_dat", \
                     "sb_graph", "sb_cg_file", "sb_mode", \
                     "hc_co", "hc_cutoff", "hc_residues", \
                     "hc_perco", "hc_persistence_cutoff", \
                     "hc_dat", "hc_graph", "hb_co", \
                     "hb_cutoff", "hb_ang", "hb_angle", \
                     "hb_dat", "hb_graph", "hb_perco", \
                     "hb_class", "hb_file", "hb_custom_group_1", \
                     "hb_custom_group_2", "kbp_ff", "force_field", \
                     "kbp_atomlist", "kbp_dat", "kbp_graph", \
                     "kbp_kbt", "ff_masses"], \
            }, \

         "filter_graph" : \
            {
                "FLAGS" : \
                    ["h", "f"],

                "DDFLAGS" : \
                    ["help"],

                "SDARGS" : \
                    ["d", "o", "c", "t", "p", "u", \
                     "l", "s", "w", "x", "k", "m", "n"], \

                "DDARGS" : \
                    ["input_dat", "output_dat", \
                     "output_clusters", "filter_threshold", \
                     "plot", "range_upper", "range_lower", \
                     "range_step", "weight_matrix", "x0"] \
            }, \

         "graph_analysis" : \
            {
                "FLAGS" : \
                    ["h", "c", "u", "p", "d"], \

                "DDFLAGS" : \
                    ["help", "components", "hubs", \
                     "all_paths", "write_paths"],

                "SDARGS" : \
                    ["r", "a", "k", "r1", "r2", \
                     "l", "s", "cb", "ub"], \

                "DDARGS" : \
                    ["reference", "adj_matrix", "hubs_cutoff", \
                     "source", "target", \
                     "maximum_path_length", "sort_paths", \
                     "components_pdb", "hubs_pdb"], \
            }, \
        }

    
    def __init__(self, \
                 command, \
                 **kwargs):
        """Initialize the caller.

        Parameters
        ----------
        command : `str`, {'pyinteraph', 'filter_graph', 
                  'graph_analysis'}
            Name of the command.

        **kwargs 
            ``Command`` options spelled as normally done
            for the corresponding command-line options,
            apart from the starting '-' or '--', that must be
            omitted. Moreover, all dashes in the options names 
            must be converted to underscores (names of function
            arguments cannot include dashes).

        Attributes
        ----------
        call_args : `list`
            List contaning all the arguments and options passed
            to the ``command``, including the ``command`` (first in
            the list). It is formatted as required by `subprocess`
            calls.

        returncode : `int` or `None`, default: `None`
            Return code of the call. Values different from 0
            indicates that something went wrong. Initialized
            to `None` before the command is run.
        """
        self.call_args = [command]
        self.returncode = None

        for arg, value in kwargs.items():
            if arg in self.__class__.ARGS[command]["FLAGS"]:
                # set the flag only if True is passed
                # as a value
                if value is True:
                    self.call_args.append("-" + arg)

            elif arg in self.__class__.ARGS[command]["DDFLAGS"]:
                # set the flag only if True is passed
                # as a value
                if value is True:
                    self.call_args.append("--" + arg)

            elif arg in self.__class__.ARGS[command]["SDARGS"]:
                # add a single dash before the argument
                if value is not None:
                    self.call_args.append("-" + arg)
                    self.call_args.append(str(value))       
           
            elif arg in self.__class__.ARGS[command]["DDARGS"]:
                # add a double dash before the argument and
                # replace all underscores in the option name
                # with dashes
                if value is not None:
                    self.call_args.append(\
                        "--" + arg.replace("_", "-"))
                    self.call_args.append(str(value))


    def run(self, \
            stdout = None):
        """Run the command.

        Parameters
        ----------
        stdout : file handler or `None`, default: `None`
            Where to redirect the standard
            output (and standard error, if `stderr`
            is `None`). 
            If `None` (default), the standard output
            will be directed as specified by the
            process called.

        Returns
        -------
        `None`
        """   
        self.returncode = \
            subprocess.call(self.call_args, \
                            stdout = stdout, \
                            stderr = subprocess.STDOUT)

        if stdout is not None:
            # close the opened file
            stdout.close()


    def check(self, \
              logfile = None):
        """Check whether the command call was successful.

        Parameters
        ----------
        logfile : file handler or `None`, default: `None`
            Where to write a warning in case something goes
            wrong during the call (return code different
            from 0). 

        Returns
        -------
        `None`
        """     
        warnstr = \
            "\nWARNING:Command {:s} executed in {:s} exited " \
            "with code {:d}. Please check the corresponding " \
            "output file for possible errors.\n"

        if self.returncode == 0:
            return
        
        else:
            logstr = warnstr.format(" ".join(self.call_args), \
                                    os.getcwd(), \
                                    self.returncode)
            if logfile == None:
                sys.stdout.write(logstr)
            else:
                logfile.write(logstr)
                logfile.close()


######################## STANDALONE FUNCTIONS #########################

def run_network_analysis(pcut, \
                         refpdb, \
                         k, \
                         graph, \
                         filtered_graph, \
                         pdb_hubs, \
                         pdb_cc, \
                         outfile_hubs, \
                         outfile_cc, \
                         logfile):
    """Wrapper function to run the network analysis 
    pipeline in one go.

    Parameters
    ----------
    pcut : `int`
        Persistence cut-off.

    refpdb : `str`
        Name of the PDB reference structure.

    k : `int`
        Minimum number of edges for a node to be 
        considered a hub.

    graph : `str`
        Name of the raw graph matrix (not filtered).

    filtered_graph : `str`
        Name of the filtered graph matrix. It is
        generated during graph filtering.

    pdb_hubs : `str`
        Name of the PDB file in which the reference
        structure is stored, but with the b-factor
        column replaced by the hub degree for hub
        nodes, 0 for others. It is generated during
        graph analysis.

    pdb_cc : `str`
        Name of the PDB file in which the reference
        structure is stored, but with the b-factor
        column replaced by the number of the connected
        component each node belongs to. It is generated
        during graph analysis.

    outfile_hubs : `str`
        Name of the text file in which the list of 
        hubs is stored.

    outfile_cc : `str`
        Name of the text file in which the list of 
        connected components is stored.

    Returns
    -------
    `None`
    """
    # filter the graph
    fg = PyInteraphCaller("filter_graph", \
                          d = graph, \
                          o = filtered_graph, \
                          t = pcut)
    
    fg.run()
    fg.check(logfile = open(logfile, "a"))

    # set the caller for the calculation
    # of the connected components
    ga_cc = PyInteraphCaller("graph_analysis", \
                             r = refpdb, \
                             a = filtered_graph, \
                             c = True, \
                             cb = pdb_cc)

    # capture the output of the connected components
    # calculation and store it into a file (it is the
    # list of the connected components)
    ga_cc.run(stdout = open(outfile_cc, "w"))
    ga_cc.check(logfile = open(logfile, "a"))

    # set the caller for hubs identification
    ga_hubs = PyInteraphCaller("graph_analysis", \
                               r = refpdb, \
                               a = filtered_graph, \
                               u = True, \
                               ub = pdb_hubs, \
                               k = k)

    # capture the output of the hubs identification and
    # store it into a file (it is the list of hubs)
    ga_hubs.run(stdout = open(outfile_hubs, "w"))
    ga_hubs.check(logfile = open(logfile, "a"))


def run_all(ntraj, \
            pydir, \
            perform_hc, \
            perform_hb, \
            perform_sb, \
            nsamplings, \
            lastframendx, \
            step, \
            traj, \
            refpdb, \
            refgro, \
            ff, \
            pcut, \
            k, \
            hc_res, \
            hc_dcuts, \
            hb_classes, \
            hb_dcuts, \
            hb_file, \
            sb_modes, \
            sb_dcuts, \
            sb_file, \
            logfile):
    """Wrap all the pipeline in one function (necessary for
    multiprocessing).

    Parameters
    ----------
    ntraj : `int`
        Trajectory number (0 correspond to the full 
        trajectory, all the others are resamplings).

    pydir : `str`
        Name of the analyses directory.

    perform_hc : `bool`
        Whether to perform the hydrophobic contacts analysis.

    perform_hb : `bool`
        Whether to perform the hydrogen bonds analysis.

    perform_sb : `bool`
        Whether to perform the salt bridges analysis.

    nsamplings : `int`
        Number of samplings to be performed.

    lastframendx : `int`
        Index of the last used frame of the trajectory
        (may not be the actual end of the trajectory,
        since the last frames may be eliminated to
        get even partitioning during resampling).
    
    step : `int`
        Width (in number of frames) of the sliding 
        window defining the portion of the trajectory
        excluded at each resampling.
    
    traj : `str`
        Name of the XTC trajectory.
    
    refpdb : `str`
        Name of the reference PDB structure.
    
    refgro : `str`
        Name of the reference GRO structure.
    
    ff : `str`
        Name of the force field from which atomic
        masses will be taken.
    
    pcut : `int`
        Persistence cut-off.
    
    k : `int`
        Minimum number of edges for a node to be
        considered a hub.
    
    hc_res : `list`
        List of residues types to be used in evaluating
        the hydrophobic interactions.
    
    hc_dcuts : `list` or `numpy.ndarray`
        List/Array of distance cut-offs for the
        analysis of hydrophobic interactions.
    
    hb_classes : `list`, {'mc_mc', 'sc_sc', 
                 'mc_sc', 'all'}
        Classes of hydrogen bonds to be analyzed.
    
    hb_dcuts : `list` or `numpy.ndarray`
        List/Array of distance cut-offs for the
        analysis of hydrogen bonds.
    
    hb_file : `str`
        Name of a custom acceptors-donors file
        to be used for hydrogen bonds analysis.
    
    sb_modes : `list`, {'all', 'different_charge',
                'same_charge'}
        Types of salt bridges to be analyzed.
    
    sb_dcuts : `list` or `numpy.ndarray`
        List/Array of distance cut-offs for the
        analysis of salt bridges.
    
    sb_file : `str`
        Name of a custom file for charged grouos
        to be used for the analysis of salt bridges.
    
    logfile : `str`
        Name of the file where to store logs.
    
    Returns
    -------
    None
    """      
    if ntraj == 0:
        # the first trajectory in the list is always the full one
        fulltrajdir = os.path.join(pydir, "fulltraj")
        os.mkdir(fulltrajdir)
        os.chdir(fulltrajdir)

        # create a new MDAnalysis Universe object from the topology
        # and the trajectory (since we are moving from one directory
        # to the other we cannot use the one created in the __main__,
        # because the file references would not be the same)
        u = mda.Universe(refpdb, traj)

        # set the currently processed trajectory to the full
        # one and the current directory to that storing the
        # analyses for the full trajectory
        curr_traj = traj
        curr_dir = fulltrajdir
        
    else:
        # if you have more than one trajectory, those from index
        # one to the end are those resampled
        curr_resampl = ntraj-1

        # make a dedicated directory for the resampling round
        # and move to it
        resampldir = \
            os.path.join(pydir, \
                "resampling{:d}".format(curr_resampl))
                    
        os.mkdir(resampldir)
        os.chdir(resampldir)
            
        # create a new MDAnalysis Universe object from the topology
        # and the trajectory (since we are moving from one directory
        # to the other we cannot use the one created in the __main__,
        # because the file references would not be the same)
        u = mda.Universe(refpdb, traj)

        # define a name for the resampled sub-trajectory
        subtraj = "subtraj{:d}.xtc".format(curr_resampl)        

        # exclude a (nframes//nsamplings)% (i.e. 10% if 100 frames
        # and nsamplings = 10) of the frames every time, starting
        # from the last (nframes//nsamplings)% of the frames
        end_first_part = lastframendx-(step*(curr_resampl+1))
        begin_second_part = lastframendx-(step*(curr_resampl))
            
        with mda.Writer(subtraj, u.atoms.n_atoms) as Wxtc:
            # there is no "first part" if the excluded frames
            # are those at the beginning of the trajectory
            if curr_resampl != (nsamplings-1):
                for ts in u.trajectory[:end_first_part]:
                    Wxtc.write(u.atoms)

            # there is no "second part" if the excluded frames
            # are those at the end of the trajectory
            if curr_resampl != 0:
                for ts in u.trajectory[begin_second_part:lastframendx]:
                    Wxtc.write(u.atoms)

        # use absolute paths to avoid relative paths issues
        # set the currently processed trajectory to the current
        # resampled trajectory and the current directory to
        # that storing the analyses for the resampled trajectory
        curr_traj = os.path.abspath(subtraj)
        curr_dir = resampldir
        
    # set names for the graph_analysis outputs
    pdb_hubs = "hubs.pdb"
    pdb_cc = "cc.pdb"
    outfile_hubs = "hubs.out"
    outfile_cc = "cc.out"


    ################## HYDROPHOBIC CONTACTS ANALYSIS ##################

    if perform_hc:
              
        # Create the directory for the hydrophobic contacts analysis
        # and move to it
        contactdir = os.path.join(curr_dir, "contact")
        os.mkdir(contactdir)
        os.chdir(contactdir)

        hc_graph = "contact_graph.dat"
        hc_filter_graph = "contact_filter.dat"

        for hc_dcut in hc_dcuts:
            hcdcutdir = os.path.join(contactdir, str(hc_dcut))
            os.mkdir(hcdcutdir)
            os.chdir(hcdcutdir)

            pyin = PyInteraphCaller("pyinteraph", \
                                    t = traj, \
                                    s = refgro, \
                                    r = refpdb, \
                                    f = True, \
                                    hc_graph = hc_graph, \
                                    ff_masses = ff, \
                                    hc_co = hc_dcut, \
                                    hc_res = hc_res)

            # capture the output of the newtwork generation
            pyin.run(stdout = open("contact.out", "w"))
            pyin.check(logfile = open(logfile, "a"))

            # run the network analysis pipeline
            run_network_analysis(pcut = pcut, \
                                 refpdb = refpdb, \
                                 k = k, \
                                 graph = hc_graph, \
                                 filtered_graph = hc_filter_graph, \
                                 pdb_hubs = pdb_hubs, \
                                 pdb_cc = pdb_cc, \
                                 outfile_hubs = outfile_hubs, \
                                 outfile_cc = outfile_cc, \
                                 logfile = logfile)

            # go back to the `contact` directory
            os.chdir(contactdir)      


    ######################### HYDROGEN BONDS ##########################

    # check if hydrogen bonds analyses were requested
    if perform_hb:               
        # Create the directory for the hydrogen bonds analysis
        # and move to it
        hbondsdir = os.path.join(curr_dir, "h_bonds")
        os.mkdir(hbondsdir)
        os.chdir(hbondsdir)

        hb_graph = "hb_graph.dat"
        hb_filter_graph = "hb_filter.dat"

        for hb_class in hb_classes:
            hbclassdir = os.path.join(hbondsdir, hb_class)
            os.mkdir(hbclassdir)
            os.chdir(hbclassdir)
            for hb_dcut in hb_dcuts:
                hbdcutdir = os.path.join(hbclassdir, str(hb_dcut))
                os.mkdir(hbdcutdir)
                os.chdir(hbdcutdir)

                pyin = PyInteraphCaller("pyinteraph", \
                                        t = traj, \
                                        s = refgro, \
                                        r = refpdb, \
                                        y = True, \
                                        hb_graph = hb_graph, \
                                        ff_masses = ff, \
                                        hb_co = hb_dcut, \
                                        hb_ad_file = hb_file, \
                                        hb_class = hb_class)

                # capture the output of the newtwork generation
                pyin.run(stdout = open("hb.out", "w"))
                pyin.check(logfile = open(logfile, "a"))

                # run the network analysis pipeline
                run_network_analysis(pcut = pcut, \
                                     refpdb = refpdb, \
                                     k = k, \
                                     graph = hb_graph, \
                                     filtered_graph = hb_filter_graph, \
                                     pdb_hubs = pdb_hubs, \
                                     pdb_cc = pdb_cc, \
                                     outfile_hubs = outfile_hubs, \
                                     outfile_cc = outfile_cc, \
                                     logfile = logfile)

                # go back to the directory with the analyses for
                # the current class of hydrogen bonds
                os.chdir(hbclassdir)

            # go back to the `h_bonds` directory
            os.chdir(hbondsdir)


    ########################## SALT BRIDGES ###########################

    # check if salt bridges analyses were requested
    if perform_sb:               
        # Create the directory for the hydrogen bonds analysis
        # and move to it
        sbridgesdir = os.path.join(curr_dir, "s_bridges")
        os.mkdir(sbridgesdir)
        os.chdir(sbridgesdir)    

        sb_graph = "sb_graph.dat"
        sb_filter_graph = "sb_filter_graph.dat"

        for sb_mode in sb_modes:
            sbmodedir = os.path.join(sbridgesdir, sb_mode)
            os.mkdir(sbmodedir)
            os.chdir(sbmodedir)
            for sb_dcut in sb_dcuts:
                sbdcutdir = os.path.join(sbmodedir, str(sb_dcut))
                os.mkdir(sbdcutdir)
                os.chdir(sbdcutdir)

                pyin = PyInteraphCaller("pyinteraph", \
                                        t = traj, \
                                        s = refgro, \
                                        r = refpdb, \
                                        b = True, \
                                        sb_mode = sb_mode, \
                                        sb_graph = sb_graph, \
                                        ff_masses = ff, \
                                        sb_co = sb_dcut, \
                                        sb_cg_file = sb_file)

                # capture the output of the newtwork generation
                pyin.run(stdout = open("sb.out", "w"))
                pyin.check(logfile = open(logfile, "a"))

                # run the network analysis pipeline
                run_network_analysis(pcut = pcut, \
                                     refpdb = refpdb, \
                                     k = k, \
                                     graph = sb_graph, \
                                     filtered_graph = sb_filter_graph, \
                                     pdb_hubs = pdb_hubs, \
                                     pdb_cc = pdb_cc, \
                                     outfile_hubs = outfile_hubs, \
                                     outfile_cc = outfile_cc, \
                                     logfile = logfile)

                # go back to the directory storing the analyses for
                # the current salt bridges mode
                os.chdir(sbmodedir)

            # go back to the `s_bridges` directory
            os.chdir(sbridgesdir)


def parse_config(configfile, \
                 sec2fields):
    """Parse the configuration file.

    Parameters
    ----------
    configfile : `str`
        Name of the config file.

    sec2fields : `dict`
        Dictionary mapping the name of each
        section to the name of its fields, which
        are in turn mapped to a tuple containing
        both the expected data type and the default
        value for that field.

    Returns
    -------
    configvalues : `dict`
        Dict mapping the name of each section to
        the name of its fields, which are in turn
        mapped to their corresponding value.
    """
    with open(configfile, "r") as f:
        configvalues = {}
        curr_sec = None
        for line in f:
            if re.match(r"^\s*$", line) \
            or line.startswith("#"):
                # ignore empty and comment lines
                continue

            elif line.startswith("["):
                # a new section has started
                sec = line.strip("[]\n")
                configvalues[sec] = {}
                curr_sec = sec

            else:
                opt, val = line.split("=")
                opt = opt.strip(" ")
                dtype, default = sec2fields[curr_sec][opt]
                val = val.strip(" \n")
                if val == "":
                    # if an empty string was parsed, assign 
                    # the default value to that field
                    val = default
                if dtype == list:
                    # if the data type is a list, convert the
                    # string into a list
                    val = \
                        [v.strip(" ") for v in val.split(",")]

                # cannot pass them directly to bool() or both
                # "True" and "False" will be evaluated to True!
                if val == "True":
                    val = True
                elif val == "False":
                    val = False
                
                configvalues[curr_sec][opt] = dtype(val)

    return configvalues


def parse_dcuts_definition(dcutsdef):
    """Parse the definition of the distance cut-offs, that can
    be either a range defined as `range(`) or `numpy.arange()`
    or a comma-separated list of numbers.

    Parameters
    ----------
    dcutsdef : str
        Distance cut-offs definition as parsed from the 
        configuration file.

    Returns
    -------
    dcuts
        Either a `numpy.ndarray` (if the range was defined
        as a Python `range()` or a `numpy.arange()`) or a
        list of floats (if it was defined as a comma-separated
        list of numbers).

    See also
    --------
    See :func:`~numpy.arange` for more information about the
    `numpy.arange()` function.
    """
    if "range" in dcutsdef:
        # parse a string defining a range as a np.arange() and
        # convert it into the corresponding range
        _dcuts = \
            [float(num) for num \
            in dcutsdef.strip("(np.range)\n").split(",")]
        dcuts = np.arange(*_dcuts)
    else:
        # parse a string defining a range as a comma-separated
        # list of values
        dcuts = [float(num.strip(" ")) for num \
                in dcutsdef.strip("\n").split(",")]

    return dcuts    



#######################################################################
###########################   MAIN   ##################################
#######################################################################

if __name__ == "__main__":

    ######################## ARGUMENT PARSER ##########################
    
    description = \
        '"""\nPipeline for network analysis with PyInteraph.\n"""'

    # create the argument parser
    parser = argparse.ArgumentParser(description = description)

    # add arguments to the parser
    parser.add_argument("-f", \
                        type = str, \
                        metavar = "CONFIG_FILE", \
                        required = True, \
                        help = "Configuration file")

    args = parser.parse_args()


    ################### CONFIGURATION FILE PARSER #####################
    
    sec2fields =   \
        {"ENVIRONMENT" : \
            {"PYINENV" : (str, None), \
             "NPROC" : (int, 1)}, \
         
         "ANALYSES" : \
            {"HC" : (bool, False), \
             "HB" : (bool, False), \
             "SB" : (bool, False)}, \
         
         "INPUT_FILES" : \
            {"TRAJ" : (str, None), \
             "PDB" : (str, None)}, \
            
         "GENERAL_OPTIONS" : \
            {"FORCE_FIELD" : (str, None), \
             "PCUT" : (int, 0), \
             "K" : (int, 1), \
             "NSAMPLINGS" : (int, 0)}, \
           
         "HC_OPTIONS" : \
            {"HC_DCUTS" : (str, None), \
             "HC_RES" : (list, [])}, \

          "HB_OPTIONS" : \
           {"HB_CLASSES" : (list, []), \
            "HB_DCUTS" : (str, None), \
            "HB_FILE" : (str, None)}, \

          "SB_OPTIONS" : \
           {"SB_MODES" : (list, []), \
            "SB_DCUTS" : (str, None), \
            "SB_FILE" : (str, None)}
        }

    options = parse_config(args.f, sec2fields)

    ########################## LOGGER SETUP ###########################
    
    loglevel = logging.INFO
    logfile = os.path.join(os.getcwd(), "pyinknife.log")
    logging.basicConfig(filename = logfile, \
                        level = loglevel, \
                        filemode = "w")


    ###################### RETRIEVE THE OPTIONS #######################

    pyin_actenv = options["ENVIRONMENT"]["PYINENV"]
    nprocesses = options["ENVIRONMENT"]["NPROC"]

    perform_hc = options["ANALYSES"]["HC"]
    perform_hb = options["ANALYSES"]["HB"]
    perform_sb = options["ANALYSES"]["SB"]

    refpdb = os.path.abspath(options["INPUT_FILES"]["PDB"])
    traj = os.path.abspath(options["INPUT_FILES"]["TRAJ"])

    pcut = options["GENERAL_OPTIONS"]["PCUT"]
    ff = options["GENERAL_OPTIONS"]["FORCE_FIELD"]
    k = options["GENERAL_OPTIONS"]["K"]
    nsamplings = options["GENERAL_OPTIONS"]["NSAMPLINGS"]

    hc_res = options["HC_OPTIONS"]["HC_RES"]
    hc_dcuts = \
        parse_dcuts_definition(options["HC_OPTIONS"]["HC_DCUTS"])

    hb_classes = options["HB_OPTIONS"]["HB_CLASSES"]
    hb_dcuts = \
        parse_dcuts_definition(options["HB_OPTIONS"]["HB_DCUTS"])
    hb_file = options["HB_OPTIONS"]["HB_FILE"]
    if hb_file is not None:
        hb_file = os.path.abspath(hb_file)

    sb_modes = options["SB_OPTIONS"]["SB_MODES"]
    sb_dcuts = \
        parse_dcuts_definition(options["SB_OPTIONS"]["SB_DCUTS"])
    sb_file = options["SB_OPTIONS"]["SB_FILE"]
    if sb_file is not None:
        sb_file = os.path.abspath(sb_file)

    # set the path to the activation file for the 
    # PyInteraph environment
    logstr = "Loaded pyinteraph environment from {:s}"
    logging.info(logstr.format(pyin_actenv))

    # Activate the PyInteraph environment
    execfile(pyin_actenv, dict(__file__= pyin_actenv))

    # Log about python and MDAnalysis versions loaded (to be
    # sure they are not those of the virtual environment)
    logstr = "Using {:s} version {:s}"
    logging.info(logstr.format("python", sys.version))
    logging.info(logstr.format("MDAnalysis", mda.__version__))

    # Create a new directory for the analyses into
    # the current working directory and move to it
    namedir = "pyinteraph_analyses"
    pydir = os.path.join(os.getcwd(), namedir)
    os.mkdir(pydir)
    os.chdir(pydir)

    # create a MDAnalysis Universe object from the topology
    # and the trajectory
    u = mda.Universe(refpdb, traj)
    logstr = \
        "MDAnalysis.Universe object created from:\n" \
        "topology: {:s}\ntrajectory {:s}"
    logging.info(logstr.format(refpdb, traj))

    # write a .gro reference structure
    refgro = os.path.join(pydir, "model0.gro")
    with mda.Writer(refgro, u.atoms.n_atoms) as Wgro:
        Wgro.write(u.atoms)
    logstr = \
        "{:s} file created from MDAnalysis.Universe object"
    logging.info(logstr.format(refgro))

    
    ntrajs = 1
    if nsamplings is not None:
        ntrajs += nsamplings
        # get the number of frames in the trajectory
        nframes = len(u.trajectory)
        logging.info("Trajectory has {:d} frames".format(nframes))
        # log the number of resamplings
        logging.info(\
            "{:d} resamplings will be performed".format(nsamplings))
        # esclude last frames if the division between the total number
        # of frames and the number of resamplings has a remainder
        lastframendx = nframes - (nframes % nsamplings)
        step = lastframendx / nsamplings

        logstr = \
            "The last {:d} frames have been excluded from the " \
            "trajectory to have resampled trajectories of the " \
            "same length"
        logging.warning(logstr.format(nframes-lastframendx))
    
    # run all the analysis pipeline in parallel on multiple
    # processes
    pool = multiprocessing.Pool(nprocesses)
    partfunc = \
        functools.partial(run_all, \
                          pydir = pydir, \
                          perform_hc = perform_hc, \
                          perform_hb = perform_hb, \
                          perform_sb = perform_sb, \
                          nsamplings = nsamplings, \
                          lastframendx = lastframendx, \
                          step = step, \
                          traj = traj, \
                          refpdb = refpdb, \
                          refgro = refgro, \
                          ff = ff, \
                          pcut = pcut, \
                          k = k, \
                          hc_res = hc_res, \
                          hc_dcuts = hc_dcuts, \
                          hb_classes = hb_classes, \
                          hb_dcuts = hb_dcuts, \
                          hb_file = hb_file, \
                          sb_modes = sb_modes, \
                          sb_dcuts = sb_dcuts, \
                          sb_file = sb_file, \
                          logfile = logfile)


    pool.map(partfunc, range(ntrajs))
    pool.close()
    pool.join()