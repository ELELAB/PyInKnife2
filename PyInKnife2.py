#!/usr/bin/env python
# -*- coding: utf-8 -*-

#    Pipeline for graph analysis with PyInteraph, a software suite 
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
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

import argparse
import sys
import os
import os.path
import subprocess
import multiprocessing
import ConfigParser
import logging
import functools

import numpy as np
import MDAnalysis as mda

__author__ = ["Valentina Sora", "Juan Salamanca Viloria", "Elena Papaleo"]
__credits__ = ["Valentina Sora", "Juan Salamanca Viloria", "Elena Papaleo"]
__license__ = "GPL"
__version__ = "0.0.1"
__date__ = "30/4/2019"
__maintainer__ = "Valentina Sora"
__email__ = "vaso@cancer.dk"
__status__ = "Development"

class PyInteraphCaller:
    
    """Class implementing a wrapper for calling PyInteraph
    tools.
    For each option provided for the command, a homonymous
    attribute is set for the corresponding instance.
    So far, the following PyInteraph commands are supported:
    -   pyinteraph
    -   filter_graph
    -   graph_analysis

    Parameters
    ----------
    command : string
        The command usually run on the command line.

    **kwargs : 
        Command options spelled as normally done
        for the corresponding command-line options,
        apart from the starting '-' or '--'. Moreover,
        all dashes in the options names must be 
        converted to underscores (names of function
        arguments cannot include dashes). 

    Attributes
    ----------
    ARGS : dict
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
                     "hb_class", "hb_ad_file", "hb_custom_group_1", \
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
        command : str
            Name of the command.

        **kwargs : dict
            Command options spelled as normally done
            for the corresponding command-line options,
            apart from the starting '-' or '--', and
            all internal dashes in the option names 
            must be substituted by underscores.
        """

        self._call_args = [command]

        for arg, value in kwargs.items():
            setattr(self, "%s" % arg, "%s" % value)

            if arg in self.__class__.ARGS[command]["FLAGS"]:
                # set the flag only if the True is passed
                # as a value
                if value is True:
                    self._call_args.append("-" + arg)

            elif arg in self.__class__.ARGS[command]["DDFLAGS"]:
                # set the flag only if the True is passed
                # as a value
                if value is True:
                    self._call_args.append("--" + arg)

            elif arg in self.__class__.ARGS[command]["SDARGS"]:
                # add a single dash before the argument
                self._call_args.append("-" + arg)
                self._call_args.append(str(value))       
           
            elif arg in self.__class__.ARGS[command]["DDARGS"]:
                # add a double dash before the argument and
                # replace all underscores in the option name
                # with dashes
                self._call_args.append(\
                    "--" + arg.replace("_", "-"))
                self._call_args.append(str(value))


    def run(self, \
            stdout = None, \
            stderr = None):

        """Run the command.

        Parameters
        ----------
        stdout : file handler or None, default: None
            Where to redirect the standard
            output (and standard error, if stderr
            is None). 
            If None (default), the standard output
            will be directed as specified by the
            process called.

        stderr : file handler or None, default: None
            Where to redirect the standard error.
            If None (default), the standard error
            will be redirected either to the standard
            output (if stdout is not None) or as
            specified by the process called (if stdout
            is None).

        Returns
        -------
        subprocess.CompletedProcess
            An instance of the subprocess.CompletedProcess
            class.
        """
        
        if stdout is not None:
            return subprocess.call(self._call_args, \
                                   stdout = stdout, \
                                   stderr = subprocess.STDOUT)
        else:
            return subprocess.call(self._call_args)



######################## STANDALONE FUNCTIONS #########################

def run_network_analysis(pcut, \
                         refpdb, \
                         k, \
                         hc_graph, \
                         filtered_graph, \
                         pdb_hubs, \
                         pdb_cc, \
                         outfile_hubs, \
                         outfile_cc):

    """Wrapper function to run the network analysis 
    pipeline in one go (necessary for multiprocessing).

    Parameters
    ----------
    pcut : int
        Persistence cut-off.

    refpdb : str
        Name of the PDB reference structure.

    k : int
        Minimum number of edges for a node to be 
        considered a hub.

    hc_graph : str
        Name of the raw graph matrix (not filtered).

    filtered_graph : str
        Name of the filtered graph matrix. It is
        generated during graph filtering.

    pdb_hubs : str
        Name of the PDB file in which the reference
        structure is stored, but with the b-factor
        column replaced by the hub degree for hub
        nodes, 0 for others. It is generated during
        graph analysis.

    pdb_cc : str
        Name of the PDB file in which the reference
        structure is stored, but with the b-factor
        column replaced by the number of the connected
        component each node belongs to. It is generated
        during graph analysis.

    outfile_hubs : str
        Name of the text file in which the list of 
        hubs is stored.

    outfile_cc : str
        Name of the text file in which the list of 
        connected components is stored.

    Returns
    -------
    None
    """

    # filter the graph
    fg = PyInteraphCaller("filter_graph", \
                          d = hc_graph, \
                          o = filtered_graph, \
                          t = pcut)
    fg.run()

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
    ga_cc.run(stdout = open(outfile_cc, "w"), \
              stderr = open("cc.err", "w"))

    # set the caller for hubs identification
    ga_hubs = PyInteraphCaller("graph_analysis", \
                               r = refpdb, \
                               a = filtered_graph, \
                               u = True, \
                               ub = pdb_hubs, \
                               k = k)

    # capture the output of the hubs identification and
    # store it into a file (it is the list of hubs)
    ga_hubs.run(stdout = open(outfile_hubs, "w"), \
                stderr = open("hubs.err", "w"))


def run_cm_pipeline(dcut, \
                    wdir, \
                    traj, \
                    refgro, \
                    refpdb, \
                    ff, \
                    hc_res, \
                    pcut, \
                    k, \
                    hc_graph, \
                    filtered_graph, \
                    pdb_hubs, \
                    pdb_cc, \
                    outfile_contact, \
                    outfile_hubs, \
                    outfile_cc):
 
    """Wrapper function to run the full pipeline for
    the contact map (building and analysis) in one go 
    (necessary for multiprocessing).

    Parameters
    ----------
    dcut : int or float
        Distance cut-off.

    wdir : str
        Path to the working directory in which the
        pipeline will be run.

    traj : str
        Name of the XTC trajectory.

    refgro : str
        Name of the GRO reference structure.

    refpdb : str
        Name of the PDB reference structure.

    ff : str
        Force field name from which atomic masses 
        will be taken.

    hc_res : list
        List of three-letters residue names to include
        in the contact map calculation.

    pcut : int
        Persistence cut-off.

    k : int
        Minimum number of edges for a node to be 
        considered a hub.

    hc_graph : str
        Name of the raw graph matrix (not filtered).

    filtered_graph : str
        Name of the filtered grap matrix. It is
        generated during graph filtering.

    pdb_hubs : str
        Name of the PDB file in which the reference
        structure is stored, but with the b-factor
        column replaced by the hub degree for hub
        nodes, 0 for others. It is generated during
        graph analysis.

    pdb_cc : str
        Name of the PDB file in which the reference
        structure is stored, but with the b-factor
        column replaced by the number of the connected
        component each node belongs to. It is generated
        during graph analysis.

    outfile_contact : str
        Name of the text file in which the output of
        the graph generation will be stored.

    outfile_hubs : str
        Name of the text file in which the list of 
        hubs is stored.

    outfile_cc : str
        Name of the text file in which the list of 
        connected components is stored.

    Returns
    -------
    None
    """

    logstr = \
        "Building and analyzing contact map for dcut {:f}"
    logging.info(logstr.format(dcut)) 

    # make sure you are in the correct directory
    os.chdir(wdir)

    # Set the path to the directory that will 
    # contain the analyses performed with this
    # distance cutoff
    dcutdir = os.path.join(wdir, str(dcut))
            
    # Create the directory and move to it
    os.mkdir(dcutdir)
    os.chdir(dcutdir)

    # Contact map generation
    pyin_cm = PyInteraphCaller("pyinteraph", \
                                v = True, \
                                t = traj, \
                                s = refgro, \
                                r = refpdb, \
                                f = True, \
                                hc_co = dcut, \
                                ff_masses = ff, \
                                hc_residues = ",".join(hc_res), \
                                hc_graph = hc_graph)

    # capture the output of the contact map generation
    out_cm = pyin_cm.run(stdout = open(outfile_contact, "w"), \
                         stderr = open("cm.err", "w"))
    
    # run the network analysis pipeline   
    run_network_analysis(pcut = pcut, \
                         refpdb = refpdb, \
                         k = k, \
                         hc_graph = hc_graph, \
                         filtered_graph = filtered_graph, \
                         pdb_hubs = pdb_hubs, \
                         pdb_cc = pdb_cc, \
                         outfile_hubs = outfile_hubs, \
                         outfile_cc = outfile_cc)


def run_hb_pipeline(hb_analysis, \
                    wdir, \
                    traj, \
                    refgro, \
                    refpdb, \
                    ff, \
                    pcut, \
                    k, \
                    hb_graph, \
                    filtered_graph, \
                    pdb_hubs, \
                    pdb_cc, \
                    outfile_hb, \
                    outfile_hubs, \
                    outfile_cc, \
                    hb_ad_file):

    """Wrapper function to run the full pipeline for the 
    hydrogen bonds network in one go (necessary for 
    multiprocessing).

    Parameters
    ----------
    hb_analysis : str, allowed: 'sc-sc', 'mc-mc', 'mc-sc', 'all'
        Type of hydrogen bonds analysis to perform. Allowed
        values are those accepted by the pyinteraph command.

    wdir : str
        Path of the working directory in which the
        pipeline will be run.

    traj : str
        Name of the .xtc trajectory.

    refgro : str
        Name of the .gro reference structure.

    refpdb : str
        Name of the PDB reference structure.

    ff : str
        Force field name from which atomic masses will be taken.

    pcut : int
        Persistence cut-off.

    k : int
        Minimum number of edges for a node to be 
        considered a hub.

    hb_graph : str
        Name of the raw graph matrix (not filtered).

    filtered_graph : str
        Name of the filtered graph matrix. It is
        generated during graph filtering.

    pdb_hubs : str
        Name of the PDB file in which the reference
        structure is stored, but with the b-factor
        column replaced by the hub degree for hub
        nodes, 0 for others. It is generated during
        graph analysis.

    pdb_cc : str
        Name of the PDB file in which the reference
        structure is stored, but with the b-factor
        column replaced by the number of the connected
        component each node belongs to. It is generated
        during graph analysis.

    outfile_hb : str
        Name of the text file in which the output of
        the graph generation will be stored.

    outfile_hubs : str
        Name of the text file in which the list of 
        hubs is stored.

    outfile_cc : str
        Name of the text file in which the list of 
        connected components is stored.

    Returns
    -------
    None
    """

    logstr = \
        "Building and analyzing {:s} hydrogen bonds network"
    logging.info(logstr.format(hb_analysis))
      
    # create a subdirectory for the specific analysis 
    # and move to it
    hbandir = os.path.join(wdir, hb_analysis.replace("-", "_"))
    os.mkdir(hbandir)
    os.chdir(hbandir)

    if hb_ad_file is not None:
        pyin_hb = PyInteraphCaller("pyinteraph", \
                                   t = traj, \
                                   s = refgro, \
                                   r = refpdb, \
                                   y = True, \
                                   hb_graph = hb_graph, \
                                   ff_masses = ff, \
                                   hb_class = hb_analysis, \
                                   hb_ad_file = hb_ad_file)
    else:
        pyin_hb = PyInteraphCaller("pyinteraph", \
                                   t = traj, \
                                   s = refgro, \
                                   r = refpdb, \
                                   y = True, \
                                   hb_graph = hb_graph, \
                                   ff_masses = ff, \
                                   hb_class = hb_analysis)
    
    # capture the output of the newtwork generation
    out_hb = pyin_hb.run(stdout = open(outfile_hb, "w"), \
                         stderr = open("hb.err", "w"))

    # run the network analysis pipeline
    run_network_analysis(pcut = pcut, \
                         refpdb = refpdb, \
                         k = k, \
                         hc_graph = hb_graph, \
                         filtered_graph = filtered_graph, \
                         pdb_hubs = pdb_hubs, \
                         pdb_cc = pdb_cc, \
                         outfile_hubs = outfile_hubs, \
                         outfile_cc = outfile_cc)


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
    
    config = ConfigParser.SafeConfigParser()
    config.read(args.f)
    
    # parse [Environment] section
    pyin_actenv = config.get("Environment", "PyInteraphEnv")

    # parse [Tasks] section
    perform_cm = config.getboolean("Tasks", "ContactMap")
    hbonds = config.get("Tasks", "HydrogenBonds").split(",")
    if hbonds[0].strip(" ") == "None":
        # values are assumed to be strings when get() is used
        hbonds = None
    perform_kbp = config.getboolean("Tasks", "Kbp")
    perform_jackknife = config.getboolean("Tasks", "Jackknife")
    
    # parse [Files] section
    traj = config.get("Files", "Traj")
    refpdb = config.get("Files", "RefPdb")
    hb_ad_file = config.get("Files", "AcceptorsDonors")
    if hb_ad_file.strip(" ") == "None":
        # values are assumed to be strings when get() is used
        hb_ad_file = None
    else:
        # get absolute path to avoid relative paths issues
        hb_ad_file = os.path.abspath(hb_ad_file)
    
    # parse [Options] section
    dcutsraw = config.get("Options", "Dcuts")
    if "range" in dcutsraw:
        # parse a string defining a range as either Python 
        # or numpy range and convert it into the corresponding range
        dcutsraw = \
            [float(num) for num in dcutsraw.strip("(numpy.range)").split(",")]
        dcuts = np.arange(*dcutsraw)
    else:
        # parse a string defining a range as a comma-separated
        # list of values
        dcuts = [float(num) for num in dcutsraw.split(",")]

    nsamplings = config.getint("Options", "Nsamplings")
    pcut = config.getint("Options", "Pcut")
    ff = config.get("Options", "ForceField")
    k = config.getint("Options", "K")
    hc_res = config.get("Options", "HcRes").split(",")
    nprocesses = config.getint("Options", "Nprocesses")


    ########################## LOGGER SETUP ###########################
    
    loglevel = logging.INFO
    logging.basicConfig(filename = "pyinknife.log", \
                        level = loglevel, \
                        filemode = "w")


    ######################## DIRECTORY SETUP ##########################

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

    # Convert o absolute paths to avoid relative paths issues
    refpdb = os.path.abspath(refpdb)
    traj = os.path.abspath(traj)

    # Create a new directory for the analyses into
    # the current working directory and move to it
    pydir = os.path.join(os.getcwd(), "pyinteraph_analyses")
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

    fulltrajdir = os.path.join(pydir, "fulltraj")
    os.mkdir(fulltrajdir)
    os.chdir(fulltrajdir)

    # create a list of trajectories and working directories
    # to be used to perform the requested analyses for both
    # the full trajectory and possibly the resampled ones
    trajs = [traj]
    wdirs = [fulltrajdir]


    ######################## JACKKNIFE RESAMPLING #########################

    if perform_jackknife:        
        # move to the main analyses directory
        os.chdir(pydir)
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

        for curr_resampl in range(nsamplings):
            logstr = "Performing resampling {:d}..."
            logging.info(logstr.format(curr_resampl))

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

            logstr = "Creating sub-trajectory {:s}..."
            logging.info(logstr.format(subtraj))          

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

                    logstr = \
                        "...frames 0:{:d} written to the sub-trajectory"
                    logging.info(logstr.format(end_first_part-1))

                # there is no "second part" if the excluded frames
                # are those at the end of the trajectory
                if curr_resampl != 0:
                    for ts in u.trajectory[begin_second_part:lastframendx]:
                        Wxtc.write(u.atoms)

                    logstr = \
                        "...frames {:d}:{:d} written to the sub-trajectory"
                    logging.info(logstr.format(\
                        begin_second_part, lastframendx-1))

            # append both the sub-trajectory and the directory in which
            # it is to the corresponding lists
            trajs.append(os.path.join(resampldir, subtraj))
            wdirs.append(resampldir)


    ############################ CONTACT MAP ############################

    if perform_cm:
        logging.info(\
            "########## STARTING CONTACT MAP ANALYSIS...##########")

        # Inform the user about which residues will be used while
        # building the contact map (should correspond to those
        # defined in the configuration file)
        logstr = \
             "Residue types considered while building " \
             "the contact map: {:s}"
        logging.info(logstr.format(",".join(hc_res)))

        # directories in which to perform the contact map analysis
        wdirs_cm = [os.path.join(_dir, "contact") for _dir in wdirs]

        for contactdir, traj in zip(wdirs_cm, trajs):
            logstr = "### ...for trajectory {:s} ###"
            logging.info(logstr.format(traj))   
            
            # Create the directory for the contact map analysis
            # and move to it
            os.mkdir(contactdir)
            os.chdir(contactdir)

            # parallelize the contact map analysis
            pool = multiprocessing.Pool(nprocesses)
            partfunc = \
                functools.partial(run_cm_pipeline,
                                  wdir = contactdir, \
                                  traj = traj, \
                                  refgro = refgro, \
                                  refpdb = refpdb, \
                                  ff = ff, \
                                  hc_res = hc_res, \
                                  pcut = pcut, \
                                  k = k, \
                                  hc_graph = "contact.dat", \
                                  filtered_graph = "filtered_graph.dat", \
                                  pdb_hubs = "hubs.pdb", \
                                  pdb_cc = "cc.pdb", \
                                  outfile_contact = "contact.out", \
                                  outfile_hubs = "hubs.out", \
                                  outfile_cc = "cc.out")

            pool.map(partfunc, dcuts)
            pool.close()
            pool.join()

    # go back to the analyses directory
    os.chdir(pydir)


    ########################### HYDROGEN BONDS ############################

    # check if hydrogen bonds analyses were requested
    if hbonds is not None:
        logging.info(\
            "########## STARTING HYDROGEN BONDS ANALYSIS...##########")

        logstr = "Using {:s} file for acceptors and donors definition"
        if hb_ad_file is None:
            logging.info(logstr.format("default"))
        else:
            logging.info(logstr.format(hb_ad_file))
        
        # directories in which to perform the contact map analysis
        wdirs_hb = [os.path.join(_dir, "h_bonds") for _dir in wdirs]  

        for hbondsdir, traj in zip(wdirs_hb, trajs):
            logstr = "### ...for trajectory {:s} ###"
            logging.info(logstr.format(traj))   
            
            # Create the directory for the hydrogen bonds analysis
            # and move to it
            os.mkdir(hbondsdir)
            os.chdir(hbondsdir)    

            # run the requested hydrogen bonds analyses in parallel
            pool = multiprocessing.Pool(nprocesses)
            partfunc = \
                functools.partial(run_hb_pipeline, \
                                  wdir = hbondsdir, \
                                  traj = traj, \
                                  refgro =  refgro, \
                                  refpdb = refpdb, \
                                  ff = ff, \
                                  pcut = pcut, \
                                  k = k, \
                                  hb_graph = "contact.dat", \
                                  filtered_graph = "filtered_graph.dat", \
                                  pdb_hubs = "hubs.pdb", \
                                  pdb_cc = "cc.pdb", \
                                  outfile_hb = "hb.out", \
                                  outfile_hubs = "hubs.out", \
                                  outfile_cc = "cc.out", \
                                  hb_ad_file = hb_ad_file)

            pool.map(partfunc, hbonds)
            pool.close()
            pool.join()

    # go back to the analyses directory
    os.chdir(pydir)


    ###################### KNOWLEDGE-BASED POTENTIAL ######################

    if perform_kbp:
        
        logging.info(\
            "########## STARTING KNOWLEDGE-BASED POTENTIAL" \
            "CALCULATION...##########")
        
        # move to the main analyses directory
        os.chdir(pydir)
        # create a subdirectory for knowledge-based
        # potential analysis and move to it
        kbpdir = os.path.join(pydir, "kbp")
        os.mkdir(kbpdir)
        os.chdir(kbpdir)

        kbp = PyInteraphCaller("pyinteraph", \
                               s = "../" + refgro, \
                               t = "../../" + traj, \
                               r = "../../" + refpdb, \
                               p = True, \
                               kbp_graph = "kbp.dat")

        kbp.run(stdout = "kbp.out", \
                stderr = "kbp.err")

    # go back to the analyses directory
    os.chdir(pydir)


