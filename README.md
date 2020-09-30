# PyInKnife 2

## Overview

PyInKnife2 is a set of Python scripts to run an updated version of the PyInKnife pipeline [^salamancaviloria2017] for Protein Structure Network (PSN) analysis using PyInteraph [^tiberti2015].

## Background

### Protein Structure Networks (PSNs)

Protein Structure Networks (PSNs) are a way to analyze a protein structure or ensemble of conformations through the mathematical concept of network. 

A network (or graph) is defined as a set of **nodes** connected by **edges** (or arcs) that represent some sort of relationship between each pair of nodes. A network can be weighted, meaning that specific weights are associated to the edges, or unweighted.

In PSNs, nodes are usually represented by the protein monomeric unit (i.e., a residue) and edges may represent different types of physical interactions between these units (such as salt bridges and hydrogen bonds) or may be computed according to some definition of "contact". PSNs are usually weighted networks.

If a PSN is calculated over an ensemble of conformations, a single network is computed for each conformation, and then all networks are merged into a single network. A way of assigning edges and edge weights in this aggregated network can be to average over the weight of an edge over the whole ensemble, or to compute the persistence of the edge in the ensemble (i.e., fraction of conformations in which the edge is present).

### PyInteraph

PyInteraph [^tiberti2014]  is a software for Protein Structure Network creation and analysis, where a PSN has residues as nodes and edges can be defined in different ways, for example representing a particular class of non-covalent interactions between residues (hydrophobic contacts, salt bridges or hydrogen bonds). 

### PyInKnife

PyInKnife [^salamancaviloria2015] is a pipeline built on top of PyInteraph to assess the robustness of the PSN analyses performed by PyInteraph on a given ensemble by running them not only on the full ensemble but also on subsets of it and comparing the networks from such resampling to the network built from the whole ensemble.

The default (and, so far, only) resampling method used is jackknife resampling.

## Requirements

The user must have `python` v3.7 or higher installed, together with the PyInteraph software (whose latest version can be found [here](https://github.com/ELELAB/pyinteraph2/tree/feature_py3)).

Required Python dependencies, if not already present, will be installed along with PyInKnife.

## Installation

To install PyInKnife, download and unzip this folder, enter the folder and run the following command:

`python3.7 setup.py install`

Upon successful installation, you should have three executable (`pyinknife_run`, `pyinknife_aggregate` and `pyinknife_plot`) available to perform the various steps of data collection and analysis.

## Usage

### pyinknife_run

This is the executable responsible for running the PyInKnife pipeline.

#### Command line

`pyinknife.py [-h] -f TRJ -s TOP [-r REF] -c CONFIGFILE [-d RUNDIR] [-n NPROC] [-ncaa [NONCANONICAL_RESIDUES [NONCANONICAL_RESIDUES ...]]]`

#### Options

`$INSTALLDIR` indicates the directory where you installed PyInKnife.

| Option                             | Meaning                                                      |
| ---------------------------------- | ------------------------------------------------------------ |
| `-h`, `--help`                     | Show the help message and exit.                              |
| `-f ` `--trj`                      | Trajectory.                                                  |
| `-s`, `--top`                      | Topology.                                                    |
| `-r`, `--ref`                      | Reference structure.                                         |
| `-c`, `--configfile`               | Configuration file to be used to run the pipeline. Default is `$INSTALLDIR`/PyInKnife2/PyInKnife/config/run.yaml. |
| `-d`, `--rundir`                   | Directory where to run the pipeline. Default is the current working directory. |
| `-n`, `--nproc`                    | Number of processes to be started in parallel. Default is 1 process(es). |
| `-ncaa`, `--noncanonical-residues` | Noncanonical residues present in your system.                |

#### Input files

##### Configuration file

A YAML file containing the script configuration (please see the `run.yaml` file in the `config` directory for an example of configuration file).

Please be aware that the following options, if passed to the corresponding command in the configuration file, will be removed during the pre-processing of the configuration, since they are either set internally or elsewhere in the configuration file:

| Command          | Ignored options                                              |
| ---------------- | ------------------------------------------------------------ |
| `pyinteraph`     | `-h`, `--help`, `-s`, `--top`, `-t`, `--trj` , `-r`, `--ref`, `--sb-mode`, `--hb-class`, `-hc-co`, `--hc-cutoff`, `-sb-co`, `--sb-cutoff`, `-hb-co`, `--hb-cutoff` |
| `filter_graph`   | `-d`, `--input-dat`, `-t`, `--filter-threshold`              |
| `graph_analysis` | `-r`, `--reference`, `-a`, `--adj-matrix`                    |

`--hc-perco` (`--hc-persistence-cutoff`), `--sb-perco` (`--sb-persistence-cutoff`), `--hb-perco` can be passed as `pyinteraph` options if you directly want a filtered network. Persistence cut-offs specified for `filter_graph` will be used to further filter this initial network.

#### Outputs

If you run with the provided configuration file in a generic test folder called "analysis", you will obtain the following directory tree representing your results:

```
analysis
|----| fulltrj
|----|----| hb
|----|----|----| sc-sc
|----|----|----|----| 3.5
|----|----|----|----|----| 50.0
|----|----|----|----|----|----| ccs
|----|----|----|----|----|----| hubs
|----|----| hc
|----|----|----| 5.0
|----|----|----|----| 50.0
|----|----|----|----|----| ccs
|----|----|----|----|----| hubs
|----|----|----| 5.5
|----|----|----|----| 50.0
|----|----|----|----|----| ccs
|----|----|----|----|----| hubs
|----|----| sb
|----|----|----| different_charge
|----|----|----|----| 4.5
|----|----|----|----|----| 50.0
|----|----|----|----|----|----| ccs
|----|----|----|----|----|----| hubs

|----| resampling0
(same tree as fulltrj)
|----| resampling1
(same tree as fulltrj)
|----| resampling2
(same tree as fulltrj)
```

The top level (`fulltrj`, `resampling*`) represents the trajectory on which the PSNs were constructed, therefore either the full trajectory or a subset of it.

The second level only appears in the hydrogen bonds (`hb`) and salt bridges (`sb`) analyses, and indicates how hydrogen bonds/salt bridges were defined. In this case, a salt bridge was defined only between residues having opposite charges (`different_charge`). In principle, one could define it as a more general interaction between all residues having a charge (`all`) or only between residues having the same charge (`same_charge`), even if it would not be a proper "salt bridge" anymore. in this example, only hydrogen bonds formed between side chains were considered (`sc-sc`).

The third level indicates the distance cut-off used to define the interaction of interest (expressed in Ångstroms), while the fourth one represents the persistence cut-off used to filter out transient interactions from the network.

The last level represents the analysis performed on the (filtered) network. So far, PyInKnife uses PyInteraph to find connected components (`ccs`) and hubs (`hubs`).

### pyinknife_aggregate

This executable is used to aggregate the results of the PSNs generated from the different trajectories.

#### Command line

`pyinknife_aggregate [-h] [-c CONFIGFILE] [-ca CONFIGFILE_AGGREGATE] [-d RUNDIR] [-od OUTDIR] [--firstccs FIRSTCCS]`

#### Options

| Option                          | Meaning                                                      |
| ------------------------------- | ------------------------------------------------------------ |
| `-h`, `--help`                  | Show the help message and exit.                              |
| `-c`, `--configfile`            | Configuration file used to run the pipeline. Default is `$INSTALLDIR`/PyInKnife2/PyInKnife/config/run.yaml. |
| `-ca`, `--configfile-aggregate` | Configuration file for data aggregation. Default is `$INSTALLDIR`/PyInKnife2/PyInKnife/config/aggregate.yaml. |
| `-d`, `--rundir`                | Directory where the pipeline was run. Default is the current working directory. |
| `-od`, `--outdir`               | Directory where to save the output files. Default is the current working directory. |
| `--firstccs`                    | First # most populated connected components to be considered. Default is 5. |

#### Outputs

If used with the default configuration files (provided that the pipeline was run with the default configuration file), `pyinknife_aggregate` produces the following outputs:

```
hb_sc-sc_3.5_50.0_ccs.csv
hb_sc-sc_3.5_50.0_hubs.csv
hc_5.0_50.0_ccs.csv
hc_5.0_50.0_hubs.csv
hc_5.5_50.0_ccs.csv
hc_5.5_50.0_hubs.csv
sb_different_charge_4.5_50.0_ccs.csv
sb_different_charge_4.5_50.0_hubs.csv
```

Each of them is a CSV file containing a dataframe where either the number of nodes in the most populated connected components or the degree distribution of the hubs is reported for each trajectory (full trajectories and resamplings), together with some across-trajectories statistics (mean, median, min, max, standard deviation and standard error of the mean).

### pyinknife_plot

This executable is used to plot the aggregated results.

#### Command line

`pyinknife_plot [-h] [-c CONFIGFILE] [-ca CONFIGFILE_AGGREGATE] [-cp CONFIGFILE_PLOT] [-p PLOT] [-d RUNDIR] [-od OUTDIR]`

#### Options

| Option                          | Description                                                  |
| ------------------------------- | ------------------------------------------------------------ |
| `-h`, `--help`                  | Show the help message and exit.                              |
| `-c`, `--configfile`            | Configuration file used to run the pipeline. Default is `$INSTALLDIR`/PyInKnife2/PyInKnife/config/plot/run.yaml. |
| `-ca`, `--configfile-aggregate` | Configuration file used for data aggregation. Default is `$INSTALLDIR`/PyInKnife2/PyInKnife/config/plot/aggregate.yaml. |
| `-cp`, `--configfile-plot`      | Configuration file used for plotting. Default depends on what you want to plot. |
| `-p`, `--plot`                  | What to plot. Default is hubs.                               |
| `-d`, `--rundir`                | Directory where the aggregate outputs are saved. Default is the current working directory. |
| `-od`, `--outdir`               | Directory where to save the output plots. Default is the current working directory. |

#### Outputs

`pyinknife_plot` will try to generate a plot for each file present in the directory containing the aggregated output files. Make sure you only have those files in the directory ot the script may either crash or produce erroneous plots. Aestehics for the plot are defined in the plot configuration file. Examples of configuration files for bar plots used to visualize aggregated data for hubs and connected components can be found in the `config` folder.

## References

[^salamancaviloria2015]: Viloria, Juan Salamanca, et al. "An optimal distance cutoff for contact-based Protein Structure Networks using side-chain centers of mass." *Scientific reports* 7.1 (2017): 1-11.
[^tiberti2014]: Tiberti, Matteo, et al. "PyInteraph: a framework for the analysis of interaction networks in structural ensembles of proteins." *Journal of chemical information and modeling* 54.5 (2014): 1537-1551.