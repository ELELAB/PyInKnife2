# PyInKnife 2

## Overview

PyInKnife2 is a set of Python scripts to run an updated version of the PyInKnife pipeline [^salamancaviloria2017] for Protein Structure Network (PSN) analysis using PyInteraph [^tiberti2015].

## Background

### Protein Structure Network (PSN)



## Requirements

The following Python requirements must be met:

* `python` v3.7 or higher

The following Python packages must be installed:

* `dask` and `dask distributed` v2.21.0 or higher
* `MDAnalysis` v0.20.0 or higher
* `numpy` v1.17.0 or higher
* `PyYAML` v5.3.1 or higher

## Installation

The scripts require no installation.

## Usage

### pyinknife.py

This is the script responsible for running the PyInKnife pipeline.

#### Command line

`pyinknife.py [-h] -f TRJ -s TOP [-r REF] -c CONFIGFILE [-d RUNDIR] [-n NPROC] [-ncaa [NONCANONICAL_RESIDUES [NONCANONICAL_RESIDUES ...]]]`

#### Options

| Option                             | Meaning                                                      |
| ---------------------------------- | ------------------------------------------------------------ |
| `-h`, `--help`                     | Show the help message and exit.                              |
| `-f ` `--trj`                      | Trajectory.                                                  |
| `-s`, `--top`                      | Topology.                                                    |
| `-r`, `--ref`                      | Reference structure.                                         |
| `-c`, `--configfile`               | Configuration file.                                          |
| `-d`, `--rundir`                   | Directory where to run the pipeline. Default is the current working directory. |
| `-n`, `--nproc`                    | Number of processes to be started in parallel. Default is 1 process(es). |
| `-ncaa`, `--noncanonical-residues` | Noncanonical residues present in your system.                |

#### Input files

##### Configuration file

A YAML file containing the script configuration (please see the `config.yaml` file for an example of configuration file).

#### Outputs

## References
