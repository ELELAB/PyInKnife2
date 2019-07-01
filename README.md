# PyInKnife2.py

Updated version of the **PyInKnife** pipeline (the older version can be found [here](https://github.com/ELELAB/PyInKnife)) to perform network generation from different types of noncovalent interactions and analyses with **PyInteraph** on a given trajectory (see the [relevant paper](https://www.nature.com/articles/s41598-017-01498-6) for further details) with the possibility of performing a jackknife resampling on the trajectory and of generating and analyzing the networks also from the resampled trajectories. The PyInteraph source code and relevant documentation can be found [here](https://github.com/ELELAB/pyinteraph).



## Requirements

`python > 2.7`

`numpy`

`MDAnalysis`  

`pyinteraph`



## Usage

### The configuration file

All settings and parameters necessary for running the pipeline are contained in a **configuration file** (i.e. `config_file.cfg`) formatted as follows:

```
[ENVIRONMENT]
PYINENV = /usr/local/envs/pyinteraph/bin/activate_this.py
NPROC = 4

[ANALYSES]
HC = True
HB = True
SB = False

[INPUT_FILES]
TRAJ = ../../thesis/psn_tests/traj_bclxl_cap_frames0-100.xtc
PDB = ../../thesis/psn_tests/model0_cap.pdb

[GENERAL_OPTIONS]
FORCE_FIELD = charmm27
PCUT = 20
K = 3
NSAMPLINGS = 10

[HC_OPTIONS]
HC_DCUTS = range(5.0,6.0,0.5)
HC_RES = ALA,ILE,LEU,VAL,PHE,TRP,TYR,MET,PRO,HIS,LYS,ARG,GLU,ASP,ASN,GLN,SER,THR,CYS,CSN,SP2

[HB_OPTIONS]
HB_CLASSES = mc-mc, sc-sc
HB_DCUTS = 1.0,1.2
HB_FILE = hydrogen_bonds_mod.ini

[SB_OPTIONS]
SB_MODES = different_charge, all
SB_DCUTS = 4.0, 4.5
SB_FILE =
```

The `[ENVIRONMENT]` section defines:

* `PYINENV`, the path to the Python script activating the PyInteraph virtual environment on your machine.
* `NPROC` is the number of processes the script will use to parallelize the calculations. Be aware that the number of processes spawned DOES NOT CORRESPOND to the number of CPU cores used. This script use the Python `multiprocessing` module for parallelization, and for this reasion the load distribution is handled directly by the OS of your machine. Please refer to the official documentation of the multiprocessing module [here](https://docs.python.org/3/library/multiprocessing.html?highlight=multiprocessing#module-multiprocessing) for further details.



The `[ANALYSES]` section defines which analyses will be performed:

* `HC` accepts either `True` or `False` and defines whether (a) network(s) based on hydrophobic contacts will be generated.
* `HB` accepts either `True` or `False` and defines whether (a) network(s) based on hydrogen bonds will be generated.
* `SB` accepts either `True` or `False` and defines whether (a) network(s) based on salt bridges will be generated.



The `[INPUT_FILES]` section defines the input files:

* `TRAJ` defines the path to the trajectory. The script assumes it has already been pre-processed (PBC artifacts removed, system centered in the box, no solvent, etc.). It can be either a relative or an absolute path.
* `PDB` defines the path to the reference PDB structure (must match the system contained in the trajectory). It can be either a relative or an absolute path.



The `[GENERAL_OPTIONS]` defines some settings that will be valid to all analyses:

* `FORCE_FIELD` is the name of the force field from the PyInteraph internal database from which the atom masses will be taken (used during the generation of the networks).
* `PCUT` is an integer defining the persistence cut-off to be used during the network filtering (the same will apply to all networks generated).
* `K` is an integer defining the minimum number of edges (=interactions) a node (=a residue) must have to be reported as hub.
* `NSAMPLINGS` is an integer defining how many samplings will be performed during the jackknife resampling procedure. Leave this field empty if you do not want to perform any resampling.



The `[HC_OPTIONS]` section defines some specific settings for the networks generated from hydrophobic contacts:

- `HC_DCUTS` defines the range of distance cut-offs to be used for the generation of the network. It can be defined either as a comma-separated list of values or with a Python-like definition of `range()` (for ranges of integers) or `numpy.arange()` (for ranges of floats).  It will be ignored if `HC` is `False`.
- `HC_RES` is a comma-separated list of residues that are considered while building the network.  It will be ignored if `HC` is `False`.



The `[HB_OPTIONS]` section defines some specific settings for the networks generated from hydrogen bonds:

* `HB_CLASSES` accepts any combination of `mc-mc`, `sc-sc`, `mc-sc`, `all` which defines which type(s) of hydrogen bonds network(s) will be generated (with main chain - main chain hydrogen bonds, with side chains - side chains hydrogen bonds, with main chain - side chains hydrogen bonds or with all hydrogen bonds). It will be ignored if `HB` is `False`.

* `HB_DCUTS` defines the range of distance cut-offs to be used for the generation of the network. It can be defined either as a comma-separated list of values or with a Python-like definition of `range()` (for ranges of integers) or `numpy.arange()` (for ranges of floats).  It will be ignored if `HB` is `False`.

* `HB_FILE` the path to a .ini files containing the definition of hydrogen bonds donors and acceptors that will be used by PyInteraph during the hydrogen bonds network(s) generation. Please refer to the PyInteraph documentation for further information on the format and usage of this file. If left empty, the default file contained in the PyInteraph installation directory will be used.  It will be ignored if `HB` is `False`.

  

The `[SB_OPTIONS]` section defines some specific settings for the networks generated from salt bridges:

* `SB_MODES` accepts any combination of `different_charge`,  `all` which defines which type(s) of salt bridges network(s) will be generated (defining a salt bridge only between residues with atom groups with different charges or between all possible pairs of charged groups contained in the .ini file for charged groups, see `SB_FILE`). It will be ignored if `SB` is `False`.

- `SB_DCUTS` defines the range of distance cut-offs to be used for the generation of the network. It can be defined either as a comma-separated list of values or with a Python-like definition of `range()` (for ranges of integers) or `numpy.arange()` (for ranges of floats).  It will be ignored if `SB` is `False`.
- `SB_FILE` the path to a .ini files containing the definition of charged groups of atoms that will be used by PyInteraph during the salt bridges network(s) generation. Please refer to the PyInteraph documentation for further information on the format and usage of this file. Must be provided if `SB` is `True`,  it will be ignored if `SB` is `False`.



### How to execute

`python PyInKnife2.py -f config_file.cfg` 



