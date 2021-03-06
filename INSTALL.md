# Installing PyInKnife2

PyInKnife2 is a Python package.

In order to install PyInKnife2 you will need **Python 3.7 or higher**, **several open-source Python packages** and **PyInteraph2**.

Python is usually available by default in any Linux installation, while the Python packages will be installed automatically during the installation process. Please see the `requirements.txt` file or the `environment.yaml` file for the complete list of Python packages required.

PyInteraph2 can be found [here](https://github.com/ELELAB/pyinteraph2). PyInteraph2 needs to be installed before proceeding with the installation of PyInKnife2. Detailed installation instructions regarding how to install both PyInteraph2 and PyInKnife2 are provided below.

## Installation instructions

Here, we provide instructions for installing PyInteraph2 and PyInKnife2 in a simple Python environment or in a `conda` environment.

### Installing in a Python environment with `pip` and `setuptools`

This is a simple installation procedure that requires a Python >= 3.7 distribution to be available. It guides you in the installation of the PyInKnife2 package in a virtual environment, meaning an instance of Python isolated from the rest of your system.

This is not strictly necessary and PyInKnife2 may be installed system-wide in a similar fashion, following steps 4 to 9.

As you will see from the instructions, PyInteraph2 is installed before PyInKnife2.

#### Step 1 - Install `virtualenv`

Check if the `virtualenv` Python package is installed in your system.

If the `virtualenv` command is available to you, it is installed and you can skip this step. 

If it is not available and you need to install it, it is usually available as a package in your distribution. For instance, on Debian-based systems (such as Debian or Ubuntu) this boils down to installing the `python-virtualenv` package:

```shell
sudo apt install python-virtualenv
```

If this is not possible for you, you may still install the `virtualenv` package for just your local user, using `pip`:

```shell
pip install --user virtualenv
```

If the installation is successful, the `virtualenv` command will be available.

#### Step 2 - Create the virtual environment

Create your Python 3.7 virtual environment, in a directory of your choice (in this case, it will be `./pyinknife2-env`):

```shell
virtualenv -p /usr/bin/python3.7 pyinknife2-env
```

You might need to replace the argument of option `-p` according to the location of your Python installation in the system.

#### Step 3 - Activate the environment

Activate the environment:

```shell
source pyinknife2-env/bin/activate
```

#### Step 4 - Get PyInteraph2

Clone the PyInteraph2 source code from its GitHub repository and enter the local copy of the repository.

```shell
git clone https://github.com/ELELAB/pyinteraph2.git
cd pyinteraph2
```

If `git` is not available to you, you can download the repository content as a ZIP file from the PyInteraph2 GitHub repository web page.

#### Step 5 - Install the required Python packages for PyInteraph2

Install all the required Python packages as specified in the `requirements.txt` file in the `pyinteraph2` directory via `pip`:

```shell
pip install -r requirements.txt
```

#### Step 6 - Install PyInteraph2

You can now install PyInteraph2:

```shell
python setup.py install
cd ..
```

PyInteraph2 is now installed and you can proceed to the installation of PyInKnife2.

#### Step 7 - Get PyInKnife2

Clone the PyInKnife2 source code from its GitHub repository and enter the local copy of the repository.

```shell
git clone https://github.com/ELELAB/PyInKnife2.git
cd PyInKnife2
```

If `git` is not available to you, you can download the repository content as a ZIP file from the PyInKnife2 GitHub repository web page.

#### Step 8 - Install the required Python packages for PyInKnife2

Install all the required Python packages as specified in the `requirements.txt` file in the `PyInKnife2` directory via `pip`:

```shell
pip install -r requirements.txt
```

#### Step 9 - Install PyInKnife2

You can now install PyInKnife2:

```shell
python setup.py install
```

That's it! The main PyInKnife2 executables (`pyinknife_run`, `pyinknife_aggregate`, `pyinknife_plot`) should now be available.

Every time you need to run PyInKnife2 after opening a new shell, just run step 3 beforehand.

### Installing with `conda`

#### Step 1 - Get PyInKnife2

Clone the PyInKnife2 source code from its GitHub repository and enter the local copy of the repository.

```shell
git clone https://github.com/ELELAB/PyInKnife2.git
```

If `git` is not available to you, you can download the repository content as a ZIP file from the PyInKnife2 GitHub repository web page.

#### Step 2 - Install `conda`

We have provided an `environment.yaml` file that can be used together with `conda` to automatically install an environment contaning all of PyInteraph2 and PyInKnife2's requirements.

In order to use it you need to have `conda` installed in your system. Go [here](https://docs.conda.io/en/latest/miniconda.html) for further instructions. Installing `miniconda` rather than full `anaconda` is advised.

Once `conda` is installed on your system, you can use it to create a virtual environment similarly to what you would do using the `virtualenv` package, as previously detailed.

#### Step 3 - Create the `conda` environment

Create your `conda` environment starting from the provided environment file:

```shell
conda env create --file PyInKnife2/environment.yaml --prefix ./pyinknife2-env
```

In this case we are asking conda to create the environment locally (`--prefix`). This is not strictly necessary.

Note that if your home directory is size-constrained the installation might fail for lack of space. In that case, you need to select another directory with more space available in which conda will download packages. This is done by running:

```shell
conda config --add pkgs_dirs /my/other/directory/conda/pkgs
```

and changing `/my/other/directory/conda/pkgs` to a directory of your choice.

#### Step 4 - Activate the environment

Activate the `conda` environment by running the command line that conda suggests for this purpose at the end of the previous step. It is usually something like:

```shell
conda activate ./pyinknife2-env
```

#### Step 5 - Get PyInteraph2

Clone the PyInteraph2 source code from its GitHub repository and enter the local copy of the repository.

```shell
git clone https://github.com/ELELAB/pyinteraph2.git
cd pyinteraph2
```

If `git` is not available to you, you can download the repository content as a ZIP file from the PyInteraph2 GitHub repository web page.

#### Step 6 - Install PyInteraph2

You can now install PyInteraph2:

```shell
python setup.py install
```

PyInteraph2 is now installed and you can proceed to the installation of PyInKnife2.

#### Step 7 - Install PyInKnife2

You can now install PyInKnife2:

```shell
cd ..
cd PyInKnife2
python setup.py install
```

That's it! The main PyInKnife2 executables (`pyinknife_run`, `pyinknife_aggregate`, `pyinknife_plot`) should now be available.

Every time you need to run PyInKnife2 after opening a new shell, just run step 4 beforehand.