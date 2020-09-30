#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

from setuptools import setup

name = "PyInKnife"
url = "https://github.com/ELELAB/PyInKnife2"
author = "Valentina Sora, Juan Salamanca Viloria, Matteo Tiberti, Elena Papaleo"
author_email = "sora.valentina1@gmail.com"
version = "2.0.0"
description = "PyInKnife pipeline for Protein Structure Network analysis with PyInteraph"
package_data = {"PyInKnife" : ["examples/*", "config/*"]}
package_dir = {"PyInKnife" : "PyInKnife"}
packages = ["PyInKnife"]
entry_points = {"console_scripts" : \
                    ["pyinknife_run = PyInKnife.pyinknife_run:main", \
                     "pyinknife_aggregate = PyInKnife.pyinknife_aggregate:main", \
                     "pyinknife_plot = PyInKnife.pyinknife_plot:main"], \
               }
# MDAnalysis 1.0.0 will be required till PyInteraph supports the hydrogen
# bonds analysis in the newer versions
install_requires = ["matplotlib", \
                    "MDAnalysis==1.0.0", \
                    "numpy", \
                    "pandas", \
                    "pyyaml", \
                    "seaborn"]

setup(name = name, \
      url = url, \
      author = author, \
      author_email = author_email, \
      version = version, \
      description = description, \
      package_data = package_data, \
      package_dir = package_dir, \
      packages = packages, \
      entry_points = entry_points, \
      install_requires = install_requires)