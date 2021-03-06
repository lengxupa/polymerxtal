PolymerXtal
==============================
[//]: # (Badges)
[![Travis Build Status](https://travis-ci.com/lengxupa/polymerXtal.svg?branch=master)](https://travis-ci.com/lengxupa/polymerXtal)
[![codecov](https://codecov.io/gh/lengxupa/polymerXtal/branch/master/graph/badge.svg)](https://codecov.io/gh/lengxupa/polymerXtal/branch/master)


In this work, we present PolymerXtal, a software designed to build and analyze molecular-level polymer crystal structures. PolymerXtal provides a standardized process to generate polymer crystal structure based on monomer, tacticity, helicity, chiriality and unit cell information and analyze the crystallinity in polymer systems with given atom trajectories. These features have allowed PolymerXtal to lead further investigations of semi-crystalline polymers where the birthplace of the important physics in play and promote future research endeavors in the area of crystalline polymer simulations.

This repository is currently under development. To do a developmental install, download this repository and type

`pip install -e .`

in the repository directory.

This package requires the following:
  - numpy
  - matplotlib
  - ovito
  - openbabel
  - scipy

Features should be developed on branches. To create and switch to a branch, use the command

`git checkout -b new_branch_name`

To switch to an existing branch, use

`git checkout new_branch_name`

To submit your feature to be incorporated to the master branch, you should submit a `Pull Request`. The repository maintainers will review your pull request before accepting your changes.

### Copyright

Copyright (c) 2020, Tongtong Shen


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.3.
