### Introduction

This project provides code for the associated paper "Proximal-based recursive implementation for model-free data-driven fault diagnosis" by J. Noom, O. Soloviev and M. Verhaegen, submitted to Automatica. 


### Installation

The code was written in Matlab R2021a on Windows 10 with a 64-bit operating system.

The Matlab scripts make use of the following toolboxes and packages:
- Control System Toolbox
- Signal Processing Toolbox
- SerDes Toolbox
- CVX (https://github.com/cvxr/CVX; only required for running `mbfd.m`)


### Usage
For batch-wise model-based fault diagnosis, run the file `mbfd.m`.

For batch-wise model-free data-driven fault diagnosis, run the file `mfddfd.m`.

For online recursive model-free data-driven fault diagnosis, run the file `rddfd.m`.


### Citations

For citing the methodology, please refer to the corresponding paper.

For citing the software, please use the DOI below.

[![DOI](https://zenodo.org/badge/738073814.svg)](https://zenodo.org/doi/10.5281/zenodo.10454000)