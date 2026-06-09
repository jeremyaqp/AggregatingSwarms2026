# Paper reference


Repository for the data and code used in :

Jeremy Fersula, Nicolas Bredeche, & Olivier Dauchot. (2026). **To be updated**

Simulations were made using the LAMMPS software, see https://github.com/lammps/lammps.

## Abstract

Morphological computing, the use of the physical design of a robot to ease the realization of a given task has been proven to be a relevant concept in the context of swarm robotics. Here we demonstrate both experimentally and numerically, that the success of such a strategy may heavily rely on the type of policy adopted by the robots, as well as on the details of the physical design. To do so, we consider a swarm of robots, composed of Kilobots embedded in an exoskeleton, the design of which controls the propensity of the robots to align or anti-align with the direction of the external force they experience. We find experimentally that the contrast that was observed between the two morphologies in the success rate of a simple phototactic task, where the robots were programmed to stop when entering a light region, becomes dramatic, if the robots are not allowed to stop, and can only slow down. Building on a faithful physical model of the self-aligning dynamics of the robots, we perform numerical simulations and demonstrate on one hand that a precise tuning of the self-aligning strength around a sweet spot is required to achieve an efficient phototactic behavior, on the other hand that exploring a range of self-alignment strength allows for a rich expressivity of collective behaviors. 
Repository for the data and code used in "Leveraging Design for Collective Phototaxis in Morphological Swarm Robotics".

# Software dependencies

* Python 3
* Git
* LAMMPS https://github.com/lammps/lammps
* C++ compiler toolchain for LAMMPS compilation (make, g++)
* MPI if using parallel computing in LAMMPS

## Tested Versions

* Ubuntu 22.04.4 LTS
* gcc version 11.4.0
* Python 3.10.12

# Installation Guide

Install the Python3 packages necessary with the provided `requirements.txt` file :

```
python3 -m pip install -r requirements.txt
```

## Figures

Extract experimental and simulated date from the provided archive

```
tar -xvf data.tar.gz Data/
```

The figures from the papers can be generated using the scripts provided in `Figures/src`.

## Simulations

Install LAMMPS from the repository :

```
git clone https://github.com/lammps/lammps
```

Revert LAMMPS to the exact commit used in the development of the plugin :

```
# From the lammps cloned repository
git reset --hard d618b0ffc05dfd86915c0d148c0a72fba995eba4 
```

Copy the files from `Simulations` to `lammps/src` and enable BROWNIAN package 

```
cp ./Simulations/LammpsPlugin/* *lammps_location*/src
# From the lammps cloned repository
make yes-brownian
make serial # compile serial lammps
make mpi # compile parallel lammps
```

These command will build two binaries of lammps : `lmp_serial` and `lmp_mpi`. You can run the simulations with the following commands :

```
# Serial execution
lmp_serial -in ./Simulations/in.phototaxis64.run
# Parallel execution (8 processes)
mpiexec -n 8 lmp_mpi -in ./Simulations/in.phototaxis64.run
```

The resulting file has extension `.cfg`. This file can either be opened with typical vizualiation tools such as Ovito (https://www.ovito.org/), or post-processed to `.csv` using the provided script `cfg_to_csv.py`.

# Installation and running time

Provided the exact gcc version is installed, typical installation time should not exceed 1 hour for a non-lammps expert.

Running time of simulations heavily depends on hardware, can be expected to be a few seconds per seed if using parallel computing on a "good computer".


# License

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
