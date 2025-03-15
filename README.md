# Penning-trap-resonance

The code in this repository is a simulation of a Penning trap. There are also tools for analyzing the results and a report written in LaTeX.

## This repository contains the following code: 

### C++:
- Particle.cpp: Is an object that stores four attributes associated with physical particles that are of interest while simulating a Penning trap. All of the variables are defined using comments in the code. The armadillo library is used in the code.

- PenningTrap.cpp: This is an object that stores seven attributes associated with a Penning trap, all of which are described using comments in the code. The object has (excluding the constructor), due to compatibility reasons, seven public methods, but the user will need only three of them. These are: (i) solveRK4 - finds time evolution using RungeKutta4 method, (ii) solveFR - finds time evolution using ForwardEuler method, (iii) particlesInsideTrap - returns the number of particles still trapped at the end of the simulation. Some of the methods take in variables, but all are explained in the code. There are also multiple places in the code with ideas for optimization. The armadillo library is used in the code.

- main_two_particles.cpp: This is the main file meant to run a simulation of one or two particles in a Penning trap. This code contains multiple variables that can change the properties of both the particles and the Penning trap. All such variables do have comments explaining them and defining units. The code creates text files containing the time, position, and velocity information of the particles. This code needs to be linked with Particle.o and PenningTrap.o. The armadillo library is used in the code, which, while linking using g++ requires including -larmadillo flag while linking.

- main_100_particles.cpp: This is the main file, similar to the previously mentioned, except for these three points: (i) It's designed to simulate more particles (it's set to 100 by default). (ii) The particles are assigned random initial values. (iii) This code creates files containing the number of particles as a function of the oscillation frequency of the electric field. Building this program is identical to the main_two_particles.cpp.


### Python:
- error_estimation.py: Contains all the calculations and plotting of all the error-related values. This code requires some specific text files that are not contained in this repository. These can be obtained from the C++ code.

- plot_results.py: This is a plotting tool. There are six different sections, each creating and storing one or more plots. This code requires text files that are not contained in this repository but may be created using code in this repository.


### LaTeX:
- main.tex: This is the source code for the report I handed in. Building this code probably requires figures that are not stored in this repository but may be created using code in this repository.

- ref.tex: This is a file that really should have held my references, but unfortunately,I didn't manage to get it working. The file contains code from the course page, and I have included it just as a proof that I have not forgotten. I was just not able to get it working. 
