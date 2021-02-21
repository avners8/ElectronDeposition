# ElectronDeposition
Avner Shultzman, The Weizmann Institute, E-mail: avner.shultzman@weizmann.ac.il

Efficienct optimization of electron deposition

## Introduction
In this work we aim to optimize the dose deposited by multiple electron beams.
The first step of the project is to efficiently compute the dose using an analytical formalism and some approximations.
Secondly, use the pencil beam metgod to approximate the 2D or 3D profile of the beam using the obtained 1D profile.
And finally, perform optimization on the deposition according to some figure-of-merit and some constraints.

## The Kinetic Model and the M1 Approximation
The first step of the project, is to efficiently compute the deposited dose from an electron beam. The standard methods in use today, is Monte-Carlo simulation, which simulate the transport of electrons in the material through a stochastic algorithm. However, this method requires a significant amount of computation power and time, and therefore is not applicable in some medical applications.
In this project, we follow the derivation made in [1], to analytically compute the dosage in the material. The motion of electrons in the material can be described by the fluence \psi, which satisfies the kinetic equation.
To efficiently compute the fluence we use the moment model. According to this model, we approximate the fluence by its first moments. In practice, this model requires a closure since the system of equations is not uniquely satisfied. In this work we use the <img src="https://latex.codecogs.com/gif.latex?M_{N} t " /> closure which minimizes the entropy of the obtained fluence. We implemented the <img src="https://latex.codecogs.com/gif.latex?M_{1} t " /> approximation which compute first two moments and compute the next one according to this closure.
The fluence is computed with an iterative algorithm until convergence.

## The Pencil Beam Method
We used the pencil beam method in order to expand the 1D profile we obtained to a 2D profile, following the work in [2]. The idea is to expand the 1D profile to a certain width and with an angular distribution. We opt to compute the profil using the pencil beam approach and not with a kinetic model approximation, since the cumpute time are much smaller and therefore suitable for optimization. Furthermore, since we are working with electron beams we can consider focused beams. The implamantation of a focused beam is done by combining pencil beams with a small width which are all pointing to the same focal point. As expected, these beams spread less in the material. 

## Optimization
The final step is to use constraints optimization algorithms in order to maximize a figure of merit. For example, we want to maximize the amount of dose in a small area in the material, while minimizing the dose in all other regions. The optimization tool enables us to modify the initial location and orientation of the beams. The use of efficient optimization techniques is essential, since we are optimizing over a large parameter space, and therefore we cannot afford an exhaustive search.


## References
1. Teddy Pichard, "A numerical approach for a system of transport equations in the field of radiotherapy", Global Science Press, 2017
2. Y. Glinec, J. Faure, and V. Malka, "Radiotherapy with laser-plasma accelerators: Monte Carlo simulation of dose deposited by an experimental quasimonoenergetic electron beam", Med Phys, 2006
