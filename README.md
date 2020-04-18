# Point Gravity (Python)
### Author: Charlie Hagedorn
### Scribe: John G. Lee

## Introduction

This python implementation of PointGravity is based on the work of Dr. Charlie 
Hagedorn which can be found [here](https://github.com/4kbt/PointGravity)!

This package allows one to calculate forces and torques from extended bodies 
due to a gravitational or Yukawa interaction. It does so in the most simple 
(although not necessarily most efficient) manner imaginable. We simply
discretize our shapes with a 3-dimensional array of points. We can then 
compute the force or torque between each unique pair of points and sum up the 
total for the entire body. Gravity depends only linearly on the source masses, 
so many complicated shapes can be modeled by just summing up individual simple 
components.

We can also simply visualize the point mass arrays using a 3-d plot.

The second portion of this package is the implementation of a multipole 
analysis for the same calculations. The gravitational potential can be 
decomposed into interactions of multipole moments allowing for accurate and 
fast calculations from relatively few low order moments. We compute several 
low order moments of basic shapes or estimate the moments from a point-mass 
array. We can rotate these moments to large orders *cite Gumerov* and 
translations coming soon. This allows us to compute the interactions in an 
entirely different and often useful perspective. This portion is based largely 
on the private FORTRAN code of Prof. Adelberger, MULTIN.


## To Do
- [X] Multipole rotations
- [X] Multipole torques
- [] Multipole forces
- [] Multipole class?
- [] Multipole translation
- [X] More tests
- [] pip package
- [] More Doc Strings
- [] Pull request to Charlie's Octave version?
- [X] Outer Multipoles from point-mass
- [] Example visualization