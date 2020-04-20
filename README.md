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

## References
https://github.com/4kbt/PointGravity
[A Sub-Millimeter Parallel-Plate Test of Gravity][1]
[Translation of multipoles for a 1/r potential][2]
[Interaction potential between extended bodies][3]
[Recursive computation of spherical harmonic rotation coefficients of large degree][4]
[Multipole calculation of gravitational forces][5]


[1]: https://digital.lib.washington.edu/researchworks/handle/1773/34135
[2]: https://journals.aps.org/prd/abstract/10.1103/PhysRevD.55.7970
[3]: https://journals.aps.org/prd/abstract/10.1103/PhysRevD.60.107501
[4]: https://arxiv.org/abs/1403.7698
[5]: https://journals.aps.org/prd/abstract/10.1103/PhysRevD.95.124059


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