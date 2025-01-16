# :small_blue_diamond: Deformation cones SageMath :small_blue_diamond:

This repository contains SageMath code to produce the deformation cone of a polytope based either in its Minkowski weights or its support function, and to construct the summands corresponding to the cone's rays.
The code was written using SageMath version 10.5.beta2 with Python 3.12.4

## Deformation cone

Un poquito de explicación teórica and a detailed explanation of the theoretical background can be found in the thesis.

## Minkowski cone

The file called `Minkowski_cone.py' provides a function `Minkowski_cone` that has a polytope as input and outputs a tuple that contains a cone, a list of edges and a list of vertices. The lists are necessary to reconstruct a weak summand from a point in the cone, because the ith coordinate in the space of the cone, represents the value the Minkowski weight takes on the ith edge on the list. The file also contains a function called `build_summand_mink' that accepts as input the cone and two lists provided by `Minkowski_cone' and the index of a ray in the cone and outputs the weak summand of the polytope generated using the Minkowski weight corresponding to a ray generator of the chosen ray. This means the weak summands produced by this function are indecomposable summands of the original polytope.

## Support function cone

gvggh

## Example

For example, if you would like to obtain an indecomposable weak summand of a cube, you could do the following

```
P = polytopes.cube()
C = Minkowski_cone(P)
weak_summand = build_summand_mink(*C,0)

´´´
To obtain the cone parametrized by support functions you would need to do the following,


