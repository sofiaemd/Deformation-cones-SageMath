# :small_blue_diamond: Deformation cones SageMath :small_blue_diamond:

This repository contains SageMath code to produce the deformation cone of a polytope based either in its Minkowski weights or its support function, and to construct the summands corresponding to the cone's rays.
The code was written using SageMath version 10.5.beta2 with Python 3.12.4

## Deformation cone

The deformation cone of a polytope is the set of weak Minkowski summands up to translation, this can be parametrized by 1-Minkowski weights or by the values the support function takes in the rays of the normal fan. The rays of this cone correspond to indecomposable weak summands of the original polytope. For more theoretical background refer to the thesis.

## Minkowski cone

The file called `Minkowski_cone.sage` provides a function `Minkowski_cone` that has a polytope as input and outputs a tuple that contains a cone, a list of edges and a list of vertices. The lists are necessary to reconstruct a weak summand from a point in the cone, because the ith coordinate in the space of the cone, represents the value the Minkowski weight takes on the ith edge on the list. The file also contains a function called `build_summand_mink` that accepts as input the cone and two lists provided by `Minkowski_cone` and the index of a ray in the cone and outputs the weak summand of the polytope generated using the Minkowski weight corresponding to a ray generator of the chosen ray.

## Support function cone

The file called `Support_functions_cone.sage` provides the function `Support_functions_cone` that accepts a polytope as input and outputs a tuple that contains a cone, and the rays of the normal fan of the polytope in corresponding order. The function also has the option of turning `wall_crossings=False` to use only the inequalities provided by Batyrev's criterion of `batyrev=False` to use only wall-crossing inequalities.
The script also provides a function called `build_summand_supp` that accepts the cone, the list of rays and the number of the ray to build and outputs the weak summand of the polytope generated using the values for the support function given by a generator of the chosen ray.

## Example

To use the functions, it is first necessary to define a polytope
```
from Minkowski_cone import Minkowski_cone, build_summand_mink
from Support_functions_cone import Support_functions_cone, build_summand_supp
P = polytopes.cube()
```
Then it is possible to build the cone of non-negative Minkowski weights and build a summand with
```
C = Minkowski_cone(P)
weak_summand = build_summand_mink(*C,0)
```
or the cone of support functions with
```
C = Support_function(P)
weak_summand = build_summand_supp(*C,0)
```
Since the cube is a simple polytope, we can choose to build the cone only using wall-crossing inequalities with
```
C = Support_function(P, batyrev=False)
```
or build the cone using only the inequalities provided by Batyrev's criterion with
```
C = Support_function(P, wall_crossings=False)
```


