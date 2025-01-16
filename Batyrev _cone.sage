# -*- coding: utf-8 -*-
r"""
Module for the computations of the deformation cones of polytopes

EXAMPLES::

REFERENCES:

    - G. Pineda-Villavicencio
      Polytopes and Graphs, Cambridge 2024.

AUTHORS:

- Sofia Errazuriz and Federico Castillo: Initial version
"""


##############################################################################
#     Copyright (C) 2023 Sofia Errazuriz <sofiaerrazurizm@uc.cl>
#                       Federico Castillo <efecastillo.math@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
##############################################################################

def Batyrev_cone(P):
    """
    INPUT: Simple polytope
    OUTPUT: Tuple containing the cone parametrizing support function of weak summands and the rays of the normal fan of the polytope.
    """
    normal_fan = NormalFan(-P)
    rays_list=[list(r) for r in normal_fan.rays()]
    ineqs=[]
    for col in normal_fan.primitive_collections():
        ineq_col = [0 if i-1 not in col else 1 for i in range(normal_fan.nrays()+1)]
        sum_rays = sum([rays[i]for i in col])
        if list(sum_rays) == [0]*len(sum_rays):
            ineqs.append(ineq_col)
        else:
            cone_containing_sum = normal_fan.cone_containing(sum_rays)
            matrix_rays = matrix([list(r) for r in cone_containing_sum.rays()])
            coeffs = matrix_rays.solve_left(matrix(sum_rays))
            for i in range(matrix_rays.nrows()):
                ineq_col[rays_list.index(list(matrix_rays[i]))+1]-=list(coeffs)[0][i]
            ineqs.append(ineq_col)
    
    equations_translation=[]
    first_cone=normal_fan.cones(codim=0)[1]
    for r in first_cone.rays():
        equation=[0]*(len(rays)+1)
        equation[rays.index(r)+1]=1
        equations_translation.append(equation)

    return (Polyhedron(ieqs=ineqs, eqns = equations_translation),rays)

def build_summand_supp(C,rays_normal_fan_p, index_ray_to_build):
    '''
    INPUT: Support functions cone, rays of the normal fan of P, number of ray to build
    OUTPUT: The polyhedron built with support function of the corresponding ray of the cone
    '''
    ray_to_build = C.rays()[index_ray_to_build]
    ineq_polytope = []
    for i in range(len(rays_normal_fan_p)):
        ineq_polytope.append([ray_to_build[i]]+list(-rays_normal_fan_p[i]))
    return Polyhedron(ieqs=ineq_polytope)