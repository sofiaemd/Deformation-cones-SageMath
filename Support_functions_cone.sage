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

from random import randint

def rays_are_wall(rays,wall_set):
    """
    INPUT: rays of a cone and set of walls in a fan
    OUTPUT: True if the set of rays corresponds to a wall in the fan, False otherwise.
    """
    for wall in wall_set:
        if set(rays).issubset(wall):
            return True
    return False


def pos_is_primitive(wall_rays,new_rays, coefficients,fan_rays, primitive):
    """
    INPUT: rays of a wall, rays of the adjacent cones, coefficients in the corresponding wall-crossing inequality, primitive collections of the fan.
    OUTPUT: True if the set of rays with positive coefficients correspond to a primitive collection, False otherwise.
    """
    coefs_pos =[fan_rays.index(r) for r in new_rays]
    for i in range(len(coefficients)-2):
        if coefficients[i]>0:
            coefs_pos.append(fan_rays.index(wall_rays[i]))
    if frozenset(coefs_pos) in primitive:
        return True
    return False
    

def Wall_crossings_cone(P, batyrev=True):
    """
    INPUT: Polytope
    OUTPUT: Tuple containing the cone parametrizing support function of weak summands and the rays of the normal fan of the polytope.
    """
    normal_fan_p = NormalFan(-P)
    wall_set_p=set([frozenset(w.rays()) for w in normal_fan_p.cones(codim=1)])
    simplicial_fan = normal_fan_p.subdivide(make_simplicial=True)
    if P.is_simple() and batyrev:
        primitive_collections = set(normal_fan_p.primitive_collections())
    
    fan_lattice = simplicial_fan.cone_lattice()
    fan_rays = normal_fan_p.rays()
    inequalities_wall = []
    equations_wall =[]
    for wall in simplicial_fan.cones(codim=1):
        wall_rays = wall.rays()
        new_rays = []
        for c in fan_lattice.upper_covers(wall):
            new_rays += set(c.rays())-set(wall_rays)
        coefficients = matrix([list(r) for r in wall_rays]+[list(r) for r in new_rays]).kernel().basis()[0]
        if coefficients[-1]<0:
            coefficients = -1*coefficients
        
        wall_crossing_ineq=[0]*(len(fan_rays)+1)
        for i in range(len(wall_rays)):
            wall_crossing_ineq[fan_rays.index(wall_rays[i])+1]=coefficients[i]
        for i in (-1,-2):
            wall_crossing_ineq[fan_rays.index(new_rays[i])+1]=coefficients[i]

        if P.is_simple():
            if not batyrev:
                inequalities_wall.append(wall_crossing_ineq)
            elif pos_is_primitive(wall_rays,new_rays,coefficients,fan_rays,primitive_collections):
                inequalities_wall.append(wall_crossing_ineq)
        elif rays_are_wall(wall.rays(), wall_set_p):
            inequalities_wall.append(wall_crossing_ineq)
        else:
            equations_wall.append(wall_crossing_ineq)

    equations_translation=[]
    first_cone=normal_fan_p.cones(codim=0)[1]
    for r in first_cone.rays():
        equation=[0]*(len(fan_rays)+1)
        equation[fan_rays.index(r)+1]=1
        equations_translation.append(equation)
    
    return (Polyhedron(ieqs=inequalities_wall, eqns=equations_translation+equations_wall),fan_rays)

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

def Support_function_cone(P,wall_crossings=True, batyrev=True):
    if type(wall_crossings)!= bool:
        raise TypeError('wall_crossings must be a Boolean')
    if type(batyrev)!= bool:
        raise TypeError('batyrev must be a Boolean')
    if wall_crossings:
        return Wall_crossings_cone(P,batyrev=batyrev)
    elif batyrev:
        if not P.is_simple():
            raise TypeError("To create the cone only with Batyrev's inequalities, the polytope must be simple")
        return Batyrev_cone(P)
    else:
        raise TypeError('wall_crossings and batyrev cannot be False at the same time')

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