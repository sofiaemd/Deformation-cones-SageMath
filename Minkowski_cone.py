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
#     Copyright (C) 2023 Sofia Errazuruz <sofiaerrazurizm@uc.cl>
#                       Federico Castillo <efecastillo.math@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
##############################################################################

from sage.geometry.polyhedron.plot import cyclic_sort_vertices_2d
from sage.matrix.constructor import matrix


def Minkowski_cone(P):
    """
    INPUT: Polytope
    OUTPUT: Cone parametrizing weak Minkowski summands.
    We follow the 1-Minkowski weights approach.
    """

    edges_p = [{f.vertices()[0], f.vertices()[1]} for f in list(P.faces(1))]
    nEdges = len(edges_p)

    # The inequalities indictaing each edge is nonnegative
    nnegativity = matrix(nEdges, 1).augment(matrix.identity(nEdges))
    edge_balancing = []
    for face in P.faces(2):
        edges_face = face.as_polyhedron().faces(1)
        ordered_vertices = cyclic_sort_vertices_2d(face.vertices())
        edge_vectors = [
            vector(ordered_vertices[i]) - vector(ordered_vertices[i + 1])
            for i in range(-1, len(ordered_vertices) - 1)
        ]
        ordered_edges = [
            {ordered_vertices[i], ordered_vertices[i + 1]}
            for i in range(-1, len(ordered_vertices) - 1)
        ]
        indices_edges = [edges_p.index(a) for a in ordered_edges]

        for i in range(P.dimension()):
            edge_balancing.append(
                [
                    edge_vectors[indices_edges.index(j - 1)][i]
                    if j - 1 in indices_edges
                    else 0
                    for j in range(len(edges_p) + 1)
                ]
            )

    return (Polyhedron(ieqs=nnegativity, eqns=edge_balancing), edges_p, list(P.vertices()))

def build_summand_mink(C,edges_p,vertices_p, build_ray):
    '''
    INPUT: Minkowski Cone, edges of P in order, vertices of P, number of ray to build
    OUTPUT: The polyhedron built with 1-weights of the corresponding ray of the cone
    '''
    vertices_p.sort(key=lambda v: sum(vector(v)))
    C_ray=list(C.rays()[build_ray])
    maximo=max(C_ray)
    C_ray=[(1/maximo)*float(c) for c in C_ray]
    dict_edges=dict()
    dict_edges[vertices_p[0]]=vector(vertices_p[0])
    for i in range(1,len(vertices_p)):
        v=vertices_p[i]
        for w in reversed(vertices_p[:i:]):
            if {v,w} in edges_p:
                dict_edges[v]= dict_edges[w] + C_ray[edges_p.index({v,w})]*(vector(v)-vector(w))
                break
    return Polyhedron(vertices=dict_edges.values())
