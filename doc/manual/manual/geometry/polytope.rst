.. _sec-zonotope:

Polytope
========

  Main author: `Damien Massé <https://labsticc.fr/fr/annuaire/masse-damien>`_

The Polytope class represent convex polytopes (or, more generally, 
convex polyhedra).
Polytopes are internally represented as intersections of linear constraints
(and a bounding box), in a class name CollectFacets.
A double-description algorithm is used to compute generators for operations
requiring them (union, projections...). Due to the imprecise nature of
floating-point computations, the algorithm is designed such that the
convex hull of the "generators" encloses the polytope. However, no guarantee
is given that each vertice of the polytope is associated to a generator.
As a result,
going back to the hyperplane representation may return a larger polytope.

Polytope
--------

.. doxygenclass:: codac2::Polytope
  :project: codac
  :members:

CollectFacets
-------------

.. doxygenclass:: codac2::CollectFacets
  :project: codac
  :members:

