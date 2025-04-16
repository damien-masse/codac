.. _sec-geom-utils:

Geometric utilities
===================

  Main author: `Simon Rohou <https://www.simon-rohou.fr/research/>`_

Alignments
----------

The following functions provide reliable outputs for alignement evaluations.

.. doxygenfunction:: codac2::orientation(const IntervalVector&,const IntervalVector&,const IntervalVector&)
  :project: codac

.. doxygenfunction:: codac2::aligned(const IntervalVector&,const IntervalVector&,const IntervalVector&)
  :project: codac

Convex hull
-----------

.. doxygenfunction:: codac2::convex_hull(std::vector<IntervalVector>)
  :project: codac