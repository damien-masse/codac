.. _sec-geom-index:

Geometry
========

Codac provides a set of utility functions for basic 2d geometric calculations.
Because computations are based on interval arithmetic, all these functions provide reliable results (unless the function name indicates otherwise). Therefore, the outputs of these functions are sets, or interval enumerations such as ``BoolInterval`` or ``OrientationInterval`` types.


.. toctree::

   utils.rst
   Segment.rst

..   Geometrical tools <http://codac.io>
..   Polygons and convex polygons <http://codac.io>
..   Polyhedron <http://codac.io>


Related interval enumerations
-----------------------------

We recall the ``BoolInterval`` enumeration:

.. doxygenenum:: codac2::BoolInterval
  :project: codac

In addition, the ``OrientationInterval`` enumeration is dedicated to the geometrical functions.

.. doxygenenum:: codac2::OrientationInterval
  :project: codac
