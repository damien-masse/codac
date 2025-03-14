.. _sec-graphics-functions:

2D drawing functions
====================

This page lists all the 2D drawing functions available in Codac. Note that the functions described here on a
:ref:`Figure2D <subsec-graphics-2d-figures-figure2d>` object also work on the :ref:`DefaultView <subsec-graphics-2d-figures-defaultview>`.

For example : 

.. tabs::

  .. code-tab:: py
    
    # On a figure
    fig.draw_box([[2.2,2.5],[2.2,2.5]])

    # On the DefaultView
    DefaultView.draw_box([[2.2,2.5],[2.2,2.5]])

  .. code-tab:: c++

    // On a figure
    fig.draw_box({{2.2,2.5},{2.2,2.5}});

    // On the DefaultView
    DefaultView::draw_box({{2.2,2.5},{2.2,2.5}});

In addition to the arguments described here for every function, an optionnal argument of type StyleProperties can be added to choose the edge
and fill color (by default black edge no fill). This object can be deduced from one or two Color objects, see :ref:`this page <sec-graphics-colors>` 
for more information.

.. tabs::

  .. code-tab:: py
    
    fig.draw_box([[2.2,2.5],[2.2,2.5]],Color.red()) # Red edge
    fig.draw_box([[2.2,2.5],[2.2,2.5]],[Color.red()]) # Red edge
    fig.draw_box([[2.2,2.5],[2.2,2.5]],[Color.red(),Color.blue()]) # Red edge and blue fill

  .. code-tab:: c++

    fig.draw_box({{2.2,2.5},{2.2,2.5}},Color.red()); // Red edge
    fig.draw_box({{2.2,2.5},{2.2,2.5}},{Color.red()}); // Red edge
    fig.draw_box({{2.2,2.5},{2.2,2.5}},{Color.red(),Color.blue()}); // Red edge and blue fill


draw_point
----------

Draws a point on the figure. This function takes **1** argument:
  - A Vector : the center of the point.

.. tabs::

  .. code-tab:: py
    
    fig.draw_point([2,2]) # Point with center (2,2)

  .. code-tab:: c++

    fig.draw_point({2,2}); // Point with center (2,2)

draw_box
--------

Draws a box on the figure. This function takes **1** argument:
  - An IntervalVector : the box to draw.

.. tabs::

  .. code-tab:: py
    
    fig.draw_box([[-1,1],[-1,1]]) # Box [-1,1]x[-1,1]

  .. code-tab:: c++

    fig.draw_box({{-1,1},{-1,1}}); // Box [-1,1]x[-1,1]

draw_circle
-----------

Draws a circle on the figure. This function takes **2** arguments:
  - A Vector : the center of the circle.
  - A double : the radius of the circle.

.. tabs::

  .. code-tab:: py
    
    fig.draw_circle([1,1],0.5) # Circle with center (1,1) and radius 0.5

  .. code-tab:: c++

    fig.draw_circle({1,1},0.5); // Circle with center (1,1) and radius 0.5

draw_ring
---------

Draws a ring on the figure. This function takes **2** arguments:
  - A Vector : the center of the ring.
  - An Interval : the inner and outer radius of the ring.

.. tabs::

  .. code-tab:: py
    
    fig.draw_ring([1,1],[4,6]) # Ring with center (1,1) and inner radius 4 and outer radius 6

  .. code-tab:: c++

    fig.draw_ring({1,1},{4,6}); // Ring with center (1,1) and inner radius 4 and outer radius 6

draw_line
---------

Draws a line on the figure. This function takes **2** arguments:
  - A Vector : the start of the line.
  - A Vector : the end of the line.

.. tabs::

  .. code-tab:: py
    
    fig.draw_line([1,1],[3,3]) # Line from (1,1) to (3,3)

  .. code-tab:: c++

    fig.draw_line({1,1},{3,3}); // Line from (1,1) to (3,3)

draw_arrow
----------

Draws an arrow on the figure. This function takes **3** arguments:
  - A Vector : the start of the arrow.
  - A Vector : the end of the arrow.
  - A double : the size of the arrow head.

.. tabs::

  .. code-tab:: py
    
    fig.draw_arrow([3,1],[2.2,2], 0.2) # Arrow from (3,1) to (2.2,2) with arrow head size 0.2

  .. code-tab:: c++

    fig.draw_arrow({3,1},{2.2,2}, 0.2); // Arrow from (3,1) to (2.2,2) with arrow head size 0.2

draw_polyline
-------------

Draws a polyline on the figure. This function takes **1** argument:
  - A vector of Vector : the polyline to draw.

This function can take an optional argument:
  - A double : the size of the arrow head (if not specified no arrow head will be displayed).

.. tabs::

  .. code-tab:: py
    
    fig.draw_polyline([[2,-0.5],[4,0.5],[3,1.5],[4,2.5],[3,3]]) # Polyline from (2,-0.5) to (4,0.5) to (3,1.5) to (4,2.5) to (3,3)
    fig.draw_polyline([[-0.8,0],[0,1.5]], 0.2) # Polyline from (-0.8,0) to (0,1.5) with arrow head size 0.2


  .. code-tab:: c++

    fig.draw_polyline({{2,-0.5},{4,0.5},{3,1.5},{4,2.5},{3,3}}); // Polyline from (2,-0.5) to (4,0.5) to (3,1.5) to (4,2.5) to (3,3)
    fig.draw_polyline({{-0.8,0},{0,1.5}}, 0.2); // Polyline from (-0.8,0) to (0,1.5) with arrow head size 0.2

draw_polygone
-------------

Draws a polygone on the figure. This function takes **1** argument:
  - A vector of Vector : the polygone to draw.

.. tabs::

  .. code-tab:: py
    
    fig.draw_polygone([[2,4.5],[4,4.5],[4.2,3.5],[3.5,3]]) # Polygone with vertices (2,4.5), (4,4.5), (4.2,3.5) and (3.5,3)

  .. code-tab:: c++

    fig.draw_polygone({{2,4.5},{4,4.5},{4.2,3.5},{3.5,3}}); // Polygone with vertices (2,4.5), (4,4.5), (4.2,3.5) and (3.5,3)

draw_pie
--------

Draws a pie on the figure. This function takes **3** arguments:
  - A Vector : the center of the pie.
  - An Interval : the inner and outer radius of the pie.
  - An Interval : the start and end angles (in radian) of the pie.


.. tabs::

  .. code-tab:: py
    
    fig.draw_pie([2,2],[1.5,2.5],[0,PI/2]) # Pie with center (2,2), inner radius 1.5, outer radius 2.5, start angle 0 and end angle PI/2

  .. code-tab:: c++

    fig.draw_pie({2,2},{1.5,2.5},{0,PI/2}); // Pie with center (2,2), inner radius 1.5, outer radius 2.5, start angle 0 and end angle PI/2

draw_ellipse
------------

Draws an ellipse on the figure. This function takes **3** argument:
  - A Vector of size 2 : the center of the ellipse.
  - A Vector of size 2 : the lengths of the semi-axes of the ellipse.
  - A double : the rotation angle of the ellipse (in radian).

.. tabs::

  .. code-tab:: py
    
    fig.draw_ellipse([1,1],[0.5,2], 0.2)

  .. code-tab:: c++

    fig.draw_ellipse({1,1},{0.5,2}, 0.2);

draw_ellipsoid
--------------

Draws an n-dimensionnal Ellipsoid projected on the figure. This function takes **1** argument:
  - An Ellipsoid : the ellipsoid to draw.

The plane used for the projection is the one defined by the axes of the figure, see :ref:`subsec-graphics-2d-figures-figure-properties`.

.. tabs::

  .. code-tab:: py
    
    mu = Vector([1., 0.]) # Center of the ellipsoid
    G = Matrix([[0.05, 0.0],
              [0.,   0.05]])  # Covariance matrix
    e = Ellipsoid(mu, G)
    fig.draw_ellipsoid(e)

  .. code-tab:: c++

    Ellipsoid e(
        {1., 0.}, // Center of the ellipsoid
        {{0.05, 0.0}, // Covariance matrix
         {0.,   0.05}}
    );
    fig.draw_ellipsoid(e);

draw_trajectory
---------------

Draws a trajectory on the figure. This function takes **1** argument:
  - An AnalyticTraj **OR** a SampledTraj: the trajectory to draw.

Trajectories can be drawn with a ColorMap instead of the classic StyleProperties. This can be done by passing a ColorMap object as the second argument.

.. tabs::

  .. code-tab:: py
    
    a = 0.8
    t=ScalarVar()
    f=AnalyticFunction([t], [a*sqrt(t)*cos(t),a*sqrt(t)*sin(t)])
    traj=AnalyticTraj(f, [0,100])
    fig.draw_trajectory(traj, ColorMap.rainbow())

  .. code-tab:: c++

    double a=0.5;
    ScalarVar t;
    AnalyticFunction f ({t},{a*sqrt(t)*cos(t),a*sqrt(t)*sin(t)});
    AnalyticTraj traj (f,{0,100});
    fig4.draw_trajectory(traj,ColorMap::rainbow());
