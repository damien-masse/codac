.. _sec-graphics-functions:

2D drawing functions
====================

This page lists all the 2D drawing functions available in Codac. Note that the functions described here on a
:ref:`Figure2D <subsec-graphics-figures-figure2d>` object also work on the :ref:`DefaultView <subsec-graphics-figures-defaultview>`.

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

.. _subsec-graphics-functions-draw-point:

draw_point
----------

Draw a point on the figure. This function takes **1** argument:
  - A Vector of size 2 : the center of the point.

.. tabs::

  .. code-tab:: py
    
    fig.draw_point([2,2]) # Implicit conversion from list to Vector

  .. code-tab:: c++

    fig.draw_point({2,2}); // Implicit conversion from std::initializer_list to Vector