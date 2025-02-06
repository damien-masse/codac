.. _sec-graphics-colors:

Colors and color maps
=====================

By default, the drawn shapes will have a black edge and no fill. Colors can be passed as an additional
argument to change it. For example:

.. tabs::

  .. code-tab:: py

    fig.draw_box([[2.2,2.5],[2.2,2.5]]) # Default style
    fig.draw_box([[2.2,2.5],[2.2,2.5]],Color.red()) # Red edge

  .. code-tab:: c++

    fig.draw_box({{2.2,2.5},{2.2,2.5}}); // Default style
    fig.draw_box({{2.2,2.5},{2.2,2.5}},Color::red()); // Red edge