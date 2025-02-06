.. _sec-graphics-2d-figures:

The 2D Figure classes
=====================

This page describes the classes used in Codac for 2D visualization.

.. _subsec-graphics-2d-figures-graphical-outputs:

Graphical outputs
-----------------

Two graphical outputs are currently supported in Codac: :ref:`VIBes <sec-graphics-vibes>` and :ref:`IPE <sec-graphics-ipe>`. VIBes is used for real-time 
visualization while IPE creates a file that can be edited by the IPE editor. These outputs are referenced by the enumeration GraphicOutput:

.. tabs::

  .. code-tab:: py

    GraphicOutput.VIBES # for VIBes
    GraphicOutput.IPE # for IPE
    GraphicOutput.VIBES | GraphicOutput.IPE # for both

  .. code-tab:: c++

    GraphicOutput::VIBES  // for VIBes
    GraphicOutput::IPE  // for IPE
    GraphicOutput::VIBES | GraphicOutput::IPE // for both

Note that for the VIBes output to work, the VIBes viewer must be launched before the program is run.

.. _subsec-graphics-2d-figures-figure2d:

Figure2D
--------

The basic class for 2D visualization is Figure2D. It is used to create a figure that can be displayed in VIBes or saved in an xml file for IPE.
The constructor takes two arguments: the name of the figure and the graphical output. A boolean can be added to specify if the figure is to be used
DefaultView (see :ref:`subsec-graphics-2d-figures-defaultview`).

.. tabs::

  .. code-tab:: py

    fig = Figure2D("My figure", GraphicOutput.VIBES | GraphicOutput.IPE)

  .. code-tab:: c++

    Figure2D fig ("My Figure",GraphicOutput::VIBES|GraphicOutput::IPE);

.. _subsec-graphics-2d-figures-defaultview:

DefaultView
-----------

A DefaultView using only VIBes as graphical output is available. This figure is the one used by the function `draw_while_paving` by default. 
Any Figure2D object can be used as DefaultView with the set method:

.. tabs::

  .. code-tab:: py

    fig = Figure2D("My figure", GraphicOutput.VIBES | GraphicOutput.IPE)
    fig.is_default() # is False
    DefaultView.set(fig)
    fig.is_default() # is True

  .. code-tab:: c++

    std::shared_ptr<codac2::Figure2D> fig = std::make_shared<Figure2D>("My Figure",GraphicOutput::VIBES|GraphicOutput::IPE);
    fig->is_default() // is false
    DefaultView::set(fig);
    fig->is_default() // is true

Note that in C++ the figure must be a shared pointer in order to be passed to the `set` method.

Equivalently, a Figure2D can be used as DefaultView by setting the flag `set_as_default` to true in the constructor:

.. tabs::

  .. code-tab:: py

    fig = Figure2D("My figure", GraphicOutput.VIBES | GraphicOutput.IPE, True)

  .. code-tab:: c++

    Figure2D fig ("My Figure",GraphicOutput::VIBES|GraphicOutput::IPE,true);


.. _subsec-graphics-2d-figures-figure-properties:

Figure properties
-----------------

Once created, the properties of a Figure2D object can be modified using the following methods:

.. tabs::

  .. code-tab:: py

    fig.set_window_properties([50,50],[500,500]) # set the window position and size
    fig.set_axes(axis(0,[-10,10]), axis(1,[-10,10])) # set the x-axis index to 0 and its range to [-10,10], same for y with index 1 

  .. code-tab:: c++

    fig.set_window_properties({50,50},{500,500}); // set the window position and size
    fig.set_axes(axis(0,{-10,10}), axis(1,{-10,10})); // set the x-axis index to 0 and its range to [-10,10], same for y with index 1 

The same methods can be applied on the DefaultView object.

.. tabs::

  .. code-tab:: py

    DefaultView.set_window_properties([50,50],[500,500]) # set the window position and size
    DefaultView.set_axes(axis(0,[-10,10]), axis(1,[-10,10])) # set the x-axis index to 0 and its range to [-10,10], same for y with index 1 

  .. code-tab:: c++

    DefaultView::set_window_properties({50,50},{500,500}); // set the window position and size
    DefaultView::set_axes(axis(0,{-10,10}), axis(1,{-10,10})); // set the x-axis index to 0 and its range to [-10,10], same for y with index 1 

Many properties have an associated getter :

- size() for the figure size (type is Index)
- axes() to get a vector of the figure axes (type is FigureAxis)
- i() and j() to get respectively the figure x-axis and y-axis index (type is Index)
- pos() and window_size() to get the window position and size (type is Vector)
- scaled_unit() to get the scaling ratio from window to figure dimensions (type is double)