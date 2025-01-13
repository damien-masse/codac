.. _sec-graphics-figures:

The Figure classes
==================

This page describes the classes used in Codac for 2D visualization.

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

Figure2D
--------




.. tabs::

  .. code-tab:: py

    from codac import *

    # .. next questions will be here

  .. code-tab:: c++

    #include <codac.h>

    using namespace std;
    using namespace codac;

    int main()
    {
      // .. next questions will be here
    }

DefaultView
-----------

A DefaultView using only VIBes as graphical output is available. This figure is for example used in `draw_while_paving` by default. 
Any Figure2D object can be used as DefaultView with the set method:

.. tabs::

  .. code-tab:: py

    figure = Figure2D("My figure", GraphicOutput.VIBES | GraphicOutput.IPE)
    DefaultView.set(figure)

  .. code-tab:: c++

    Figure2D figure = std::make_shared<Figure2D>("My Figure",GraphicOutput::VIBES|GraphicOutput::IPE);
    DefaultView::set(figure);