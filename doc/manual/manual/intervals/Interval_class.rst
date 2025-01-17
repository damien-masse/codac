.. _sec-intervals-class:

The Interval class
==================

``Interval x(double lb, double ub)`` will define the interval :math:`[x]`.
It is made of its lower and upper bounds :math:`[x^{-},x^{+}]`.

.. tabs::
  
  .. code-tab:: py

    x = Interval(0,10)                          # [0,10]
    x = Interval(1,oo)                          # [1,∞]
    x = Interval(-oo,-10)                       # [-∞,-10]

  .. code-tab:: c++

    Interval x(0,10);                           // [0,10]
    Interval x(1,oo);                           // [1,∞]
    Interval x(-oo,-10);                        // [-∞,-10]

Some pre-defined values are also at hand:

.. tabs::
  
  .. code-tab:: py

    x = Interval()                              # [-∞,∞] (default value)
    x = Interval.empty()                        # ∅
    x = Interval.pi()                           # [π]
    x = Interval.two_pi()                       # [2π]
    x = Interval.half_pi()                      # [π/2]

  .. code-tab:: c++

    Interval x;                                 // [-∞,∞] (default value)
    Interval x = Interval::empty();             // ∅
    Interval x = Interval::pi();                // [π]
    Interval x = Interval::two_pi();            // [2π]
    Interval x = Interval::half_pi();           // [π/2]

.. _sec-intervals-class-pi:

Note that the constant :math:`[\pi]` is a reliable enclosure of the :math:`\pi` value, that cannot be exactly represented in a computer with a single floating-point value.

.. tabs::
  
  .. code-tab:: py

    x = Interval.pi()                           # [π]
    # x = [3.141592653589793, 3.141592653589794]

  .. code-tab:: c++

    Interval x = Interval::pi();                // [π]
    // x = [3.141592653589793, 3.141592653589794]