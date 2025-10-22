.. _sec-functions-parallelepiped-eval:

Parallelepiped inclusion function
=================================

  Main authors: `Maël Godard <https://godardma.github.io>`_


Use case
--------

Consider a function :math:`\mathbf{f}:\mathbb{R}^n\to\mathbb{R}^m`. In the case where :math:`0<n<m`,
a parallelepiped inclusion function is available in the library.

A parallelepiped inclusion function of :math:`\mathbf{f}` is noted :math:`\langle\mathbf{f}\rangle`.
For a :math:`n`-dimensional box :math:`[\mathbf{x}]`, :math:`\langle\mathbf{f}\rangle([\mathbf{x}])` is a parallelepiped which enclose :math:`\mathbf{f}([\mathbf{x}])`.


Examples
--------

:math:`\mathbb{R}^1\to\mathbb{R}^2`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Define the parabolic function :math:`\mathbf{f}:\mathbb{R}\to\mathbb{R}^2` such that

.. math:: 
  \mathbf{f}(x)=
  \left(
  \begin{array}{c}
  x\\
  x^{2}
  \end{array}
  \right)

It can be defined in Codac as follows:

.. tabs::

  .. group-tab:: Python

    .. literalinclude:: src.py
      :language: py
      :start-after: [1-beg]
      :end-before: [1-end]
      :dedent: 0

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [1-beg]
      :end-before: [1-end]
      :dedent: 4

  .. group-tab:: Matlab

    .. literalinclude:: src.m
      :language: matlab
      :start-after: [1-beg]
      :end-before: [1-end]
      :dedent: 0

The performances of the parallelepiped inclusion function can be visualized as follows:

.. tabs::

  .. group-tab:: Python

    .. literalinclude:: src.py
      :language: py
      :start-after: [2-beg]
      :end-before: [2-end]
      :dedent: 0

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [2-beg]
      :end-before: [2-end]
      :dedent: 4

  .. group-tab:: Matlab

    .. literalinclude:: src.m
      :language: matlab
      :start-after: [2-beg]
      :end-before: [2-end]
      :dedent: 0

The resulting plot is shown below. In blue is the interval evaluation of Codac (see :ref:`subsec-functions-analytic-class-evaluations`).
In green is the evaluation using parallelepipeds.

.. image:: img/2d_parabola.png
   :alt: 2D Parabola with parallelepiped inclusion
   :align: center
   :width: 400px

:math:`\mathbb{R}^2\to\mathbb{R}^3`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Define the parabolic function :math:`\mathbf{f}:\mathbb{R}^2\to\mathbb{R}^3` such that

.. math:: 
  \mathbf{f}(\mathbf{x})=
  \left(
  \begin{array}{c}
  \mathbf{x}_{1}\\
  \mathbf{x}_{2}\\
  \mathbf{x}_{1}^{2}+\mathbf{x}_{2}^{2}
  \end{array}
  \right)

It can be defined in Codac as follows:

.. tabs::

  .. group-tab:: Python

    .. literalinclude:: src.py
      :language: py
      :start-after: [3-beg]
      :end-before: [3-end]
      :dedent: 0

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [3-beg]
      :end-before: [3-end]
      :dedent: 4

  .. group-tab:: Matlab

    .. literalinclude:: src.m
      :language: matlab
      :start-after: [3-beg]
      :end-before: [3-end]
      :dedent: 0

The performances of the parallelepiped inclusion function can be visualized as follows:

.. tabs::

  .. group-tab:: Python

    .. literalinclude:: src.py
      :language: py
      :start-after: [4-beg]
      :end-before: [4-end]
      :dedent: 0

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [4-beg]
      :end-before: [4-end]
      :dedent: 4

  .. group-tab:: Matlab

    .. literalinclude:: src.m
      :language: matlab
      :start-after: [4-beg]
      :end-before: [4-end]
      :dedent: 0

The resulting figures are shown below. On the left is the result of the interval evaluation of Codac (layer "box") (see :ref:`subsec-functions-analytic-class-evaluations`).
On the right is the evaluation using parallelepipeds (layer "parallelepiped").

.. image:: img/3d_parabola.png
   :alt: 3D Parabola with boxes and parallelepipeds
   :align: center
   :width: 600px

Related work
------------

This method comes from `this article <https://www.sciencedirect.com/science/article/pii/S0888613X25002154?via%3Dihub>`_