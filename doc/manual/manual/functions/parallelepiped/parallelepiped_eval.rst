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

:math:`\mathbb{R}^2\to\mathbb{R}^3`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Related work
------------

This method comes from `this article <https://www.sciencedirect.com/science/article/pii/S0888613X25002154?via%3Dihub>`_