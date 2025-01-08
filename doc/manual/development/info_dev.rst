Information for developers
==========================

Building the Manual with Sphinx
-------------------------------

To build this manual using Sphinx, follow these steps:

1. **Requirements**: Install `Sphinx <https://www.sphinx-doc.org/>`_ together with the extensions used by Codac, if they are not already installed.

   .. code-block:: bash

      pip install sphinx furo sphinx-math-dollar

2. **Building the manual**: After configuring the CMake project, the manual can be build with:

   .. code-block:: bash

      make manual

   The generated website will be locally available in ``./build/doc/manual``.

To contribute and extend this manual, please consult the Sphinx documentation:  
https://www.sphinx-doc.org