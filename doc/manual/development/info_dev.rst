.. _sec-dev-info:

Information for developers
==========================

Building the Manual with Sphinx
-------------------------------

To build this manual using Sphinx, follow these steps:

1. **Requirements**: Install `Sphinx <https://www.sphinx-doc.org/>`_ together with the extensions used by Codac, if they are not already installed.

   .. code-block:: bash

      pip install sphinx sphinx_rtd_theme furo sphinx-math-dollar sphinx_tabs

2. **Building the manual**: After configuring the CMake project, the manual can be build with:

   .. code-block:: bash

      make manual

   The generated website will be locally available in ``./build/doc/manual``.

To contribute and extend this manual, please consult the Sphinx documentation:  
https://www.sphinx-doc.org

.. _sec-dev-info-binding:

Building a local Python binding for Codac
-----------------------------------------

You can compile Codac's Python binding binaries locally on your machine to take advantage of the library's latest features.
If you simply want to use the latest Codac release in Python, you can download the binaries directly from `PyPi <https://pypi.org/project/codac/>`_ or :ref:`follow the standard installation procedure <sec-install-py>`.

1. **Ensure the following prerequisites are met**:

   - the prerequisites for the :ref:`C++ installation of Codac <sec-install-cpp-prerequisites>`.
   - a recent `Doxygen <https://www.doxygen.nl>`_ version. On Linux systems, Debian packages are available:

   .. code-block:: bash

      sudo apt-get install -y doxygen

   - a supported version of Python (>=3.6).

2. **Configure IBEX prior to compiling Codac**:
   
   We recall that IBEX sources can be obtained with:

   .. code-block:: bash
      
      git clone -b master https://github.com/lebarsfa/ibex-lib.git $HOME/ibex-lib

   You will need to compile both IBEX and Codac using the ``-fPIC`` options. This can be done with the following CMake configuration:

   .. code-block:: bash
      
      cd $HOME/ibex-lib/build
      cmake -DCMAKE_CXX_FLAGS="-fPIC" -DCMAKE_C_FLAGS="-fPIC" -DCMAKE_INSTALL_PREFIX=$HOME/ibex-lib/build_install -DCMAKE_BUILD_TYPE=Release ..
      make
      make install

3. **Compile Codac with Python binding**:

   We recall that Codac sources can be obtained with:
   
   .. code-block:: bash
      
      git clone https://github.com/codac-team/codac $HOME/codac

   In addition to the ``-fPIC`` options, you will have to configure ``WITH_PYTHON=ON``. Note that the ``git submodule`` commands will automatically get the `pybind11 <https://pybind11.readthedocs.io>`_ files required for the binding.
   
   .. code-block:: bash
      
      cd $HOME/codac/build
      # Get automatically pybind11 and eigen submodules:
      git submodule init ; git submodule update
      # Configure CMake
      cmake -DCMAKE_CXX_FLAGS="-fPIC" -DCMAKE_C_FLAGS="-fPIC" -DWITH_PYTHON=ON -DCMAKE_INSTALL_PREFIX=$HOME/codac/build_install -DCMAKE_PREFIX_PATH=$HOME/ibex-lib/build_install -DCMAKE_BUILD_TYPE=Release ..
      make

4. **Configure your Python environment**:

   Finally, you need to configure your system so that Python can find access to your Codac binding binaries:

   .. code-block:: bash
      
      cd $HOME/codac/build/python/python_package
      python setup.py develop --user

5. **Verify the installation** (optional):

   To ensure that the installation has worked properly, the unit tests of the library can be run:

   .. code-block:: bash
      
      python -m unittest discover codac.tests

6. **Try an example** (optional):

   You may want to try Codac in Python by running one of the proposed examples. After the installation, you can run the following commands:

   .. code-block:: bash

      cd $HOME/codac/examples/03_sivia
      python main.py

   Note that before executing the example, you will have to launch the VIBes viewer.
   You should obtain a graphical output corresponding to a set inversion.


.. admonition:: About the use of Doxygen
   
   Doxygen software extracts C++ documentation from header files into XML format. We then convert this data into docstring format before embedding it into the binding binaries. In this way, the writing of the documentation is centralized in a single location in the C++ header files.


.. admonition:: For admins
   
   The command for uploading the generated wheels on PyPi is:

   .. code-block:: bash
      
      python -m twine upload --repository pypi *