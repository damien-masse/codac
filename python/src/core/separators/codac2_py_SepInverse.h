/** 
 *  Codac binding (core)
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <codac2_SepInverse.h>
#include "codac2_Sep.h"
#include "codac2_py_Sep.h"
#include "codac2_py_SepInverse_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py):
#include "codac2_py_SepCtcPair_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py):

using namespace std;
using namespace codac2;
namespace py = pybind11;
using namespace pybind11::literals;

template<typename T>
void export_SepInverse(py::module& m, const std::string& export_name, py::class_<SepBase,pySep>& pysep)
{
  using D = typename T::Domain;
  py::class_<SepInverse<D>> exported(m, export_name.c_str(), pysep, SEPINVERSE_MAIN);

  exported
    .def(py::init<const AnalyticFunction<T>&,const D&,bool>(),
      SEPINVERSE_Y_SEPINVERSE_CONST_ANALYTICFUNCTION_TYPENAME_VALUETYPE_Y_TYPE_REF_CONST_Y_REF_BOOL,
      "f"_a, "y"_a, "with_centered_form"_a = true);

  if constexpr(std::is_same_v<T,IntervalVector>) // separators only associated with interval vectors
  {
    exported
    .def(py::init<const AnalyticFunction<T>&,const pySep&,bool>(),
      SEPINVERSE_Y_SEPINVERSE_CONST_ANALYTICFUNCTION_TYPENAME_VALUETYPE_Y_TYPE_REF_CONST_S_REF_BOOL,
      "f"_a, "s"_a, "with_centered_form"_a = true);
  }

  exported

    .def("separate", &SepInverse<D>::separate,
      BOXPAIR_SEPCTCPAIR_SEPARATE_CONST_INTERVALVECTOR_REF_CONST,
      "x"_a)
  ;
}