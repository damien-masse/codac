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
#include <codac2_geometry.h>
#include "codac2_py_geometry_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py):

using namespace std;
using namespace codac2;
namespace py = pybind11;
using namespace pybind11::literals;

void export_geometry(py::module& m)
{
  py::enum_<OrientationInterval>(m, "OrientationInterval")

    .value("EMPTY", OrientationInterval::EMPTY)
    .value("COLINEAR", OrientationInterval::COLINEAR)
    .value("CLOCKWISE", OrientationInterval::CLOCKWISE)
    .value("COUNTERCLOCKWISE", OrientationInterval::COUNTERCLOCKWISE)
    .value("UNKNOWN", OrientationInterval::UNKNOWN)

    .def(py::self | py::self,
      ORIENTATIONINTERVAL_OPERATOROR_ORIENTATIONINTERVAL_ORIENTATIONINTERVAL)

    .def("__repr__", [](const OrientationInterval& x) {
          std::ostringstream s;
          s << x;
          return string(s.str()); 
        },
      OSTREAM_REF_OPERATOROUT_OSTREAM_REF_CONST_ORIENTATIONINTERVAL_REF)
  ;

  m.def("orientation", &orientation,
    ORIENTATIONINTERVAL_ORIENTATION_CONST_INTERVALVECTOR_REF_CONST_INTERVALVECTOR_REF_CONST_INTERVALVECTOR_REF,
    "p1"_a, "p2"_a, "p3"_a);

  m.def("aligned", &aligned,
    BOOLINTERVAL_ALIGNED_CONST_INTERVALVECTOR_REF_CONST_INTERVALVECTOR_REF_CONST_INTERVALVECTOR_REF,
    "p1"_a, "p2"_a, "p3"_a);

  m.def("convex_hull", &convex_hull,
    VECTOR_INTERVALVECTOR_CONVEX_HULL_VECTOR_INTERVALVECTOR,
    "pts"_a);
}