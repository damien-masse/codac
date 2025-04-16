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
#include <codac2_Edge.h>
#include "codac2_py_Edge_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py):

using namespace std;
using namespace codac2;
namespace py = pybind11;
using namespace pybind11::literals;

void export_Edge(py::module& m)
{
  py::class_<Edge> exported(m, "Edge", EDGE_MAIN);
  exported

    .def(py::init<const std::array<IntervalVector,2>&>(),
      EDGE_EDGE_CONST_ARRAY_INTERVALVECTOR2_REF,
      "x"_a)

    .def(py::init<const IntervalVector&,const IntervalVector&>(),
      EDGE_EDGE_CONST_INTERVALVECTOR_REF_CONST_INTERVALVECTOR_REF,
      "x1"_a, "x2"_a)

    .def("box", &Edge::box,
      INTERVALVECTOR_EDGE_BOX_CONST)

    .def("intersects", &Edge::intersects,
      BOOLINTERVAL_EDGE_INTERSECTS_CONST_EDGE_REF_CONST,
      "e"_a)

    .def("contains", &Edge::contains,
      BOOLINTERVAL_EDGE_CONTAINS_CONST_INTERVALVECTOR_REF_CONST,
      "p"_a)

    .def(py::self == py::self,
      BOOL_EDGE_OPERATOREQ_CONST_EDGE_REF_CONST,
      "e"_a)

    .def("__and__",
        [](const Edge& e1, const Edge& e2) -> IntervalVector
        {
          return e1 & e2;
        },
      INTERVALVECTOR_OPERATORAND_CONST_EDGE_REF_CONST_EDGE_REF,
      "e2"_a);

  ;

  m.def("proj_intersection", &proj_intersection,
    INTERVALVECTOR_PROJ_INTERSECTION_CONST_EDGE_REF_CONST_EDGE_REF,
    "e1"_a, "e2"_a);
  
  py::implicitly_convertible<py::list,Edge>();
}