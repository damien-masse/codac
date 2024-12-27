/** 
 *  Codac binding (core)
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include <sstream>
#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include "codac2_AnalyticExprWrapper.h"
#include "codac2_analytic_operations.h"
#include "codac2_py_matlab.h"
#include "codac2_py_AnalyticExprWrapper_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py)

using namespace codac2;
namespace py = pybind11;
using namespace pybind11::literals;

inline void export_ScalarExpr(py::module& m)
{
  py::class_<ScalarExpr>
    exported(m, "ScalarExpr");

  exported

    .def(py::init<int>())
    .def(py::init<double>())
    .def(py::init<Interval>())
    .def(py::init<ScalarExpr>())
    .def(py::init<ScalarVar>())

    .def("__pos__",  [](const ScalarExpr& e1) { return e1; }, py::is_operator())
    .def(py::self + py::self)
    .def("__add__",  [](const ScalarExpr& e1, const Interval& e2)   { return e1+e2; }, py::is_operator())
    .def("__radd__", [](const ScalarExpr& e1, const Interval& e2)   { return e2+e1; }, py::is_operator())

    .def(- py::self)
    .def(py::self - py::self)
    .def("__sub__",  [](const ScalarExpr& e1, const Interval& e2)   { return e1-e2; }, py::is_operator())
    .def("__rsub__", [](const ScalarExpr& e1, const Interval& e2)   { return e2-e1; }, py::is_operator())

    .def(py::self * py::self)
    .def("__mul__",  [](const ScalarExpr& e1, const Interval& e2)   { return e1*e2; }, py::is_operator())
    .def("__rmul__", [](const ScalarExpr& e1, const Interval& e2)   { return e2*e1; }, py::is_operator())
    .def("__mul__",  [](const ScalarExpr& e1, const VectorExpr& e2) { return e1*e2; }, py::is_operator())

    .def("__truediv__",  [](const ScalarExpr& e1, const ScalarExpr& e2) { return e1/e2; }, py::is_operator())
    .def("__truediv__",  [](const ScalarExpr& e1, const Interval& e2)   { return e1/e2; }, py::is_operator())
    .def("__rtruediv__", [](const ScalarExpr& e1, const Interval& e2)   { return e2/e1; }, py::is_operator())

  ;

  py::implicitly_convertible<int,ScalarExpr>();
  py::implicitly_convertible<double,ScalarExpr>();
  py::implicitly_convertible<Interval,ScalarExpr>();
  py::implicitly_convertible<ScalarVar,ScalarExpr>();
}

inline void export_VectorExpr(py::module& m)
{
  py::class_<VectorExpr>
    exported(m, "VectorExpr");

  exported

    .def(py::init<Vector>())
    .def(py::init<IntervalVector>())
    .def(py::init<VectorExpr>())
    .def(py::init<VectorVar>())

    .def("__getitem__", [](const VectorExpr& x, Index_type index)
        {
          matlab::test_integer(index);
          return x[matlab::input_index(index)];
        }, py::return_value_policy::reference_internal,
      ANALYTICEXPRWRAPPER_SCALAROPVALUE_ANALYTICEXPRWRAPPER_T_OPERATORCOMPO_INDEX_CONST)

    .def("__pos__",  [](const VectorExpr& e1)                           { return e1; }, py::is_operator())
    .def(py::self + py::self)
    .def("__add__",  [](const VectorExpr& e1, const IntervalVector& e2) { return e1+e2; }, py::is_operator())
    .def("__radd__", [](const VectorExpr& e1, const IntervalVector& e2) { return e1+e2; }, py::is_operator())

    .def(- py::self)
    .def(py::self - py::self)
    .def("__sub__",  [](const VectorExpr& e1, const IntervalVector& e2) { return e1-e2; }, py::is_operator())
    .def("__rsub__", [](const VectorExpr& e1, const IntervalVector& e2) { return e2-e1; }, py::is_operator())

    .def("__mul__",  [](const VectorExpr& e1, const ScalarExpr& e2)     { return e1*e2; }, py::is_operator())
    .def("__rmul__", [](const VectorExpr& e1, const ScalarExpr& e2)     { return e2*e1; }, py::is_operator())
    .def("__mul__",  [](const VectorExpr& e1, const Interval& e2)       { return e1*e2; }, py::is_operator())
    .def("__rmul__", [](const VectorExpr& e1, const Interval& e2)       { return e2*e1; }, py::is_operator())

    .def("__truediv__",  [](const VectorExpr& e1, const ScalarExpr& e2) { return e1/e2; }, py::is_operator())
    .def("__truediv__",  [](const VectorExpr& e1, const Interval& e2)   { return e1/e2; }, py::is_operator())

  ;

  py::implicitly_convertible<Vector,VectorExpr>();
  py::implicitly_convertible<IntervalVector,VectorExpr>();
  py::implicitly_convertible<VectorVar,VectorExpr>();
}