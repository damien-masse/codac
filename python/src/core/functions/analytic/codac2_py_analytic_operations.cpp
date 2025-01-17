/** 
 *  Codac binding (core)
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <sstream>
#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <codac2_analytic_variables.h>
#include <codac2_analytic_operations.h>
#include <codac2_AnalyticFunction.h>
#include "codac2_py_analytic_operations_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py)
#include "codac2_py_AnalyticExprWrapper.h"

using namespace codac2;
namespace py = pybind11;
using namespace pybind11::literals;


void export_expression_operations(py::module& m)
{
  m 

  .def("pow", [](const ScalarExpr& e1, const ScalarExpr& e2) { return pow(e1,e2); },
    SCALAREXPR_POW_CONST_SCALAREXPR_REF_CONST_SCALAREXPR_REF,
    "x1"_a, "x2"_a)

  .def("sqr", [](const ScalarExpr& e1) { return sqr(e1); },
    SCALAREXPR_SQR_CONST_SCALAREXPR_REF,
    "x1"_a)

  .def("sqrt", [](const ScalarExpr& e1) { return sqrt(e1); },
    SCALAREXPR_SQRT_CONST_SCALAREXPR_REF,
    "x1"_a)

  .def("exp", [](const ScalarExpr& e1) { return exp(e1); },
    SCALAREXPR_EXP_CONST_SCALAREXPR_REF,
    "x1"_a)

  .def("log", [](const ScalarExpr& e1) { return log(e1); },
    SCALAREXPR_LOG_CONST_SCALAREXPR_REF,
    "x1"_a)

  .def("cos", [](const ScalarExpr& e1) { return cos(e1); },
    SCALAREXPR_COS_CONST_SCALAREXPR_REF,
    "x1"_a)

  .def("sin", [](const ScalarExpr& e1) { return sin(e1); },
    SCALAREXPR_SIN_CONST_SCALAREXPR_REF,
    "x1"_a)

  .def("tan", [](const ScalarExpr& e1) { return tan(e1); },
    SCALAREXPR_TAN_CONST_SCALAREXPR_REF,
    "x1"_a)

  .def("acos", [](const ScalarExpr& e1) { return acos(e1); },
    SCALAREXPR_ACOS_CONST_SCALAREXPR_REF,
    "x1"_a)

  .def("asin", [](const ScalarExpr& e1) { return asin(e1); },
    SCALAREXPR_ASIN_CONST_SCALAREXPR_REF,
    "x1"_a)

  .def("atan", [](const ScalarExpr& e1) { return atan(e1); },
    SCALAREXPR_ATAN_CONST_SCALAREXPR_REF,
    "x1"_a)

  .def("atan2", [](const ScalarExpr& e1, const ScalarExpr& e2) { return atan2(e1,e2); },
    SCALAREXPR_ATAN2_CONST_SCALAREXPR_REF_CONST_SCALAREXPR_REF,
    "y"_a, "x"_a)

  .def("cosh", [](const ScalarExpr& e1) { return cosh(e1); },
    SCALAREXPR_COSH_CONST_SCALAREXPR_REF,
    "x1"_a)

  .def("sinh", [](const ScalarExpr& e1) { return sinh(e1); },
    SCALAREXPR_SINH_CONST_SCALAREXPR_REF,
    "x1"_a)

  .def("tanh", [](const ScalarExpr& e1) { return tanh(e1); },
    SCALAREXPR_TANH_CONST_SCALAREXPR_REF,
    "x1"_a)

  .def("abs", [](const ScalarExpr& e1) { return abs(e1); },
    SCALAREXPR_ABS_CONST_SCALAREXPR_REF,
    "x1"_a)

  ;

  const char* EVAL_VEC_DOC = VECTOREXPR_VEC_CONST_ANALYTICEXPRWRAPPER_X_REF_VARIADIC;
  using S = const ScalarExpr&;

  m .def("vec", [](S e1) { return vec(e1); }, EVAL_VEC_DOC, "x1"_a)
    .def("vec", [](S e1, S e2) { return vec(e1,e2); }, EVAL_VEC_DOC, "x1"_a,"x2"_a)
    .def("vec", [](S e1, S e2, S e3) { return vec(e1,e2,e3); }, EVAL_VEC_DOC, "x1"_a,"x2"_a,"x3"_a)
    .def("vec", [](S e1, S e2, S e3, S e4) { return vec(e1,e2,e3,e4); }, EVAL_VEC_DOC, "x1"_a,"x2"_a,"x3"_a,"x4"_a)
    .def("vec", [](S e1, S e2, S e3, S e4, S e5) { return vec(e1,e2,e3,e4,e5); }, EVAL_VEC_DOC, "x1"_a,"x2"_a,"x3"_a,"x4"_a,"x5"_a)
    .def("vec", [](S e1, S e2, S e3, S e4, S e5, S e6) { return vec(e1,e2,e3,e4,e5,e6); }, EVAL_VEC_DOC, "x1"_a,"x2"_a,"x3"_a,"x4"_a,"x5"_a,"x6"_a)
    .def("vec", [](S e1, S e2, S e3, S e4, S e5, S e6, S e7) { return vec(e1,e2,e3,e4,e5,e6,e7); }, EVAL_VEC_DOC, "x1"_a,"x2"_a,"x3"_a,"x4"_a,"x5"_a,"x6"_a,"x7"_a)
    .def("vec", [](S e1, S e2, S e3, S e4, S e5, S e6, S e7, S e8) { return vec(e1,e2,e3,e4,e5,e6,e7,e8); }, EVAL_VEC_DOC, "x1"_a,"x2"_a,"x3"_a,"x4"_a,"x5"_a,"x6"_a,"x7"_a,"x8"_a)
    .def("vec", [](S e1, S e2, S e3, S e4, S e5, S e6, S e7, S e8, S e9) { return vec(e1,e2,e3,e4,e5,e6,e7,e8,e9); }, EVAL_VEC_DOC, "x1"_a,"x2"_a,"x3"_a,"x4"_a,"x5"_a,"x6"_a,"x7"_a,"x8"_a,"x9"_a)
    .def("vec", [](S e1, S e2, S e3, S e4, S e5, S e6, S e7, S e8, S e9, S e10) { return vec(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10); }, EVAL_VEC_DOC, "x1"_a,"x2"_a,"x3"_a,"x4"_a,"x5"_a,"x6"_a,"x7"_a,"x8"_a,"x9"_a,"x10"_a)
  ;
}