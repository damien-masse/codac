/** 
 *  TrajectoryBase Python binding
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include <sstream>
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <codac2_Matrix.h>
#include "codac2_py_matlab.h"

using namespace std;
using namespace codac2;
namespace py = pybind11;
using namespace pybind11::literals;

template<typename S>
void export_TrajectoryBase(py::class_<S>& pyclass)
{
  pyclass

    .def("size", [](const S& x)
        {
          return x.size();
        },
      VIRTUAL_INDEX_TRAJECTORYBASE_T_SIZE_CONST)

    .def("is_empty", [](const S& x)
        {
          return x.is_empty();
        },
      VIRTUAL_BOOL_TRAJECTORYBASE_T_IS_EMPTY_CONST)

    .def("tdomain", [](const S& x)
        {
          return x.tdomain();
        },
      VIRTUAL_INTERVAL_TRAJECTORYBASE_T_TDOMAIN_CONST)

    .def("truncate_tdomain", [](S& x, const Interval& new_tdomain)
        {
          x.truncate_tdomain(new_tdomain);
        },
      VIRTUAL_VOID_TRAJECTORYBASE_T_TRUNCATE_TDOMAIN_CONST_INTERVAL_REF,
      "new_tdomain"_a)

    .def("codomain", [](const S& x)
        {
          return x.codomain();
        },
      VIRTUAL_WRAPPER_T_DOMAIN_TRAJECTORYBASE_T_CODOMAIN_CONST)

    .def("__call__", [](const S& x, double t)
        {
          return x(t);
        },
      VIRTUAL_T_TRAJECTORYBASE_T_OPERATORCALL_DOUBLE_CONST,
      "t"_a)

    .def("__call__", [](const S& x, const Interval& t)
        {
          return x(t);
        },
      VIRTUAL_WRAPPER_T_DOMAIN_TRAJECTORYBASE_T_OPERATORCALL_CONST_INTERVAL_REF_CONST,
      "t"_a)

    .def("nan_value", [](const S& x)
        {
          return x.nan_value();
        },
      AUTO_TRAJECTORYBASE_T_NAN_VALUE_CONST)

    .def("sampled", [](const S& x, double dt)
        {
          return x.sampled(dt);
        },
      VIRTUAL_SAMPLEDTRAJECTORY_T_TRAJECTORYBASE_T_SAMPLED_DOUBLE_CONST,
      "dt"_a)

    .def("primitive", [](const S& x, const typename S::Type::Scalar& y0, double dt)
        {
          return x.primitive(y0,dt);
        },
      SAMPLEDTRAJECTORY_T_TRAJECTORYBASE_T_PRIMITIVE_CONST_T_REF_DOUBLE_CONST,
      "y0"_a, "dt"_a)

  ;
}