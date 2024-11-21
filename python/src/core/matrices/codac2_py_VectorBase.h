/** 
 *  VectorBase binding
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

#include "codac2_py_MatrixBase.h"

using namespace std;
using namespace codac2;
namespace py = pybind11;
using namespace pybind11::literals;

template<typename S,typename M,typename T>
void export_VectorBase(py::module& m, py::class_<S>& pyclass)
{
  //export_MatrixBase<S,T,true>(m, pyclass);
  // ^ We do not "inherit" from export_MatrixBase here, in order to
  // avoid double inheritance from IntervalVector or IntervalMatrix.
  // Inheritance is compensated in Vector binding.

  pyclass

    .def("__getitem__", [](const S& x, size_t_type index) -> const T&
        {
          matlab::test_integer(index);
          return x[matlab::input_index(index)];
        }, py::return_value_policy::reference_internal,
      VECTORBASE_EIGENADDONS_CONST_SCALAR_REF_OPERATORCOMPO_SIZET_CONST)

    .def("__setitem__", [](S& x, size_t_type index, const T& a)
        {
          matlab::test_integer(index);
          x[matlab::input_index(index)] = a;
        },
      VECTORBASE_EIGENADDONS_SCALAR_REF_OPERATORCOMPO_SIZET)

    .def("subvector", [](const S& x, size_t_type start_id, size_t_type end_id) -> S
        {
          matlab::test_integer(start_id, end_id);
          return x.subvector(matlab::input_index(start_id), matlab::input_index(end_id));
        },
      VECTORBASE_EIGENADDONS_AUTO_SUBVECTOR_SIZET_SIZET_CONST,
      "start_id"_a, "end_id"_a)

    .def("resize", [](S& x, size_t_type n)
        {
          matlab::test_integer(n);
          x.resize(n);
        },
      DOC_TO_BE_DEFINED,
      "n"_a)

    .def("resize_save_values", [](S& x, size_t_type n)
        {
          matlab::test_integer(n);
          x.resize_save_values(n);
        },
      VECTORBASE_EIGENADDONS_VOID_RESIZE_SAVE_VALUES_SIZET,
      "n"_a)

    .def("put", [](S& x, size_t_type start_id, const S& x1)
        {
          matlab::test_integer(start_id);
          x.put(matlab::input_index(start_id), x1);
        },
      VECTORBASE_EIGENADDONS_VOID_PUT_SIZET_CONST_MATRIX_SCALARRC_REF,
      "start_id"_a, "x"_a)
    
    .def_static("zeros", [](size_t_type n)
        {
          matlab::test_integer(n);
          return S::zeros(n);
        },
      VECTORBASE_EIGENADDONS_STATIC_MATRIX_SCALARRC_ZEROS_SIZET,
      "n"_a)
    
    .def_static("ones", [](size_t_type n)
        {
          matlab::test_integer(n);
          return S::ones(n);
        },
      VECTORBASE_EIGENADDONS_STATIC_MATRIX_SCALARRC_ONES_SIZET,
      "n"_a)

    .def("__repr__", [](const S& x)
        {
          std::ostringstream s;
          s << x;
          return string(s.str()); 
        },
      DOC_TO_BE_DEFINED)

    .def("__iter__", [](const S &x)
        {
          return py::make_iterator(x.begin(), x.end());
        },
      py::keep_alive<0, 1>() /*  keep object alive while iterator exists */)
  ;
}