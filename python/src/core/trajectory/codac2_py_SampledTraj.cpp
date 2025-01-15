/** 
 *  SampledTraj binding
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <sstream>
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <codac2_SampledTraj.h>
#include <codac2_SampledTraj_operations.h>
#include "codac2_py_SampledTraj_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py):
#include "codac2_py_TrajBase_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py):
#include "codac2_py_SampledTraj_operations_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py):
#include "codac2_py_TrajBase.h"

using namespace std;
using namespace codac2;
namespace py = pybind11;
using namespace pybind11::literals;

template<typename T>
void _export_SampledTraj(py::module& m, const string& class_name)
{
  py::class_<SampledTraj<T>> exported_class(m, class_name.c_str(), SAMPLEDTRAJ_MAIN);
  export_TrajBase<SampledTraj<T>>(exported_class);

  exported_class

    .def(py::init<>(),
      SAMPLEDTRAJ_T_SAMPLEDTRAJ)
    
  ;

  if constexpr(std::is_same_v<T,double>)
  {
    exported_class
      .def(py::init<const std::list<double>&,const std::list<T>&>(),
        SAMPLEDTRAJ_T_SAMPLEDTRAJ_CONST_LIST_DOUBLE_REF_CONST_LIST_T_REF,
        "l_t"_a, "l_x"_a)
    ;
  }

  else if constexpr(std::is_same_v<T,Vector>)
  {
    exported_class
      .def(py::init(
          [](const py::array_t<double>& l_t, const py::array_t<double>& l_x)
            {
              if constexpr(std::is_same_v<T,Vector>)
              {
                assert_release(l_t.size() > 0 && l_x.size() > 0);

                list<Vector> l_vx;

                py::buffer_info info = l_x.request();
                auto ptr = static_cast<double*>(info.ptr);

                size_t nb_values = 1;
                for(const auto& r : info.shape)
                  nb_values *= r;
                
                assert(nb_values % l_t.size() == 0);
                size_t n = nb_values / l_t.size();

                for(size_t i = 0 ; i < nb_values ; i+=n)
                {
                  Vector x(n);
                  for(size_t k = 0 ; k < n ; k++)
                    x[k] = *ptr++;
                  l_vx.push_back(x);
                }

                return std::make_unique<SampledTraj<T>>(l_t.cast<std::list<double>>(), l_vx);
              }

              else
                return std::make_unique<SampledTraj<T>>();
            }),
        SAMPLEDTRAJ_T_SAMPLEDTRAJ_CONST_LIST_DOUBLE_REF_CONST_LIST_T_REF,
        "l_t"_a, "l_x"_a)
    ;
  }

  exported_class

    .def(py::init<const std::map<double,T>&>(),
      SAMPLEDTRAJ_T_SAMPLEDTRAJ_CONST_MAP_DOUBLET_REF,
      "m"_a)

    .def("nb_samples", &SampledTraj<T>::nb_samples,
      SIZET_SAMPLEDTRAJ_T_NB_SAMPLES_CONST)

    .def("sampled", [](const SampledTraj<T>& x, double dt, bool keep_original_values)
        {
          return x.sampled(dt,keep_original_values);
        },
      SAMPLEDTRAJ_T_SAMPLEDTRAJ_T_SAMPLED_DOUBLE_BOOL_CONST,
      "dt"_a, "keep_original_values"_a)

    .def("__getitem__", [](const SampledTraj<T>& x, Index_type index) -> const T&
        {
          matlab::test_integer(index);
          return x.at(matlab::input_index(index));
        }, py::return_value_policy::reference_internal)

    .def("__setitem__", [](SampledTraj<T>& x, Index_type index, const T& a)
        {
          matlab::test_integer(index);
          x[matlab::input_index(index)] = a;
        })

    .def("__repr__", [](const SampledTraj<T>& x) {
          std::ostringstream stream;
          stream << x;
          return string(stream.str()); 
        },
      OSTREAM_REF_OPERATOROUT_OSTREAM_REF_CONST_SAMPLEDTRAJ_T_REF)
  ;
}

void export_SampledTraj(py::module& m)
{
  _export_SampledTraj<double>(m, "SampledTraj");
  _export_SampledTraj<Vector>(m, "SampledVectorTraj");

  m

    .def("sqr", (SampledTraj<double> (*)(const SampledTraj<double>&)) &codac2::sqr,
      py::is_operator(),
      SAMPLEDTRAJ_T_SQR_CONST_SAMPLEDTRAJ_T_REF,
      "x1"_a)

    .def("sqrt", (SampledTraj<double> (*)(const SampledTraj<double>&)) &codac2::sqrt,
      py::is_operator(),
      SAMPLEDTRAJ_T_SQRT_CONST_SAMPLEDTRAJ_T_REF,
      "x1"_a)

    .def("pow", (SampledTraj<double> (*)(const SampledTraj<double>&,int)) &codac2::pow,
      py::is_operator(),
      SAMPLEDTRAJ_T_POW_CONST_SAMPLEDTRAJ_T_REF_INT,
      "x1"_a, "x2"_a)

    .def("pow", (SampledTraj<double> (*)(const SampledTraj<double>&,double)) &codac2::pow,
      py::is_operator(),
      SAMPLEDTRAJ_T_POW_CONST_SAMPLEDTRAJ_T_REF_DOUBLE,
      "x1"_a, "x2"_a)

    .def("root", (SampledTraj<double> (*)(const SampledTraj<double>&,int)) &codac2::root,
      py::is_operator(),
      SAMPLEDTRAJ_T_ROOT_CONST_SAMPLEDTRAJ_T_REF_INT,
      "x1"_a, "x2"_a)

    .def("exp", (SampledTraj<double> (*)(const SampledTraj<double>&)) &codac2::exp,
      py::is_operator(),
      SAMPLEDTRAJ_T_EXP_CONST_SAMPLEDTRAJ_T_REF,
      "x1"_a)

    .def("log", (SampledTraj<double> (*)(const SampledTraj<double>&)) &codac2::log,
      py::is_operator(),
      SAMPLEDTRAJ_T_LOG_CONST_SAMPLEDTRAJ_T_REF,
      "x1"_a)

    .def("cos", (SampledTraj<double> (*)(const SampledTraj<double>&)) &codac2::cos,
      py::is_operator(),
      SAMPLEDTRAJ_T_COS_CONST_SAMPLEDTRAJ_T_REF,
      "x1"_a)

    .def("sin", (SampledTraj<double> (*)(const SampledTraj<double>&)) &codac2::sin,
      py::is_operator(),
      SAMPLEDTRAJ_T_SIN_CONST_SAMPLEDTRAJ_T_REF,
      "x1"_a)

    .def("tan", (SampledTraj<double> (*)(const SampledTraj<double>&)) &codac2::tan,
      py::is_operator(),
      SAMPLEDTRAJ_T_TAN_CONST_SAMPLEDTRAJ_T_REF,
      "x1"_a)

    .def("acos", (SampledTraj<double> (*)(const SampledTraj<double>&)) &codac2::acos,
      py::is_operator(),
      SAMPLEDTRAJ_T_ACOS_CONST_SAMPLEDTRAJ_T_REF,
      "x1"_a)

    .def("asin", (SampledTraj<double> (*)(const SampledTraj<double>&)) &codac2::asin,
      py::is_operator(),
      SAMPLEDTRAJ_T_ASIN_CONST_SAMPLEDTRAJ_T_REF,
      "x1"_a)

    .def("atan", (SampledTraj<double> (*)(const SampledTraj<double>&)) &codac2::atan,
      py::is_operator(),
      SAMPLEDTRAJ_T_ATAN_CONST_SAMPLEDTRAJ_T_REF,
      "x1"_a)

#if 0
  20
  inline SampledTraj<T> atan2(const SampledTraj<T>& x1, const SampledTraj<T>& x2)
  inline SampledTraj<T> atan2(const SampledTraj<T>& x1, const T& x2)
  inline SampledTraj<T> atan2(const T& x1, const SampledTraj<T>& x2)
  inline SampledTraj<T> cosh(const SampledTraj<T>& x1)
  inline SampledTraj<T> sinh(const SampledTraj<T>& x1)
  inline SampledTraj<T> tanh(const SampledTraj<T>& x1)
  inline SampledTraj<T> acosh(const SampledTraj<T>& x1)
  inline SampledTraj<T> asinh(const SampledTraj<T>& x1)
  inline SampledTraj<T> atanh(const SampledTraj<T>& x1)
  inline SampledTraj<T> abs(const SampledTraj<T>& x1)
  inline SampledTraj<T> min(const SampledTraj<T>& x1, const SampledTraj<T>& x2)
  inline SampledTraj<T> min(const SampledTraj<T>& x1, const T& x2)
  inline SampledTraj<T> min(const T& x1, const SampledTraj<T>& x2)
  inline SampledTraj<T> max(const SampledTraj<T>& x1, const SampledTraj<T>& x2)
  inline SampledTraj<T> max(const SampledTraj<T>& x1, const T& x2)
  inline SampledTraj<T> max(const T& x1, const SampledTraj<T>& x2)
  inline SampledTraj<T> sign(const SampledTraj<T>& x1)
  inline SampledTraj<T> integer(const SampledTraj<T>& x1)
  inline SampledTraj<T> floor(const SampledTraj<T>& x1)
  inline SampledTraj<T> ceil(const SampledTraj<T>& x1)

  #endif
  ;

}