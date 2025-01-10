/** 
 *  SampledTrajectory Python binding
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <limits>
#include <sstream>
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>
#include <codac2_SampledTrajectory.h>
#include "codac2_py_SampledTrajectory_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py):
#include "codac2_py_TrajectoryBase_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py):
#include "codac2_py_TrajectoryBase.h"
#include "codac2_Trajectory_operator.h"
#include "codac2_py_doc.h"

using namespace std;
using namespace codac2;
namespace py = pybind11;
using namespace pybind11::literals;

template<typename T>
void _export_SampledTrajectory(py::module& m, const string& class_name)
{
  py::class_<SampledTrajectory<T>> exported_class(m, class_name.c_str(), SAMPLEDTRAJECTORY_MAIN);
  export_TrajectoryBase<SampledTrajectory<T>>(exported_class);

  exported_class

    .def(py::init<>(),
      SAMPLEDTRAJECTORY_T_SAMPLEDTRAJECTORY)
    
  ;

  if constexpr(std::is_same_v<T,double>)
  {
    exported_class
      .def(py::init<const std::list<double>&,const std::list<T>&>(),
        SAMPLEDTRAJECTORY_T_SAMPLEDTRAJECTORY_CONST_LIST_DOUBLE_REF_CONST_LIST_T_REF,
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

                return std::make_unique<SampledTrajectory<T>>(l_t.cast<std::list<double>>(), l_vx);
              }

              else
                return std::make_unique<SampledTrajectory<T>>();
            }),
        SAMPLEDTRAJECTORY_T_SAMPLEDTRAJECTORY_CONST_LIST_DOUBLE_REF_CONST_LIST_T_REF,
        "l_t"_a, "l_x"_a)
    ;
  }

  exported_class

    .def(py::init<const std::map<double,T>&>(),
      SAMPLEDTRAJECTORY_T_SAMPLEDTRAJECTORY_CONST_MAP_DOUBLET_REF,
      "m"_a)

    .def("nb_samples", &SampledTrajectory<T>::nb_samples,
      SIZET_SAMPLEDTRAJECTORY_T_NB_SAMPLES_CONST)

    .def("sampled", [](const SampledTrajectory<T>& x, double dt, bool keep_original_values)
        {
          return x.sampled(dt,keep_original_values);
        },
      SAMPLEDTRAJECTORY_T_SAMPLEDTRAJECTORY_T_SAMPLED_DOUBLE_BOOL_CONST,
      "dt"_a, "keep_original_values"_a)

    .def("__getitem__", [](const SampledTrajectory<T>& x, Index_type index) -> const T&
        {
          matlab::test_integer(index);
          return x.at(matlab::input_index(index));
        }, py::return_value_policy::reference_internal)

    .def("__setitem__", [](SampledTrajectory<T>& x, Index_type index, const T& a)
        {
          matlab::test_integer(index);
          x[matlab::input_index(index)] = a;
        })

    .def("__repr__", [](const SampledTrajectory<T>& x) {
          std::ostringstream stream;
          stream << x;
          return string(stream.str()); 
        },
      OSTREAM_REF_OPERATOROUT_OSTREAM_REF_CONST_SAMPLEDTRAJECTORY_T_REF)
  ;
}

void export_SampledTrajectory(py::module& m)
{
  _export_SampledTrajectory<double>(m, "SampledTrajectory_double");
  _export_SampledTrajectory<Vector>(m, "SampledTrajectory_Vector");
}