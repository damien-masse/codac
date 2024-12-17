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
#include <codac2_SampledTrajectory.h>
#include "codac2_py_SampledTrajectory_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py):
#include "codac2_py_TrajectoryBase_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py):
#include "codac2_py_TrajectoryBase.h"
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

  ;
}

void export_SampledTrajectory(py::module& m)
{
  _export_SampledTrajectory<double>(m, "SampledTrajectory_double");
  _export_SampledTrajectory<Vector>(m, "SampledTrajectory_Vector");
}