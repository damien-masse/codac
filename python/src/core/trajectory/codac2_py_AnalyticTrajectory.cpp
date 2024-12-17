/** 
 *  AnalyticTrajectory Python binding
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
#include <codac2_AnalyticTrajectory.h>
#include <codac2_SampledTrajectory.h>
#include "codac2_py_AnalyticTrajectory_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py):
#include "codac2_py_TrajectoryBase_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py):
#include "codac2_py_TrajectoryBase.h"
#include "codac2_py_doc.h"

using namespace std;
using namespace codac2;
namespace py = pybind11;
using namespace pybind11::literals;

template<typename O>
void _export_AnalyticTrajectory(py::module& m, const string& class_name)
{
  py::class_<AnalyticTrajectory<O>> exported_class(m, class_name.c_str(), ANALYTICTRAJECTORY_MAIN);
  export_TrajectoryBase<AnalyticTrajectory<O>>(exported_class);

  exported_class

    .def(py::init<const AnalyticFunction<O>&,const Interval&>(),
      ANALYTICTRAJECTORY_OS_ANALYTICTRAJECTORY_CONST_ANALYTICFUNCTION_O_REF_CONST_INTERVAL_REF,
      "f"_a, "tdomain"_a)

  ;
}

void export_AnalyticTrajectory(py::module& m)
{
  _export_AnalyticTrajectory<ScalarOpValue>(m, "AnalyticTrajectory_Scalar");
  _export_AnalyticTrajectory<VectorOpValue>(m, "AnalyticTrajectory_Vector");
}