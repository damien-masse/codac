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
#include <codac2_pave.h>
#include "codac2_py_pave_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py):
#include "codac2_py_cast.h"

using namespace std;
using namespace codac2;
namespace py = pybind11;
using namespace pybind11::literals;

void export_pave(py::module& m)
{
  m.def("pave", (PavingOut (*)(const IntervalVector&,const CtcBase<IntervalVector>&,double))&codac2::pave,
    PAVINGOUT_PAVE_CONST_INTERVALVECTOR_REF_CONST_CTCBASE_INTERVALVECTOR_REF_DOUBLE,
    "x"_a, "c"_a, "eps"_a);

  m.def("pave", (PavingInOut (*)(const IntervalVector&,const SepBase&,double))&codac2::pave,
    PAVINGINOUT_PAVE_CONST_INTERVALVECTOR_REF_CONST_SEPBASE_REF_DOUBLE,
    "x"_a, "s"_a, "eps"_a);

  m.def("sivia",
      [](const IntervalVector& x, const py::object& f, const py::object& y, double eps)
      {
        if(!is_instance<AnalyticFunction<ScalarType>>(f)
          && !is_instance<AnalyticFunction<VectorType>>(f)
          && !is_instance<AnalyticFunction<VectorType>>(f))
        {
          assert_release("sivia: invalid function type");
        }

        if(is_instance<AnalyticFunction<ScalarType>>(f))
          return sivia(x, cast<AnalyticFunction<ScalarType>>(f), y.cast<Interval>(), eps);

        else if(is_instance<AnalyticFunction<VectorType>>(f))
          return sivia(x, cast<AnalyticFunction<VectorType>>(f), y.cast<IntervalVector>(), eps);

        else
          return sivia(x, cast<AnalyticFunction<MatrixType>>(f), y.cast<IntervalMatrix>(), eps);
      },
    PAVINGINOUT_SIVIA_CONST_INTERVALVECTOR_REF_CONST_ANALYTICFUNCTION_Y_REF_CONST_TYPENAME_Y_DOMAIN_REF_DOUBLE,
    "x"_a, "f"_a, "y"_a, "eps"_a);
}