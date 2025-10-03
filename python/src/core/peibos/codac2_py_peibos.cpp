/** 
 *  Codac binding (unsupported)
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Maël Godard
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <codac2_template_tools.h>
 #include <codac2_py_cast.h>
#include <codac2_Vector.h>
#include <codac2_AnalyticFunction.h>
#include <codac2_Parallelepiped.h>
#include <codac2_peibos.h>
#include "codac2_py_peibos_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py):

using namespace std;
using namespace codac2;
namespace py = pybind11;
using namespace pybind11::literals;

void export_peibos(py::module& m)
{
  m.def("parallelepiped_inclusion", 
    [](const py::object& f, const IntervalVector& X)
    {
      return parallelepiped_inclusion(cast<AnalyticFunction<VectorType>>(f), X);
    },
    PARALLELEPIPED_PARALLELEPIPED_INCLUSION_CONST_ANALYTICFUNCTION_VECTORTYPE_REF_CONST_INTERVALVECTOR_REF,
    "f"_a, "X"_a);

  m.def("parallelepiped_inclusion", 
    [](const Vector& z, const IntervalMatrix& JJf, const IntervalMatrix& JJf_punc, const py::object& psi_0, const OctaSym& symmetry, const IntervalVector& X)
    {
      return parallelepiped_inclusion(z, JJf, JJf_punc, cast<AnalyticFunction<VectorType>>(psi_0), symmetry, X);
    },
    PARALLELEPIPED_PARALLELEPIPED_INCLUSION_CONST_VECTOR_REF_CONST_INTERVALMATRIX_REF_CONST_INTERVALMATRIX_REF_CONST_ANALYTICFUNCTION_VECTORTYPE_REF_CONST_OCTASYM_REF_CONST_INTERVALVECTOR_REF,
    "z"_a, "JJf"_a, "JJf_punc"_a, "psi_0"_a, "symmetry"_a, "X"_a);


  m.def("PEIBOS", 
    [](const py::object& f, const py::object& psi_0, const vector<OctaSym>& symmetries, double epsilon, bool verbose = false)
    {
      return PEIBOS(cast<AnalyticFunction<VectorType>>(f), cast<AnalyticFunction<VectorType>>(psi_0), symmetries, epsilon, verbose);
    },
    VECTOR_PARALLELEPIPED_PEIBOS_CONST_ANALYTICFUNCTION_VECTORTYPE_REF_CONST_ANALYTICFUNCTION_VECTORTYPE_REF_CONST_VECTOR_OCTASYM_REF_DOUBLE_BOOL,
    "f"_a, "psi_0"_a, "symmetries"_a, "epsilon"_a, "verbose"_a = false);

  m.def("PEIBOS", 
    [](const py::object& f, const py::object& psi_0, const vector<OctaSym>& symmetries, double epsilon, const Vector& offset, bool verbose = false)
    {
      return PEIBOS(cast<AnalyticFunction<VectorType>>(f), cast<AnalyticFunction<VectorType>>(psi_0), symmetries, epsilon, offset, verbose);
    },
    VECTOR_PARALLELEPIPED_PEIBOS_CONST_ANALYTICFUNCTION_VECTORTYPE_REF_CONST_ANALYTICFUNCTION_VECTORTYPE_REF_CONST_VECTOR_OCTASYM_REF_DOUBLE_CONST_VECTOR_REF_BOOL,
    "f"_a, "psi_0"_a, "symmetries"_a, "epsilon"_a, "offset"_a, "verbose"_a = false);
}