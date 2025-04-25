/** 
 *  Codac binding (core)
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Damien
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <codac2_IntFullPivLU.h>
#include "codac2_py_IntFullPivLU_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py):

using namespace std;
using namespace codac2;
namespace py = pybind11;
using namespace pybind11::literals;

void export_IntFullPivLU(py::module& m)
{
  py::class_<IntFullPivLU> exported(m, "IntFullPivLU", INTFULLPIVLU_MAIN);
  exported

    .def(py::init<const Matrix &>(),INTFULLPIVLU_INTFULLPIVLU_CONST_MATRIX_REF,
      "M"_a)

    .def(py::init<const IntervalMatrix &>(),INTFULLPIVLU_INTFULLPIVLU_CONST_INTERVALMATRIX_REF,
      "M"_a)

    .def("isInjective",&IntFullPivLU::isInjective,BOOLINTERVAL_INTFULLPIVLU_ISINJECTIVE_CONST)

    .def("isInvertible",&IntFullPivLU::isInvertible,BOOLINTERVAL_INTFULLPIVLU_ISINVERTIBLE_CONST)

    .def("isSurjective",&IntFullPivLU::isSurjective,BOOLINTERVAL_INTFULLPIVLU_ISSURJECTIVE_CONST)

    .def("determinant",&IntFullPivLU::determinant,INTERVAL_INTFULLPIVLU_DETERMINANT_CONST)

    .def("rank",&IntFullPivLU::rank,INTERVAL_INTFULLPIVLU_RANK_CONST)

    .def("dimensionOfKernel",&IntFullPivLU::dimensionOfKernel,INTERVAL_INTFULLPIVLU_DIMENSIONOFKERNEL_CONST)

    .def("kernel",&IntFullPivLU::kernel,INTERVALMATRIX_INTFULLPIVLU_KERNEL_CONST)

    .def("image",&IntFullPivLU::image<IntervalMatrix>,DERIVED_INTFULLPIVLU_IMAGE_CONST_EIGEN_MATRIXBASE_DERIVED_REF_CONST,"M"_a)

    .def("image",[](const IntFullPivLU &S, const Matrix& M) { return S.image(M); },DERIVED_INTFULLPIVLU_IMAGE_CONST_EIGEN_MATRIXBASE_DERIVED_REF_CONST,"M"_a)

    .def("solve",&IntFullPivLU::solve,INTERVALMATRIX_INTFULLPIVLU_SOLVE_CONST_INTERVALMATRIX_REF_CONST,"rhs"_a)

    .def("ReconstructedMatrix",&IntFullPivLU::ReconstructedMatrix,INTERVALMATRIX_INTFULLPIVLU_RECONSTRUCTEDMATRIX_CONST)

    .def("maxPivot",&IntFullPivLU::maxPivot,DOUBLE_INTFULLPIVLU_MAXPIVOT_CONST)

    .def("permutationP",&IntFullPivLU::permutationP,CONST_EIGEN_FULLPIVLU_MATRIX_PERMUTATIONPTYPE_REF_INTFULLPIVLU_PERMUTATIONP_CONST)

    .def("permutationQ",&IntFullPivLU::permutationQ,CONST_EIGEN_FULLPIVLU_MATRIX_PERMUTATIONQTYPE_REF_INTFULLPIVLU_PERMUTATIONQ_CONST)

    .def("matrixLU",&IntFullPivLU::matrixLU,CONST_INTERVALMATRIX_REF_INTFULLPIVLU_MATRIXLU_CONST)
     
  ;
}
