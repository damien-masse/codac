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
#include <codac2_template_tools.h>
#include <codac2_SlicedTube.h>
#include <codac2_CtcDeriv.h>
#include "codac2_py_Ctc.h"
#include "codac2_py_CtcDeriv_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py):

using namespace std;
using namespace codac2;
namespace py = pybind11;
using namespace pybind11::literals;

void export_CtcDeriv(py::module& m)
{
  py::class_<CtcDeriv> exported(m, "CtcDeriv", CTCDERIV_MAIN);
  exported

    .def(py::init<>(),
      CTCDERIV_CTCDERIV)

    .def("contract", (void (CtcDeriv::*)(Slice<Interval>&,const Slice<Interval>&) const) &CtcDeriv::contract,
      VOID_CTCDERIV_CONTRACT_SLICE_INTERVAL_REF_CONST_SLICE_INTERVAL_REF_CONST,
      "x"_a, "v"_a)
    .def("contract", (void (CtcDeriv::*)(Slice<IntervalVector>&,const Slice<IntervalVector>&) const) &CtcDeriv::contract,
      VOID_CTCDERIV_CONTRACT_SLICE_INTERVALVECTOR_REF_CONST_SLICE_INTERVALVECTOR_REF_CONST,
      "x"_a, "v"_a)
    .def("contract", (void (CtcDeriv::*)(SlicedTube<Interval>&,const SlicedTube<Interval>&) const) &CtcDeriv::contract,
      VOID_CTCDERIV_CONTRACT_SLICEDTUBE_T_REF_CONST_SLICEDTUBE_T_REF_CONST,
      "x"_a, "v"_a)
    .def("contract", (void (CtcDeriv::*)(SlicedTube<IntervalVector>&,const SlicedTube<IntervalVector>&) const) &CtcDeriv::contract,
      VOID_CTCDERIV_CONTRACT_SLICEDTUBE_T_REF_CONST_SLICEDTUBE_T_REF_CONST,
      "x"_a, "v"_a)

  ;
}