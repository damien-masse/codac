/** 
 *  Codac binding (graphics)
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <codac2_Interval.h>
#include <codac2_IntervalVector.h>
#include <codac2_drawwhilepaving.h>
#include "codac2_py_drawwhilepaving_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py):

using namespace std;
using namespace codac2;
namespace py = pybind11;
using namespace pybind11::literals;

void export_drawwhilepaving(py::module& m)
{
  m.def("cartesian_drawing", &cartesian_drawing);
  m.def("polar_drawing", &polar_drawing);
  
  m.def("draw_while_paving", (void (*)(const IntervalVector& x0, const CtcBase<IntervalVector>& c, double eps,
      const std::function<void(Figure2D&,const IntervalVector&,const StyleProperties&)>& draw_box,
      std::shared_ptr<Figure2D> fig)) &draw_while_paving,
    VOID_DRAW_WHILE_PAVING_CONST_INTERVALVECTOR_REF_CONST_C_REF_DOUBLE_CONST_FUNCTION_VOID_FIGURE2D_REFCONST_INTERVALVECTOR_REFCONST_STYLEPROPERTIES_REF__REF_SHARED_PTR_FIGURE2D,
    "x0"_a, "c"_a, "eps"_a, "draw_box"_a, "fig"_a=nullptr);

  m.def("draw_while_paving_polar", (void (*)(const IntervalVector& x0, const CtcBase<IntervalVector>& c, double eps,
      std::shared_ptr<Figure2D> fig)) &draw_while_paving,
    VOID_DRAW_WHILE_PAVING_CONST_INTERVALVECTOR_REF_CONST_C_REF_DOUBLE_SHARED_PTR_FIGURE2D,
    "x0"_a, "c"_a, "eps"_a, "fig"_a=nullptr);

  m.def("draw_while_paving", (void (*)(const IntervalVector& x0, const SepBase& s, double eps,
      const std::function<void(Figure2D&,const IntervalVector&,const StyleProperties&)>& draw_box,
      std::shared_ptr<Figure2D> fig)) &draw_while_paving,
    VOID_DRAW_WHILE_PAVING_CONST_INTERVALVECTOR_REF_CONST_S_REF_DOUBLE_CONST_FUNCTION_VOID_FIGURE2D_REFCONST_INTERVALVECTOR_REFCONST_STYLEPROPERTIES_REF__REF_SHARED_PTR_FIGURE2D,
    "x0"_a, "s"_a, "eps"_a, "draw_box"_a, "fig"_a=nullptr);

  m.def("draw_while_paving_polar", (void (*)(const IntervalVector& x0, const SepBase& s, double eps,
      std::shared_ptr<Figure2D> fig)) &draw_while_paving,
    VOID_DRAW_WHILE_PAVING_CONST_INTERVALVECTOR_REF_CONST_S_REF_DOUBLE_SHARED_PTR_FIGURE2D,
    "x0"_a, "s"_a, "eps"_a, "fig"_a=nullptr);
}