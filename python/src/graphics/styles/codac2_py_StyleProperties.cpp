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
#include <codac2_StyleProperties.h>
#include "codac2_py_StyleProperties_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py):

using namespace std;
using namespace codac2;
namespace py = pybind11;
using namespace pybind11::literals;

void export_StyleProperties(py::module& m)
{
  py::enum_<LineStyle>(m, "LineStyle")
    .value("SOLID", LineStyle::SOLID)
    .value("DASHED", LineStyle::DASHED)
    .value("DOTTED", LineStyle::DOTTED)
    .value("DASH_DOTTED", LineStyle::DASH_DOTTED)
    .value("DASH_DOT_DOTTED", LineStyle::DASH_DOT_DOTTED)
  ;

  py::class_<StyleProperties> exported(m, "StyleProperties", STYLEPROPERTIES_MAIN);
  exported
  
    .def(py::init<>(),
      STYLEPROPERTIES_STYLEPROPERTIES)
  
    .def(py::init<const Color&>(),
      STYLEPROPERTIES_STYLEPROPERTIES_CONST_COLOR_REF,
      "color"_a)

    .def(py::init<const Color&, const std::string&>(),
      STYLEPROPERTIES_STYLEPROPERTIES_CONST_COLOR_REF_CONST_STRING_REF,
      "color"_a, "layer"_a)

    .def(py::init<const Color&, const LineStyle&>(),
      STYLEPROPERTIES_STYLEPROPERTIES_CONST_COLOR_REF_CONST_LINESTYLE_REF,
      "color"_a, "line_style"_a)

    .def(py::init<const Color&, const LineStyle&, const std::string&>(),
      STYLEPROPERTIES_STYLEPROPERTIES_CONST_COLOR_REF_CONST_LINESTYLE_REF_CONST_STRING_REF,
      "color"_a, "line_style"_a, "layer"_a)

    .def(py::init(
        [](const std::vector<Color>& v)
        {
          if(v.size() == 1)
            return std::make_unique<StyleProperties>(v[0]);
          else if(v.size() == 2)
            return std::make_unique<StyleProperties,std::initializer_list<Color>>({ v[0], v[1] });
          else
          {
            throw invalid_argument("StyleProperties must be built from one (edge) or two (edge/fill) colors.");
            return std::make_unique<StyleProperties>();
          }
        }),
      STYLEPROPERTIES_STYLEPROPERTIES_INITIALIZER_LIST_COLOR,
      "v"_a)

    .def(py::init(
        [](const std::vector<Color>& v, const std::string& layer)
        {
          if(v.size() == 1)
            return std::make_unique<StyleProperties>(v[0],layer);
          else if(v.size() == 2)
            return std::make_unique<StyleProperties,std::initializer_list<Color>>({ v[0], v[1] },layer);
          else
          {
            throw invalid_argument("StyleProperties must be built from one (edge) or two (edge/fill) colors.");
            return std::make_unique<StyleProperties>();
          }
        }),
      STYLEPROPERTIES_STYLEPROPERTIES_INITIALIZER_LIST_COLOR_CONST_STRING_REF,
      "v"_a, "layer"_a)

    .def(py::init(
        [](const std::vector<Color>& v, const LineStyle& line_style)
        {
          if(v.size() == 1)
            return std::make_unique<StyleProperties>(v[0],line_style);
          else if(v.size() == 2)
            return std::make_unique<StyleProperties,std::initializer_list<Color>>({ v[0], v[1] },line_style);
          else
          {
            throw invalid_argument("StyleProperties must be built from one (edge) or two (edge/fill) colors.");
            return std::make_unique<StyleProperties>();
          }
        }),
      STYLEPROPERTIES_STYLEPROPERTIES_INITIALIZER_LIST_COLOR_CONST_LINESTYLE_REF,
      "v"_a, "line_style"_a)

    .def(py::init(
        [](const std::vector<Color>& v, const LineStyle& line_style, const std::string& layer)
        {
          if(v.size() == 1)
            return std::make_unique<StyleProperties>(v[0],line_style,layer);
          else if(v.size() == 2)
            return std::make_unique<StyleProperties,std::initializer_list<Color>>({ v[0], v[1] },line_style,layer);
          else
          {
            throw invalid_argument("StyleProperties must be built from one (edge) or two (edge/fill) colors.");
            return std::make_unique<StyleProperties>();
          }
        }),
      STYLEPROPERTIES_STYLEPROPERTIES_INITIALIZER_LIST_COLOR_CONST_LINESTYLE_REF_CONST_STRING_REF,
      "v"_a, "line_style"_a, "layer"_a)

    .def_static("inside", &StyleProperties::inside,
      STATIC_STYLEPROPERTIES_STYLEPROPERTIES_INSIDE)

    .def_static("outside", &StyleProperties::outside,
      STATIC_STYLEPROPERTIES_STYLEPROPERTIES_OUTSIDE)

    .def_static("boundary", &StyleProperties::boundary,
      STATIC_STYLEPROPERTIES_STYLEPROPERTIES_BOUNDARY)

  ;
  
  py::implicitly_convertible<py::list,StyleProperties>();
  py::implicitly_convertible<Color,StyleProperties>();
}