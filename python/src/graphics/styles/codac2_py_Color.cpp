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
#include <codac2_Color.h>
#include <codac2_StyleProperties.h>
#include "codac2_py_Color_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py):

using namespace std;
using namespace codac2;
namespace py = pybind11;
using namespace pybind11::literals;


void export_Color(py::module& m)
{

  py::enum_<Color::Model>(m, "Model")
    .value("RGB", Color::Model::RGB)
    .value("HSV", Color::Model::HSV)
  ;

  py::class_<Color> exported_color(m, "Color", COLOR_MAIN);
  exported_color
    
    .def_readwrite("m", &Color::m,
      MODEL_COLOR_M)

    .def(py::init<>(),COLOR_COLOR)

    .def(py::init<const std::array<float,3>&,Color::Model>(),
      COLOR_COLOR_CONST_ARRAY_FLOAT3_REF_MODEL,
      "xyz"_a, "m_"_a=Color::Model::RGB)

    .def(py::init<const std::array<float,4>&,Color::Model>(),
      COLOR_COLOR_CONST_ARRAY_FLOAT4_REF_MODEL,
      "xyza"_a, "m_"_a=Color::Model::RGB)

    .def(py::init<const std::initializer_list<float>&,Color::Model>(),
      COLOR_COLOR_CONST_INITIALIZER_LIST_FLOAT_MODEL,
      "xyza"_a, "m_"_a=Color::Model::RGB)

    .def(py::init<const std::string&>(),
      COLOR_COLOR_CONST_STRING_REF,
      "hex_str"_a)


    // Html color

    .def("hex_str", &Color::hex_str,
      STRING_COLOR_HEX_STR_CONST)

    // Conversions

    .def("rgb", &Color::rgb,
      COLOR_COLOR_RGB_CONST)

    .def("hsv", &Color::hsv,
      COLOR_COLOR_HSV_CONST)

    // Overload flux operator

    .def("__str__", [](const Color& c) {
      std::ostringstream oss;
      oss << c;
      return oss.str();
    })

    // Predefined colors

    .def_static("none", &Color::none,
      STATIC_COLOR_COLOR_NONE)

    .def_static("black", &Color::black,
      STATIC_COLOR_COLOR_BLACK_FLOAT,
      "alpha"_a=1.)

    .def_static("white", &Color::white,
      STATIC_COLOR_COLOR_WHITE_FLOAT,
      "alpha"_a=1.)

    .def_static("green", &Color::green,
      STATIC_COLOR_COLOR_GREEN_FLOAT,
      "alpha"_a=1.)

    .def_static("blue", &Color::blue,
      STATIC_COLOR_COLOR_BLUE_FLOAT,
      "alpha"_a=1.)

    .def_static("cyan", &Color::cyan,
      STATIC_COLOR_COLOR_CYAN_FLOAT,
      "alpha"_a=1.)

    .def_static("yellow", &Color::yellow,
      STATIC_COLOR_COLOR_YELLOW_FLOAT,
      "alpha"_a=1.)

    .def_static("red", &Color::red,
      STATIC_COLOR_COLOR_RED_FLOAT,
      "alpha"_a=1.)

    .def_static("dark_gray", &Color::dark_gray,
      STATIC_COLOR_COLOR_DARK_GRAY_FLOAT,
      "alpha"_a=1.)
  ;
}