/** 
 *  \file codac2_StyleProperties.h
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include <string>
#include <set>
#include "codac2_Color.h"

namespace codac2
{

  /**
   * \struct StyleProperties
   * \brief Style properties structure, to specify the style of a shape
   * 
   * This class is used to specify the style of a shape, including the stroke color, fill color, line style and layer.
   */
  struct StyleProperties
  {
    Color stroke_color = Color::black();
    Color fill_color = Color::none();
    std::string line_style = "-";
    double line_width = 0.;
    std::string layer = "alpha";

    /**
     * \brief Default constructor
     * 
     * Black stroke color, transparent fill color, solid line style and layer "alpha"
     */
    StyleProperties();

    /**
     * \brief Constructor from a stroke color, and two optional parameters
     * 
     * \param stroke_color Stroke color
     * \param param1 Optional parameter, can be layer name or line style
     * \param param2 Optional parameter, can be layer name or line style
     * \param param3 Optional parameter, can be layer name or line style
     */
    StyleProperties(const Color& stroke_color, const std::string& param1 = "", const std::string& param2 = "", const std::string& param3 = "");

    /**
     * \brief Constructor from a stroke and eventually a fill color, and two optional parameters
     * 
     * \param colors list of colors (stroke color and eventually fill color)
     * \param param1 Optional parameter, can be layer name or line style
     * \param param2 Optional parameter, can be layer name or line style
     * \param param3 Optional parameter, can be layer name or line style
     */
    StyleProperties(std::initializer_list<Color> colors, const std::string& param1 = "", const std::string& param2 = "", const std::string& param3 = "");

    /**
     * \brief Parse a parameter and update the style properties accordingly
     * 
     * \param param Parameter to parse, can be a line style, line width or layer name
     * \return true if the parameter was successfully parsed, false otherwise
     */
    bool parse_parameter(const std::string& param);

    /**
     * \brief Update the style properties from three parameters
     * 
     * \param param1 First parameter, can be a line style, line width or layer name
     * \param param2 Second parameter, can be a line style, line width or layer name
     * \param param3 Third parameter, can be a line style, line width or layer name
     * \return true if the parameters were successfully parsed, false otherwise
     */
    bool update_properties(const std::string& param1, const std::string& param2, const std::string& param3);

    /**
     * \brief Predefined "inside" style, dark gray edge and green fill
     */
    static StyleProperties inside();

    /**
     * \brief Predefined "outside" style, dark gray edge and light blue fill
     */
    static StyleProperties outside();

    /**
     * \brief Predefined "boundary" style, dark gray edge and yellow fill
     */
    static StyleProperties boundary();

    /**
     * \brief Set of available line styles
     */
    static std::set<std::string> available_line_styles()
    {
      std::set<std::string> line_styles={"-", "--", "..", "-.", "-.."};
      return line_styles;
    };

    /**
     * \brief Set of available characters for line width
     */
    static std::set<char> available_linewdith_chars()
    {
      std::set<char> available_chars = {'0','1','2','3','4','5','6','7','8','9','.'};
      return available_chars;
    };

  };
}