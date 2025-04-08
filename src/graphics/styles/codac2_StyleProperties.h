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
#include "codac2_Color.h"

namespace codac2
{
  /**
   * \enum LineStyle
   * \brief Line style enumeration, for dashes and dots
   * 
   * This enum is used to specify the line style of a shape, either solid, dashed, dotted, dash-dotted or dash-dot-dotted.
   */
  enum LineStyle
  {
    SOLID,
    DASHED,
    DOTTED,
    DASH_DOTTED,
    DASH_DOT_DOTTED
  };

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
    LineStyle line_style = LineStyle::SOLID;
    std::string layer = "alpha";

    /**
     * \brief Default constructor
     * 
     * Black stroke color, transparent fill color, solid line style and layer "alpha"
     */
    StyleProperties();

    /**
     * \brief Constructor from a stroke color
     * 
     * \param stroke_color Stroke color
     */
    StyleProperties(const Color& stroke_color);

    /**
     * \brief Constructor from a stroke color and a layer name
     * 
     * \param stroke_color Stroke color
     * \param layer Layer name
     */
    StyleProperties(const Color& stroke_color, const std::string& layer);

    /**
     * \brief Constructor from a stroke color and a line style
     * 
     * \param stroke_color Stroke color
     * \param line_style Line style
     */
    StyleProperties(const Color& stroke_color, const LineStyle& line_style);

    /**
     * \brief Constructor from a stroke color, a line style and a layer name
     * 
     * \param stroke_color Stroke color
     * \param line_style Line style
     * \param layer Layer name
     */
    StyleProperties(const Color& stroke_color, const LineStyle& line_style, const std::string& layer);

    /**
     * \brief Constructor from a stroke and eventually a fill color
     * 
     * \param colors list of colors (stroke color and eventually fill color)
     */
    StyleProperties(std::initializer_list<Color> colors);

    /**
     * \brief Constructor from a stroke and eventually a fill color and a layer name
     * 
     * \param colors list of colors (stroke color and eventually fill color)
     * \param layer Layer name
     */
    StyleProperties(std::initializer_list<Color> colors, const std::string& layer);

    /**
     * \brief Constructor from a stroke and eventually a fill color and a line style
     * 
     * \param colors list of colors (stroke color and eventually fill color)
     * \param line_style Line style
     */
    StyleProperties(std::initializer_list<Color> colors, const LineStyle& line_style);

    /**
     * \brief Constructor from a stroke and eventually a fill color, a line style and a layer name
     * 
     * \param colors list of colors (stroke color and eventually fill color)
     * \param line_style Line style
     * \param layer Layer name
     */
    StyleProperties(std::initializer_list<Color> colors, const LineStyle& line_style, const std::string& layer);

    /**
     * \brief Predefined "inside" style, dark gray edge and green fill
     */
    static StyleProperties inside();

    /**
     * \brief Predefined "outside" style, dark gray edge and cyan fill
     */
    static StyleProperties outside();

    /**
     * \brief Predefined "boundary" style, dark gray edge and yellow fill
     */
    static StyleProperties boundary();
  };
}