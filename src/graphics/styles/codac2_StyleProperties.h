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
  enum LineStyle
  {
    SOLID,
    DASHED,
    DOTTED,
    DASH_DOTTED,
    DASH_DOT_DOTTED
  };

  struct StyleProperties
  {
    Color stroke_color = Color::black();
    Color fill_color = Color::none();
    LineStyle line_style = LineStyle::SOLID;
    std::string layer = "alpha";

    StyleProperties();
    StyleProperties(const Color& stroke_color);
    StyleProperties(const Color& stroke_color, const std::string& layer);
    StyleProperties(const Color& stroke_color, const LineStyle& line_style);
    StyleProperties(const Color& stroke_color, const LineStyle& line_style, const std::string& layer);
    StyleProperties(std::initializer_list<Color> colors);
    StyleProperties(std::initializer_list<Color> colors, const std::string& layer);
    StyleProperties(std::initializer_list<Color> colors, const LineStyle& line_style);
    StyleProperties(std::initializer_list<Color> colors, const LineStyle& line_style, const std::string& layer);

    static StyleProperties inside();
    static StyleProperties outside();
    static StyleProperties boundary();
  };
}