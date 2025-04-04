/** 
 *  codac2_StyleProperties.cpp
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include "codac2_StyleProperties.h"

using namespace std;
using namespace codac2;

StyleProperties::StyleProperties()
{ }

StyleProperties::StyleProperties(const Color& stroke_color_)
  : stroke_color(stroke_color_), fill_color(Color::none())
{ }

StyleProperties::StyleProperties(const Color& stroke_color_, const string& layer_)
  : stroke_color(stroke_color_), fill_color(Color::none()), layer(layer_)
{ }

StyleProperties::StyleProperties(const Color& stroke_color_, const LineStyle& line_style_)
  : stroke_color(stroke_color_), fill_color(Color::none()), line_style(line_style_)
{ }

StyleProperties::StyleProperties(const Color& stroke_color_, const LineStyle& line_style_, const string& layer_)
  : stroke_color(stroke_color_), fill_color(Color::none()), line_style(line_style_), layer(layer_)
{ }

StyleProperties::StyleProperties(std::initializer_list<Color> colors)
  : stroke_color(*colors.begin())
{
  assert(colors.size() <= 2);
  if (colors.size() == 1)
    fill_color = Color::none();
  else
    fill_color = *std::prev(colors.end());
}

StyleProperties::StyleProperties(std::initializer_list<Color> colors, const string& layer_)
  : stroke_color(*colors.begin()), layer(layer_)
{
  assert(colors.size() <= 2);
  if (colors.size() == 1)
    fill_color = Color::none();
  else
    fill_color = *std::prev(colors.end());
}

StyleProperties::StyleProperties(std::initializer_list<Color> colors, const LineStyle& line_style_)
  : stroke_color(*colors.begin()), line_style(line_style_)
{
  assert(colors.size() <= 2);
  if (colors.size() == 1)
    fill_color = Color::none();
  else
    fill_color = *std::prev(colors.end());
}

StyleProperties::StyleProperties(std::initializer_list<Color> colors, const LineStyle& line_style_, const string& layer_)
  : stroke_color(*colors.begin()), line_style(line_style_), layer(layer_)
{
  assert(colors.size() <= 2);
  if (colors.size() == 1)
    fill_color = Color::none();
  else
    fill_color = *std::prev(colors.end());
}

StyleProperties StyleProperties::inside()
{
  StyleProperties s;
  s.stroke_color = Color::dark_gray();
  s.fill_color = Color::green();
  s.layer = "inside";
  return s;
}

StyleProperties StyleProperties::outside()
{
  StyleProperties s;
  s.stroke_color = Color::dark_gray();
  s.fill_color = Color::cyan();
  s.layer = "outside";
  return s;
}

StyleProperties StyleProperties::boundary()
{
  StyleProperties s;
  s.stroke_color = Color::dark_gray();
  s.fill_color = Color::yellow();
  s.layer = "boundary";
  return s;
}