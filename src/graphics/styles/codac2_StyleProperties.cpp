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

StyleProperties::StyleProperties(const Color& stroke_color_, const std::string& param1, const std::string& param2, const std::string& param3)
  : stroke_color(stroke_color_), fill_color(Color::none())
{ 
  bool successful_parsing = update_properties(param1, param2, param3);
  assert_release(successful_parsing && "Warning: StyleProperties could not parse line width");
}

StyleProperties::StyleProperties(std::initializer_list<Color> colors, const std::string& param1, const std::string& param2, const std::string& param3)
  : stroke_color(*colors.begin())
{
  assert(colors.size() <= 2);
  if (colors.size() == 1)
    fill_color = Color::none();
  else
    fill_color = *std::prev(colors.end());

  bool successful_parsing = update_properties(param1, param2, param3);
  assert_release(successful_parsing && "Warning: StyleProperties could not parse line width");
}

bool is_valid_number(const string& str) 
{
  for (char c : str)
      if (! StyleProperties::available_linewdith_chars().contains(c)) 
          return false;
  return true;
}

bool StyleProperties::parse_parameter(const std::string& param)
{
  if (param != "")
  {
    if (StyleProperties::available_line_styles().contains(param))
      line_style = param;

    else if (is_valid_number(param))
    {
      try
      {
        line_width = std::stod(param);
        return true;
      }
      catch (const std::invalid_argument& e)
      {
        assert_release(false && "Invalid line width value");
      }
    }

    else
      layer = param;
  }
}

bool StyleProperties::update_properties(const std::string& param1, const std::string& param2, const std::string& param3)
{
  return parse_parameter(param1) && 
         parse_parameter(param2) && 
         parse_parameter(param3);
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
  s.fill_color = Color::light_blue();
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