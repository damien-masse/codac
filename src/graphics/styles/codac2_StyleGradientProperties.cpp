/** 
 *  codac2_StyleGradientProperties.cpp
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Simon Rohou, Maël Godard
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include "codac2_StyleGradientProperties.h"

using namespace std;
using namespace codac2;

StyleGradientProperties::StyleGradientProperties()
{ }

StyleGradientProperties::StyleGradientProperties(const ColorMap& cmap_, const std::string& param1, const std::string& param2, const std::string& param3)
  : cmap(cmap_)
{ 
  parse_parameter(param1); parse_parameter(param2); parse_parameter(param3);
}

StyleGradientProperties::StyleGradientProperties(const std::string& param1, const std::string& param2, const std::string& param3)
{
  parse_parameter(param1); parse_parameter(param2); parse_parameter(param3);
}

void StyleGradientProperties::parse_parameter(const std::string& param)
{
  if (param != "")
  {
    if (StyleGradientProperties::available_line_styles().contains(param))
      line_style = param;

    else if (param[0]>='0' && param[0]<='9')
    {
      try
      {
        line_width = std::stod(param);
      }
      catch (const std::invalid_argument& e)
      {
        assert_release_constexpr(false && "Invalid line width value");
      }
    }

    else
      layer = param;
  }
}