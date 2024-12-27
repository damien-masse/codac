/** 
 *  \file codac2_math.h
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

//#include <numbers>
#include <3rd/gaol/gaol_port.h>

namespace codac2
{
  // As long as GAOL defines pi in common namespace,
  // we use it as default pi value instead of the pi
  // from the standard library.
  //constexpr double pi = std::numbers::pi; // Need C++20
}