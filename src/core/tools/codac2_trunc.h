/** 
 *  \file codac2_trunc.h
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Simon Rohou
 *  \copyright  Copyright 2025 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include "codac2_Interval.h"

namespace codac2
{
  double trunc(double x);
  Interval trunc(const Interval& x);
  
  double untrunc(double x);
  Interval untrunc(const Interval& x);
}