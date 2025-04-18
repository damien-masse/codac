/** 
 *  \file codac2_Segment.h
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include <vector>
#include <utility>
#include "codac2_Vector.h"
#include "codac2_BoolInterval.h"
#include "codac2_IntervalVector.h"

namespace codac2
{
  class Segment : public std::array<IntervalVector,2>
  {
    public:

      Segment(const std::array<IntervalVector,2>& x);
      Segment(const IntervalVector& x1, const IntervalVector& x2);

      IntervalVector box() const;

      BoolInterval intersects(const Segment& e) const;
      BoolInterval contains(const IntervalVector& p) const;

      bool operator==(const Segment& p) const;
  };

  IntervalVector operator&(const Segment& e1, const Segment& e2);
  IntervalVector proj_intersection(const Segment& e1, const Segment& e2);
  BoolInterval colinear(const Segment& e1, const Segment& e2);

  std::ostream& operator<<(std::ostream& str, const Segment& e);
}