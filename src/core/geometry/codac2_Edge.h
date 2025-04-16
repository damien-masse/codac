/** 
 *  \file codac2_Edge.h
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
  class Edge : public std::array<IntervalVector,2>
  {
    public:

      Edge(const std::array<IntervalVector,2>& x);
      Edge(const IntervalVector& x1, const IntervalVector& x2);

      IntervalVector box() const;

      BoolInterval intersects(const Edge& e) const;
      BoolInterval contains(const IntervalVector& p) const;

      bool operator==(const Edge& p) const;
  };

  IntervalVector operator&(const Edge& e1, const Edge& e2);
  IntervalVector proj_intersection(const Edge& e1, const Edge& e2);
}