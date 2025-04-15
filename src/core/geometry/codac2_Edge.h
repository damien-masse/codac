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
  class Edge : public std::array<Vector,2>
  {
    public:

      Edge(const std::array<Vector,2>& x);
      Edge(const Vector& x1, const Vector& x2);

      IntervalVector box() const;

      BoolInterval intersects(const Edge& e) const;
      BoolInterval contains(const Vector& p) const;
  };

  IntervalVector proj_intersection(const Edge& e1, const Edge& e2);
}