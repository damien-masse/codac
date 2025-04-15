/** 
 *  codac2_Edge.cpp
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <iostream>
#include "codac2_Polygon.h"
#include "codac2_geometry.h"
#include "codac2_hull.h"

using namespace std;
using namespace codac2;

namespace codac2
{
  Edge::Edge(const std::array<Vector,2>& x)
    : std::array<Vector,2>(x)
  {
    assert_release(x[0].size() == 2 && x[1].size() == 2);
  }

  Edge::Edge(const Vector& x1, const Vector& x2)
    : Edge(std::array<Vector,2>({ x1, x2 }))
  { }

  IntervalVector Edge::box() const
  {
    return hull((*this)[0], (*this)[1]);
  }

  BoolInterval Edge::intersects(const Edge& e) const
  {
    const auto& x1 = (*this)[0];
    const auto& x2 = (*this)[1];
    const auto& e1 = e[0];
    const auto& e2 = e[1];

    auto o1 = orientation(x1,x2,e1);
    auto o2 = orientation(x1,x2,e2);
    auto o3 = orientation(e1,e2,x1);
    auto o4 = orientation(e1,e2,x2);

    if(o1 == OrientationInterval::EMPTY || o2 == OrientationInterval::EMPTY || 
        o3 == OrientationInterval::EMPTY || o4 == OrientationInterval::EMPTY)
      return BoolInterval::EMPTY;

    else if(o1 == OrientationInterval::UNKNOWN || o2 == OrientationInterval::UNKNOWN || 
        o3 == OrientationInterval::UNKNOWN || o4 == OrientationInterval::UNKNOWN)
      return BoolInterval::UNKNOWN;

    else if(o1 != o2 && o3 != o4) // general case
      return BoolInterval::TRUE;

    else // special cases of colinearity
    {
      if(o1 == OrientationInterval::COLINEAR && box().contains(e1))
        return BoolInterval::TRUE;

      else if(o2 == OrientationInterval::COLINEAR && box().contains(e2))
        return BoolInterval::TRUE;

      else if(o3 == OrientationInterval::COLINEAR && e.box().contains(x1))
        return BoolInterval::TRUE;

      else if(o4 == OrientationInterval::COLINEAR && e.box().contains(x2))
        return BoolInterval::TRUE;

      else
        return BoolInterval::FALSE;
    }
  }

  BoolInterval Edge::contains(const Vector& p) const
  {
    if(box().contains(p))
    {
      switch(orientation((*this)[0], (*this)[1], p))
      {
        case OrientationInterval::COLINEAR:
          return BoolInterval::TRUE;
        case OrientationInterval::EMPTY:
          return BoolInterval::EMPTY;
        case OrientationInterval::UNKNOWN:
          return BoolInterval::UNKNOWN;
        default:
          return BoolInterval::FALSE;
      }
    }

    return BoolInterval::FALSE;
  }
}