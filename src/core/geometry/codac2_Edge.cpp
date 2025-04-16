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
  Edge::Edge(const std::array<IntervalVector,2>& x)
    : std::array<IntervalVector,2>(x)
  {
    assert_release(x[0].size() == 2 && x[1].size() == 2);
  }

  Edge::Edge(const IntervalVector& x1, const IntervalVector& x2)
    : Edge(std::array<IntervalVector,2>({ x1, x2 }))
  { }

  IntervalVector Edge::box() const
  {
    return hull((*this)[0],(*this)[1]);
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
      if(o1 == OrientationInterval::COLINEAR && box().is_superset(e1))
        return BoolInterval::TRUE;

      else if(o2 == OrientationInterval::COLINEAR && box().is_superset(e2))
        return BoolInterval::TRUE;

      else if(o3 == OrientationInterval::COLINEAR && e.box().is_superset(x1))
        return BoolInterval::TRUE;

      else if(o4 == OrientationInterval::COLINEAR && e.box().is_superset(x2))
        return BoolInterval::TRUE;

      else
        return BoolInterval::FALSE;
    }
  }

  BoolInterval Edge::contains(const IntervalVector& p) const
  {
    if(box().is_superset(p))
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

  bool Edge::operator==(const Edge& p) const
  {
    return ((*this)[0] == p[0] && (*this)[1] == p[1])
      || ((*this)[1] == p[0] && (*this)[0] == p[1]);
  }
  
  IntervalVector operator&(const Edge& e1, const Edge& e2)
  {
    return e1.box() & e2.box() & proj_intersection(e1,e2);
  }
  
  IntervalVector proj_intersection(const Edge& e1, const Edge& e2)
  {
    Interval x1 = e1[0][0], y1 = e1[0][1];
    Interval x2 = e1[1][0], y2 = e1[1][1];
    Interval x3 = e2[0][0], y3 = e2[0][1];
    Interval x4 = e2[1][0], y4 = e2[1][1];

    Interval c = ((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4));

    if(c == 0.) // e1 and e2 are parallel
    {
      // They are then either colinear
      //   => infinite intersection points
      // Or not
      //   => no intersection point => empty output

      if((aligned(e1[0],e1[1],e2[0]) & BoolInterval::TRUE) == BoolInterval::TRUE)
        return IntervalVector(2);
      else
        return IntervalVector::empty(2);
    }

    Interval a = (x1*y2-y1*x2);
    Interval b = (x3*y4-y3*x4);

    return {
      (a*(x3-x4)-(x1-x2)*b)/c,
      (a*(y3-y4)-(y1-y2)*b)/c
    };
  }
}