/** 
 *  codac2_geometry.cpp
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include "codac2_geometry.h"

using namespace std;
using namespace codac2;

namespace codac2
{
  OrientationInterval orientation(const IntervalVector& p1, const IntervalVector& p2, const IntervalVector& p3)
  {
    assert_release(p1.size() == 2 && p2.size() == 2 && p3.size() == 2);

    Interval o = (p2[1]-p1[1])*(p3[0]-p2[0])-(p2[0]-p1[0])*(p3[1]-p2[1]); 

    if(o.is_empty())
      return OrientationInterval::EMPTY;
    else if(o == 0.)
      return OrientationInterval::COLINEAR;
    else if(o.contains(0.))
      return OrientationInterval::UNKNOWN;
    else if(o.lb() > 0.)
      return OrientationInterval::CLOCKWISE;
    else 
      return OrientationInterval::COUNTERCLOCKWISE;
  }

  BoolInterval aligned(const IntervalVector& p1, const IntervalVector& p2, const IntervalVector& p3)
  {
    assert_release(p1.size() == 2 && p2.size() == 2 && p3.size() == 2);

    Interval cross_product = (p2[0]-p1[0])*(p3[1]-p1[1])-(p2[1]-p1[1])*(p3[0]-p1[0]);

    if(cross_product == 0.)
      return BoolInterval::TRUE;

    else if(!cross_product.contains(0.))
      return BoolInterval::FALSE;

    else
      return BoolInterval::UNKNOWN;
  }
}