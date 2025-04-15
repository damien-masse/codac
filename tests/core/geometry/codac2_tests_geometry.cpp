/** 
 *  Codac tests
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <catch2/catch_test_macros.hpp>
#include <codac2_geometry.h>

using namespace std;
using namespace codac2;

TEST_CASE("orientation")
{
  // Align points
  CHECK(orientation({0.,0.},{5.,0.},{8.,0.}) == OrientationInterval::COLINEAR);
  CHECK(orientation({0.,0.},{5.,0.},{-2.,0.}) == OrientationInterval::COLINEAR);
  CHECK(orientation({0.,0.},{-2.,0.},{5.,0.}) == OrientationInterval::COLINEAR);
  CHECK(orientation({0.,0.},{0.,0.},{5.,0.}) == OrientationInterval::COLINEAR);
  CHECK(orientation({0.,0.},{0.,0.},{0.,0.}) == OrientationInterval::COLINEAR);

  // Clockwise points
  CHECK(orientation({0.,0.},{0.,-5.},{-2.,-5.}) == OrientationInterval::CLOCKWISE);
  CHECK(orientation({-999.,0.1},{0.,0.},{-999.,-0.1}) == OrientationInterval::CLOCKWISE);

  // Anticlockwise points
  CHECK(orientation({0.,0.},{0.,-5.},{2.,-5.}) == OrientationInterval::COUNTERCLOCKWISE);
  CHECK(orientation({999.,0.1},{0.,0.},{999.,-0.1}) == OrientationInterval::COUNTERCLOCKWISE);
}

TEST_CASE("Alignments")
{
  IntervalVector p1({0.,0.});
  IntervalVector p2({0.,1.});
  IntervalVector p3({0.,10.});
  IntervalVector p4({1.,10.});
  IntervalVector p5({0.,9.+1e-8*Interval(-1.,1.)});
  IntervalVector p6({0.+1e-8*Interval(-1.,1.),9.});
  IntervalVector p7({0.+1e-8*Interval(-1.,1.),9.+1e-8*Interval(-1.,1.)});

  CHECK(aligned(p1,p2,p3) == BoolInterval::TRUE);
  CHECK(aligned(p1,p2,p4) == BoolInterval::FALSE);
  CHECK(aligned(p1,p1,p4) == BoolInterval::TRUE);
  CHECK(aligned(p1,p2,p5) == BoolInterval::TRUE);
  CHECK(aligned(p1,p2,p6) == BoolInterval::UNKNOWN);
  CHECK(aligned(p1,p2,p7) == BoolInterval::UNKNOWN);
  CHECK(aligned({0,0},{1,1},{10,10}) == BoolInterval::TRUE);
  CHECK(aligned({{-oo,0},{-oo,0}},{1,1},{10,10}) == BoolInterval::UNKNOWN);
  CHECK(aligned(IntervalVector(2),{1,1},{10,10}) == BoolInterval::UNKNOWN);
}