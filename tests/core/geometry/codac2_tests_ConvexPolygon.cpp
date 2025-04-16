/** 
 *  Codac tests
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <iostream>
#include <catch2/catch_test_macros.hpp>
#include <codac2_ConvexPolygon.h>

using namespace std;
using namespace codac2;

TEST_CASE("ConvexPolygon")
{
  ConvexPolygon p1({{0,4},{4,8},{7,0}});
  ConvexPolygon p2({{-1,2},{-1,3},{5,3},{5,2}});

  CHECK((p1 & p2) == ConvexPolygon({
      {3.5,2},{5,2},{5,3},{1.75,3}
    }));
}