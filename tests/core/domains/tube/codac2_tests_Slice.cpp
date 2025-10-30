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
#include <codac2_Slice.h>
#include <codac2_TDomain.h>
#include <codac2_SlicedTube.h>

using namespace std;
using namespace codac2;

TEST_CASE("Slice")
{
  SECTION("Slice inversion")
  {
    auto tdomain = create_tdomain({0,1});
    SlicedTube<Interval> x(tdomain, {0,10});
    auto sx = x.first_slice();
    x.set({2,3},0);
    x.set({5,6},1);
    CHECK(sx->invert({4,6}, 0) == Interval::empty());
    CHECK(sx->invert({2.5,6}, 0) == 0.);
    CHECK(sx->invert({0,1}, 1) == Interval::empty());
    CHECK(sx->invert({2.5,6}, 1) == 1.);
    CHECK(sx->invert({2.5,6}, {0.2,0.5}) == Interval(0.2,0.5));
    CHECK(sx->invert({2.5,6}, 0.2) == Interval(0.2));
  }

  SECTION("Slice inversion, unbounded derivative")
  {
    auto tdomain = create_tdomain({0,1});
    SlicedTube<Interval> x(tdomain, {0,10});
    SlicedTube<Interval> v(tdomain, Interval());
    auto sx = x.first_slice();
    auto sv = v.first_slice();
    CHECK(sx->invert(5., *sv, sx->t0_tf()) == Interval(0,1));
    CHECK(sx->invert(15., *sv, sx->t0_tf()) == Interval::empty());
  }
}