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
#include <codac2_Figure2D.h>
#include <codac2_ConvexPolygon.h>
#include <codac2_Approx.h>

using namespace std;
using namespace codac2;

TEST_CASE("geometry - manual")
{
  {
    // [geometry-1-beg]
    OrientationInterval a = orientation({0,0},{1,1},{10,10});
    // a == OrientationInterval::COLINEAR
    OrientationInterval b = orientation({0,0},{1.1,1.1},{10,10});
    // b == OrientationInterval::UNKNOWN (due to floating point uncertainty)

    IntervalVector p1({0,0}), p2({1,1}), p3({-1,2});
    OrientationInterval c = orientation(p1,p2,p3);
    // c == OrientationInterval::CLOCKWISE
    OrientationInterval d = orientation(p3,p2,p1);
    // d == OrientationInterval::COUNTERCLOCKWISE
    // [geometry-1-end]
    
    CHECK(a == OrientationInterval::COLINEAR);
    CHECK(b == OrientationInterval::UNKNOWN);
    CHECK(c == OrientationInterval::CLOCKWISE);
    CHECK(d == OrientationInterval::COUNTERCLOCKWISE);
  }

  {
    // [geometry-2-beg]
    BoolInterval a = aligned({0,0},{1,1},{10,10});
    // a == BoolInterval::TRUE
    BoolInterval b = aligned({0,0},{1.1,1.1},{10,10});
    // b == BoolInterval::UNKNOWN (due to floating point uncertainty)

    IntervalVector p1({0,0}), p2({0,1}), p3({1,10});
    BoolInterval c = aligned(p1,p2,p3);
    // c == BoolInterval::FALSE
    // [geometry-2-end]
    
    CHECK(a == BoolInterval::TRUE);
    CHECK(b == BoolInterval::UNKNOWN);
    CHECK(c == BoolInterval::FALSE);
  }

  {
    //DefaultView::set_axes(axis(0,{-1,4.5}), axis(1,{0,5}));

    // [geometry-3-beg]
    vector<IntervalVector> v {
      {0,3},{1,1},{2,2},{3.5,3.5},
      {0.5,1.2},{1,2},{3,1},{3,3},
      {1,3.5},{2.5,4},{0.1,2.7},{3.2,3.6}
    };

    vector<IntervalVector> hull = convex_hull(v);
    // hull == { {1,1}, {3,1}, {3.5,3.5}, {2.5,4}, {1,3.5}, {0,3}, {0.5,1.2} }
    // [geometry-3-end]
    CHECK(hull[0] == IntervalVector({1,1}));
    CHECK(hull[1] == IntervalVector({3,1}));
    CHECK(hull[2] == IntervalVector({3.5,3.5}));
    CHECK(hull[3] == IntervalVector({2.5,4}));
    CHECK(hull[4] == IntervalVector({1,3.5}));
    CHECK(hull[5] == IntervalVector({0,3}));
    CHECK(hull[6] == IntervalVector({0.5,1.2}));

    //DefaultView::draw_polygone(ConvexPolygon(v), {Color::blue(),Color::blue(0.1)});
    //for(const auto& vi : v)
    //  DefaultView::draw_point(vi.mid());
  }

  {
    // [segment-1-beg]
    Segment u1({{0,0},{2,2}}), v1({{2,4},{0,6}});
    auto p1 = u1 & v1; // the two segments do not intersect
    // p1 == IntervalVector::empty(2)

    Segment u2({{4,0},{0,4}}), v2({{2,0},{4,2}});
    auto p2 = u2 & v2;
    // p2 == IntervalVector({3,1})

    Segment u3({{1,1},{4,4}}), v3({{2,2},{5,5}});
    auto p3 = u3 & v3; // the two segments are colinear
    // p3 == IntervalVector({{2,4},{2,4}})
    // [segment-1-end]
    CHECK(p1 == IntervalVector::empty(2));
    CHECK(p2 == IntervalVector({3,1}));
    CHECK(p3 == IntervalVector({{2,4},{2,4}}));
  }

  {
    // [segment-2-beg]
    Segment u1({{0,0},{2,2}}), v1({{2,4},{0,6}});
    auto p1 = proj_intersection(u1,v1);
    // p1 == IntervalVector({3,3})

    Segment u2({{4,0},{0,4}}), v2({{2,0},{4,2}});
    auto p2 = proj_intersection(u2,v2);
    // p2 == IntervalVector({3,1})

    Segment u3({{1,1},{4,4}}), v3({{2,2},{5,5}});
    auto p3 = proj_intersection(u3,v3);
    // p3 == IntervalVector({{-oo,oo},{-oo,oo}})
    // [segment-2-end]
    CHECK(p1 == IntervalVector({3,3}));
    CHECK(p2 == IntervalVector({3,1}));
    CHECK(p3 == IntervalVector(2));
  }

  {
    // [segment-3-beg]
    IntervalVector p1({0,0}), p2({1,1}), p3({1.1,1.1}), p4({10,10});

    BoolInterval a = colinear(Segment(p1,p2),Segment(p1,p4));
    // a == BoolInterval::TRUE

    BoolInterval b = colinear(Segment(p1,p3),Segment(p1,p4));
    // b == BoolInterval::UNKNOWN (due to floating point uncertainty)
    // [segment-3-end]
    CHECK(a == BoolInterval::TRUE);
    CHECK(b == BoolInterval::UNKNOWN);
  }
}