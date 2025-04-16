/** 
 *  Codac tests
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <catch2/catch_test_macros.hpp>
#include <codac2_Edge.h>
#include <codac2_Approx.h>

using namespace std;
using namespace codac2;

Interval hull(double x)
{
  if(x == -oo) return { -oo, next_float(-oo) };
  if(x == oo)  return { previous_float(oo), oo };
  else return { x };
}

TEST_CASE("Edge")
{
  SECTION("intersects")
  {
    CHECK(Edge({{0,0},{10,0}}).intersects(Edge({{4,0},{6,0}})) == BoolInterval::TRUE);
    CHECK(Edge({{0,0},{10,0}}).intersects(Edge({{5,0},{15,0}})) == BoolInterval::TRUE);
    CHECK(Edge({{0,0},{10,0}}).intersects(Edge({{5,-5},{5,10}})) == BoolInterval::TRUE);
    CHECK(Edge({{0,0},{10,0}}).intersects(Edge({{15,-5},{15,10}})) == BoolInterval::FALSE);
    CHECK(Edge({{0,0},{10,10}}).intersects(Edge({{1,1},{9,9}})) == BoolInterval::TRUE);
    CHECK(Edge({{1,1},{9,9}}).intersects(Edge({{0,0},{10,10}})) == BoolInterval::TRUE);

    CHECK(Edge(IntervalVector({hull(-oo),2}),Vector({0,2})).intersects(Edge(Vector({0,0}),Vector({0,1}))) == BoolInterval::FALSE);
    CHECK(Edge(IntervalVector({hull(-oo),2}),Vector({0,2})).intersects(Edge(Vector({0,0}),Vector({0,2}))) == BoolInterval::TRUE);
  }

  SECTION("contains bis")
  {
    Vector p1({0,0.}), p2({1.5,3.}), p3({0,2.});
    Edge e1(Vector({-10,-10.}), Vector({10,10.}));
    Edge e2(Vector({1.,1.}), Vector({10,10.}));
    Edge e3(Vector({0,0.}), Vector({2.,4.}));

    CHECK(e1.contains(p1) != BoolInterval::FALSE);
    CHECK(e2.contains(p1) == BoolInterval::FALSE);
    CHECK(e3.contains(p2) != BoolInterval::FALSE);
    CHECK(e3.contains(p3) == BoolInterval::FALSE);
  }

  SECTION("proj_intersection")
  {
    CHECK(proj_intersection({{0,0},{1,1}}, {{3,3},{4,4}})
      == IntervalVector(2));

    CHECK(proj_intersection({{0,0},{1,0}}, {{0,2},{1,2}})
      == IntervalVector::empty(2));

    CHECK(proj_intersection({{2,0},{6,4}},{{6,5},{5,6}})
      == Approx(IntervalVector({6.5,4.5})));

    CHECK(proj_intersection({{-0.7,0.7},{0,1}}, {{-1,0},{-0.7,-0.7}})
      == Approx(IntervalVector({-1.2069,0.4828}), 1e-4));

    CHECK(proj_intersection({{1,2},{2,3}}, {{8,4},{7,4}})
      == IntervalVector({3,4}));

    CHECK(proj_intersection({{4,11},{8,11}}, {{12,12},{12,16}})
      == IntervalVector({12,11}));

    CHECK(proj_intersection({{1,2},{1,2}}, {{8,4},{7,4}})
      == IntervalVector(2));

    CHECK(proj_intersection({{0,1},{0,2}}, {{2,1},{2,2}})
      == IntervalVector::empty(2));

    CHECK(proj_intersection({{-4,4},{0,0}}, {{0,0},{4,4}})
      == IntervalVector({0,0}));

    CHECK(proj_intersection({{-4,4},{0,0}}, {{4,4},{0,0}})
      == IntervalVector({0,0}));

    CHECK(proj_intersection({{0,4},{4,0}}, {{4,4},{0,0}})
      == IntervalVector({2,2}));
  }

  SECTION("intersection")
  {
    CHECK((Edge({0,0},{0,0}) & Edge({0,0},{5,0})) == IntervalVector({0,0})); // degenerate line, horizontal edge line
    CHECK((Edge({0,0},{0,0}) & Edge({0,0},{0,5})) == IntervalVector({0,0})); // degenerate line, vertical edge line

    CHECK((Edge({0,0},{0,1}) & Edge({0,0},{0,5.0})) == IntervalVector({0,{0,1.}})); // vertical edge line (colinear)
    CHECK((Edge({0,0},{0,1}) & Edge({0,0},{0,0.5})) == IntervalVector({0,{0,0.5}})); // vertical edge line (colinear)
    CHECK((Edge({0,0},{1,0}) & Edge({0,0},{5.0,0})) == IntervalVector({{0,1.},0})); // horizontal edge line (colinear)
    CHECK((Edge({0,0},{1,0}) & Edge({0,0},{0.5,0})) == IntervalVector({{0,0.5},0})); // horizontal edge line (colinear)

    CHECK((Edge({7,4},{7,5}) & Edge({7,6},{7,8})) == IntervalVector::empty(2)); // vertical edge line (colinear), no intersection
    CHECK((Edge({4,7},{5,7}) & Edge({6,7},{8,7})) == IntervalVector::empty(2)); // horizontal edge line (colinear), no intersection

    CHECK((Edge({0,0},{3,0}) & Edge({1,2},{1,-999})) == IntervalVector({1,0})); // perpendicular intersection
    CHECK((Edge({0.5,-1},{0.5,5}) & Edge({-3,2},{5,2})) == IntervalVector({0.5,2})); // perpendicular intersection

    CHECK((Edge({0,0},{3,0}) & Edge({1,-10},{1,-999})) == IntervalVector::empty(2)); // perpendicular lines, no intersection
    CHECK((Edge({0.5,-1},{0.5,5}) & Edge({-3,-2},{5,-2})) == IntervalVector::empty(2)); // perpendicular lines, no intersection

    CHECK((Edge({8,4},{9,2}) & Edge({7,3},{9,3})) == IntervalVector({8.5,3})); // perpendicular oblique lines
    CHECK((Edge({8,4},{9,2}) & Edge({8.5,8},{8.5,0})) == IntervalVector({8.5,3})); // perpendicular oblique lines

    CHECK((Edge({8,4},{9,2}) & Edge({7.,3},{7.5,3})) == IntervalVector::empty(2)); // secant oblique lines, no intersection
    CHECK((Edge({8,4},{9,2}) & Edge({8.5,8},{8.5,7})) == IntervalVector::empty(2)); // secant oblique lines, no intersection
    CHECK((Edge({6,-1},{8,1}) & Edge({7.5,0},{9,0})) == IntervalVector::empty(2)); // secant oblique lines, no intersection
    CHECK((Edge({6,-1},{8,1}) & Edge({6.5,0.5},{6.5,2})) == IntervalVector::empty(2)); // secant oblique lines, no intersection
  
    // Other tests
    CHECK((Edge({8,14},{6,13}) & Edge({1,1},{1,14})) == IntervalVector::empty(2));
    CHECK((Edge({hull(-oo),0},{1,0}) & Edge({0,1},{0,-1})) == IntervalVector::zero(2));
    CHECK((Edge({0,0},{1,2}) & Edge({0.5,0},{0.5,2})) == IntervalVector({0.5,1}));
    CHECK((Edge({-1,-1},{0,0}) & Edge({9,3},{0,0})) == IntervalVector({0,0}));
    CHECK((Edge({0,0},{10,0}) & Edge({2,0},{8,0})) == IntervalVector({{2,8},{0}})); // colinear case
    CHECK((Edge({0,0},{10,10}) & Edge({2,2},{8,8})) == IntervalVector({{2,8},{2,8}})); // colinear case
    CHECK((Edge({0,0},{10,0}) & Edge({2,1},{8,1})) == IntervalVector::empty(2)); // parallel case

    CHECK((Edge({2,0},{6,4}) & Edge({6,5},{5,6})) == IntervalVector::empty(2));
    CHECK(proj_intersection(Edge({2,0},{6,4}), Edge({6,5},{5,6})) == IntervalVector({6.5,4.5}));
  }
}