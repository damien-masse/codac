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

    CHECK(Edge(Vector({next_float(-oo),2}),Vector({0,2})).intersects(Edge(Vector({0,0}),Vector({0,1}))) == BoolInterval::FALSE);
    CHECK(Edge(Vector({next_float(-oo),2}),Vector({0,2})).intersects(Edge(Vector({0,0}),Vector({0,2}))) == BoolInterval::TRUE);
  }

  SECTION("intersects bis")
  {
    Vector p1({0.,0.}), p2({1.5,3.}), p3({0.,2.});
    Edge e1(Vector({-10.,-10.}), Vector({10.,10.}));
    Edge e2(Vector({1.,1.}), Vector({10.,10.}));
    Edge e3(Vector({0.,0.}), Vector({2.,4.}));

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

/* Former tests from Codac 1
  SECTION("Intersection, line/line")
  {
    CHECK((ThickEdge(ThickPoint(0.,0.), ThickPoint(0.,0.)) & ThickEdge(ThickPoint(0.,0.), ThickPoint(5.,0.))) == ThickPoint(0.,0.)); // degenerate line, horizontal edge line
    CHECK((ThickEdge(ThickPoint(0.,0.), ThickPoint(0.,0.)) & ThickEdge(ThickPoint(0.,0.), ThickPoint(0.,5.))) == ThickPoint(0.,0.)); // degenerate line, vertical edge line

    CHECK((ThickEdge(ThickPoint(0.,0.), ThickPoint(0.,1.)) & ThickEdge(ThickPoint(0.,0.), ThickPoint(0.,5.))) == ThickPoint(0.,Interval(0.,1.))); // vertical edge line (colinear)
    CHECK((ThickEdge(ThickPoint(0.,0.), ThickPoint(0.,1.)) & ThickEdge(ThickPoint(0.,0.), ThickPoint(0.,0.5))) == ThickPoint(0.,Interval(0.,0.5))); // vertical edge line (colinear)
    CHECK((ThickEdge(ThickPoint(0.,0.), ThickPoint(1.,0.)) & ThickEdge(ThickPoint(0.,0.), ThickPoint(5.,0.))) == ThickPoint(Interval(0.,1.),0.)); // horizontal edge line (colinear)
    CHECK((ThickEdge(ThickPoint(0.,0.), ThickPoint(1.,0.)) & ThickEdge(ThickPoint(0.,0.), ThickPoint(0.5,0.))) == ThickPoint(Interval(0.,0.5),0.)); // horizontal edge line (colinear)

    CHECK((ThickEdge(ThickPoint(7.,4.), ThickPoint(7.,5.)) & ThickEdge(ThickPoint(7.,6.), ThickPoint(7.,8.))) == ThickPoint(Interval::EMPTY_SET,Interval::EMPTY_SET)); // vertical edge line (colinear), no intersection
    CHECK((ThickEdge(ThickPoint(4.,7.), ThickPoint(5.,7.)) & ThickEdge(ThickPoint(6.,7.), ThickPoint(8.,7.))) == ThickPoint(Interval::EMPTY_SET,Interval::EMPTY_SET)); // horizontal edge line (colinear), no intersection

    CHECK((ThickEdge(ThickPoint(0.,0.), ThickPoint(3.,0.)) & ThickEdge(ThickPoint(1.,2.), ThickPoint(1.,-999.))) == ThickPoint(1.,0.)); // perpendicular intersection
    CHECK((ThickEdge(ThickPoint(0.5,-1.), ThickPoint(0.5,5.)) & ThickEdge(ThickPoint(-3.,2.), ThickPoint(5.,2.))) == ThickPoint(0.5,2.)); // perpendicular intersection

    CHECK((ThickEdge(ThickPoint(0.,0.), ThickPoint(3.,0.)) & ThickEdge(ThickPoint(1.,-10.), ThickPoint(1.,-999.))) == ThickPoint(Interval::EMPTY_SET,Interval::EMPTY_SET)); // perpendicular lines, no intersection
    CHECK((ThickEdge(ThickPoint(0.5,-1.), ThickPoint(0.5,5.)) & ThickEdge(ThickPoint(-3.,-2.), ThickPoint(5.,-2.))) == ThickPoint(Interval::EMPTY_SET,Interval::EMPTY_SET)); // perpendicular lines, no intersection

    CHECK((ThickEdge(ThickPoint(8.,4.), ThickPoint(9.,2.)) & ThickEdge(ThickPoint(7.,3.), ThickPoint(9.,3.))) == ThickPoint(8.5, 3.)); // perpendicular oblique lines
    CHECK((ThickEdge(ThickPoint(8.,4.), ThickPoint(9.,2.)) & ThickEdge(ThickPoint(8.5,8.), ThickPoint(8.5,0.))) == ThickPoint(8.5, 3.)); // perpendicular oblique lines

    CHECK((ThickEdge(ThickPoint(8.,4.), ThickPoint(9.,2.)) & ThickEdge(ThickPoint(7.,3.), ThickPoint(7.5,3.))) == ThickPoint(Interval::EMPTY_SET,Interval::EMPTY_SET)); // secant oblique lines, no intersection
    CHECK((ThickEdge(ThickPoint(8.,4.), ThickPoint(9.,2.)) & ThickEdge(ThickPoint(8.5,8.), ThickPoint(8.5,7.))) == ThickPoint(Interval::EMPTY_SET,Interval::EMPTY_SET)); // secant oblique lines, no intersection
    CHECK((ThickEdge(ThickPoint(6.,-1.), ThickPoint(8.,1.)) & ThickEdge(ThickPoint(7.5,0.), ThickPoint(9.,0.))) == ThickPoint(Interval::EMPTY_SET,Interval::EMPTY_SET)); // secant oblique lines, no intersection
    CHECK((ThickEdge(ThickPoint(6.,-1.), ThickPoint(8.,1.)) & ThickEdge(ThickPoint(6.5,0.5), ThickPoint(6.5,2.))) == ThickPoint(Interval::EMPTY_SET,Interval::EMPTY_SET)); // secant oblique lines, no intersection
  
    // Other tests
    CHECK((ThickEdge(ThickPoint(8.,14.), ThickPoint(6.,13.)) & ThickEdge(ThickPoint(1.,1.), ThickPoint(1.,14.))) == ThickPoint(Interval::EMPTY_SET,Interval::EMPTY_SET));
  }
*/
}