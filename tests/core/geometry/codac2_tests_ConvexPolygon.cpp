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
#include <codac2_Figure2D.h>
#include <codac2_Approx.h>

using namespace std;
using namespace codac2;

TEST_CASE("ConvexPolygon - degenerate cases")
{
  {
    ConvexPolygon p1({{0,4}});
    CHECK(p1.edges().size() == 1);
    CHECK(p1.edges()[0][0] == IntervalVector({0,4}));
    CHECK(p1.edges()[0][1] == IntervalVector({0,4}));
  }

  {
    ConvexPolygon p1({{0,4},{2,8}});
    CHECK(p1.edges().size() == 2);
    CHECK(p1.edges()[0][0] == IntervalVector({2,8}));
    CHECK(p1.edges()[0][1] == IntervalVector({0,4}));
    CHECK(p1.edges()[1][0] == IntervalVector({0,4}));
    CHECK(p1.edges()[1][1] == IntervalVector({2,8}));
    CHECK(p1 == ConvexPolygon({{2,8},{0,4}}));
  }
}

TEST_CASE("ConvexPolygon - intersection")
{
  {
    ConvexPolygon p1({{0,4},{4,8},{7,0}});
    ConvexPolygon p2({{-1,2},{-1,3},{5,3},{5,2}});

    CHECK((p1 & p2) == ConvexPolygon({
        {3.5,2},{5,2},{5,3},{1.75,3}
      }));
  }
  
  {
    ConvexPolygon p1({{3,2},{1,6},{6,5}});
    ConvexPolygon p2(IntervalVector({{0,4},{0,4}}));

    CHECK((p1 & p2) == ConvexPolygon({
        {3,2},{2,4},{4,4},{4,3}
      }));
  }

  {
    ConvexPolygon p1({{1,2},{3,4},{5,1},{2,1}});
    ConvexPolygon p2(IntervalVector({{2,6},{2,6}}));

    CHECK(Approx<Polygon>(p1 & p2) == ConvexPolygon({
        {2,2},{2,3},{3,4},{4.+1./3.,2.}
      }));
  }

  { // Polygons intersections, test 3 (big box)
    ConvexPolygon p1({{1,2},{3,4},{5,1},{2,1}});
    ConvexPolygon p2(IntervalVector({{-10,10},{-10,10}}));

    CHECK((p1 & p2) == p1); // same polygon
  }

  { // Polygons intersections, test 4 (inner box)
    ConvexPolygon p1({{1,2},{3,4},{5,1},{2,1}});
    ConvexPolygon p2(IntervalVector({{2.8,3},{2.8,3}}));

    CHECK((p1 & p2) == p2); // same box
  }

  { // Polygons intersections, test 5
    ConvexPolygon p1({{2,1},{3,1},{4,2},{4,3},{3,4},{2,4},{1,3},{1,2}});
    ConvexPolygon p2(IntervalVector({{1,4},{1,4}}));

    CHECK((p1 & p2) == p1); // same polygon
  }

  { // Polygons intersections, test 6 (shifted polygon points declaration)
    ConvexPolygon p1({{3,4},{2,4},{1,3},{1,2},{2,1},{3,1},{4,2},{4,3}});
    ConvexPolygon p2(IntervalVector({{1,4},{1,4}}));

    CHECK((p1 & p2) == p1); // same polygon
  }

  { // Polygons intersections, test 7 (degenerate case)
    ConvexPolygon p1({{4000,200}});
    ConvexPolygon p2(IntervalVector({4000,200}));

    CHECK((p1 & p2) == p1); // same polygon
  }

  { // Polygons intersections, test 8
    ConvexPolygon p1({{1,1},{2,4},{7,5},{6,2}});
    ConvexPolygon p2(IntervalVector({{2,6},{1,5}}));

    auto q = p1 & p2;
    CHECK(Approx<Polygon>(q) == Polygon({{2,1.2},{6,2},{6,4.8},{2,4}}));
    CHECK(q.edges().size() == 4);
    CHECK(q.unsorted_vertices().size() == 4);
  }

  { // Polygons intersections, test 9
    ConvexPolygon p1({{1,1},{2,4},{7,5},{6,2}});
    ConvexPolygon p2(IntervalVector({{3,5},{1,5}}));

    auto q = p1 & p2;
    CHECK(Approx<Polygon>(q) == Polygon({{3,1.4},{5,1.8},{5,4.6},{3,4.2}}));
    CHECK(q.edges().size() == 4);
    CHECK(q.unsorted_vertices().size() == 4);
  }

  { // Polygons intersections, test 10 (degenerated box)
    ConvexPolygon p1({{1,1},{2,4},{7,5},{6,2}});
    ConvexPolygon p2(IntervalVector({{4},{1,5}}));

    auto q = p1 & p2;
    CHECK(Approx<Polygon>(q) == Polygon({{4,4.4},{4,1.6}}));
    CHECK(q.unsorted_vertices().size() == 2);
  }

  { // Polygons intersections, test 10 (degenerated polygon)
    ConvexPolygon p1({{1,1},{2,4},{7,5},{6,2}});
    ConvexPolygon p2({{4,1},{4,5}});

    auto q = p1 & p2;
    
    CHECK(Approx<Polygon>(q) == Polygon({{4,4.4},{4,1.6}}));
    CHECK(q.unsorted_vertices().size() == 2);
  }
}