/** 
 *  codac2_ConvexPolygon.cpp
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <iostream>
#include "codac2_ConvexPolygon.h"
#include "codac2_geometry.h"
#include "codac2_template_tools.h"

using namespace std;
using namespace codac2;

namespace codac2
{
  ConvexPolygon::ConvexPolygon(std::initializer_list<Vector> vertices)
    : ConvexPolygon(std::vector<Vector>(vertices))
  { }

  ConvexPolygon::ConvexPolygon(const std::vector<Vector>& vertices)
    : ConvexPolygon(vectorVector_to_vectorIntervalVector(vertices))
  { }

  ConvexPolygon::ConvexPolygon(const std::vector<IntervalVector>& vertices)
    : Polygon(convex_hull(vertices))
  { }

  ConvexPolygon::ConvexPolygon(std::initializer_list<Edge> edges)
    : ConvexPolygon(std::vector<Edge>(edges))
  { }

  ConvexPolygon::ConvexPolygon(const vector<Edge>& edges)
    : Polygon([edges] {
      auto l = Polygon(edges).unsorted_vertices();
      return convex_hull({ l.begin(), l.end() });
    }())
  { }

  ConvexPolygon::ConvexPolygon(const IntervalVector& x)
    : Polygon(x)
  { }

  ConvexPolygon operator&(const ConvexPolygon& p1, const ConvexPolygon& p2)
  {
    vector<IntervalVector> inter_vertices;

    for(const auto& v1 : p1.unsorted_vertices())
      if((p2.contains(v1) & BoolInterval::TRUE) == BoolInterval::TRUE)
        inter_vertices.push_back(v1);

    for(const auto& v2 : p2.unsorted_vertices())
      if((p1.contains(v2) & BoolInterval::TRUE) == BoolInterval::TRUE)
        inter_vertices.push_back(v2);

    for(const auto& e1 : p1.edges())
      for(const auto& e2 : p2.edges())
      {
        auto x = e1 & e2;
        if(!x.is_empty())
          inter_vertices.push_back(x.mid());
      }
    
    return ConvexPolygon(inter_vertices);
  }
}