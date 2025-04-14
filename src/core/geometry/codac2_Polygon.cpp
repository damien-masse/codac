/** 
 *  codac2_Polygon.cpp
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <iostream>
#include "codac2_Polygon.h"
#include "codac2_geometry.h"

using namespace std;
using namespace codac2;

namespace codac2
{
  Polygon::Polygon(initializer_list<Vector> vertices)
    : Polygon(std::vector<Vector>(vertices))
  { }

  Polygon::Polygon(const std::vector<Vector>& vertices)
    : _edges(
      [vertices]
      {
        assert_release(vertices.size() > 1);
        vector<Edge> edges;

        Index i = 0;
        for(const auto& vi : vertices)
        {
          assert_release(vi.size() == 2);

          if(i == 0) edges.push_back({ vi, vi });
          else       edges.push_back({ edges[i-1][1], vi });

          i++;
        }

        edges[0][0] = edges[edges.size()-1][1]; // closing the polygon
        return edges;
      }())
  { }

  Polygon::Polygon(initializer_list<Edge> edges)
    : _edges(edges)
  { }

  Polygon::Polygon(const vector<Edge>& edges)
    : _edges(edges)
  { }
  
  const vector<Edge>& Polygon::edges() const
  {
    return _edges;
  }
  
  list<Vector> Polygon::unsorted_vertices() const
  {
    list<Vector> l;
    for(const auto& ei : _edges)
      l.push_back(ei[0]);

    // Removing duplicates
    l.sort([](const Vector& a, const Vector& b) {
        return a[0] < b[0] || (a[0] == b[0] && a[1] < b[1]);
      });
    l.unique();
    return l;
  }

  BoolInterval Polygon::contains(const Vector& p) const
  {
    assert_release(p.size() == 2);

    for(const auto& edge_k : _edges)
    {
      auto b = edge_k.contains(p);
      if(b == BoolInterval::TRUE || b == BoolInterval::UNKNOWN)
        return b;
    }

    bool retry;
    double eps = 0.;
    Edge transect(Vector(2),Vector(2));
    // ^ selected transect (horizontal ray) for crossing the polygon.
    // Odd number of crossing => point is inside
    // Even number of crossing => point is outside

    // Some limit cases must be considered, for instance when the
    // transect crosses exactly one vertex of the polygon, and/or 
    // when the point to be tested is this very vertex. Below, we
    // generate a convenient transect avoiding such configuration.
    do
    {
      retry = false;

      // Horizontal ray candidate:
      Edge try_transect { Vector({next_float(-oo),p[1] + eps}), p };

      // The ray may pass through the vertices, we must double counting
      for(const auto& pi : unsorted_vertices())
      {
        auto b = try_transect.contains(pi);
        if(b == BoolInterval::UNKNOWN)
          return BoolInterval::UNKNOWN;
        else if(b == BoolInterval::TRUE)
        {
          eps = Interval(0,1).rand();
          retry = true;
          break;
        }
      }

      if(!retry)
        transect = try_transect;

    } while(retry);

    // Now the number i of crossing can be computed.
    Index i = 0;

    for(const auto& edge_k : _edges)
    {
      if(transect.box().is_strict_superset(edge_k.box()))
        continue; // case of a ray passing over a colinear edge

      switch(transect.intersects(edge_k))
      {
        case BoolInterval::TRUE:
          i++;
          break;

        case BoolInterval::FALSE:
          // no intersection
          break;

        case BoolInterval::UNKNOWN:
          return BoolInterval::UNKNOWN;

        case BoolInterval::EMPTY:
        default:
          assert(false && "BoolInterval::EMPTY should not happen");
      }
    }
    
    return (i%2 == 0) ? BoolInterval::FALSE : BoolInterval::TRUE;
  }
}