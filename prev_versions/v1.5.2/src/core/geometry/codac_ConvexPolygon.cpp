/** 
 *  ConvexPolygon class
 * ----------------------------------------------------------------------------
 *  \date       2018
 *  \author     Simon Rohou
 *  \copyright  Copyright 2021 Codac Team
 *  \license    This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 */

#include <iostream>
#include <iomanip>
#include "codac_IntervalMatrix.h"
#include "codac_ConvexPolygon.h"
#include "codac_GrahamScan.h"
#include "codac_VIBesFig.h"
#include "codac_polygon_arithmetic.h"

using namespace std;
using namespace ibex;

namespace codac
{
  // Definition

  ConvexPolygon::ConvexPolygon()
    : ConvexPolygon(IntervalVector(2,Interval(-9999.,9999.)))
  {

  }

  ConvexPolygon::ConvexPolygon(const ConvexPolygon& p)
    : Polygon(p.m_v_floating_pts)
  {
    // Already convex
  }

  ConvexPolygon::ConvexPolygon(const IntervalVector& box)
    : Polygon()
  {
    assert(box.size() == 2);
    assert(!box.is_empty());

    ThickPoint::push(box, m_v_floating_pts);
    m_v_floating_pts = GrahamScan::convex_hull(m_v_floating_pts);
  }

  ConvexPolygon::ConvexPolygon(const vector<ThickPoint>& v_thick_pts)
    : Polygon()
  {
    for(const auto& thick_pt : v_thick_pts)
      if(thick_pt.does_not_exist())
        return; // one undefined point means undefined polygon

    Vector center(2,0.);
    for(const auto& thick_pt : v_thick_pts)
      center += thick_pt.box().mid();
    center *= 1./v_thick_pts.size();

    vector<Vector> v_pts;
    for(const auto& thick_pt : v_thick_pts)
    {
      if(thick_pt.does_not_exist()) // empty polygon
      {
        m_v_floating_pts.clear();
        return;
      }

      if(thick_pt.box().max_diam() < 1e-10)
      {
        Vector outside(2); // attention: this is not reliable
        outside[0] = center[0] < thick_pt.box()[0].lb() ? thick_pt.box()[0].ub() : thick_pt.box()[0].lb();
        outside[1] = center[1] < thick_pt.box()[1].lb() ? thick_pt.box()[1].ub() : thick_pt.box()[1].lb();
        v_pts.push_back(outside);
      }

      else
      {
        vector<Vector> v_pts_ = thick_pt.bounds_pts();
        for(const auto& pt : v_pts_)
          v_pts.push_back(pt);
      }
    }

    m_v_floating_pts = GrahamScan::convex_hull(v_pts);
  }

  ConvexPolygon::ConvexPolygon(const vector<Vector>& v_floating_pts, bool convex_and_convention_order)
    : Polygon(v_floating_pts)
  {
    if(!convex_and_convention_order)
      m_v_floating_pts = GrahamScan::convex_hull(m_v_floating_pts);
  }


  // Tests

  const BoolInterval ConvexPolygon::is_subset(const ConvexPolygon& p) const
  {
    BoolInterval is_subset = YES;

    for(const auto& pt : vertices())
    {
      is_subset = is_subset && p.contains(ThickPoint(pt));
      if(is_subset == NO)
        return NO;
    }

    return is_subset;
  }
  
  bool ConvexPolygon::is_unbounded() const
  {
    return false; // todo
  }

  const BoolInterval ConvexPolygon::encloses(const ThickPoint& p) const
  {
    cout << "encloses(): deprecated. Use contains() instead." << endl;
    return contains(p);
  }

  const BoolInterval ConvexPolygon::intersects(const ConvexPolygon& p) const
  {
    // todo: optimize this
    ConvexPolygon inter = *this & p;
    return inter.is_empty() ? BoolInterval::NO : BoolInterval::MAYBE;
  }

  const BoolInterval ConvexPolygon::contains(const ThickPoint& p) const
  {
    if(p.does_not_exist() || is_empty())
      return NO;

    if(!p.box().intersects(box()))
      return NO; // fast test

    // Using the ray tracing method:
    //   A ray is defined from p to the right ; if it crosses
    //   'a' times one of the edges, and (a & 1), then p is inside

    vector<ThickEdge> v_edges = edges();
    int a = 0; // the crossing number counter
    size_t n = v_edges.size();
    const ThickEdge ray(p, ThickPoint(box()[0].ub()+1., p[1])); // horizontal edge to the right
    ThickPoint prev_e = v_edges[n-1] & ray;

    // Loop through all edges of the polygon, looking for intersections
    for(size_t i = 0 ; i < n ; i++)
    {
      if(!v_edges[i].box()[1].intersects(p[1])) // no feasible intersection
        continue;

      // Intersecting point
      const ThickPoint e = v_edges[i] & ray;
      if(!e.does_not_exist()) // intersection to the left of p, not considered
      {
        if(e[0].intersects(p[0]))
          return MAYBE; // uncertainty

        if(e[1].intersects(box()[1].lb()) || e[1].intersects(box()[1].ub()))
          continue; // the ray is horizontally tangent to the polygon

        if(prev_e[0].intersects(e[0]))
          continue; // point already dealt

        if(e[0].lb() > p[0].ub()) // intersection on the right
          a++;
      }

      prev_e = e;
    }

    return (a & 1) ? YES : NO;
  }


  // Setting values

  const ConvexPolygon& ConvexPolygon::inflate(double rad)
  {
    // todo
    return *this;
  }

  const Interval surface(const IntervalVector& p1, const IntervalVector& p2, const IntervalVector& p3)
  {
    // Shoelace formula
    return 0.5*abs(p1[0]*p2[1] + p2[0]*p3[1] + p3[0]*p1[1] - p2[0]*p1[1] - p3[0]*p2[1] - p1[0]*p3[1]);
  }

  const ConvexPolygon& ConvexPolygon::simplify(size_t max_edges)
  {
    assert(max_edges >= 3);
    IntervalVector box_limit = box();
    box_limit.inflate(box_limit.max_diam()*2.); // avoid wrong optimizations in case of almost parallel edges

    while(m_v_floating_pts.size() > max_edges)
    {
      size_t n = m_v_floating_pts.size();

      // Finding shortest edge, to be removed
      double min_surf = 0.;
      size_t min_i = 0;
      ThickPoint min_inter;

      for(size_t i = 0 ; i < n ; i++)
      {
        const ThickEdge e1 = ThickEdge(ThickPoint(m_v_floating_pts[(i-1+n)%n]), ThickPoint(m_v_floating_pts[i]));
        const ThickEdge e2 = ThickEdge(ThickPoint(m_v_floating_pts[(i+1)%n]), ThickPoint(m_v_floating_pts[(i+2)%n]));

        if(ThickEdge::parallel(e1, e2) == NO)
        {
          // Computing new vertex: intersection of neighbor edges
          ThickPoint inter = ThickEdge::proj_intersection(e1, e2);

          assert(!inter.is_unbounded());
          assert(e1 != e2);
          assert(!e1.does_not_exist() && !e2.does_not_exist());
          if(inter.does_not_exist())
          {
            cout << "inter.does_not_exist():" << endl;
            cout << "inter=" << inter << ", e1=" << e1 << ", e2=" << e2 << endl;
            cout << box_limit << endl;
            *this = ConvexPolygon(box_limit);
            return *this;
          }
          assert(!inter.does_not_exist());

          if(!box_limit.contains(inter.mid()))
            continue;

          if(GrahamScan::orientation(m_v_floating_pts[i], inter.box(), m_v_floating_pts[(i+1)%n]) == OrientationInterval::CLOCKWISE)
            continue;

          double surf = surface(m_v_floating_pts[i], inter.box(), m_v_floating_pts[(i+1)%n]).ub();

          if(min_i == 0 || surf < min_surf) // keeping the simplification that has less impact
          {
            min_inter = inter;
            min_surf = surf;
            min_i = i;
          }
        }
      }

      if(min_inter.does_not_exist()) // if not able to simplify (ex: square case)
        return *this;

      // Updating one of the vertices, removing the other one
      m_v_floating_pts[min_i] = min_inter.mid(); // todo: attention: reliability is lost here
      m_v_floating_pts.erase(m_v_floating_pts.begin() + ((min_i+1)%n));
    }
    
    assert(m_v_floating_pts.size() <= max_edges);
    return *this;
  }

  const ConvexPolygon& ConvexPolygon::rotate(const Interval& theta, const IntervalVector& center)
  {
    IntervalMatrix rot(3,3);
    rot[0][0] = cos(theta); rot[0][1] = -sin(theta); rot[0][2] = 0.;
    rot[1][0] = sin(theta); rot[1][1] = cos(theta);  rot[1][2] = 0.;
    rot[2][0] = 0.;         rot[2][1] = 0.;          rot[2][2] = 1.;

    IntervalMatrix tra(3,3);
    tra[0][0] = 1.; tra[0][1] = 0.; tra[0][2] = center[0];
    tra[1][0] = 0.; tra[1][1] = 1.; tra[1][2] = center[1];
    tra[2][0] = 0.; tra[2][1] = 0.; tra[2][2] = 1.;

    IntervalMatrix rtra(3,3);
    rtra[0][0] = 1.; rtra[0][1] = 0.; rtra[0][2] = -center[0];
    rtra[1][0] = 0.; rtra[1][1] = 1.; rtra[1][2] = -center[1];
    rtra[2][0] = 0.; rtra[2][1] = 0.; rtra[2][2] = 1.;

    vector<ThickPoint> v_thick_pts(m_v_floating_pts.size());
    for(size_t i = 0 ; i < m_v_floating_pts.size() ; i++)
    {
      IntervalVector pt(3);
      pt[0] = m_v_floating_pts[i][0];
      pt[1] = m_v_floating_pts[i][1];
      pt[2] = 1.;
      v_thick_pts[i] = ThickPoint((tra * rot * rtra * pt).subvector(0,1));
    }

    *this = ConvexPolygon(v_thick_pts);
    return *this;
  }


  // Operators

  const IntervalVector ConvexPolygon::fast_intersection(const IntervalVector& x) const
  {
    // todo: keep this? is it faster than the generic one?

    assert(x.size() == 2);

    if(is_empty() || x.is_empty())
      return IntervalVector(2, Interval::EMPTY_SET);

    IntervalVector reduced_x = x & box();

    if(reduced_x.is_empty())
      return reduced_x;

    IntervalVector inter(2, Interval::EMPTY_SET);

    vector<ThickEdge> v_edges = edges();
    for(const auto& edge : v_edges)
      inter |= edge & reduced_x;

    vector<ThickPoint> v_x_vertices;
    ThickPoint::push(reduced_x, v_x_vertices);

    for(const auto& vertex : v_x_vertices)
      if(contains(vertex) != NO)
        inter |= vertex.box();

    return inter;
  }

  const ConvexPolygon& ConvexPolygon::operator&=(const ConvexPolygon& x)
  {
    // todo: optimize this?
    *this = *this & x;
    return *this;
  }

  const ConvexPolygon& ConvexPolygon::operator|=(const ConvexPolygon& x)
  {
    // todo: optimize this?
    *this = *this | x;
    return *this;
  }
}