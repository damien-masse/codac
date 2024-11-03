/** 
 *  Arithmetic operations on convex polygons
 * ----------------------------------------------------------------------------
 *  \date       2020
 *  \author     Simon Rohou
 *  \bugfixes   Bertrand Neveu
 *  \copyright  Copyright 2021 Codac Team
 *  \license    This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 */

#include "codac_polygon_arithmetic.h"

using namespace std;
using namespace ibex;

namespace codac
{
  const ConvexPolygon operator+(const ConvexPolygon& x)
  {
    return x;
  }

  const ConvexPolygon operator+(const ConvexPolygon& x, const IntervalVector& v)
  {
    return operator+(v, x);
  }

  const ConvexPolygon operator+(const IntervalVector& v, const ConvexPolygon& x)
  {
    assert(x.size() == v.size());

    vector<ThickPoint> v_result_thick_pts = ThickPoint::to_ThickPoints(x.vertices());
    for(auto& pt : v_result_thick_pts)
      pt = ThickPoint(pt.box() + v);
    // ^ The operation may transform a degenerate point-box into a large box

    return ConvexPolygon(v_result_thick_pts);
  }

  const ConvexPolygon operator-(const ConvexPolygon& x)
  {
    vector<Vector> v_result_pts = x.vertices();
    for(auto& pt : v_result_pts)
      pt *= -1.;
    return ConvexPolygon(v_result_pts, true);
  }

  const ConvexPolygon operator-(const ConvexPolygon& x, const IntervalVector& v)
  {
    return operator+(x, -v);
  }

  const ConvexPolygon operator-(const IntervalVector& v, const ConvexPolygon& x)
  {
    return operator+(-x, v);
  }

  const ConvexPolygon operator*(const IntervalMatrix& m, const ConvexPolygon& x)
  {
    assert(x.size() == m.nb_cols() && x.size() == m.nb_rows());

    vector<ThickPoint> v_result_thick_pts = ThickPoint::to_ThickPoints(x.vertices());
    for(auto& pt : v_result_thick_pts)
      pt = ThickPoint(m * pt.box());

    return ConvexPolygon(v_result_thick_pts);
  }
  
  vector<ThickPoint> inter_thickpoints(const ConvexPolygon& p1, const ConvexPolygon& p2)
  {
    vector<ThickPoint> v_pts;

    // Add all vertices of p1 that are inside p2
    for(const auto& pt_ : p1.vertices())
    {
      ThickPoint pt(pt_);
      if(p2.contains(pt) != NO)
        v_pts.push_back(pt);
    }

    // Add all vertices of p2 that are inside p1
    for(const auto& pt_ : p2.vertices())
    {
      ThickPoint pt(pt_);
      if(p1.contains(pt) != NO)
        v_pts.push_back(pt);
    }

    // Add all intersection points
    for(const auto& e1 : p1.edges())
      for(const auto& e2 : p2.edges())
      {
        const ThickPoint intersection_pt = e1 & e2;

        if(!intersection_pt.does_not_exist())
        {
          // If edges are possibly parallel:
          if(ThickEdge::parallel(e1, e2) != NO)
          {
            if(e1.contains(e2.p1()) != NO)
              v_pts.push_back(e2.p1());

            if(e1.contains(e2.p2()) != NO)
              v_pts.push_back(e2.p2());

            if(e2.contains(e1.p1()) != NO)
              v_pts.push_back(e1.p1());

            if(e2.contains(e1.p2()) != NO)
              v_pts.push_back(e1.p2());
          }

          else
            v_pts.push_back(intersection_pt);
        }
      }

    return v_pts;
  }

  const ConvexPolygon operator&(const ConvexPolygon& p1, const ConvexPolygon& p2)
  {
    return ConvexPolygon(inter_thickpoints(p1,p2));
  }

  const ConvexPolygon operator&(const IntervalVector& p1, const ConvexPolygon& p2)
  {
    assert(p1.size() == 2 && "other dimensions not supported");
    return operator&(p2, ConvexPolygon(p1));
  }

  const ConvexPolygon operator&(const ConvexPolygon& p1, const IntervalVector& p2)
  {
    assert(p2.size() == 2 && "other dimensions not supported");
    return operator&(p1, ConvexPolygon(p2));
  }
  
  const ConvexPolygon operator|(const ConvexPolygon& p1, const ConvexPolygon& p2)
  {
    vector<ThickPoint> v_pts;
    for(const auto& pt_ : p1.vertices())
      v_pts.push_back(ThickPoint(pt_));
    for(const auto& pt_ : p2.vertices())
      v_pts.push_back(ThickPoint(pt_));
    return ConvexPolygon(v_pts);
  }

  const ConvexPolygon operator|(const IntervalVector& p1, const ConvexPolygon& p2)
  {
    assert(p1.size() == 2 && "other dimensions not supported");
    return operator|(p2, ConvexPolygon(p1));
  }

  const ConvexPolygon operator|(const ConvexPolygon& p1, const IntervalVector& p2)
  {
    assert(p2.size() == 2 && "other dimensions not supported");
    return operator|(p1, ConvexPolygon(p2));
  }
}