/** 
 *  \file codac2_ConvexPolygon.h
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include "codac2_Polygon.h"

namespace codac2
{
  class ConvexPolygon : public Polygon
  {
    public:

      ConvexPolygon(std::initializer_list<Vector> vertices);
      ConvexPolygon(const std::vector<Vector>& vertices);
      explicit ConvexPolygon(const std::vector<IntervalVector>& vertices);
      ConvexPolygon(std::initializer_list<Edge> edges);
      ConvexPolygon(const std::vector<Edge>& edges);
      explicit ConvexPolygon(const IntervalVector& x);

      /**
       * \brief Provides an empty convex polygon
       * 
       * \return an empty set
       */
      static ConvexPolygon empty();

    protected:

      ConvexPolygon();
  };

  // The intersection of two convex polygons is a convex polygon
  ConvexPolygon operator&(const ConvexPolygon& p1, const ConvexPolygon& p2);
}