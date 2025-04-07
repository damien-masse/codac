/** 
 *  \file codac2_Figure2DInterface.h
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include <vector>
#include <string>
#include "codac2_Vector.h"
#include "codac2_IntervalVector.h"
#include "codac2_StyleProperties.h"
#include "codac2_Ellipsoid.h"
#include "codac2_SampledTraj.h"
#include "codac2_AnalyticTraj.h"

namespace codac2
{
  /*
   * This interface contains the “low-level” methods that need to be overloaded
   * to implement graphic primitives in the output view. The Figure2D class may
   * contain methods not listed here: these methods are more “high-level” and do
   * not require a specific implementation related to the graphic tool.
   * For example, Figure2D::draw_line calls Figure2D::draw_polyline (a line is a
   * special case of a two-point polyline), so only draw_polyline needs to be
   * implemented in this interface.
  */
  class Figure2DInterface
  {
    public:

      // Geometric shapes
      virtual void draw_point(const Vector& c, const StyleProperties& s = StyleProperties()) = 0;
      virtual void draw_box(const IntervalVector& x, const StyleProperties& s = StyleProperties()) = 0;
      virtual void draw_circle(const Vector& c, double r, const StyleProperties& s = StyleProperties()) = 0;
      virtual void draw_ring(const Vector& c, const Interval& r, const StyleProperties& s = StyleProperties()) = 0;
      virtual void draw_polyline(const std::vector<Vector>& x, float tip_length, const StyleProperties& s = StyleProperties()) = 0;
      virtual void draw_polygone(const std::vector<Vector>& x, const StyleProperties& s = StyleProperties()) = 0;
      virtual void draw_pie(const Vector& c, const Interval& r, const Interval& theta, const StyleProperties& s = StyleProperties()) = 0;
      virtual void draw_ellipse(const Vector& c, const Vector& ab, double theta, const StyleProperties& s = StyleProperties()) = 0;

      // Robots
      virtual void draw_tank(const Vector& x, float size, const StyleProperties& s = StyleProperties()) = 0;
      virtual void draw_AUV(const Vector& x, float size, const StyleProperties& s = StyleProperties()) = 0;
      virtual void draw_motor_boat(const Vector& x, float size, const StyleProperties& s = StyleProperties()) = 0;
  };
}