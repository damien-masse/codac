/** 
 *  \file codac2_geometry.h
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include "codac2_BoolInterval.h"
#include "codac2_IntervalVector.h"

namespace codac2
{
  enum class OrientationInterval
  {
    EMPTY = 0x00,
    COLINEAR = 0x01,
    CLOCKWISE = 0x02,
    COUNTERCLOCKWISE = 0x04,
    UNKNOWN = 0x01 | 0x02 | 0x04
  };

  constexpr OrientationInterval operator|(OrientationInterval a, OrientationInterval b)
  { return static_cast<OrientationInterval>(static_cast<int>(a) | static_cast<int>(b)); }

  /**
   * \brief Streams out a OrientationInterval
   * 
   * \param os the stream to be updated
   * \param x the orientation interval to stream out
   * \return a reference to the updated stream
   */
  inline std::ostream& operator<<(std::ostream& os, const OrientationInterval& x)
  {
    if(x == OrientationInterval::EMPTY)
      os << "[ empty ]";
    else if(x == OrientationInterval::COLINEAR)
      os << "[ col ]";
    else if(x == OrientationInterval::CLOCKWISE)
      os << "[ cw ]";
    else if(x == OrientationInterval::COUNTERCLOCKWISE)
      os << "[ ccw ]";
    else if(x == (OrientationInterval::COLINEAR | OrientationInterval::CLOCKWISE))
      os << "[ col, cw ]";
    else if(x == (OrientationInterval::COLINEAR | OrientationInterval::COUNTERCLOCKWISE))
      os << "[ col, ccw ]";
    else if(x == OrientationInterval::UNKNOWN)
      os << "[ col, cw, ccw ]";
    return os;
  }

  /**
   * \brief Computes the orientation of an ordered triplet of 2D points.
   *
   * Determines whether the sequence of points (``p1``,``p2``,``p3``) makes a
   * counterclockwise turn, a clockwise turn, or if the points are colinear.
   * Depending on floating point uncertainties, the test may not be able to conclude
   * (a ``UNKNOWN`` value would then be returned).
   *
   * \param p1 First point (2d ``Vector``) of the triplet.
   * \param p2 Second point (2d ``Vector``) of the triplet.
   * \param p3 Third point (2d ``Vector``) of the triplet.
   * 
   * \return an orientation of type ``OrientationInterval``
   */
  OrientationInterval orientation(const IntervalVector& p1, const IntervalVector& p2, const IntervalVector& p3);

  /**
   * \brief Checks whether three 2D points are aligned (colinear).
   *
   * Determines if the points lie on the same straight line using an 
   * orientation test (cross product). Depending on floating point uncertainties,
   * the test may not be able to conclude (a ``UNKNOWN`` value would then be returned).
   *
   * \param p1 First point (2d ``Vector``) of the triplet.
   * \param p2 Second point (2d ``Vector``) of the triplet.
   * \param p3 Third point (2d ``Vector``) of the triplet.
   * 
   * \return a ``BooleanInterval``
   */
  BoolInterval aligned(const IntervalVector& p1, const IntervalVector& p2, const IntervalVector& p3);
}