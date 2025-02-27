/** 
 *  \file codac2_cart_prod.h
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2023 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include "codac2_IntervalVector.h"

namespace codac2
{
  inline IntervalVector to_IntervalVector(const Interval& x)
  {
    IntervalVector a(1);
    a[0] = x;
    return a;
  }

  inline IntervalVector to_IntervalVector(const IntervalVector& x)
  {
    return x;
  }

  inline IntervalVector to_IntervalVector(const Vector& x)
  {
    return x.template cast<Interval>();
  }

  inline IntervalVector cart_prod()
  {
    return IntervalVector(0);
  }

  template<typename... X>
    requires ((is_interval_based_v<X>) || ...)
      && ((!is_ctc_v<X>) && ...) && ((!is_sep_v<X>) && ...)
  inline IntervalVector cart_prod(const X&... x)
  {
    std::vector<IntervalVector> v_x;
    ((v_x.push_back(to_IntervalVector(x))), ...);

    Index n = 0;
    for(const auto& xi : v_x)
      n += xi.size();
    IntervalVector x_(n);

    Index i = 0;
    for(const auto& xi : v_x)
    {
      x_.put(i, xi);
      i += xi.size();
    }
    return x_;
  }

  template<typename... X>
    requires ((!is_interval_based_v<X>) && ...)
      && ((!is_ctc_v<X>) && ...) && ((!is_sep_v<X>) && ...)
  inline Vector cart_prod(const X&... x)
  {
    return cart_prod(to_IntervalVector(x)...).mid();
  }
}