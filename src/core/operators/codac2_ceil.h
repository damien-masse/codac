/** 
 *  \file codac2_ceil.h
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include "codac2_Interval.h"
#include "codac2_AnalyticType.h"
#include "codac2_AnalyticExprWrapper.h"

namespace codac2
{
  struct CeilOp
  {
    static Interval fwd(const Interval& x1);
    static void bwd(const Interval& y, Interval& x1);
  };

  // Inline functions

  inline Interval CeilOp::fwd(const Interval& x1)
  {
    return ceil(x1);
  }

  inline void CeilOp::bwd(const Interval& y, Interval& x1)
  {
    // The content of this function comes from the IBEX library.
    // See ibex::Interval (IBEX lib, main author: Gilles Chabert)
    //   https://ibex-lib.readthedocs.io

    if(y.is_empty())
      x1.set_empty();

    else
    {
      double r = std::floor(y.ub());
      double l = std::ceil(y.lb());

      if(r < l)
        x1.set_empty();

      else
        x1 &= Interval(l,r) + Interval(-1,0);
    }
  }
}