/** 
 *  \file codac2_min.h
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
  struct MinOp
  {
    static Interval fwd(const Interval& x1, const Interval& x2);
    static void bwd(const Interval& y, Interval& x1, Interval& x2);
  };

  // Inline functions

  inline Interval MinOp::fwd(const Interval& x1, const Interval& x2)
  {
    return min(x1,x2);
  }

  inline void MinOp::bwd(const Interval& y, Interval& x1, Interval& x2)
  {
    // The content of this function comes from the IBEX library.
    // See ibex::Interval (IBEX lib, main author: Gilles Chabert)
    //   https://ibex-lib.readthedocs.io

    Interval mx1 = -x1;
    Interval mx2 = -x2;

    MaxOp::bwd(-y,mx1,mx2);

    if(mx1.is_empty() || mx2.is_empty())
    {
      x1.set_empty();
      x2.set_empty();
    }

    else
    {
      x1 = -mx1;
      x2 = -mx2;
    }
  }
}