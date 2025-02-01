/** 
 *  \file codac2_atan.h
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
  struct AtanOp
  {
    static Interval fwd(const Interval& x1);
    static ScalarType fwd_natural(const ScalarType& x1);
    static ScalarType fwd_centered(const ScalarType& x1);
    static void bwd(const Interval& y, Interval& x1);
  };

  // Analytic operator

  inline ScalarExpr
  atan(const ScalarExpr& x1)
  {
    return { std::make_shared<AnalyticOperationExpr<AtanOp,ScalarType,ScalarType>>(x1) };
  }
  
  // Inline functions

  inline Interval AtanOp::fwd(const Interval& x1)
  {
    return atan(x1);
  }

  inline ScalarType AtanOp::fwd_natural(const ScalarType& x1)
  {
    return {
      fwd(x1.a),
      x1.def_domain
    };
  }

  inline ScalarType AtanOp::fwd_centered(const ScalarType& x1)
  {
    IntervalMatrix d(1,x1.da.size());
    for(Index i = 0 ; i < d.size() ; i++)
      d(0,i) = x1.da(0,i)/(1.+sqr(x1.a));

    return {
      fwd(x1.m),
      fwd(x1.a),
      d,
      x1.def_domain
    };
  }

  inline void AtanOp::bwd(const Interval& y, Interval& x1)
  {
    bwd_atan(y, x1);
  }
}