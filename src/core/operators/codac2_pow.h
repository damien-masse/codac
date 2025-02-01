/** 
 *  \file codac2_pow.h
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
  struct PowOp
  {
    static Interval fwd(const Interval& x1, const Interval& x2);
    static ScalarType fwd_natural(const ScalarType& x1, const ScalarType& x2);
    static ScalarType fwd_centered(const ScalarType& x1, const ScalarType& x2);
    static void bwd(const Interval& y, Interval& x1, Interval& x2);
  };

  // Analytic operator

  inline ScalarExpr
  pow(const ScalarExpr& x1, const ScalarExpr& x2)
  {
    return { std::make_shared<AnalyticOperationExpr<PowOp,ScalarType,ScalarType,ScalarType>>(x1,x2) };
  }

  // Inline functions

  inline Interval PowOp::fwd(const Interval& x1, const Interval& x2)
  {
    return pow(x1,x2);
  }

  inline ScalarType PowOp::fwd_natural(const ScalarType& x1, const ScalarType& x2)
  {
    return {
      fwd(x1.a, x2.a),
      x1.def_domain && x2.def_domain
    };
  }

  inline ScalarType PowOp::fwd_centered(const ScalarType& x1, const ScalarType& x2)
  {
    IntervalMatrix d(1,x1.da.size());
    for(Index i = 0 ; i < d.size() ; i++)
      d(0,i) = x2.a*x1.da(0,i)*pow(x1.a,x2.a-1.);

    return {
      fwd(x1.m, x2.m),
      fwd(x1.a, x2.a),
      d,
      x1.def_domain && x2.def_domain
    };
  }

  inline void PowOp::bwd(const Interval& y, Interval& x1, Interval& x2)
  {
    bwd_pow(y, x1, x2);
  }
}