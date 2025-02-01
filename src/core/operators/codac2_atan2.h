/** 
 *  \file codac2_atan2.h
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
  struct Atan2Op
  {
    static Interval fwd(const Interval& x1, const Interval& x2);
    static ScalarType fwd_natural(const ScalarType& x1, const ScalarType& x2);
    static ScalarType fwd_centered(const ScalarType& x1, const ScalarType& x2);
    static void bwd(const Interval& y, Interval& x1, Interval& x2);
  };

  // Analytic operator

  inline ScalarExpr
  atan2(const ScalarExpr& x1, const ScalarExpr& x2)
  {
    return { std::make_shared<AnalyticOperationExpr<Atan2Op,ScalarType,ScalarType,ScalarType>>(x1,x2) };
  }
  
  // Inline functions

  inline Interval Atan2Op::fwd(const Interval& x1, const Interval& x2)
  {
    return atan2(x1,x2);
  }

  inline ScalarType Atan2Op::fwd_natural(const ScalarType& x1, const ScalarType& x2)
  {
    return {
      fwd(x1.a, x2.a),
      x1.def_domain && x2.def_domain
      && !(x1.a == 0. && x2.a == 0.) /* def domain of the derivative of atan2 */
    };
  }

  inline ScalarType Atan2Op::fwd_centered(const ScalarType& x1, const ScalarType& x2)
  {
    assert(x1.da.rows() == 1);
    assert(x1.da.rows() == x2.da.rows() && x1.da.cols() == x2.da.cols());

    IntervalMatrix d(1,x1.da.cols());
    for(Index i = 0 ; i < d.size() ; i++)
      d(0,i) = (-x1.a*x2.da(0,i)/(sqr(x2.a)+sqr(x1.a)))+(x2.a*x1.da(0,i)/(sqr(x2.a)+sqr(x1.a)));

    return {
      fwd(x1.m, x2.m),
      fwd(x1.a, x2.a),
      d,
      x1.def_domain && x2.def_domain
      && !(x1.a == 0. && x2.a == 0.) /* def domain of the derivative of atan2 */
    };
  }

  inline void Atan2Op::bwd(const Interval& y, Interval& x1, Interval& x2)
  {
    bwd_atan2(y, x1, x2);
  }
}