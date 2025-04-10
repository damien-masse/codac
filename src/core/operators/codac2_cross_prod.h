/** 
 *  \file codac2_cross_prod.h
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou, Gilles Chabert
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include "codac2_Interval.h"
#include "codac2_AnalyticType.h"
#include "codac2_AnalyticExprWrapper.h"

namespace codac2
{
  struct CrossProdOp
  {
    static IntervalVector fwd(const IntervalVector& x1, const IntervalVector& x2);
    static VectorType fwd_natural(const VectorType& x1, const VectorType& x2);
    static VectorType fwd_centered(const VectorType& x1, const VectorType& x2);
    static void bwd(const IntervalVector& y, IntervalVector& x1, IntervalVector& x2);
  };

  // Analytic operator
  // The following function can be used to build analytic expressions.

  inline VectorExpr
  cross_prod(const VectorExpr& x1, const VectorExpr& x2)
  {
    return { std::make_shared<AnalyticOperationExpr<CrossProdOp,VectorType,VectorType,VectorType>>(x1,x2) };
  }
  
  // Inline functions

  inline IntervalVector CrossProdOp::fwd(const IntervalVector& x1, const IntervalVector& x2)
  {
    return {
      x1[1]*x2[2] - x1[2]*x2[1],
      x1[2]*x2[0] - x1[0]*x2[2],
      x1[0]*x2[1] - x1[1]*x2[0]
    };
  }

  inline VectorType CrossProdOp::fwd_natural(const VectorType& x1, const VectorType& x2)
  {
    return {
      fwd(x1.a,x2.a),
      x1.def_domain && x2.def_domain
    };
  }

  inline VectorType CrossProdOp::fwd_centered(const VectorType& x1, const VectorType& x2)
  {
    return {
      fwd(x1.m,x2.m),
      fwd(x1.a,x2.a),
      IntervalMatrix(0,0), // not supported yet for auto diff
      x1.def_domain && x2.def_domain
    };
  }

  inline void CrossProdOp::bwd([[maybe_unused]] const IntervalVector& y, [[maybe_unused]] IntervalVector& x1, [[maybe_unused]] IntervalVector& x2)
  {
    // todo
  }
}