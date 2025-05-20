/** 
 *  \file codac2_extend.h
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Simon Rohou, Damien Massé
 *  \copyright  Copyright 2025 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include "codac2_Interval.h"
#include "codac2_AnalyticType.h"
#include "codac2_AnalyticExprWrapper.h"

namespace codac2
{
  // extend(f(x),g(x)) = f(x) if x \in D_f, otherwise g(x)
  struct ExtendOp
  {
    template<typename X1,typename X2>
    static std::string str(const X1& x1, const X2& x2)
    {
      return "extend(" + x1->str() + "," + x2->str() + ")";
    }

    template<typename X1, typename X2>
    static std::pair<Index,Index> output_shape(const X1& s1, const X2& s2)
    {
      auto shape1 = s1->output_shape();
      auto shape2 = s2->output_shape();
      assert(shape1 == shape2);
      return shape1;
    }

    /* ScalarType version */
    static Interval fwd(const Interval& x1, const Interval& x2);
    static ScalarType fwd_natural(const ScalarType& x1, const ScalarType& x2);
    static ScalarType fwd_centered(const ScalarType& x1, const ScalarType& x2);
    static void bwd(const Interval& y, Interval& x1, Interval& x2);

    /* VectorType version */
    static IntervalVector fwd(const IntervalVector& x1, const IntervalVector& x2);
    static VectorType fwd_natural(const VectorType& x1, const VectorType& x2);
    static VectorType fwd_centered(const VectorType& x1, const VectorType& x2);
    static void bwd(const IntervalVector& y, IntervalVector& x1, IntervalVector& x2);

    /* MatrixType version */
    static IntervalMatrix fwd(const IntervalMatrix& x1, const IntervalMatrix& x2);
    static MatrixType fwd_natural(const MatrixType& x1, const MatrixType& x2);
    static MatrixType fwd_centered(const MatrixType& x1, const MatrixType& x2);
    static void bwd(const IntervalMatrix& y, IntervalMatrix& x1, IntervalMatrix& x2);
  };

  // Analytic operator
  // The following function can be used to build analytic expressions.

  inline ScalarExpr
  extend(const ScalarExpr& x1, const ScalarExpr& x2)
  {
    return { std::make_shared<AnalyticOperationExpr<ExtendOp,ScalarType,ScalarType,ScalarType>>(x1,x2) };
  }
  
  inline VectorExpr
  extend(const VectorExpr& x1, const VectorExpr& x2)
  {
    return { std::make_shared<AnalyticOperationExpr<ExtendOp,VectorType,VectorType,VectorType>>(x1,x2) };
  }
  
  inline MatrixExpr
  extend(const MatrixExpr& x1, const MatrixExpr& x2)
  {
    return { std::make_shared<AnalyticOperationExpr<ExtendOp,MatrixType,MatrixType,MatrixType>>(x1,x2) };
  }
  
  // Inline functions
  // 1 ScalarType

  inline Interval ExtendOp::fwd(const Interval& x1, const Interval& x2)
  {
    return x1 | x2;
  }

  inline ScalarType ExtendOp::fwd_natural(const ScalarType& x1, const ScalarType& x2)
  {
    
    return {
      (x1.def_domain ? x1.a : fwd(x1.a,x2.a)),
      x1.def_domain || x2.def_domain
    };
  }

  inline ScalarType ExtendOp::fwd_centered(const ScalarType& x1, const ScalarType& x2)
  {
    return {
      (x1.def_domain ? x1.m : fwd(x1.m,x2.m)),
      (x1.def_domain ? x1.a : fwd(x1.a,x2.a)),
      (x1.def_domain ? x1.da : (x1.da | x2.da)),
      x1.def_domain || x2.def_domain
    };
  }

  inline void ExtendOp::bwd(const Interval& y, Interval& x1, Interval& x2)
  {
    // any value in x1 must be in y
    x1 &= y;
    if (x1.is_empty()) 
    {
       // only if x1 is empty, values of x2 must be in y
       x2 &= y;
    }
  }

  // 2 VectorType

  inline IntervalVector ExtendOp::fwd(const IntervalVector& x1, const IntervalVector& x2)
  {
    return x1 | x2;
  }

  inline VectorType ExtendOp::fwd_natural(const VectorType& x1, const VectorType& x2)
  {
    
    return {
      (x1.def_domain ? x1.a : fwd(x1.a,x2.a)),
      x1.def_domain || x2.def_domain
    };
  }

  inline VectorType ExtendOp::fwd_centered(const VectorType& x1, const VectorType& x2)
  {
    return {
      (x1.def_domain ? x1.m : fwd(x1.m,x2.m)),
      (x1.def_domain ? x1.a : fwd(x1.a,x2.a)),
      (x1.def_domain ? x1.da : (x1.da | x2.da)),
      x1.def_domain || x2.def_domain
    };
  }

  inline void ExtendOp::bwd(const IntervalVector& y, IntervalVector& x1, IntervalVector& x2)
  {
    // any value in x1 must be in y
    x1 &= y;
    if (x1.is_empty()) 
    {
       // only if x1 is empty, values of x2 must be in y
       x2 &= y;
    }
  }

  // 2 MatrixType

  inline IntervalMatrix ExtendOp::fwd(const IntervalMatrix& x1, const IntervalMatrix& x2)
  {
    return x1 | x2;
  }

  inline MatrixType ExtendOp::fwd_natural(const MatrixType& x1, const MatrixType& x2)
  {
    
    return {
      (x1.def_domain ? x1.a : fwd(x1.a,x2.a)),
      x1.def_domain || x2.def_domain
    };
  }

  inline MatrixType ExtendOp::fwd_centered(const MatrixType& x1, const MatrixType& x2)
  {
    return {
      (x1.def_domain ? x1.m : fwd(x1.m,x2.m)),
      (x1.def_domain ? x1.a : fwd(x1.a,x2.a)),
      (x1.def_domain ? x1.da : (x1.da | x2.da)),
      x1.def_domain || x2.def_domain
    };
  }

  inline void ExtendOp::bwd(const IntervalMatrix& y, IntervalMatrix& x1, IntervalMatrix& x2)
  {
    // any value in x1 must be in y
    x1 &= y;
    if (x1.is_empty()) 
    {
       // only if x1 is empty, values of x2 must be in y
       x2 &= y;
    }
  }

}
