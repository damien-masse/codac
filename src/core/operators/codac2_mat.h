/** 
 *  \file codac2_matrix.h
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include "codac2_IntervalVector.h"
#include "codac2_IntervalMatrix.h"
#include "codac2_AnalyticType.h"
#include "codac2_AnalyticExprWrapper.h"

namespace codac2
{
  struct MatrixOp
  {
    static inline void set_col_i(IntervalMatrix& m, const IntervalVector& x, Index i)
    {
      assert(i >= 0 && i < m.cols());
      assert_release(x.size() == m.rows());
      m.col(i) = x;
    }

    template<typename X1,typename... X>
      requires ((X1::ColsAtCompileTime == 1) && ((X::ColsAtCompileTime == 1) && ...))
    static inline IntervalMatrix fwd(const X1& x1, const X&... x)
    {
      IntervalMatrix m(x1.size(), 1+sizeof...(X));
      Index i = 0;
      MatrixOp::set_col_i(m, x1, i++); // first column x1
      (MatrixOp::set_col_i(m, x, i++), ...); // other columns
      return m;
    }

    template<typename... X>
      requires (std::is_base_of_v<VectorType,X> && ...)
    static inline MatrixType fwd_natural(const X&... x)
    {
      return {
        MatrixOp::fwd(x.a...),
        ((true && x.def_domain), ...)
      };
    }

    template<typename... X>
      requires (std::is_base_of_v<VectorType,X> && ...)
    static inline MatrixType fwd_centered(const X&... x)
    {
      throw std::runtime_error("MatrixOp not fully implemented yet");
      return {
        MatrixOp::fwd(x.m...),
        MatrixOp::fwd(x.a...),
        IntervalMatrix(0,0), // not supported yet for matrices
        ((true && x.def_domain), ...)
      };
    }

    template<typename... X>
      requires (std::is_base_of_v<IntervalVector,X> && ...)
    static inline void bwd(const IntervalMatrix& y, X&... x)
    {
      throw std::runtime_error("MatrixOp not fully implemented yet");
      Index i = 0;
      ((x &= y.col(i++)), ...);
    }
  };

  // Analytic operator
  // The following function can be used to build analytic expressions.

  // Generic variadic case, cannot handle const values (int, double) for now

  template<typename... X>
  inline MatrixExpr
  mat(const std::shared_ptr<AnalyticExpr<X>>&... x)
  {
    return { std::make_shared<AnalyticOperationExpr<MatrixOp,MatrixType,X...>>(
      AnalyticOperationExpr<MatrixOp,MatrixType,X...>(x...)) };
  }
}