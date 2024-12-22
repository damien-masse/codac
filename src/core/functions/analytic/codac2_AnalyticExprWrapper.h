/** 
 *  \file codac2_AnalyticExprWrapper.h
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include "codac2_AnalyticExpr.h"

namespace codac2
{
  template<typename T>
  struct AnalyticExprWrapper : public std::shared_ptr<AnalyticExpr<T>>
  {
    AnalyticExprWrapper(const std::shared_ptr<AnalyticExpr<T>>& e)
      : std::shared_ptr<AnalyticExpr<T>>(e)
    { }

    explicit AnalyticExprWrapper(const typename T::Domain& e)
      : std::shared_ptr<AnalyticExpr<T>>(const_value(e))
    { }

    template<typename T_=T>
      requires std::is_same_v<T_,VectorOpValue>
    inline AnalyticExprWrapper<ScalarOpValue> operator[](Index i) const
    {
      return { std::make_shared<AnalyticOperationExpr<ComponentOp,ScalarOpValue,VectorOpValue>>(*this,i) };
    }
  };

  using ScalarExpr = AnalyticExprWrapper<ScalarOpValue>;
  using VectorExpr = AnalyticExprWrapper<VectorOpValue>;
  using MatrixExpr = AnalyticExprWrapper<MatrixOpValue>;

  template<class X>
  concept IsScalarExprOrVar = (std::is_base_of_v<VarBase,X> || std::is_base_of_v<ScalarExpr,X>);

  template<class X>
  concept IsVectorExprOrVar = (std::is_base_of_v<VarBase,X> || std::is_base_of_v<VectorExpr,X>);

  template<class X>
  concept IsMatrixExprOrVar = (std::is_base_of_v<VarBase,X> || std::is_base_of_v<MatrixExpr,X>);
}