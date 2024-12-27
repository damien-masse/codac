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
#include "codac2_ValueType.h"
#include "codac2_analytic_variables.h"

namespace codac2
{
  template<typename T>
  struct AnalyticExprWrapper;

  template<typename T>
  AnalyticExprWrapper<typename ValueType<T>::Type> const_value(const T& x);

  template<typename T>
  struct AnalyticExprWrapper : public std::shared_ptr<AnalyticExpr<T>>
  {
    AnalyticExprWrapper(const std::shared_ptr<AnalyticExpr<T>>& e)
      : std::shared_ptr<AnalyticExpr<T>>(e)
    { }
    
    AnalyticExprWrapper(const ScalarVar& e)
      requires std::is_same_v<T,ScalarType>
      : std::shared_ptr<AnalyticExpr<T>>({ std::dynamic_pointer_cast<AnalyticExpr<T>>(e.copy()) })
    { }
    
    AnalyticExprWrapper(const VectorVar& e)
      requires std::is_same_v<T,VectorType>
      : std::shared_ptr<AnalyticExpr<T>>({ std::dynamic_pointer_cast<AnalyticExpr<T>>(e.copy()) })
    { }
    
    template<typename V>
    AnalyticExprWrapper(const V& e)
      requires std::is_same_v<typename ValueType<V>::Type,T>
      : std::shared_ptr<AnalyticExpr<T>>(const_value(e))
    { }
    
    inline AnalyticExprWrapper<ScalarType> operator[](Index i) const
      requires std::is_same_v<T,VectorType>
    {
      return { std::make_shared<AnalyticOperationExpr<ComponentOp,ScalarType,VectorType>>(*this,i) };
    }
  };

  using ScalarExpr = AnalyticExprWrapper<ScalarType>;
  using VectorExpr = AnalyticExprWrapper<VectorType>;
  using MatrixExpr = AnalyticExprWrapper<MatrixType>;

  template<class X>
  concept IsScalarExprOrVar = (std::is_base_of_v<VarBase,X> || std::is_base_of_v<ScalarExpr,X>);

  template<class X>
  concept IsVectorExprOrVar = (std::is_base_of_v<VarBase,X> || std::is_base_of_v<VectorExpr,X>);

  template<class X>
  concept IsMatrixExprOrVar = (std::is_base_of_v<VarBase,X> || std::is_base_of_v<MatrixExpr,X>);
}