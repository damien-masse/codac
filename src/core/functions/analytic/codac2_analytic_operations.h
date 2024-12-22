/** 
 *  \file codac2_analytic_operations.h
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include <memory>
#include "codac2_AnalyticExpr.h"
#include "codac2_analytic_variables.h"
#include "codac2_analytic_constants.h"
#include "codac2_AnalyticExprWrapper.h"

namespace codac2
{
  using ScalarExpr = AnalyticExprWrapper<ScalarOpValue>;
  using VectorExpr = AnalyticExprWrapper<VectorOpValue>;
  using MatrixExpr = AnalyticExprWrapper<MatrixOpValue>;

  template<class X>
  concept IsScalarExprOrVar = (std::is_base_of_v<VarBase,X> || std::is_base_of_v<ScalarExpr,X>);

  template<class X>
  concept IsVectorExprOrVar = (std::is_base_of_v<VarBase,X> || std::is_base_of_v<VectorExpr,X>);

  template<class X>
  concept IsMatrixExprOrVar = (std::is_base_of_v<VarBase,X> || std::is_base_of_v<MatrixExpr,X>);

  // operator+ (unary case)
    
    inline ScalarExpr
    operator+(const ScalarExpr& x1)
    {
      return { std::make_shared<AnalyticOperationExpr<AddOp,ScalarOpValue,ScalarOpValue>>(x1) };
    }

    inline VectorExpr
    operator+(const VectorExpr& x1)
    {
      return { std::make_shared<AnalyticOperationExpr<AddOp,VectorOpValue,VectorOpValue>>(x1) };
    }

    inline MatrixExpr
    operator+(const MatrixExpr& x1)
    {
      return { std::make_shared<AnalyticOperationExpr<AddOp,MatrixOpValue,MatrixOpValue>>(x1) };
    }

  // operator+
    
    inline ScalarExpr
    operator+(ScalarExpr x1, ScalarExpr x2)
    {
      return { std::make_shared<AnalyticOperationExpr<AddOp,ScalarOpValue,ScalarOpValue,ScalarOpValue>>(x1,x2) };
    }

    template<typename X1>
      requires (!IsScalarExprOrVar<X1>)
    inline ScalarExpr
    operator+(const X1& x1, ScalarExpr x2)
    {
      return operator+(const_value(x1),x2);
    }

    template<typename X2>
      requires (!IsScalarExprOrVar<X2>)
    inline ScalarExpr
    operator+(ScalarExpr x1, const X2& x2)
    {
      return operator+(x1, const_value(x2));
    }

    inline VectorExpr
    operator+(VectorExpr x1, VectorExpr x2)
    {
      return { std::make_shared<AnalyticOperationExpr<AddOp,VectorOpValue,VectorOpValue,VectorOpValue>>(x1,x2) };
    }

    template<typename X1>
      requires (!IsVectorExprOrVar<X1>)
    inline VectorExpr
    operator+(const X1& x1, VectorExpr x2)
    {
      return operator+(const_value(x1),x2);
    }

    template<typename X2>
      requires (!IsVectorExprOrVar<X2>)
    inline VectorExpr
    operator+(VectorExpr x1, const X2& x2)
    {
      return operator+(x1, const_value(x2));
    }

    inline MatrixExpr
    operator+(MatrixExpr x1, MatrixExpr x2)
    {
      return { std::make_shared<AnalyticOperationExpr<AddOp,MatrixOpValue,MatrixOpValue,MatrixOpValue>>(x1,x2) };
    }

    template<typename X1>
      requires (!IsMatrixExprOrVar<X1>)
    inline MatrixExpr
    operator+(const X1& x1, MatrixExpr x2)
    {
      return operator+(const_value(x1),x2);
    }

    template<typename X2>
      requires (!IsMatrixExprOrVar<X2>)
    inline MatrixExpr
    operator+(MatrixExpr x1, const X2& x2)
    {
      return operator+(x1, const_value(x2));
    }

  // operator- (unary case)
    
    inline ScalarExpr
    operator-(const ScalarExpr& x1)
    {
      return { std::make_shared<AnalyticOperationExpr<SubOp,ScalarOpValue,ScalarOpValue>>(x1) };
    }

    inline VectorExpr
    operator-(const VectorExpr& x1)
    {
      return { std::make_shared<AnalyticOperationExpr<SubOp,VectorOpValue,VectorOpValue>>(x1) };
    }

    inline MatrixExpr
    operator-(const MatrixExpr& x1)
    {
      return { std::make_shared<AnalyticOperationExpr<SubOp,MatrixOpValue,MatrixOpValue>>(x1) };
    }

  // operator-
    
    inline ScalarExpr
    operator-(const ScalarExpr& x1, const ScalarExpr& x2)
    {
      return { std::make_shared<AnalyticOperationExpr<SubOp,ScalarOpValue,ScalarOpValue,ScalarOpValue>>(x1,x2) };
    }

    template<typename X1>
      requires (!IsScalarExprOrVar<X1>)
    inline ScalarExpr
    operator-(const X1& x1, const ScalarExpr& x2)
    {
      return operator-(const_value(x1),x2);
    }

    template<typename X2>
      requires (!IsScalarExprOrVar<X2>)
    inline ScalarExpr
    operator-(const ScalarExpr& x1, const X2& x2)
    {
      return operator-(x1, const_value(x2));
    }

    inline VectorExpr
    operator-(const VectorExpr& x1, const VectorExpr& x2)
    {
      return { std::make_shared<AnalyticOperationExpr<SubOp,VectorOpValue,VectorOpValue,VectorOpValue>>(x1,x2) };
    }

    template<typename X1>
      requires (!IsVectorExprOrVar<X1>)
    inline VectorExpr
    operator-(const X1& x1, const VectorExpr& x2)
    {
      return operator-(const_value(x1),x2);
    }

    template<typename X2>
      requires (!IsVectorExprOrVar<X2>)
    inline VectorExpr
    operator-(const VectorExpr& x1, const X2& x2)
    {
      return operator-(x1, const_value(x2));
    }

    inline MatrixExpr
    operator-(const MatrixExpr& x1, const MatrixExpr& x2)
    {
      return { std::make_shared<AnalyticOperationExpr<SubOp,MatrixOpValue,MatrixOpValue,MatrixOpValue>>(x1,x2) };
    }

    template<typename X1>
      requires (!IsMatrixExprOrVar<X1>)
    inline MatrixExpr
    operator-(const X1& x1, const MatrixExpr& x2)
    {
      return operator-(const_value(x1),x2);
    }

    template<typename X2>
      requires (!IsMatrixExprOrVar<X2>)
    inline MatrixExpr
    operator-(const MatrixExpr& x1, const X2& x2)
    {
      return operator-(x1, const_value(x2));
    }

  // operator*
    
    inline ScalarExpr
    operator*(const ScalarExpr& x1, const ScalarExpr& x2)
    {
      return { std::make_shared<AnalyticOperationExpr<MulOp,ScalarOpValue,ScalarOpValue,ScalarOpValue>>(x1,x2) };
    }
    
    inline VectorExpr
    operator*(const ScalarExpr& x1, const VectorExpr& x2)
    {
      return { std::make_shared<AnalyticOperationExpr<MulOp,VectorOpValue,ScalarOpValue,VectorOpValue>>(x1,x2) };
    }
    
    inline VectorExpr
    operator*(const VectorExpr& x1, const ScalarExpr& x2)
    {
      return { std::make_shared<AnalyticOperationExpr<MulOp,VectorOpValue,VectorOpValue,ScalarOpValue>>(x1,x2) };
    }

    inline ScalarExpr
    operator*(const Interval& x1, const ScalarExpr& x2)
    {
      return operator*(const_value(x1),x2);
    }

    inline ScalarExpr
    operator*(const ScalarExpr& x1, const Interval& x2)
    {
      return operator*(x1,const_value(x2));
    }

    inline VectorExpr
    operator*(const Interval& x1, const VectorExpr& x2)
    {
      return operator*(const_value(x1),x2);
    }

    inline VectorExpr
    operator*(const VectorExpr& x1, const Interval& x2)
    {
      return operator*(x1,const_value(x2));
    }

    inline VectorExpr
    operator*(const ScalarExpr& x1, const IntervalVector& x2)
    {
      return operator*(x1,const_value(x2));
    }

    inline VectorExpr
    operator*(const IntervalVector& x1, const ScalarExpr& x2)
    {
      return operator*(const_value(x1),x2);
    }

    inline VectorExpr
    operator*(const MatrixExpr& x1, const VectorExpr& x2)
    {
      return { std::make_shared<AnalyticOperationExpr<MulOp,VectorOpValue,MatrixOpValue,VectorOpValue>>(x1,x2) };
    }

    inline VectorExpr
    operator*(const MatrixExpr& x1, const IntervalVector& x2)
    {
      return { std::make_shared<AnalyticOperationExpr<MulOp,VectorOpValue,MatrixOpValue,VectorOpValue>>(x1,const_value(x2)) };
    }

    inline VectorExpr
    operator*(const IntervalMatrix& x1, const VectorExpr& x2)
    {
      return { std::make_shared<AnalyticOperationExpr<MulOp,VectorOpValue,MatrixOpValue,VectorOpValue>>(const_value(x1),x2) };
    }

  // operator/
    
    inline ScalarExpr
    operator/(const ScalarExpr& x1, const ScalarExpr& x2)
    {
      return { std::make_shared<AnalyticOperationExpr<DivOp,ScalarOpValue,ScalarOpValue,ScalarOpValue>>(x1,x2) };
    }

    template<typename X1>
      requires (!IsScalarExprOrVar<X1>)
    inline ScalarExpr
    operator/(const X1& x1, const ScalarExpr& x2)
    {
      return operator/(const_value(x1),x2);
    }

    template<typename X2>
      requires (!IsScalarExprOrVar<X2>)
    inline ScalarExpr
    operator/(const ScalarExpr& x1, const X2& x2)
    {
      return operator/(x1, const_value(x2));
    }
    
    inline VectorExpr
    operator/(const VectorExpr& x1, const ScalarExpr& x2)
    {
      return { std::make_shared<AnalyticOperationExpr<DivOp,VectorOpValue,VectorOpValue,ScalarOpValue>>(x1,x2) };
    }

  // Other operators

    inline ScalarExpr
    pow(const ScalarExpr& x1, const ScalarExpr& x2)
    {
      return { std::make_shared<AnalyticOperationExpr<PowOp,ScalarOpValue,ScalarOpValue,ScalarOpValue>>(x1,x2) };
    }

    template<typename X1>
      requires (!IsScalarExprOrVar<X1>)
    inline ScalarExpr
    pow(const X1& x1, const ScalarExpr& x2)
    {
      return pow(const_value(x1),x2);
    }

    template<typename X2>
      requires (!IsScalarExprOrVar<X2>)
    inline ScalarExpr
    pow(const ScalarExpr& x1, const X2& x2)
    {
      return pow(x1,const_value(x2));
    }

    inline ScalarExpr
    sqr(const ScalarExpr& x1)
    {
      return { std::make_shared<AnalyticOperationExpr<SqrOp,ScalarOpValue,ScalarOpValue>>(x1) };
    }

    inline ScalarExpr 
    sqrt(const ScalarExpr& x1)
    {
      return { std::make_shared<AnalyticOperationExpr<SqrtOp,ScalarOpValue,ScalarOpValue>>(x1) };
    }

    inline ScalarExpr
    exp(const ScalarExpr& x1)
    {
      return { std::make_shared<AnalyticOperationExpr<ExpOp,ScalarOpValue,ScalarOpValue>>(x1) };
    }

    inline ScalarExpr
    log(const ScalarExpr& x1)
    {
      return { std::make_shared<AnalyticOperationExpr<LogOp,ScalarOpValue,ScalarOpValue>>(x1) };
    }

    inline ScalarExpr
    cos(const ScalarExpr& x1)
    {
      return { std::make_shared<AnalyticOperationExpr<CosOp,ScalarOpValue,ScalarOpValue>>(x1) };
    }

    inline ScalarExpr
    sin(const ScalarExpr& x1)
    {
      return { std::make_shared<AnalyticOperationExpr<SinOp,ScalarOpValue,ScalarOpValue>>(x1) };
    }

    inline ScalarExpr
    tan(const ScalarExpr& x1)
    {
      return { std::make_shared<AnalyticOperationExpr<TanOp,ScalarOpValue,ScalarOpValue>>(x1) };
    }

    inline ScalarExpr
    acos(const ScalarExpr& x1)
    {
      return { std::make_shared<AnalyticOperationExpr<AcosOp,ScalarOpValue,ScalarOpValue>>(x1) };
    }

    inline ScalarExpr
    asin(const ScalarExpr& x1)
    {
      return { std::make_shared<AnalyticOperationExpr<AsinOp,ScalarOpValue,ScalarOpValue>>(x1) };
    }

    inline ScalarExpr
    atan(const ScalarExpr& x1)
    {
      return { std::make_shared<AnalyticOperationExpr<AtanOp,ScalarOpValue,ScalarOpValue>>(x1) };
    }

    inline ScalarExpr
    atan2(const ScalarExpr& x1, const ScalarExpr& x2)
    {
      return { std::make_shared<AnalyticOperationExpr<Atan2Op,ScalarOpValue,ScalarOpValue,ScalarOpValue>>(x1,x2) };
    }

    inline ScalarExpr
    cosh(const ScalarExpr& x1)
    {
      return { std::make_shared<AnalyticOperationExpr<CoshOp,ScalarOpValue,ScalarOpValue>>(x1) };
    }

    inline ScalarExpr
    sinh(const ScalarExpr& x1)
    {
      return { std::make_shared<AnalyticOperationExpr<SinhOp,ScalarOpValue,ScalarOpValue>>(x1) };
    }

    inline ScalarExpr
    tanh(const ScalarExpr& x1)
    {
      return { std::make_shared<AnalyticOperationExpr<TanhOp,ScalarOpValue,ScalarOpValue>>(x1) };
    }

    inline ScalarExpr
    abs(const ScalarExpr& x1)
    {
      return { std::make_shared<AnalyticOperationExpr<AbsOp,ScalarOpValue,ScalarOpValue>>(x1) };
    }

  // Structural operators

    inline ScalarExpr _add_to_vec(const ScalarExpr& x)
    {
      return x;
    }

    inline ScalarExpr _add_to_vec(double x)
    {
      return const_value(x);
    }

    template<class X>
    concept IsNotScalarExpr = !std::is_base_of_v<ScalarExpr,X>;

    template<typename X1>
      requires IsNotScalarExpr<X1>
    inline VectorExpr
    vec(const X1& x1)
    {
      return { std::make_shared<AnalyticOperationExpr<VectorOp,VectorOpValue,
        ScalarOpValue>>(
          _add_to_vec(x1)) };
    }

    template<typename X1, typename X2>
      requires (IsNotScalarExpr<X1> && IsNotScalarExpr<X2>)
    inline VectorExpr
    vec(const X1& x1, const X2& x2)
    {
      return { std::make_shared<AnalyticOperationExpr<VectorOp,VectorOpValue,
        ScalarOpValue,ScalarOpValue>>(
          _add_to_vec(x1),_add_to_vec(x2)) };
    }

    template<typename X1, typename X2, typename X3>
      requires (IsNotScalarExpr<X1> && IsNotScalarExpr<X2> && IsNotScalarExpr<X3>)
    inline VectorExpr
    vec(const X1& x1, const X2& x2, const X3& x3)
    {
      return { std::make_shared<AnalyticOperationExpr<VectorOp,VectorOpValue,
        ScalarOpValue,ScalarOpValue,ScalarOpValue>>(
          _add_to_vec(x1),_add_to_vec(x2),_add_to_vec(x3)) };
    }

    template<typename X1, typename X2, typename X3, typename X4>
      requires (IsNotScalarExpr<X1> && IsNotScalarExpr<X2> && IsNotScalarExpr<X3> && IsNotScalarExpr<X4>)
    inline VectorExpr
    vec(const X1& x1, const X2& x2, const X3& x3, const X4& x4)
    {
      return { std::make_shared<AnalyticOperationExpr<VectorOp,VectorOpValue,
        ScalarOpValue,ScalarOpValue,ScalarOpValue,ScalarOpValue>>(
          _add_to_vec(x1),_add_to_vec(x2),_add_to_vec(x3),_add_to_vec(x4)) };
    }

    template<typename X1, typename X2, typename X3, typename X4, typename X5>
      requires (IsNotScalarExpr<X1> && IsNotScalarExpr<X2> && IsNotScalarExpr<X3> && IsNotScalarExpr<X4> && IsNotScalarExpr<X5>)
    inline VectorExpr
    vec(const X1& x1, const X2& x2, const X3& x3, const X4& x4, const X5& x5)
    {
      return { std::make_shared<AnalyticOperationExpr<VectorOp,VectorOpValue,
        ScalarOpValue,ScalarOpValue,ScalarOpValue,ScalarOpValue,ScalarOpValue>>(
          _add_to_vec(x1),_add_to_vec(x2),_add_to_vec(x3),_add_to_vec(x4),_add_to_vec(x5)) };
    }

    template<typename X1, typename X2, typename X3, typename X4, typename X5, typename X6>
      requires (IsNotScalarExpr<X1> && IsNotScalarExpr<X2> && IsNotScalarExpr<X3> && IsNotScalarExpr<X4> && IsNotScalarExpr<X5> && IsNotScalarExpr<X6>)
    inline VectorExpr
    vec(const X1& x1, const X2& x2, const X3& x3, const X4& x4, const X5& x5, const X6& x6)
    {
      return { std::make_shared<AnalyticOperationExpr<VectorOp,VectorOpValue,
        ScalarOpValue,ScalarOpValue,ScalarOpValue,ScalarOpValue,ScalarOpValue,ScalarOpValue>>(
          _add_to_vec(x1),_add_to_vec(x2),_add_to_vec(x3),_add_to_vec(x4),_add_to_vec(x5),_add_to_vec(x6)) };
    }

    template<typename X1, typename X2, typename X3, typename X4, typename X5, typename X6, typename X7>
      requires (IsNotScalarExpr<X1> && IsNotScalarExpr<X2> && IsNotScalarExpr<X3> && IsNotScalarExpr<X4> && IsNotScalarExpr<X5> && IsNotScalarExpr<X6> && IsNotScalarExpr<X7>)
    inline VectorExpr
    vec(const X1& x1, const X2& x2, const X3& x3, const X4& x4, const X5& x5, const X6& x6, const X7& x7)
    {
      return { std::make_shared<AnalyticOperationExpr<VectorOp,VectorOpValue,
        ScalarOpValue,ScalarOpValue,ScalarOpValue,ScalarOpValue,ScalarOpValue,ScalarOpValue,ScalarOpValue>>(
          _add_to_vec(x1),_add_to_vec(x2),_add_to_vec(x3),_add_to_vec(x4),_add_to_vec(x5),_add_to_vec(x6),_add_to_vec(x7)) };
    }

    template<typename X1, typename X2, typename X3, typename X4, typename X5, typename X6, typename X7, typename X8>
      requires (IsNotScalarExpr<X1> && IsNotScalarExpr<X2> && IsNotScalarExpr<X3> && IsNotScalarExpr<X4> && IsNotScalarExpr<X5> && IsNotScalarExpr<X6> && IsNotScalarExpr<X7> && IsNotScalarExpr<X8>)
    inline VectorExpr
    vec(const X1& x1, const X2& x2, const X3& x3, const X4& x4, const X5& x5, const X6& x6, const X7& x7, const X8& x8)
    {
      return { std::make_shared<AnalyticOperationExpr<VectorOp,VectorOpValue,
        ScalarOpValue,ScalarOpValue,ScalarOpValue,ScalarOpValue,ScalarOpValue,ScalarOpValue,ScalarOpValue,ScalarOpValue>>(
          _add_to_vec(x1),_add_to_vec(x2),_add_to_vec(x3),_add_to_vec(x4),_add_to_vec(x5),_add_to_vec(x6),_add_to_vec(x7),_add_to_vec(x8)) };
    }

    template<typename X1, typename X2, typename X3, typename X4, typename X5, typename X6, typename X7, typename X8, typename X9>
      requires (IsNotScalarExpr<X1> && IsNotScalarExpr<X2> && IsNotScalarExpr<X3> && IsNotScalarExpr<X4> && IsNotScalarExpr<X5> && IsNotScalarExpr<X6> && IsNotScalarExpr<X7> && IsNotScalarExpr<X8> && IsNotScalarExpr<X9>)
    inline VectorExpr
    vec(const X1& x1, const X2& x2, const X3& x3, const X4& x4, const X5& x5, const X6& x6, const X7& x7, const X8& x8, const X9& x9)
    {
      return { std::make_shared<AnalyticOperationExpr<VectorOp,VectorOpValue,
        ScalarOpValue,ScalarOpValue,ScalarOpValue,ScalarOpValue,ScalarOpValue,ScalarOpValue,ScalarOpValue,ScalarOpValue,ScalarOpValue>>(
          _add_to_vec(x1),_add_to_vec(x2),_add_to_vec(x3),_add_to_vec(x4),_add_to_vec(x5),_add_to_vec(x6),_add_to_vec(x7),_add_to_vec(x8),_add_to_vec(x9)) };
    }

    template<typename X1, typename X2, typename X3, typename X4, typename X5, typename X6, typename X7, typename X8, typename X9, typename X10>
      requires (IsNotScalarExpr<X1> && IsNotScalarExpr<X2> && IsNotScalarExpr<X3> && IsNotScalarExpr<X4> && IsNotScalarExpr<X5> && IsNotScalarExpr<X6> && IsNotScalarExpr<X7> && IsNotScalarExpr<X8> && IsNotScalarExpr<X9> && IsNotScalarExpr<X10>)
    inline VectorExpr
    vec(const X1& x1, const X2& x2, const X3& x3, const X4& x4, const X5& x5, const X6& x6, const X7& x7, const X8& x8, const X9& x9, const X10& x10)
    {
      return { std::make_shared<AnalyticOperationExpr<VectorOp,VectorOpValue,
        ScalarOpValue,ScalarOpValue,ScalarOpValue,ScalarOpValue,ScalarOpValue,ScalarOpValue,ScalarOpValue,ScalarOpValue,ScalarOpValue,ScalarOpValue>>(
          _add_to_vec(x1),_add_to_vec(x2),_add_to_vec(x3),_add_to_vec(x4),_add_to_vec(x5),_add_to_vec(x6),_add_to_vec(x7),_add_to_vec(x8),_add_to_vec(x9),_add_to_vec(x10)) };
    }

    // Generic variadic case, cannot handle const values (int, double) for now

    template<typename... X>
    inline VectorExpr
    vec(const AnalyticExprWrapper<X>&... x)
    {
      return { std::make_shared<AnalyticOperationExpr<VectorOp,VectorOpValue,X...>>(x...) };
    }

    inline ScalarExpr
    ind(Index i, const VectorExpr& x1)
    {
      return { std::make_shared<AnalyticOperationExpr<ComponentOp,ScalarOpValue,VectorOpValue>>(x1,i) };
    }

    template<typename... X>
    inline MatrixExpr
    mat(const std::shared_ptr<AnalyticExpr<X>>&... x)
    {
      return { std::make_shared<AnalyticOperationExpr<MatrixOp,MatrixOpValue,X...>>(
        AnalyticOperationExpr<MatrixOp,MatrixOpValue,X...>(x...)) };
    }

  // Vectorial operations

    inline ScalarExpr
    det(const MatrixExpr& x1)
    {
      return { std::make_shared<AnalyticOperationExpr<DetOp,ScalarOpValue,MatrixOpValue>>(x1) };
    }

    inline ScalarExpr
    det(const VectorExpr& x1, const VectorExpr& x2)
    {
      return { std::make_shared<AnalyticOperationExpr<DetOp,ScalarOpValue,VectorOpValue,VectorOpValue>>(x1,x2) };
    }

    inline ScalarExpr
    det(const VectorExpr& x1, const VectorExpr& x2, const VectorExpr& x3)
    {
      return { std::make_shared<AnalyticOperationExpr<DetOp,ScalarOpValue,VectorOpValue,VectorOpValue,VectorOpValue>>(x1,x2,x3) };
    }
}