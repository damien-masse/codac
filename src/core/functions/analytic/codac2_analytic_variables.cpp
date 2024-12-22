/** 
 *  codac2_ExprBase.cpp
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include "codac2_analytic_variables.h"
#include "codac2_AnalyticExprWrapper.h"

using namespace std;
using namespace codac2;


// ScalarVar class

  ScalarVar::ScalarVar()
  { }

  ScalarVar::ScalarVar(const ScalarVar& x)
    : AnalyticVarExpr<ScalarOpValue>(x)
  { }

  std::shared_ptr<VarBase> ScalarVar::arg_copy() const
  {
    return std::make_shared<ScalarVar>(*this);
  }

  std::shared_ptr<ExprBase> ScalarVar::copy() const
  {
    return std::make_shared<ScalarVar>(*this);
  }

  Index ScalarVar::size() const
  {
    return 1;
  }

  AnalyticExprWrapper<ScalarOpValue> ScalarVar::operator-() const
  {
    return { std::make_shared<AnalyticOperationExpr<SubOp,ScalarOpValue,ScalarOpValue>>(*this) };
  }


// VectorVar class

  VectorVar::VectorVar(Index n)
    : _n(n)
  {
    assert_release(n > 0);
  }

  VectorVar::VectorVar(const VectorVar& x)
    : AnalyticVarExpr<VectorOpValue>(x), _n(x._n)
  { }

  std::shared_ptr<VarBase> VectorVar::arg_copy() const
  {
    return std::make_shared<VectorVar>(*this);
  }

  std::shared_ptr<ExprBase> VectorVar::copy() const
  {
    return std::make_shared<VectorVar>(*this);
  }

  Index VectorVar::size() const
  {
    return _n;
  }

  AnalyticExprWrapper<ScalarOpValue> VectorVar::operator[](Index i) const
  {
    assert_release(i >= 0 && i < _n);
    return { std::make_shared<AnalyticOperationExpr<ComponentOp,ScalarOpValue,VectorOpValue>>(
      std::dynamic_pointer_cast<AnalyticExpr<VectorOpValue>>(this->copy()), i) };
  }

  AnalyticExprWrapper<VectorOpValue> VectorVar::subvector(Index i, Index j) const
  {
    assert_release(i >= 0 && i < _n && j >= i && j < _n);
    return { std::make_shared<AnalyticOperationExpr<SubvectorOp,VectorOpValue,VectorOpValue>>(
      std::dynamic_pointer_cast<AnalyticExpr<VectorOpValue>>(this->copy()), i, j) };
  }