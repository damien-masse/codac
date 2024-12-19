/** 
 *  \file codac2_Trajectory_operator.h
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include "codac2_analytic_operations.h"
#include "codac2_SampledTrajectory.h"

namespace codac2
{
  struct TrajectoryOp
  {
    static IntervalVector fwd(const SampledTrajectory<Vector>& x1, const Interval& x2);
    static VectorOpValue fwd(const SampledTrajectory<Vector>& x1, const ScalarOpValue& x2);
    static void bwd(const SampledTrajectory<Vector>& x1, const IntervalVector& y, Interval& x2);
  };


  template<>
  class AnalyticOperationExpr<TrajectoryOp,VectorOpValue,ScalarOpValue>
    : public AnalyticExpr<VectorOpValue>, public OperationExprBase<AnalyticExpr<ScalarOpValue>>
  {
    public:

      AnalyticOperationExpr(const SampledTrajectory<Vector>& x1, const std::shared_ptr<AnalyticExpr<ScalarOpValue>>& x2)
        : OperationExprBase<AnalyticExpr<ScalarOpValue>>(x2), _x1(x1)
      { }

      AnalyticOperationExpr(const AnalyticOperationExpr& e)
        : OperationExprBase<AnalyticExpr<ScalarOpValue>>(e), _x1(e._x1)
      { }

      std::shared_ptr<ExprBase> copy() const
      {
        return std::make_shared<AnalyticOperationExpr<TrajectoryOp,VectorOpValue,ScalarOpValue>>(*this);
      }

      void replace_expr(const ExprID& old_expr_id, const std::shared_ptr<ExprBase>& new_expr)
      {
        return OperationExprBase<AnalyticExpr<ScalarOpValue>>::replace_expr(old_expr_id, new_expr);
      }
      
      VectorOpValue fwd_eval(ValuesMap& v, Index total_input_size) const
      {
        return AnalyticExpr<VectorOpValue>::init_value(
          v, TrajectoryOp::fwd(_x1, std::get<0>(this->_x)->fwd_eval(v, total_input_size)));
      }
      
      void bwd_eval(ValuesMap& v) const
      {
        TrajectoryOp::bwd(_x1, AnalyticExpr<VectorOpValue>::value(v).a, std::get<0>(this->_x)->value(v).a);
        std::get<0>(this->_x)->bwd_eval(v);
      }

      virtual bool belongs_to_args_list(const FunctionArgsList& args) const
      {
        return std::get<0>(this->_x)->belongs_to_args_list(args);
      }

    protected:

      const SampledTrajectory<Vector> _x1;
  };

}