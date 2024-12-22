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
  template<typename T, typename D = typename Wrapper<T>::Domain>
    requires (std::is_same_v<D,Interval> || std::is_same_v<D,IntervalVector>)
  struct TrajectoryOp
  {
    static D fwd(const TrajectoryBase<T>& x1, const Interval& x2)
    {
      return x1(x2);
    }

    static typename ArgWrapper<T>::Domain fwd(const TrajectoryBase<T>& x1, const ScalarOpValue& x2)
    {
      IntervalMatrix d(x1.size(),x2.da.cols());

      return {
        fwd(x1,x2.m),
        fwd(x1,x2.a),
        d,
        x2.def_domain && x1.tdomain().is_superset(x2.m)
      };
    }

    static void bwd(
      [[maybe_unused]] const TrajectoryBase<T>& x1,
      [[maybe_unused]] const D& y,
      [[maybe_unused]] Interval& x2)
    {
      // todo
    }
  };

  template<typename T>
  class AnalyticOperationExpr<TrajectoryOp<T>,typename ArgWrapper<T>::Domain,ScalarOpValue>
    : public AnalyticExpr<typename ArgWrapper<T>::Domain>, public OperationExprBase<AnalyticExpr<ScalarOpValue>>
  {
    public:

      using O = typename ArgWrapper<T>::Domain;

      AnalyticOperationExpr(const TrajectoryBase<T>& x1, const std::shared_ptr<AnalyticExpr<ScalarOpValue>>& x2)
        : OperationExprBase<AnalyticExpr<ScalarOpValue>>(x2), _x1(x1)
      { }

      AnalyticOperationExpr(const AnalyticOperationExpr& e)
        : OperationExprBase<AnalyticExpr<ScalarOpValue>>(e), _x1(e._x1)
      { }

      std::shared_ptr<ExprBase> copy() const
      {
        return std::make_shared<AnalyticOperationExpr<TrajectoryOp<T>,O,ScalarOpValue>>(*this);
      }

      void replace_expr(const ExprID& old_expr_id, const std::shared_ptr<ExprBase>& new_expr)
      {
        return OperationExprBase<AnalyticExpr<ScalarOpValue>>::replace_expr(old_expr_id, new_expr);
      }
      
      O fwd_eval(ValuesMap& v, Index total_input_size) const
      {
        return AnalyticExpr<O>::init_value(
          v, TrajectoryOp<T>::fwd(_x1, std::get<0>(this->_x)->fwd_eval(v, total_input_size)));
      }
      
      void bwd_eval(ValuesMap& v) const
      {
        TrajectoryOp<T>::bwd(_x1, AnalyticExpr<O>::value(v).a, std::get<0>(this->_x)->value(v).a);
        std::get<0>(this->_x)->bwd_eval(v);
      }

      virtual bool belongs_to_args_list(const FunctionArgsList& args) const
      {
        return std::get<0>(this->_x)->belongs_to_args_list(args);
      }

    protected:

      const TrajectoryBase<T>& _x1;
  };

  template<typename T>
  AnalyticFunction<typename ArgWrapper<T>::Domain> TrajectoryBase<T>::as_function() const
  {
    using D = typename ArgWrapper<T>::Domain;

    ScalarVar t;
    return {{t},
      AnalyticExprWrapper<D>(std::make_shared<
        AnalyticOperationExpr
        <TrajectoryOp<T>, D, ScalarOpValue>>(*this,t))
    };
  }
}