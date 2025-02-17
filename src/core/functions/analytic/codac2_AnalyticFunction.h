/** 
 *  \file codac2_AnalyticFunction.h
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include <map>
#include "codac2_AnalyticExpr.h"
#include "codac2_Domain.h"
#include "codac2_analytic_variables.h"
#include "codac2_FunctionBase.h"
#include "codac2_template_tools.h"
#include "codac2_AnalyticExprWrapper.h"
#include "codac2_ScalarExprList.h"
#include "codac2_operators.h"

namespace codac2
{
  enum class EvalMode
  {
    NATURAL = 0x01,
    CENTERED = 0x02,
    DEFAULT = 0x03 // corresponds to (NATURAL|CENTERED)
  };

  EvalMode operator&(EvalMode a, EvalMode b);
  EvalMode operator|(EvalMode a, EvalMode b);

  template<typename T>
    requires (std::is_base_of_v<AnalyticTypeBase,T>)
  class AnalyticFunction : public FunctionBase<AnalyticExpr<T>>
  {
    public:

      // The implementation of constructors cannot be placed in the _impl.h file,
      // probably due to the use of deduction guidelines (see after this class).

      AnalyticFunction(const FunctionArgsList& args, const ScalarExprList& y)
        requires(std::is_same_v<T,VectorType>)
        : FunctionBase<AnalyticExpr<T>>(args, y)
      {
        assert_release(y->belongs_to_args_list(this->args()) && 
          "Invalid argument: variable not present in input arguments");
      }

      AnalyticFunction(const FunctionArgsList& args, const AnalyticExprWrapper<T>& y)
        : FunctionBase<AnalyticExpr<T>>(args, y)
      {
        assert_release(y->belongs_to_args_list(this->args()) && 
          "Invalid argument: variable not present in input arguments");
      }

      AnalyticFunction(const FunctionArgsList& args, const AnalyticVarExpr<T>& y)
        : AnalyticFunction(args, { std::dynamic_pointer_cast<AnalyticExpr<T>>(y.copy()) })
      { }

      AnalyticFunction(const AnalyticFunction<T>& f)
        : FunctionBase<AnalyticExpr<T>>(f)
      { }

      template<typename... X>
      AnalyticExprWrapper<T> operator()(const X&... x) const;

      template<typename... Args>
      auto real_eval(const Args&... x) const;

      template<typename... Args>
      typename T::Domain eval(const EvalMode& m, const Args&... x) const;

      template<typename... Args>
      typename T::Domain eval(const Args&... x) const;

      template<typename... Args>
      auto diff(const Args&... x) const;

      Index output_size() const;

      template<typename U>
      friend std::ostream& operator<<(std::ostream& os, [[maybe_unused]] const AnalyticFunction<U>& f);

    protected:

      template<typename Y>
      friend class CtcInverse;

      template<typename D>
      void add_value_to_arg_map(ValuesMap& v, const D& x, Index i) const;

      template<typename... Args>
      void fill_from_args(ValuesMap& v, const Args&... x) const;

      template<typename D>
      void intersect_value_from_arg_map(const ValuesMap& v, D& x, Index i) const;

      template<typename... Args>
      void intersect_from_args(const ValuesMap& v, Args&... x) const;

      template<bool NATURAL_EVAL,typename... Args>
      auto eval_(const Args&... x) const;

      template<typename... Args>
      void check_valid_inputs(const Args&... x) const;
  };

  AnalyticFunction(const FunctionArgsList&, double) -> 
    AnalyticFunction<ScalarType>;

  AnalyticFunction(const FunctionArgsList&, const Interval&) -> 
    AnalyticFunction<ScalarType>;

  AnalyticFunction(const FunctionArgsList&, std::initializer_list<int>) -> 
    AnalyticFunction<VectorType>;

  AnalyticFunction(const FunctionArgsList&, std::initializer_list<double>) -> 
    AnalyticFunction<VectorType>;

  AnalyticFunction(const FunctionArgsList&, std::initializer_list<Interval>) -> 
    AnalyticFunction<VectorType>;

  AnalyticFunction(const FunctionArgsList&, std::initializer_list<ScalarExpr>) -> 
    AnalyticFunction<VectorType>;

  AnalyticFunction(const FunctionArgsList&, std::initializer_list<ScalarVar>) -> 
    AnalyticFunction<VectorType>;

  template<typename T>
  AnalyticFunction(const FunctionArgsList&, const T&) -> 
    AnalyticFunction<typename ValueType<T>::Type>;
}

#include "codac2_AnalyticFunction_impl.h"