/** 
 *  Codac binding (core)
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include <sstream>
#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <codac2_analytic_variables.h>
#include <codac2_analytic_operations.h>
#include <codac2_AnalyticFunction.h>
#include "codac2_py_AnalyticFunction_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py)
#include "codac2_py_FunctionBase_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py)
#include "codac2_py_ExprWrapper.h"

using namespace codac2;
namespace py = pybind11;
using namespace pybind11::literals;

#define bind_mode_(exported, op_name, op, doc) \
  \
  exported \
  \
    /* Several cases of scalar inputs */ \
    .def(op_name, [](AnalyticFunction<T>& f, const EvalMode& m) { return f.op(m); }, doc) \
    .def(op_name, [](AnalyticFunction<T>& f, const EvalMode& m, I x1) { return f.op(m,x1); }, doc) \
    .def(op_name, [](AnalyticFunction<T>& f, const EvalMode& m, I x1, I x2) { return f.op(m,x1,x2); }, doc) \
    .def(op_name, [](AnalyticFunction<T>& f, const EvalMode& m, I x1, I x2, I x3) { return f.op(m,x1,x2,x3); }, doc) \
    .def(op_name, [](AnalyticFunction<T>& f, const EvalMode& m, I x1, I x2, I x3, I x4) { return f.op(m,x1,x2,x3,x4); }, doc) \
    .def(op_name, [](AnalyticFunction<T>& f, const EvalMode& m, I x1, I x2, I x3, I x4, I x5) { return f.op(m,x1,x2,x3,x4,x5); }, doc) \
    .def(op_name, [](AnalyticFunction<T>& f, const EvalMode& m, I x1, I x2, I x3, I x4, I x5, I x6) { return f.op(m,x1,x2,x3,x4,x5,x6); }, doc) \
    .def(op_name, [](AnalyticFunction<T>& f, const EvalMode& m, I x1, I x2, I x3, I x4, I x5, I x6, I x7) { return f.op(m,x1,x2,x3,x4,x5,x6,x7); }, doc) \
    .def(op_name, [](AnalyticFunction<T>& f, const EvalMode& m, I x1, I x2, I x3, I x4, I x5, I x6, I x7, I x8) { return f.op(m,x1,x2,x3,x4,x5,x6,x7,x8); }, doc) \
    .def(op_name, [](AnalyticFunction<T>& f, const EvalMode& m, I x1, I x2, I x3, I x4, I x5, I x6, I x7, I x8, I x9) { return f.op(m,x1,x2,x3,x4,x5,x6,x7,x8,x9); }, doc) \
    .def(op_name, [](AnalyticFunction<T>& f, const EvalMode& m, I x1, I x2, I x3, I x4, I x5, I x6, I x7, I x8, I x9, I x10) { return f.op(m,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10); }, doc) \
  \
    /* Several cases of vector inputs */ \
    .def(op_name, [](AnalyticFunction<T>& f, const EvalMode& m) { return f.op(m); }, doc) \
    .def(op_name, [](AnalyticFunction<T>& f, const EvalMode& m, IV x1) { return f.op(m,x1); }, doc) \
    .def(op_name, [](AnalyticFunction<T>& f, const EvalMode& m, IV x1, IV x2) { return f.op(m,x1,x2); }, doc) \
    .def(op_name, [](AnalyticFunction<T>& f, const EvalMode& m, IV x1, IV x2, IV x3) { return f.op(m,x1,x2,x3); }, doc) \
    .def(op_name, [](AnalyticFunction<T>& f, const EvalMode& m, IV x1, IV x2, IV x3, IV x4) { return f.op(m,x1,x2,x3,x4); }, doc) \
    .def(op_name, [](AnalyticFunction<T>& f, const EvalMode& m, IV x1, IV x2, IV x3, IV x4, IV x5) { return f.op(m,x1,x2,x3,x4,x5); }, doc) \
    .def(op_name, [](AnalyticFunction<T>& f, const EvalMode& m, IV x1, IV x2, IV x3, IV x4, IV x5, IV x6) { return f.op(m,x1,x2,x3,x4,x5,x6); }, doc) \
    .def(op_name, [](AnalyticFunction<T>& f, const EvalMode& m, IV x1, IV x2, IV x3, IV x4, IV x5, IV x6, IV x7) { return f.op(m,x1,x2,x3,x4,x5,x6,x7); }, doc) \
    .def(op_name, [](AnalyticFunction<T>& f, const EvalMode& m, IV x1, IV x2, IV x3, IV x4, IV x5, IV x6, IV x7, IV x8) { return f.op(m,x1,x2,x3,x4,x5,x6,x7,x8); }, doc) \
    .def(op_name, [](AnalyticFunction<T>& f, const EvalMode& m, IV x1, IV x2, IV x3, IV x4, IV x5, IV x6, IV x7, IV x8, IV x9) { return f.op(m,x1,x2,x3,x4,x5,x6,x7,x8,x9); }, doc) \
    .def(op_name, [](AnalyticFunction<T>& f, const EvalMode& m, IV x1, IV x2, IV x3, IV x4, IV x5, IV x6, IV x7, IV x8, IV x9, IV x10) { return f.op(m,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10); }, doc) \
  \
  ; \

#define bind_(exported, op_name, op, doc) \
  \
  exported \
  \
    /* Several cases of scalar inputs */ \
    .def(op_name, [](AnalyticFunction<T>& f) { return f.op(); }, doc) \
    .def(op_name, [](AnalyticFunction<T>& f, I x1) { return f.op(x1); }, doc) \
    .def(op_name, [](AnalyticFunction<T>& f, I x1, I x2) { return f.op(x1,x2); }, doc) \
    .def(op_name, [](AnalyticFunction<T>& f, I x1, I x2, I x3) { return f.op(x1,x2,x3); }, doc) \
    .def(op_name, [](AnalyticFunction<T>& f, I x1, I x2, I x3, I x4) { return f.op(x1,x2,x3,x4); }, doc) \
    .def(op_name, [](AnalyticFunction<T>& f, I x1, I x2, I x3, I x4, I x5) { return f.op(x1,x2,x3,x4,x5); }, doc) \
    .def(op_name, [](AnalyticFunction<T>& f, I x1, I x2, I x3, I x4, I x5, I x6) { return f.op(x1,x2,x3,x4,x5,x6); }, doc) \
    .def(op_name, [](AnalyticFunction<T>& f, I x1, I x2, I x3, I x4, I x5, I x6, I x7) { return f.op(x1,x2,x3,x4,x5,x6,x7); }, doc) \
    .def(op_name, [](AnalyticFunction<T>& f, I x1, I x2, I x3, I x4, I x5, I x6, I x7, I x8) { return f.op(x1,x2,x3,x4,x5,x6,x7,x8); }, doc) \
    .def(op_name, [](AnalyticFunction<T>& f, I x1, I x2, I x3, I x4, I x5, I x6, I x7, I x8, I x9) { return f.op(x1,x2,x3,x4,x5,x6,x7,x8,x9); }, doc) \
    .def(op_name, [](AnalyticFunction<T>& f, I x1, I x2, I x3, I x4, I x5, I x6, I x7, I x8, I x9, I x10) { return f.op(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10); }, doc) \
  \
    /* Several cases of vector inputs */ \
    .def(op_name, [](AnalyticFunction<T>& f) { return f.op(); }, doc) \
    .def(op_name, [](AnalyticFunction<T>& f, IV x1) { return f.op(x1); }, doc) \
    .def(op_name, [](AnalyticFunction<T>& f, IV x1, IV x2) { return f.op(x1,x2); }, doc) \
    .def(op_name, [](AnalyticFunction<T>& f, IV x1, IV x2, IV x3) { return f.op(x1,x2,x3); }, doc) \
    .def(op_name, [](AnalyticFunction<T>& f, IV x1, IV x2, IV x3, IV x4) { return f.op(x1,x2,x3,x4); }, doc) \
    .def(op_name, [](AnalyticFunction<T>& f, IV x1, IV x2, IV x3, IV x4, IV x5) { return f.op(x1,x2,x3,x4,x5); }, doc) \
    .def(op_name, [](AnalyticFunction<T>& f, IV x1, IV x2, IV x3, IV x4, IV x5, IV x6) { return f.op(x1,x2,x3,x4,x5,x6); }, doc) \
    .def(op_name, [](AnalyticFunction<T>& f, IV x1, IV x2, IV x3, IV x4, IV x5, IV x6, IV x7) { return f.op(x1,x2,x3,x4,x5,x6,x7); }, doc) \
    .def(op_name, [](AnalyticFunction<T>& f, IV x1, IV x2, IV x3, IV x4, IV x5, IV x6, IV x7, IV x8) { return f.op(x1,x2,x3,x4,x5,x6,x7,x8); }, doc) \
    .def(op_name, [](AnalyticFunction<T>& f, IV x1, IV x2, IV x3, IV x4, IV x5, IV x6, IV x7, IV x8, IV x9) { return f.op(x1,x2,x3,x4,x5,x6,x7,x8,x9); }, doc) \
    .def(op_name, [](AnalyticFunction<T>& f, IV x1, IV x2, IV x3, IV x4, IV x5, IV x6, IV x7, IV x8, IV x9, IV x10) { return f.op(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10); }, doc) \
  \
  ; \


AnalyticExprWrapper<ScalarOpValue> pyobj_to_ScalarExpr(py::object& obj)
{
  if(py::isinstance<ScalarExpr>(obj))
    return obj.cast<ScalarExpr>();
  assert(false && "unable to convert py::object to ScalarExpr");
  return { nullptr };
}

std::shared_ptr<ExprBase> pyobj_to_ExprBase(const py::object& obj)
{
  if(obj)
  {
    if(py::isinstance<AnalyticExprWrapper<ScalarOpValue>>(obj))
      return obj.cast<AnalyticExprWrapper<ScalarOpValue>>()->copy();

    else if(py::isinstance<AnalyticExprWrapper<VectorOpValue>>(obj))
      return obj.cast<AnalyticExprWrapper<VectorOpValue>>()->copy();
  }

  assert_release("invalid expression");
  return nullptr;
}

FunctionArgsList create_FunctionArgsList(const std::vector<py::object>& l)
{
  FunctionArgsList args {};
  Index i = 0;

  for(const auto& li : l)
  {
    i++;

    if(py::isinstance<ScalarVar>(li))
      args.push_back(li.cast<ScalarVar>().arg_copy());

    else if(py::isinstance<VectorVar>(li))
      args.push_back(li.cast<VectorVar>().arg_copy());

    else
      throw std::invalid_argument("Argument " + std::to_string(i) + " is invalid. Only variables are accepted.");
  }

  return args;
}

template<typename T>
void export_AnalyticFunction(py::module& m, const std::string& export_name)
{
  py::class_<AnalyticFunction<T>, std::shared_ptr<AnalyticFunction<T>> /* due to enable_shared_from_this */>
    exported(m, export_name.c_str(), ANALYTICFUNCTION_MAIN);
  exported
  
    .def(py::init(
      [](const std::vector<py::object>& l, const AnalyticExprWrapper<T>& expr)
      {
        return std::make_unique<AnalyticFunction<T>>(
          create_FunctionArgsList(l),
          std::dynamic_pointer_cast<AnalyticExpr<T>>(expr.get()->copy())
        );
      }
    ), ANALYTICFUNCTION_T_ANALYTICFUNCTION_CONST_FUNCTIONARGSLIST_REF_CONST_SHARED_PTR_ANALYTICEXPR_T_REF)
  ;

  if constexpr(std::is_same_v<T,VectorOpValue>)
  {
    exported.def(py::init(
      [](const std::vector<py::object>& l, std::vector<py::object>& v_expr)
      {
        return std::make_unique<AnalyticFunction<T>>(
          create_FunctionArgsList(l),
          [&]() -> VectorExpr {
            switch(v_expr.size())
            {
              case 1:
                return vec(
                  pyobj_to_ScalarExpr(v_expr[0])
                );
              case 2:
                return vec(
                  pyobj_to_ScalarExpr(v_expr[0]), pyobj_to_ScalarExpr(v_expr[1])
                );
              case 3:
                return vec(
                  pyobj_to_ScalarExpr(v_expr[0]), pyobj_to_ScalarExpr(v_expr[1]),
                  pyobj_to_ScalarExpr(v_expr[2])
                );
              case 4:
                return vec(
                  pyobj_to_ScalarExpr(v_expr[0]), pyobj_to_ScalarExpr(v_expr[1]),
                  pyobj_to_ScalarExpr(v_expr[2]), pyobj_to_ScalarExpr(v_expr[3])
                );
              case 5:
                return vec(
                  pyobj_to_ScalarExpr(v_expr[0]), pyobj_to_ScalarExpr(v_expr[1]),
                  pyobj_to_ScalarExpr(v_expr[2]), pyobj_to_ScalarExpr(v_expr[3]),
                  pyobj_to_ScalarExpr(v_expr[4])
                );
              case 6:
                return vec(
                  pyobj_to_ScalarExpr(v_expr[0]), pyobj_to_ScalarExpr(v_expr[1]),
                  pyobj_to_ScalarExpr(v_expr[2]), pyobj_to_ScalarExpr(v_expr[3]),
                  pyobj_to_ScalarExpr(v_expr[4]), pyobj_to_ScalarExpr(v_expr[5])
                );
              case 7:
                return vec(
                  pyobj_to_ScalarExpr(v_expr[0]), pyobj_to_ScalarExpr(v_expr[1]),
                  pyobj_to_ScalarExpr(v_expr[2]), pyobj_to_ScalarExpr(v_expr[3]),
                  pyobj_to_ScalarExpr(v_expr[4]), pyobj_to_ScalarExpr(v_expr[5]),
                  pyobj_to_ScalarExpr(v_expr[6])
                );
              case 8:
                return vec(
                  pyobj_to_ScalarExpr(v_expr[0]), pyobj_to_ScalarExpr(v_expr[1]),
                  pyobj_to_ScalarExpr(v_expr[2]), pyobj_to_ScalarExpr(v_expr[3]),
                  pyobj_to_ScalarExpr(v_expr[4]), pyobj_to_ScalarExpr(v_expr[5]),
                  pyobj_to_ScalarExpr(v_expr[6]), pyobj_to_ScalarExpr(v_expr[7])
                );
              case 9:
                return vec(
                  pyobj_to_ScalarExpr(v_expr[0]), pyobj_to_ScalarExpr(v_expr[1]),
                  pyobj_to_ScalarExpr(v_expr[2]), pyobj_to_ScalarExpr(v_expr[3]),
                  pyobj_to_ScalarExpr(v_expr[4]), pyobj_to_ScalarExpr(v_expr[5]),
                  pyobj_to_ScalarExpr(v_expr[6]), pyobj_to_ScalarExpr(v_expr[7]),
                  pyobj_to_ScalarExpr(v_expr[8])
                );
              case 10:
                return vec(
                  pyobj_to_ScalarExpr(v_expr[0]), pyobj_to_ScalarExpr(v_expr[1]),
                  pyobj_to_ScalarExpr(v_expr[2]), pyobj_to_ScalarExpr(v_expr[3]),
                  pyobj_to_ScalarExpr(v_expr[4]), pyobj_to_ScalarExpr(v_expr[5]),
                  pyobj_to_ScalarExpr(v_expr[6]), pyobj_to_ScalarExpr(v_expr[7]),
                  pyobj_to_ScalarExpr(v_expr[8]), pyobj_to_ScalarExpr(v_expr[9])
                );
              default:
                assert_release("cannot create AnalyticFunction with provided vector arguments");
                return { nullptr };
            }
          }()
        );
      }
    ), ANALYTICFUNCTION_T_ANALYTICFUNCTION_CONST_FUNCTIONARGSLIST_REF_CONST_SCALAREXPRLIST_REF);
  }

  exported

    .def("input_size", &AnalyticFunction<T>::input_size,
      INDEX_FUNCTIONBASE_E_INPUT_SIZE_CONST)

    .def("__call__", [](const AnalyticFunction<T>& f, py::list& x)
      {
        std::vector<std::shared_ptr<ExprBase>> v(x.size());
        for(size_t i = 0 ; i < x.size() ; i++)
          v[i] = pyobj_to_ExprBase(x[i]);
        return AnalyticExprWrapper<T>(std::dynamic_pointer_cast<AnalyticExpr<T>>(f(v)->copy()));
      },
      SHARED_PTR_E_FUNCTIONBASE_E_OPERATORCALL_CONST_X_REF_VARIADIC_CONST)
  ;

  using I = const Interval&; using IV = const IntervalVector&;

  bind_(exported, "real_eval", real_eval, AUTO_ANALYTICFUNCTION_T_REAL_EVAL_CONST_ARGS_REF_VARIADIC_CONST);
  bind_mode_(exported, "eval", eval, T_DOMAIN_ANALYTICFUNCTION_T_EVAL_CONST_ARGS_REF_VARIADIC_CONST);
  bind_(exported, "eval", eval, T_DOMAIN_ANALYTICFUNCTION_T_EVAL_CONST_ARGS_REF_VARIADIC_CONST);
  bind_(exported, "diff", diff, AUTO_ANALYTICFUNCTION_T_DIFF_CONST_ARGS_REF_VARIADIC_CONST);

  exported

    // Mixed scalar/vector inputs are not supported yet

    .def("__repr__", [](const AnalyticFunction<T>& f) {
          std::ostringstream stream;
          stream << f;
          return std::string(stream.str()); 
        },
      OSTREAM_REF_OPERATOROUT_OSTREAM_REF_CONST_ANALYTICFUNCTION_T_REF)
  ;
}