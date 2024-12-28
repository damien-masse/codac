/** 
 *  \file
 *  Codac binding (core)
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2021 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include <pybind11/pybind11.h>
#include <codac2_SampledTrajectory.h>
#include <codac2_AnalyticTrajectory.h>

namespace py = pybind11;
using namespace pybind11::literals;

namespace codac2
{
  // Common caster

    template<typename T>
    bool is_instance(const py::object& x)
    {
      return py::isinstance<T>(x);
    }

    template<typename T>
    const T& cast(const py::object& x)
    {
      if(!is_instance<T>(x))
      { assert_release("cast error"); }
      return x.cast<const T&>();
    }

  // Casting to trajectory types

    template<typename Tr>
      requires std::is_same_v<typename Tr::TrajType,VectorType>
    bool is_instance(const py::object& x)
    {
      const py::object& x_traj = x.attr("traj");
      return x_traj && py::isinstance<Tr>(x_traj);
    }

    template<typename Tr>
      requires std::is_same_v<typename Tr::TrajType,VectorType>
    const Tr& cast(const py::object& x)
    {
      assert(is_instance<Tr>(x));
      const py::object& x_traj = x.attr("traj");
      return x_traj.cast<const Tr&>();
    }

  // Casting to function/expression types

    template<typename T>
      requires std::is_same_v<T,ScalarExpr>
    const T& cast(const py::object& x)
    {
      if(!is_instance<T>(x))
      { assert("unable to cast the provided py::object into a ScalarExpr"); }
      return x.cast<const T&>();
    }

    template<typename T>
      requires std::is_same_v<T,std::shared_ptr<ExprBase>>
    T cast(const py::object& x)
    {
      if(!py::isinstance<ScalarExpr>(x)
        && !py::isinstance<VectorExpr>(x)
        && !py::isinstance<MatrixExpr>(x))
      { assert("unable to cast the provided py::object into a ExprBase"); }

      if(py::isinstance<ScalarExpr>(x))
        return x.cast<ScalarExpr>()->copy();

      else if(py::isinstance<VectorExpr>(x))
        return x.cast<VectorExpr>()->copy();

      else
        return x.cast<MatrixExpr>()->copy();
    }

    template<typename T>
      requires (std::is_same_v<AnalyticFunction<ScalarType>,T>
        || std::is_same_v<AnalyticFunction<VectorType>,T>
        || std::is_same_v<AnalyticFunction<MatrixType>,T>)
    bool is_instance(const py::object& x)
    {
      const py::object& f_ = x.attr("f");
      return f_ && py::isinstance<T>(f_);
    }

    template<typename T>
      requires (std::is_same_v<AnalyticFunction<ScalarType>,T>
        || std::is_same_v<AnalyticFunction<VectorType>,T>
        || std::is_same_v<AnalyticFunction<MatrixType>,T>)
    const T& cast(const py::object& x)
    {
      assert(is_instance<T>(x));
      const py::object& f_ = x.attr("f");
      return f_.cast<const T&>();
    }
}