/** 
 *  \file codac2_analytic_values.h
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include "codac2_Interval.h"
#include "codac2_Vector.h"
#include "codac2_Matrix.h"
#include "codac2_IntervalVector.h"
#include "codac2_IntervalMatrix.h"

namespace codac2
{
  struct OpValueBase
  {
    virtual ~OpValueBase() = default;
  };

  template<typename S, typename T, typename M>
  struct OpValue : public OpValueBase
  {
    using Scalar = S;
    using Domain = T;

    T m;
    T a;
    M da;
    bool def_domain;

    OpValue() = delete;

    OpValue(const T& m_, const T& a_, const M& da_, bool def_domain_)
      : m(m_), a(a_), da(da_), def_domain(def_domain_)
    { }

    OpValue<S,T,M>& operator&=(const OpValue<S,T,M>& x)
    {
      a &= x.a;
      // restore this? da &= x.da;
      def_domain &= x.def_domain;
      return *this;
    }
  };

  using ScalarOpValue = OpValue<double,Interval,IntervalMatrix>;
  using VectorOpValue = OpValue<Vector,IntervalVector,IntervalMatrix>;
  using MatrixOpValue = OpValue<Matrix,IntervalMatrix,IntervalMatrix>;

  template<typename T>
  struct ArgWrapper
  { };

  template<>
  struct ArgWrapper<int> {
    using Domain = ScalarOpValue;
  };

  template<>
  struct ArgWrapper<Index> {
    using Domain = ScalarOpValue;
  };

  template<>
  struct ArgWrapper<double> {
    using Domain = ScalarOpValue;
  };

  template<>
  struct ArgWrapper<Interval> {
    using Domain = ScalarOpValue;
  };

  template<>
  struct ArgWrapper<Vector> {
    using Domain = VectorOpValue;
  };

  template<>
  struct ArgWrapper<IntervalVector> {
    using Domain = VectorOpValue;
  };

  template<>
  struct ArgWrapper<Matrix> {
    using Domain = MatrixOpValue;
  };

  template<>
  struct ArgWrapper<IntervalMatrix> {
    using Domain = MatrixOpValue;
  };
}