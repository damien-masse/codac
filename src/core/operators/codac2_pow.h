/** 
 *  \file codac2_pow.h
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou, Damien Massé
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include "codac2_Interval.h"
#include "codac2_AnalyticType.h"
#include "codac2_AnalyticExprWrapper.h"
#include "codac2_sqr.h"

namespace codac2
{
  struct PowOp
  {
    template<typename X1,typename X2>
    static std::string str(const X1& x1, const X2& x2)
    {
      return "(" + x1->str(!x1->is_str_leaf()) + ")^(" + x2->str(!x2->is_str_leaf()) + ")";
    }

    template<typename X1, typename X2>
    static std::pair<Index,Index> output_shape([[maybe_unused]] const X1& s1, [[maybe_unused]] const X2& s2)
    {
      return {1,1};
    }
    
    static Interval fwd(const Interval& x1, const Interval& x2);
    static ScalarType fwd_natural(const ScalarType& x1, const ScalarType& x2);
    static ScalarType fwd_centered(const ScalarType& x1, const ScalarType& x2);
    static void bwd(const Interval& y, Interval& x1, int x2);
    static void bwd(const Interval& y, Interval& x1, Interval& x2);
  };

  // Analytic operator
  // The following function can be used to build analytic expressions.

  inline ScalarExpr
  pow(const ScalarExpr& x1, const ScalarExpr& x2)
  {
    auto x2_as_const_value = std::dynamic_pointer_cast<ConstValueExpr<ScalarType>>(x2->copy());
    if(x2_as_const_value && x2_as_const_value->value() == 2.)
    {
      // SqrOp may be better than the generic PowOp with a simple power of 2
      return { std::make_shared<AnalyticOperationExpr<SqrOp,ScalarType,ScalarType>>(x1) };
    }

    return { std::make_shared<AnalyticOperationExpr<PowOp,ScalarType,ScalarType,ScalarType>>(x1,x2) };
  }

  inline ScalarExpr
  operator^(const ScalarExpr& x1, const ScalarExpr& x2)
  {
    std::cout
      << "Operator '^' is discouraged in Codac. Use pow(x,y), or sqr(x) instead."
      << std::endl;
    return pow(x1,x2);
  }

  // Inline functions

  inline Interval PowOp::fwd(const Interval& x1, const Interval& x2)
  {
    return pow(x1,x2);
  }

  inline ScalarType PowOp::fwd_natural(const ScalarType& x1, const ScalarType& x2)
  {
    bool x2isint = x2.a.is_integer();
    bool x2positive = x2.a.lb()>=0.0;
    return {
      fwd(x1.a, x2.a),
      x1.def_domain && x2.def_domain 
        && (x2isint || x1.a.lb()>=0.0)
        && (x2positive || !x1.a.contains(0.0))
    };
  }

  inline ScalarType PowOp::fwd_centered(const ScalarType& x1, const ScalarType& x2)
  {
    if(centered_form_not_available_for_args(x1,x2))
      return fwd_natural(x1,x2);

    bool x2isint = x2.a.is_integer();
    bool x2positive = x2.a.lb()>=0.0;
    IntervalMatrix d(1,x1.da.size());
    
    if (x2.a.is_degenerated())
    { 
      // To avoid calling log(x1.a), as it would return emptyset if x1.a<=0
      for(Index i = 0 ; i < d.size() ; i++)
        d(0,i) = x2.a*x1.da(0,i)*pow(x1.a,x2.a-1.);
    }

    else
    {
      for(Index i = 0 ; i < d.size() ; i++)
      {
        d(0,i) = x2.a*x1.da(0,i)*pow(x1.a,x2.a-1.) 
                + x2.da(0,i)*log(x1.a)*pow(x1.a,x2.a); 
			  // Not good when x1 close to 0
      }
    }

    return {
      fwd(x1.m, x2.m),
      fwd(x1.a, x2.a),
      d,
      x1.def_domain && x2.def_domain 
        && (x2isint || x1.a.lb()>=0.0)
        && (x2positive || !x1.a.contains(0.0))
    };
  }

  inline void PowOp::bwd(const Interval& y, Interval& x1, int x2)
  {
    // The content of this function comes from the IBEX library.
    // See ibex::Interval (IBEX lib, main author: Gilles Chabert)
    //   https://ibex-lib.readthedocs.io

    if(x2 % 2 == 0)
    {
      Interval proj = root(y, x2);
      Interval pos_proj = proj & x1;
      Interval neg_proj = (-proj) & x1;
      x1 = pos_proj | neg_proj;
    }

    else
      x1 &= root(y, x2);
  }

  inline void PowOp::bwd(const Interval& y, Interval& x1, Interval& x2)
  {
    if (x2.is_degenerated() && x2 == (int)(x2.mid())) {
       if (x2.mid()==0.0) {
         if (y.contains(1.0)) return;
         x1.set_empty(); x2.set_empty();
         return;
       }
       PowOp::bwd(y, x1, (int)(x2.mid()));
       if (x1.is_empty()) x2.set_empty();
    } else {
       assert_release(x1.lb()>0.0 && "PowOp::bwd(y,x1,x2) (with x1 not positive and x2 interval) not implemented yet with Gaol");
       Interval ly = log(y);
       Interval lx = log(x1);
       Interval x2p = x2 & (ly / (lx & Interval(-oo,0.0)));
       Interval x1p = (x2p.is_empty() ? Interval::empty() :
			(x1 & exp(ly/x2p)));
       Interval x2n = x2 & (ly / (lx & Interval(0.0,+oo)));
       Interval x1n = (x2n.is_empty() ? Interval::empty() :
			(x1 & exp(ly/x2n)));
       x2 = x2p | x2n;
       x1 = x1p | x1n;
    }         
  }
}
