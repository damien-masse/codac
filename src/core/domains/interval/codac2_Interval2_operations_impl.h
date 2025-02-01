/** 
 *  \file codac2_Interval2_operations_impl.h
 * 
 *  This class reuses several functions developed for ibex::Interval2.
 *  See ibex::Interval2 (IBEX lib, main author: Gilles Chabert)
 *    https://ibex-lib.readthedocs.io
 *  
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Gilles Chabert, Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

// Inline functions

namespace codac2
{
  inline Interval2 sqr(const Interval2& x)
  {
    return gaol::sqr(x);
  }

  inline Interval2 sqrt(const Interval2& x)
  {
    return gaol::sqrt(x);
  }

  inline Interval2 pow(const Interval2& x, int n)
  {
    return gaol::pow(x,n);
  }

  inline Interval2 pow(const Interval2& x, double d)
  {
    if(d == -oo || d == oo)
      return Interval2::empty();

    else
      return gaol::pow(x,d);
  }

  inline Interval2 pow(const Interval2& x, const Interval2& y)
  {
    return gaol::pow(x,y);
  }

  inline Interval2 root(const Interval2& x, int n)
  {
    // Get the root of the positive part (gaol does
    // not consider negative values to be in the definition
    // domain of the root function)

    gaol::interval res = gaol::nth_root(x, n>=0 ? n : -n);

    if(n%2 == 1 && x.lb() < 0)
      res |= -gaol::nth_root(-x, n >= 0 ? n : -n);

    if(n < 0)
      res = 1.0/res;

    return res;
  }

  inline Interval2 exp(const Interval2& x)
  {
    return gaol::exp(x);
  }

  inline Interval2 log(const Interval2& x)
  {
    if(x.ub() <= 0) // gaol returns (-oo,-DBL_MAX) if x.ub()==0, instead of empty set
      return Interval2::empty();

    else
      return gaol::log(x);
  }

  inline Interval2 cos(const Interval2& x)
  {
    return gaol::cos(x);
  }

  inline Interval2 sin(const Interval2& x)
  {
    return gaol::sin(x);
  }

  inline Interval2 tan(const Interval2& x)
  {
    return gaol::tan(x);
  }

  inline Interval2 acos(const Interval2& x)
  {
    return gaol::acos(x);
  }

  inline Interval2 asin(const Interval2& x)
  {
    return gaol::asin(x);
  }

  inline Interval2 atan(const Interval2& x)
  {
    return gaol::atan(x);
  }

  inline Interval2 atan2(const Interval2& y, const Interval2& x)
  {
    if(y.is_empty() || x.is_empty())
      return Interval2::empty();

    // We handle the special case x=[0,0] separately
    else if(x == Interval2::zero())
    {
      if(y.lb() >= 0)
      {
        if(y.ub() == 0)
          return Interval2::empty(); // atan2(0,0) is undefined
        else
          return Interval2::half_pi();
      }

      else if(y.ub() <= 0)
        return -Interval2::half_pi();

      else
        return Interval2(-1,1)*Interval2::half_pi();
    }

    else if(x.lb() >= 0)
      return atan(y/x); // now, x.ub()>0 -> atan does not give an empty set

    else if(x.ub() <= 0)
    {
      if(y.lb() >= 0)
        return atan(y/x) + Interval2::pi(); // x.lb()<0
      else if(y.ub() < 0)
        return atan(y/x) - Interval2::pi();
      else
        return Interval2(-1,1)*Interval2::pi();
    }

    else
    {
      if(y.lb() >= 0)
        return atan(y/x.ub()) | (atan(y/x.lb()) + Interval2::pi());

      else if(y.ub() <= 0)
      {
        if(x.lb() != -oo)
        {
          if(x.ub() != oo)
            return (atan(y/x.lb())-Interval2::pi()) | atan(y/x.ub());
          else
            return (atan(y/x.lb())-Interval2::pi()) | Interval2::zero();
        }

        else
        {
          if(x.ub() != oo)
            return (-Interval2::pi()) | atan(y/x.ub());
          else
            return -Interval2::pi() | Interval2::zero();
        }
      }

      else
        return Interval2(-1,1)*Interval2::pi();
    }
  }

  inline Interval2 cosh(const Interval2& x)
  {
    if(x.is_unbounded()) 
      return Interval2(gaol::cosh(x).left(),oo);
    else
      return gaol::cosh(x);
  }

  inline Interval2 sinh(const Interval2& x)
  {
    return gaol::sinh(x);
  }

  inline Interval2 tanh(const Interval2& x)
  {
    return gaol::tanh(x);
  }

  inline Interval2 acosh(const Interval2& x)
  {
    return gaol::acosh(x);
  }

  inline Interval2 asinh(const Interval2& x)
  {
    if(x.is_empty())
      return Interval2::empty();

    else if(x.lb() >= 0)
      return gaol::asinh(x);

    else if(x.ub() <= 0)
      return -gaol::asinh(-x);

    else
      return {
        -gaol::asinh(gaol::interval(0,-x.lb())).right(),
        gaol::asinh(gaol::interval(0,x.ub())).right()
      };
  }

  inline Interval2 atanh(const Interval2& x)
  {
    return gaol::atanh(x);
  }

  inline Interval2 abs(const Interval2& x)
  {
    return gaol::abs(x);
  }

  inline Interval2 min(const Interval2& x, const Interval2& y)
  {
    return gaol::min(x,y);
  }

  inline Interval2 max(const Interval2& x, const Interval2& y)
  {
    return gaol::max(x,y);
  }

  inline Interval2 sign(const Interval2& x)
  {
    return x.ub() < 0 ? -1. : x.lb() > 0 ? 1. : Interval2(-1.,1.);
  }

  inline Interval2 integer(const Interval2& x)
  {
    return gaol::integer(x);
  }

  inline Interval2 floor(const Interval2& x)
  {
    return gaol::floor(x);
  }

  inline Interval2 ceil(const Interval2& x)
  {
    return gaol::ceil(x);
  }

#if 0
  /**
   * \brief Computes the backward (reverse) addition
   * 
   * \note Contracts \f$[x_1]\f$, \f$[x_2]\f$ as: 
   * \f$([x_1],[x_2])\cap\{(x_1,x_2)\in[x_1]\times[x_2]\ \mid\ \exists y\in[y],\ y=x_1+x_2\}\f$.
   * 
   * \param y interval value (result of the forward operation)
   * \param x1 prior value for \f$[x_1]\f$, may be contracted
   * \param x2 prior value for \f$[x_2]\f$, may be contracted
   */
  inline void bwd_add(const Interval2& y, Interval2& x1, Interval2& x2)
  {
    gaol::bwd_add(y,x1,x2);
  }

  /**
   * \brief Computes the backward (reverse) substraction
   * 
   * \note Contracts \f$[x_1]\f$, \f$[x_2]\f$ as: 
   * \f$([x_1],[x_2])\cap\{(x_1,x_2)\in[x_1]\times[x_2]\ \mid\ \exists y\in[y],\ y=x_1-x_2\}\f$.
   * 
   * \param y interval value (result of the forward operation)
   * \param x1 prior value for \f$[x_1]\f$, may be contracted
   * \param x2 prior value for \f$[x_2]\f$, may be contracted
   */
  inline void bwd_sub(const Interval2& y, Interval2& x1, Interval2& x2)
  {
    gaol::bwd_sub(y,x1,x2);
  }

  /**
   * \brief Computes the backward (reverse) multiplication
   * 
   * \note Contracts \f$[x_1]\f$, \f$[x_2]\f$ as: 
   * \f$([x_1],[x_2])\cap\{(x_1,x_2)\in[x_1]\times[x_2]\ \mid\ \exists y\in[y],\ y=x_1*x_2\}\f$.
   * 
   * \param y interval value (result of the forward operation)
   * \param x1 prior value for \f$[x_1]\f$, may be contracted
   * \param x2 prior value for \f$[x_2]\f$, may be contracted
   */
  inline void bwd_mul(const Interval2& y, Interval2& x1, Interval2& x2)
  {
    gaol::bwd_mul(y,x1,x2);
  }

  /**
   * \brief Computes the backward (reverse) division
   * 
   * \note Contracts \f$[x_1]\f$, \f$[x_2]\f$ as: 
   * \f$([x_1],[x_2])\cap\{(x_1,x_2)\in[x_1]\times[x_2]\ \mid\ \exists y\in[y],\ y=x_1/x_2\}\f$.
   * 
   * \param y interval value (result of the forward operation)
   * \param x1 prior value for \f$[x_1]\f$, may be contracted
   * \param x2 prior value for \f$[x_2]\f$, may be contracted
   */
  inline void bwd_div(const Interval2& y, Interval2& x1, Interval2& x2)
  {
    gaol::bwd_div(y,x1,x2);
  }

  /**
   * \brief Computes the backward (reverse) squared operation
   * 
   * \note Contracts \f$[x]\f$ as: 
   * \f$[x]\cap\{x\in[x]\ \mid\ \exists y\in[y],\ y=x^2\}\f$.
   * 
   * \param y interval value (result of the forward operation)
   * \param x prior value for \f$[x]\f$, may be contracted
   */
  inline void bwd_sqr(const Interval2& y, Interval2& x)
  {
    gaol::bwd_sqr(y,x);
  }

  /**
   * \brief Computes the backward (reverse) squared-root operation
   * 
   * \note Contracts \f$[x]\f$ as: 
   * \f$[x]\cap\{x\in[x]\ \mid\ \exists y\in[y],\ y=\sqrt{x}\}\f$.
   * 
   * \param y interval value (result of the forward operation)
   * \param x prior value for \f$[x]\f$, may be contracted
   */
  inline void bwd_sqrt(const Interval2& y, Interval2& x)
  {
    gaol::bwd_sqrt(y,x);
  }

  /**
   * \brief Computes the backward (reverse) power operation
   * 
   * \note Contracts \f$[x]\f$ as: 
   * \f$[x]\cap\{x\in[x]\ \mid\ \exists y\in[y],\ y=x^p\}\f$.
   * 
   * \param y interval value (result of the forward operation)
   * \param x prior value for \f$[x]\f$, may be contracted
   * \param p power integer value
   */
  inline void bwd_pow(const Interval2& y, Interval2& x, int p)
  {
    gaol::bwd_pow(y,p,x);
  }

  /**
   * \brief Computes the backward (reverse) power operation
   * 
   * \note Contracts \f$[x]\f$, \f$[p]\f$ as: 
   * \f$([x],[p])\cap\{(x,p)\in[x]\times[p]\ \mid\ \exists y\in[y],\ y=x^{p}\}\f$.
   * 
   * \param y interval value (result of the forward operation)
   * \param x prior value for \f$[x]\f$, may be contracted
   * \param p prior value for \f$[p]\f$, may be contracted
   */
  inline void bwd_pow(const Interval2& y, Interval2& x, Interval2& p)
  {
    assert(p.is_degenerated() && "bwd_power(y,x1,x2) (with x1 and x2 intervals) not implemented yet with Gaol");
    gaol::bwd_pow(y,p.mid(),x);
  }

  /**
   * \brief Computes the backward (reverse) root operation
   * 
   * \note Contracts \f$[x]\f$ as: 
   * \f$[x]\cap\{x\in[x]\ \mid\ \exists y\in[y],\ y=\sqrt[p]{x}\}\f$.
   * 
   * \param y interval value (result of the forward operation)
   * \param x prior value for \f$[x]\f$, may be contracted
   * \param p root integer value
   */
  inline void bwd_root(const Interval2& y, Interval2& x, int p)
  {
    gaol::bwd_root(y,p,x);
  }

  /**
   * \brief Computes the backward (reverse) exponential operation
   * 
   * \note Contracts \f$[x]\f$ as: 
   * \f$[x]\cap\{x\in[x]\ \mid\ \exists y\in[y],\ y=\exp(x)\}\f$.
   * 
   * \param y interval value (result of the forward operation)
   * \param x prior value for \f$[x]\f$, may be contracted
   */
  inline void bwd_exp(const Interval2& y, Interval2& x)
  {
    gaol::bwd_exp(y,x);
  }

  /**
   * \brief Computes the backward (reverse) logarithmic operation
   * 
   * \note Contracts \f$[x]\f$ as: 
   * \f$[x]\cap\{x\in[x]\ \mid\ \exists y\in[y],\ y=\log(x)\}\f$.
   * 
   * \param y interval value (result of the forward operation)
   * \param x prior value for \f$[x]\f$, may be contracted
   */
  inline void bwd_log(const Interval2& y, Interval2& x)
  {
    gaol::bwd_log(y,x);
  }

  /**
   * \brief Computes the backward (reverse) cosine operation
   * 
   * \note Contracts \f$[x]\f$ as: 
   * \f$[x]\cap\{x\in[x]\ \mid\ \exists y\in[y],\ y=\cos(x)\}\f$.
   * 
   * \param y interval value (result of the forward operation)
   * \param x prior value for \f$[x]\f$, may be contracted
   */
  inline void bwd_cos(const Interval2& y, Interval2& x)
  {
    gaol::bwd_cos(y,x);
  }

  /**
   * \brief Computes the backward (reverse) sine operation
   * 
   * \note Contracts \f$[x]\f$ as: 
   * \f$[x]\cap\{x\in[x]\ \mid\ \exists y\in[y],\ y=\sin(x)\}\f$.
   * 
   * \param y interval value (result of the forward operation)
   * \param x prior value for \f$[x]\f$, may be contracted
   */
  inline void bwd_sin(const Interval2& y, Interval2& x)
  {
    gaol::bwd_sin(y,x);
  }

  /**
   * \brief Computes the backward (reverse) tangent operation
   * 
   * \note Contracts \f$[x]\f$ as: 
   * \f$[x]\cap\{x\in[x]\ \mid\ \exists y\in[y],\ y=\tan(x)\}\f$.
   * 
   * \param y interval value (result of the forward operation)
   * \param x prior value for \f$[x]\f$, may be contracted
   */
  inline void bwd_tan(const Interval2& y, Interval2& x)
  {
    gaol::bwd_tan(y,x);
  }

  /**
   * \brief Computes the backward (reverse) arccos operation
   * 
   * \note Contracts \f$[x]\f$ as: 
   * \f$[x]\cap\{x\in[x]\ \mid\ \exists y\in[y],\ y=\acos(x)\}\f$.
   * 
   * \param y interval value (result of the forward operation)
   * \param x prior value for \f$[x]\f$, may be contracted
   */
  inline void bwd_acos(const Interval2& y, Interval2& x)
  {
    gaol::bwd_acos(y,x);
  }

  /**
   * \brief Computes the backward (reverse) arcsin operation
   * 
   * \note Contracts \f$[x]\f$ as: 
   * \f$[x]\cap\{x\in[x]\ \mid\ \exists y\in[y],\ y=\asin(x)\}\f$.
   * 
   * \param y interval value (result of the forward operation)
   * \param x prior value for \f$[x]\f$, may be contracted
   */
  inline void bwd_asin(const Interval2& y, Interval2& x)
  {
    gaol::bwd_asin(y,x);
  }

  /**
   * \brief Computes the backward (reverse) arctan operation
   * 
   * \note Contracts \f$[x]\f$ as: 
   * \f$[x]\cap\{x\in[x]\ \mid\ \exists y\in[y],\ y=\atan(x)\}\f$.
   * 
   * \param y interval value (result of the forward operation)
   * \param x prior value for \f$[x]\f$, may be contracted
   */
  inline void bwd_atan(const Interval2& y, Interval2& x)
  {
    gaol::bwd_atan(y,x);
  }

  /**
   * \brief Computes the backward (reverse) arctan2 operation
   * 
   * \note Contracts \f$[x_1]\f$, \f$[x_2]\f$ as: 
   * \f$([x_1],[x_2])\cap\{(x_1,x_2)\in[x_1]\times[x_2]\ \mid\ \exists y\in[y],\ y=\atan2(x_1,x_2)\}\f$.
   * 
   * \param y interval value (result of the forward operation)
   * \param x1 prior value for \f$[x_1]\f$, may be contracted
   * \param x2 prior value for \f$[x_2]\f$, may be contracted
   */
  inline void bwd_atan2(const Interval2& y, Interval2& x1, Interval2& x2)
  {
    gaol::bwd_atan2(y,x1,x2);
  }

  /**
   * \brief Computes the backward (reverse) hyperbolic cosine operation
   * 
   * \note Contracts \f$[x]\f$ as: 
   * \f$[x]\cap\{x\in[x]\ \mid\ \exists y\in[y],\ y=\cosh(x)\}\f$.
   * 
   * \param y interval value (result of the forward operation)
   * \param x prior value for \f$[x]\f$, may be contracted
   */
  inline void bwd_cosh(const Interval2& y, Interval2& x)
  {
    gaol::bwd_cosh(y,x);
  }

  /**
   * \brief Computes the backward (reverse) hyperbolic sine operation
   * 
   * \note Contracts \f$[x]\f$ as: 
   * \f$[x]\cap\{x\in[x]\ \mid\ \exists y\in[y],\ y=\sinh(x)\}\f$.
   * 
   * \param y interval value (result of the forward operation)
   * \param x prior value for \f$[x]\f$, may be contracted
   */
  inline void bwd_sinh(const Interval2& y, Interval2& x)
  {
    gaol::bwd_sinh(y,x);
  }

  /**
   * \brief Computes the backward (reverse) hyperbolic tangent operation
   * 
   * \note Contracts \f$[x]\f$ as: 
   * \f$[x]\cap\{x\in[x]\ \mid\ \exists y\in[y],\ y=\tanh(x)\}\f$.
   * 
   * \param y interval value (result of the forward operation)
   * \param x prior value for \f$[x]\f$, may be contracted
   */
  inline void bwd_tanh(const Interval2& y, Interval2& x)
  {
    gaol::bwd_tanh(y,x);
  }

  /**
   * \brief Computes the backward (reverse) hyperbolic arccos operation
   * 
   * \note Contracts \f$[x]\f$ as: 
   * \f$[x]\cap\{x\in[x]\ \mid\ \exists y\in[y],\ y=\acosh(x)\}\f$.
   * 
   * \param y interval value (result of the forward operation)
   * \param x prior value for \f$[x]\f$, may be contracted
   */
  inline void bwd_acosh(const Interval2& y, Interval2& x)
  {
    gaol::bwd_acosh(y,x);
  }

  /**
   * \brief Computes the backward (reverse) hyperbolic arcsin operation
   * 
   * \note Contracts \f$[x]\f$ as: 
   * \f$[x]\cap\{x\in[x]\ \mid\ \exists y\in[y],\ y=\asinh(x)\}\f$.
   * 
   * \param y interval value (result of the forward operation)
   * \param x prior value for \f$[x]\f$, may be contracted
   */
  inline void bwd_asinh(const Interval2& y, Interval2& x)
  {
    gaol::bwd_asinh(y,x);
  }

  /**
   * \brief Computes the backward (reverse) hyperbolic arctan operation
   * 
   * \note Contracts \f$[x]\f$ as: 
   * \f$[x]\cap\{x\in[x]\ \mid\ \exists y\in[y],\ y=\atanh(x)\}\f$.
   * 
   * \param y interval value (result of the forward operation)
   * \param x prior value for \f$[x]\f$, may be contracted
   */
  inline void bwd_atanh(const Interval2& y, Interval2& x)
  {
    gaol::bwd_atanh(y,x);
  }

  /**
   * \brief Computes the backward (reverse) absolute-value operation
   * 
   * \note Contracts \f$[x]\f$ as: 
   * \f$[x]\cap\{x\in[x]\ \mid\ \exists y\in[y],\ y=\mid x\mid\}\f$.
   * 
   * \param y interval value (result of the forward operation)
   * \param x prior value for \f$[x]\f$, may be contracted
   */
  inline void bwd_abs(const Interval2& y, Interval2& x)
  {
    gaol::bwd_abs(y,x);
  }

  /**
   * \brief Computes the backward (reverse) of the max operation
   * 
   * \note Contracts \f$[x_1]\f$, \f$[x_2]\f$ as: 
   * \f$([x_1],[x_2])\cap\{(x_1,x_2)\in[x_1]\times[x_2]\ \mid\ \exists y\in[y],\ y=\max(x_1,x_2)\}\f$.
   * 
   * \param y interval value (result of the forward operation)
   * \param x1 prior value for \f$[x_1]\f$, may be contracted
   * \param x2 prior value for \f$[x_2]\f$, may be contracted
   */
  inline void bwd_min(const Interval2& y, Interval2& x1, Interval2& x2)
  {
    gaol::bwd_min(y,x1,x2);
  }

  /**
   * \brief Computes the backward (reverse) of the min operation
   * 
   * \note Contracts \f$[x_1]\f$, \f$[x_2]\f$ as: 
   * \f$([x_1],[x_2])\cap\{(x_1,x_2)\in[x_1]\times[x_2]\ \mid\ \exists y\in[y],\ y=\min(x_1,x_2)\}\f$.
   * 
   * \param y interval value (result of the forward operation)
   * \param x1 prior value for \f$[x_1]\f$, may be contracted
   * \param x2 prior value for \f$[x_2]\f$, may be contracted
   */
  inline void bwd_max(const Interval2& y, Interval2& x1, Interval2& x2)
  {
    gaol::bwd_max(y,x1,x2);
  }

  /**
   * \brief Computes the backward (reverse) sign operation
   * 
   * \note Contracts \f$[x]\f$ as: 
   * \f$[x]\cap\{x\in[x]\ \mid\ \exists y\in[y],\ y=\sign(x)\}\f$.
   * 
   * \param y interval value (result of the forward operation)
   * \param x prior value for \f$[x]\f$, may be contracted
   */
  inline void bwd_sign(const Interval2& y, Interval2& x)
  {
    gaol::bwd_sign(y,x);
  }

  /**
   * \brief Computes the backward (reverse) floor operation
   * 
   * \note Contracts \f$[x]\f$ as: 
   * \f$[x]\cap\{x\in[x]\ \mid\ \exists y\in[y],\ y=\floor(x)\}\f$.
   * 
   * \param y interval value (result of the forward operation)
   * \param x prior value for \f$[x]\f$, may be contracted
   */
  inline void bwd_floor(const Interval2& y, Interval2& x)
  {
    gaol::bwd_floor(y,x);
  }

  /**
   * \brief Computes the backward (reverse) ceil operation
   * 
   * \note Contracts \f$[x]\f$ as: 
   * \f$[x]\cap\{x\in[x]\ \mid\ \exists y\in[y],\ y=\ceil(x)\}\f$.
   * 
   * \param y interval value (result of the forward operation)
   * \param x prior value for \f$[x]\f$, may be contracted
   */
  inline void bwd_ceil(const Interval2& y, Interval2& x)
  {
    gaol::bwd_ceil(y,x);
  }

  /**
   * \brief Contract \f$[x_1]\f$ and \f$[x_2]\f$ w.r.t. the fact that they are equivalent modulo the period \f$p\f$
   * 
   * \param x1 prior value for \f$[x_1]\f$, may be contracted
   * \param x2 prior value for \f$[x_2]\f$, may be contracted
   * \param p period value
   */
  inline void bwd_imod(Interval2& x1, Interval2& x2, double p)
  {
    gaol::bwd_imod(x1,x2,p);
  }

  #endif

}