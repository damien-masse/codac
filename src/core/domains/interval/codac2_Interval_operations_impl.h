/** 
 *  \file codac2_Interval_operations_impl.h
 * 
 *  This class reuses several functions developed for ibex::Interval.
 *  See ibex::Interval (IBEX lib, main author: Gilles Chabert)
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
  inline Interval sqr(const Interval& x)
  {
    return gaol::sqr(x);
  }

  inline Interval sqrt(const Interval& x)
  {
    return gaol::sqrt(x);
  }

  inline Interval pow(const Interval& x, int n)
  {
    return gaol::pow(x,n);
  }

  inline Interval pow(const Interval& x, double d)
  {
    if(d == -oo || d == oo)
      return Interval::empty();

    else
      return gaol::pow(x,d);
  }

  inline Interval pow(const Interval& x, const Interval& y)
  {
    return gaol::pow(x,y);
  }

  inline Interval root(const Interval& x, int n)
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

  inline Interval exp(const Interval& x)
  {
    return gaol::exp(x);
  }

  inline Interval log(const Interval& x)
  {
    if(x.ub() <= 0) // gaol returns (-oo,-DBL_MAX) if x.ub()==0, instead of empty set
      return Interval::empty();

    else
      return gaol::log(x);
  }

  inline Interval cos(const Interval& x)
  {
    return gaol::cos(x);
  }

  inline Interval sin(const Interval& x)
  {
    return gaol::sin(x);
  }

  inline Interval tan(const Interval& x)
  {
    return gaol::tan(x);
  }

  inline Interval acos(const Interval& x)
  {
    return gaol::acos(x);
  }

  inline Interval asin(const Interval& x)
  {
    return gaol::asin(x);
  }

  inline Interval atan(const Interval& x)
  {
    return gaol::atan(x);
  }

  inline Interval atan2(const Interval& y, const Interval& x)
  {
    if(y.is_empty() || x.is_empty())
      return Interval::empty();

    // We handle the special case x=[0,0] separately
    else if(x == Interval::zero())
    {
      if(y.lb() >= 0)
      {
        if(y.ub() == 0)
          return Interval::empty(); // atan2(0,0) is undefined
        else
          return Interval::half_pi();
      }

      else if(y.ub() <= 0)
        return -Interval::half_pi();

      else
        return Interval(-1,1)*Interval::half_pi();
    }

    else if(x.lb() >= 0)
      return atan(y/x); // now, x.ub()>0 -> atan does not give an empty set

    else if(x.ub() <= 0)
    {
      if(y.lb() >= 0)
        return atan(y/x) + Interval::pi(); // x.lb()<0
      else if(y.ub() < 0)
        return atan(y/x) - Interval::pi();
      else
        return Interval(-1,1)*Interval::pi();
    }

    else
    {
      if(y.lb() >= 0)
        return atan(y/x.ub()) | (atan(y/x.lb()) + Interval::pi());

      else if(y.ub() <= 0)
      {
        if(x.lb() != -oo)
        {
          if(x.ub() != oo)
            return (atan(y/x.lb())-Interval::pi()) | atan(y/x.ub());
          else
            return (atan(y/x.lb())-Interval::pi()) | Interval::zero();
        }

        else
        {
          if(x.ub() != oo)
            return (-Interval::pi()) | atan(y/x.ub());
          else
            return -Interval::pi() | Interval::zero();
        }
      }

      else
        return Interval(-1,1)*Interval::pi();
    }
  }

  inline Interval cosh(const Interval& x)
  {
    if(x.is_unbounded()) 
      return Interval(gaol::cosh(x).left(),oo);
    else
      return gaol::cosh(x);
  }

  inline Interval sinh(const Interval& x)
  {
    return gaol::sinh(x);
  }

  inline Interval tanh(const Interval& x)
  {
    return gaol::tanh(x);
  }

  inline Interval acosh(const Interval& x)
  {
    return gaol::acosh(x);
  }

  inline Interval asinh(const Interval& x)
  {
    if(x.is_empty())
      return Interval::empty();

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

  inline Interval atanh(const Interval& x)
  {
    return gaol::atanh(x);
  }

  inline Interval abs(const Interval& x)
  {
    return gaol::abs(x);
  }

  inline Interval min(const Interval& x, const Interval& y)
  {
    return gaol::min(x,y);
  }

  inline Interval max(const Interval& x, const Interval& y)
  {
    return gaol::max(x,y);
  }

  inline Interval sign(const Interval& x)
  {
    return x.ub() < 0 ? -1. : x.lb() > 0 ? 1. : Interval(-1.,1.);
  }

  inline Interval integer(const Interval& x)
  {
    return gaol::integer(x);
  }

  inline Interval floor(const Interval& x)
  {
    return gaol::floor(x);
  }

  inline Interval ceil(const Interval& x)
  {
    return gaol::ceil(x);
  }
}