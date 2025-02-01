/** 
 *  \file codac2_Interval2_operations.h
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

#include "codac2_Interval2.h"

namespace codac2
{
  /**
   * \brief Returns \f$[x]^2\f$
   * 
   * \param x interval value
   * \return the operation result
   */
  Interval2 sqr(const Interval2& x);

  /**
   * \brief Returns \f$\sqrt{[x]}\f$
   * 
   * \param x interval value
   * \return the operation result
   */
  Interval2 sqrt(const Interval2& x);

  /**
   * \brief Returns \f$[x]^n\f$, \f$n\in\mathbb{Z}\f$
   * 
   * \param x interval value
   * \param n integer power value
   * \return the operation result
   */
  Interval2 pow(const Interval2& x, int n);

  /**
   * \brief Returns \f$[x]^d\f$, \f$d\in\mathbb{R}\f$
   * 
   * \param x interval value
   * \param d real power value
   * \return the operation result
   */
  Interval2 pow(const Interval2& x, double d);

  /**
   * \brief Returns \f$[x]^{[y]}\f$, \f$y\in\mathbb{IR}\f$
   * 
   * \param x interval value
   * \param y interval power value
   * \return the operation result
   */
  Interval2 pow(const Interval2& x, const Interval2& y);

  /**
   * \brief Returns the n-th root: \f$\sqrt[n]{[x]}\f$
   * 
   * \param x interval value
   * \param n integer root
   * \return the operation result
   */
  Interval2 root(const Interval2& x, int n);

  /**
   * \brief Returns \f$\exp([x])\f$
   * 
   * \param x interval value
   * \return the operation result
   */
  Interval2 exp(const Interval2& x);

  /**
   * \brief Returns \f$\log([x])\f$
   * 
   * \param x interval value
   * \return the operation result
   */
  Interval2 log(const Interval2& x);

  /**
   * \brief Returns \f$\cos([x])\f$
   * 
   * \param x interval value
   * \return the operation result
   */
  Interval2 cos(const Interval2& x);

  /**
   * \brief Returns \f$\sin([x])\f$
   * 
   * \param x interval value
   * \return the operation result
   */
  Interval2 sin(const Interval2& x);

  /**
   * \brief Returns \f$\tan([x])\f$
   * 
   * \param x interval value
   * \return the operation result
   */
  Interval2 tan(const Interval2& x);

  /**
   * \brief Returns \f$\acos([x])\f$
   * 
   * \param x interval value
   * \return the operation result
   */
  Interval2 acos(const Interval2& x);

  /**
   * \brief Returns \f$\asin([x])\f$
   * 
   * \param x interval value
   * \return the operation result
   */
  Interval2 asin(const Interval2& x);

  /**
   * \brief Returns \f$\atan([x])\f$
   * 
   * \param x interval value
   * \return the operation result
   */
  Interval2 atan(const Interval2& x);

  /**
   * \brief Returns \f$\mathrm{arctan2}([y],[x])\f$
   * 
   * \param y interval value
   * \param x interval value
   * \return the operation result
   */
  Interval2 atan2(const Interval2& y, const Interval2& x);

  /**
   * \brief Returns \f$\cosh([x])\f$
   * 
   * \param x interval value
   * \return the operation result
   */
  Interval2 cosh(const Interval2& x);

  /**
   * \brief Returns \f$\sinh([x])\f$
   * 
   * \param x interval value
   * \return the operation result
   */
  Interval2 sinh(const Interval2& x);

  /**
   * \brief Returns \f$\tanh([x])\f$
   * 
   * \param x interval value
   * \return the operation result
   */
  Interval2 tanh(const Interval2& x);

  /**
   * \brief Returns \f$\acosh([x])\f$
   * 
   * \param x interval value
   * \return the operation result
   */
  Interval2 acosh(const Interval2& x);

  /**
   * \brief Returns \f$\asinh([x])\f$
   * 
   * \param x interval value
   * \return the operation result
   */
  Interval2 asinh(const Interval2& x);

  /**
   * \brief Returns \f$\atanh([x])\f$
   * 
   * \param x interval value
   * \return the operation result
   */
  Interval2 atanh(const Interval2& x);

  /**
   * \brief Returns \f$\mid[x]\mid = \left\{\mid x \mid, x\in[x]\right\}\f$
   * 
   * \param x interval value
   * \return the operation result
   */
  Interval2 abs(const Interval2& x);

  /**
   * \brief Returns \f$\min([x],[y])=\left\{\min(x,y), x\in[x], y\in[y]\right\}\f$
   * 
   * \param x interval value
   * \param y interval value
   * \return the operation result
   */
  Interval2 min(const Interval2& x, const Interval2& y);

  /**
   * \brief Returns \f$\max([x],[y])=\left\{\max(x,y), x\in[x], y\in[y]\right\}\f$
   * 
   * \param x interval value
   * \param y interval value
   * \return the operation result
   */
  Interval2 max(const Interval2& x, const Interval2& y);

  /**
   * \brief Returns \f$\sign([x])=\left[\left\{\sign(x), x\in[x]\right\}\right]\f$
   * 
   * \note By convention, \f$ 0\in[x] \Longrightarrow \sign([x])=[-1,1]\f$.
   * 
   * \param x interval value
   * \return the operation result
   */
  Interval2 sign(const Interval2& x);

  /**
   * \brief Returns the largest integer interval included in \f$[x]\f$
   * 
   * \param x interval value
   * \return the operation result
   */
  Interval2 integer(const Interval2& x);

  /**
   * \brief Returns floor of \f$[x]\f$
   * 
   * \param x interval value
   * \return the operation result
   */
  Interval2 floor(const Interval2& x);

  /**
   * \brief Returns ceil of \f$[x]\f$
   * 
   * \param x interval value
   * \return the operation result
   */
  Interval2 ceil(const Interval2& x);

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
  void bwd_add(const Interval2& y, Interval2& x1, Interval2& x2)
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
  void bwd_sub(const Interval2& y, Interval2& x1, Interval2& x2)
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
  void bwd_mul(const Interval2& y, Interval2& x1, Interval2& x2)
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
  void bwd_div(const Interval2& y, Interval2& x1, Interval2& x2)
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
  void bwd_sqr(const Interval2& y, Interval2& x)
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
  void bwd_sqrt(const Interval2& y, Interval2& x)
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
  void bwd_pow(const Interval2& y, Interval2& x, int p)
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
  void bwd_pow(const Interval2& y, Interval2& x, Interval2& p)
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
  void bwd_root(const Interval2& y, Interval2& x, int p)
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
  void bwd_exp(const Interval2& y, Interval2& x)
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
  void bwd_log(const Interval2& y, Interval2& x)
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
  void bwd_cos(const Interval2& y, Interval2& x)
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
  void bwd_sin(const Interval2& y, Interval2& x)
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
  void bwd_tan(const Interval2& y, Interval2& x)
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
  void bwd_acos(const Interval2& y, Interval2& x)
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
  void bwd_asin(const Interval2& y, Interval2& x)
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
  void bwd_atan(const Interval2& y, Interval2& x)
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
  void bwd_atan2(const Interval2& y, Interval2& x1, Interval2& x2)
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
  void bwd_cosh(const Interval2& y, Interval2& x)
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
  void bwd_sinh(const Interval2& y, Interval2& x)
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
  void bwd_tanh(const Interval2& y, Interval2& x)
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
  void bwd_acosh(const Interval2& y, Interval2& x)
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
  void bwd_asinh(const Interval2& y, Interval2& x)
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
  void bwd_atanh(const Interval2& y, Interval2& x)
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
  void bwd_abs(const Interval2& y, Interval2& x)
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
  void bwd_min(const Interval2& y, Interval2& x1, Interval2& x2)
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
  void bwd_max(const Interval2& y, Interval2& x1, Interval2& x2)
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
  void bwd_sign(const Interval2& y, Interval2& x)
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
  void bwd_floor(const Interval2& y, Interval2& x)
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
  void bwd_ceil(const Interval2& y, Interval2& x)
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
  void bwd_imod(Interval2& x1, Interval2& x2, double p)
  {
    gaol::bwd_imod(x1,x2,p);
  }
#endif
}

#include "codac2_Interval2_operations_impl.h"