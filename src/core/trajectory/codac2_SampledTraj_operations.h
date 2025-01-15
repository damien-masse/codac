/** 
 *  \file codac2_SampledTraj_operations.h
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include "codac2_SampledTraj.h"
#include "codac2_math.h"
#include "codac2_template_tools.h"

namespace codac2
{
  #define macro_unary_traj(f) \
  { \
    auto y = x1; \
    for(auto it = y.begin() ; it != y.end() ; it++) \
      it->second = f(it->second); \
    return y; \
  }; \

  #define macro_binary_traj_traj(f) \
  { \
    assert_release(x1.nb_samples() == x2.nb_samples()); \
    auto y = x1; \
    auto it_y = y.begin(); \
    auto it_x2 = x2.begin(); \
    while(it_y != y.end()) \
    { \
      assert_release(it_y->first == it_x2->first \
        && "inconsistent dates between the two trajectories"); \
      it_y->second = f((double)it_y->second,(double)it_x2->second); \
      it_y++; it_x2++; \
    } \
    return y; \
  } \

  #define macro_binary_real_traj(f) \
  { \
    auto y = x2; \
    for(auto& [ti,xi] : y) \
      xi = f((double)x1,(double)xi); \
    return y; \
  } \

  #define macro_binary_traj_real(f) \
  { \
    auto y = x1; \
    for(auto& [ti,xi] : y) \
      xi = f((double)xi,(double)x2); \
    return y; \
  } \

  /** \brief \f$x1(\cdot)\f$
    * \param x1
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  SampledTraj<T> operator+(const SampledTraj<T>& x1);

  /** \brief \f$x_1(\cdot)+x_2(\cdot)\f$
    * \param x1
    * \param x2
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  SampledTraj<T> operator+(const SampledTraj<T>& x1, const SampledTraj<T>& x2);

  /** \brief \f$x_1(\cdot)+x_2\f$
    * \param x1
    * \param x2
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  SampledTraj<T> operator+(const SampledTraj<T>& x1, const T& x2);

  /** \brief \f$x+x_2(\cdot)\f$
    * \param x1
    * \param x2
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  SampledTraj<T> operator+(const T& x1, const SampledTraj<T>& x2);

  /** \brief \f$-x_1(\cdot)\f$
    * \param x1
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  SampledTraj<T> operator-(const SampledTraj<T>& x1);

  /** \brief \f$x_1(\cdot)-x_2(\cdot)\f$
    * \param x1
    * \param x2
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  SampledTraj<T> operator-(const SampledTraj<T>& x1, const SampledTraj<T>& x2);

  /** \brief \f$x_1(\cdot)-x_2\f$
    * \param x1
    * \param x2
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  SampledTraj<T> operator-(const SampledTraj<T>& x1, const T& x2);

  /** \brief \f$x-x_2(\cdot)\f$
    * \param x1
    * \param x2
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  SampledTraj<T> operator-(const T& x1, const SampledTraj<T>& x2);

  /** \brief \f$x_1(\cdot)\cdot x_2(\cdot)\f$
    * \param x1
    * \param x2
    * \return trajectory output
    */
  SampledTraj<double> operator*(const SampledTraj<double>& x1, const SampledTraj<double>& x2);

  /** \brief \f$x_1(\cdot)\cdot x_2\f$
    * \param x1
    * \param x2
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  SampledTraj<T> operator*(const SampledTraj<T>& x1, const T& x2);

  /** \brief \f$x\cdot x_2(\cdot)\f$
    * \param x1
    * \param x2
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  SampledTraj<T> operator*(const T& x1, const SampledTraj<T>& x2);


  /** \brief \f$x_1(\cdot)/x_2(\cdot)\f$
    * \param x1
    * \param x2
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  SampledTraj<T> operator/(const SampledTraj<T>& x1, const SampledTraj<T>& x2);

  /** \brief \f$x_1(\cdot)/x_2\f$
    * \param x1
    * \param x2
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  SampledTraj<T> operator/(const SampledTraj<T>& x1, const T& x2);

  /** \brief \f$x/x_2(\cdot)\f$
    * \param x1
    * \param x2
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  SampledTraj<T> operator/(const T& x1, const SampledTraj<T>& x2);

  /** \brief \f$x^2(\cdot)\f$
    * \param x1
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  inline SampledTraj<T> sqr(const SampledTraj<T>& x1)
  {
    return pow(x1,2);
  }

  /** \brief \f$\sqrt{x_1(\cdot)}\f$
    * \param x1
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  inline SampledTraj<T> sqrt(const SampledTraj<T>& x1)
  /// \cond
    macro_unary_traj(std::sqrt)
  /// \endcond
  ;

  /** \brief \f$x^x_2(\cdot)\f$
    * \param x1
    * \param x2
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  inline SampledTraj<T> pow(const SampledTraj<T>& x1, int x2)
  /// \cond
    macro_binary_traj_real(codac2::pow)
  /// \endcond
  ;

  /** \brief \f$x^x_2(\cdot)\f$
    * \param x1
    * \param x2
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  inline SampledTraj<T> pow(const SampledTraj<T>& x1, double x2)
  /// \cond
    macro_binary_traj_real(codac2::pow)
  /// \endcond
  ;

  /** \brief \f$x^x_2(\cdot)\f$
    * \param x1
    * \param x2
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  inline SampledTraj<T> root(const SampledTraj<T>& x1, int x2)
  /// \cond
    macro_binary_traj_real(codac2::root)
  /// \endcond
  ;

  /** \brief \f$\exp(x_1(\cdot))\f$
    * \param x1
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  inline SampledTraj<T> exp(const SampledTraj<T>& x1)
  /// \cond
    macro_unary_traj(std::exp)
  /// \endcond
  ;

  /** \brief \f$\log(x_1(\cdot))\f$
    * \param x1
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  inline SampledTraj<T> log(const SampledTraj<T>& x1)
  /// \cond
    macro_unary_traj(std::log)
  /// \endcond
  ;

  /** \brief \f$\cos(x_1(\cdot))\f$
    * \param x1
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  inline SampledTraj<T> cos(const SampledTraj<T>& x1)
  /// \cond
    macro_unary_traj(std::cos)
  /// \endcond
  ;

  /** \brief \f$\sin(x_1(\cdot))\f$
    * \param x1
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  inline SampledTraj<T> sin(const SampledTraj<T>& x1)
  /// \cond
    macro_unary_traj(std::sin)
  /// \endcond
  ;

  /** \brief \f$\tan(x_1(\cdot))\f$
    * \param x1
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  inline SampledTraj<T> tan(const SampledTraj<T>& x1)
  /// \cond
    macro_unary_traj(std::tan)
  /// \endcond
  ;

  /** \brief \f$\arccos(x_1(\cdot))\f$
    * \param x1
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  inline SampledTraj<T> acos(const SampledTraj<T>& x1)
  /// \cond
    macro_unary_traj(std::acos)
  /// \endcond
  ;

  /** \brief \f$\arcsin(x_1(\cdot))\f$
    * \param x1
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  inline SampledTraj<T> asin(const SampledTraj<T>& x1)
  /// \cond
    macro_unary_traj(std::asin)
  /// \endcond
  ;

  /** \brief \f$\arctan(x_1(\cdot))\f$
    * \param x1
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  inline SampledTraj<T> atan(const SampledTraj<T>& x1)
  /// \cond
    macro_unary_traj(std::atan)
  /// \endcond
  ;

  /** \brief \f$\mathrm{arctan2}(x_1(\cdot),x_2(\cdot))\f$
    * \param x1
    * \param x2
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  inline SampledTraj<T> atan2(const SampledTraj<T>& x1, const SampledTraj<T>& x2)
  /// \cond
    macro_binary_traj_traj(std::atan2)
  /// \endcond
  ;

  /** \brief \f$\mathrm{arctan2}(x_1(\cdot),x_2)\f$
    * \param x1
    * \param x2
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  inline SampledTraj<T> atan2(const SampledTraj<T>& x1, const T& x2)
  /// \cond
    macro_binary_traj_real(std::atan2)
  /// \endcond
  ;

  /** \brief \f$\mathrm{arctan2}(x_1, x_2(\cdot))\f$
    * \param x1
    * \param x2
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  inline SampledTraj<T> atan2(const T& x1, const SampledTraj<T>& x2)
  /// \cond
    macro_binary_real_traj(std::atan2)
  /// \endcond
  ;

  /** \brief \f$\cosh(x_1(\cdot))\f$
    * \param x1
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  inline SampledTraj<T> cosh(const SampledTraj<T>& x1)
  /// \cond
    macro_unary_traj(std::cosh)
  /// \endcond
  ;

  /** \brief \f$\sinh(x_1(\cdot))\f$
    * \param x1
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  inline SampledTraj<T> sinh(const SampledTraj<T>& x1)
  /// \cond
    macro_unary_traj(std::sinh)
  /// \endcond
  ;

  /** \brief \f$\tanh(x_1(\cdot))\f$
    * \param x1
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  inline SampledTraj<T> tanh(const SampledTraj<T>& x1)
  /// \cond
    macro_unary_traj(std::tanh)
  /// \endcond
  ;

  /** \brief \f$\mathrm{arccosh}(x_1(\cdot))\f$
    * \param x1
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  inline SampledTraj<T> acosh(const SampledTraj<T>& x1)
  /// \cond
    macro_unary_traj(std::acosh)
  /// \endcond
  ;

  /** \brief \f$\mathrm{arcsinh}(x_1(\cdot))\f$
    * \param x1
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  inline SampledTraj<T> asinh(const SampledTraj<T>& x1)
  /// \cond
    macro_unary_traj(std::asinh)
  /// \endcond
  ;

  /** \brief \f$\mathrm{arctanh}(x_1(\cdot))\f$
    * \param x1
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  inline SampledTraj<T> atanh(const SampledTraj<T>& x1)
  /// \cond
    macro_unary_traj(std::atanh)
  /// \endcond
  ;

  /** \brief \f$\mid x_1(\cdot)\mid\f$
    * \param x1
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  inline SampledTraj<T> abs(const SampledTraj<T>& x1)
  /// \cond
    macro_unary_traj(std::abs)
  /// \endcond
  ;

  /** \brief \f$\min(x_1(\cdot),x_2(\cdot))\f$
    * \param x1
    * \param x2
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  inline SampledTraj<T> min(const SampledTraj<T>& x1, const SampledTraj<T>& x2)
  /// \cond
    macro_binary_traj_traj(std::min)
  /// \endcond
  ;

  /** \brief \f$\min(x_1(\cdot),x_2)\f$
    * \param x1
    * \param x2
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  inline SampledTraj<T> min(const SampledTraj<T>& x1, const T& x2)
  /// \cond
    macro_binary_traj_real(min)
  /// \endcond
  ;

  /** \brief \f$\min(x_1, x_2(\cdot))\f$
    * \param x1
    * \param x2
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  inline SampledTraj<T> min(const T& x1, const SampledTraj<T>& x2)
  /// \cond
    macro_binary_real_traj(min)
  /// \endcond
  ;

  /** \brief \f$\max(x_1(\cdot),x_2(\cdot))\f$
    * \param x1
    * \param x2
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  inline SampledTraj<T> max(const SampledTraj<T>& x1, const SampledTraj<T>& x2)
  /// \cond
    macro_binary_traj_traj(std::max)
  /// \endcond
  ;

  /** \brief \f$\max(x_1(\cdot),x_2)\f$
    * \param x1
    * \param x2
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  inline SampledTraj<T> max(const SampledTraj<T>& x1, const T& x2)
  /// \cond
    macro_binary_traj_real(max)
  /// \endcond
  ;

  /** \brief \f$\max(x_1, x_2(\cdot))\f$
    * \param x1
    * \param x2
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  inline SampledTraj<T> max(const T& x1, const SampledTraj<T>& x2)
  /// \cond
    macro_binary_real_traj(max)
  /// \endcond
  ;

  /** \brief \f$\mid x_1(\cdot)\mid\f$
    * \param x1
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  inline SampledTraj<T> sign(const SampledTraj<T>& x1)
  /// \cond
    macro_unary_traj(codac2::sign)
  /// \endcond
  ;

  /** \brief \f$\mid x_1(\cdot)\mid\f$
    * \param x1
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  inline SampledTraj<T> integer(const SampledTraj<T>& x1)
  /// \cond
    macro_unary_traj(codac2::integer)
  /// \endcond
  ;

  /** \brief \f$\mid x_1(\cdot)\mid\f$
    * \param x1
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  inline SampledTraj<T> floor(const SampledTraj<T>& x1)
  /// \cond
    macro_unary_traj(std::floor)
  /// \endcond
  ;

  /** \brief \f$\mid x_1(\cdot)\mid\f$
    * \param x1
    * \return trajectory output
    */
  template<typename T>
    requires IsRealType<T>
  inline SampledTraj<T> ceil(const SampledTraj<T>& x1)
  /// \cond
    macro_unary_traj(std::ceil)
  /// \endcond
  ;
}