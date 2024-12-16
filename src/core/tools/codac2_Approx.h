/** 
 *  \file codac2_Approx.h
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include <limits>
#include "codac2_matrices.h"
#include "codac2_Vector.h"
#include "codac2_IntervalVector.h"
#include "codac2_Row.h"
#include "codac2_IntervalRow.h"
#include "codac2_Matrix.h"
#include "codac2_IntervalMatrix.h"

namespace codac2
{
  #define DEFAULT_EPS std::numeric_limits<double>::epsilon()*10

  template<typename T>
  class Approx
  {
    public:

      explicit Approx(const T& x, double eps = DEFAULT_EPS)
        : _x(x), _eps(eps)
      { }

      template<typename OtherDerived>
      explicit Approx(const Eigen::MatrixBase<OtherDerived>& x, double eps = DEFAULT_EPS)
        : _x(x.eval()), _eps(eps)
      { }

      friend bool operator==(const T& x1, const Approx<T>& x2)
      {
        if constexpr(std::is_same_v<T,double>)
          return std::fabs(x1-x2._x) < x2._eps;

        else if(x1.size() != x2._x.size())
          return false;

        else if(x1 == x2._x)
          return true;

        else if constexpr(std::is_same_v<T,Interval>)
        {
          return (x1.lb() == x2._x.lb() || x1.lb() == Approx<double>(x2._x.lb(),x2._eps))
              && (x1.ub() == x2._x.ub() || x1.ub() == Approx<double>(x2._x.ub(),x2._eps));
        }

        else if constexpr(std::is_same_v<T,Vector>
          || std::is_same_v<T,IntervalVector>
          || std::is_same_v<T,Row>
          || std::is_same_v<T,IntervalRow>
          || std::is_same_v<T,Matrix>
          || std::is_same_v<T,IntervalMatrix>)
        {
          for(Index i = 0 ; i < x1.rows() ; i++)
            for(Index j = 0 ; j < x1.cols() ; j++)
              if(!(x1(i,j) == Approx<typename T::Scalar>(x2._x(i,j), x2._eps)))
                return false;
          return true;
        }

        else
        {
          assert_release(false && "Approx::operator== unhandled case")
          return false;
        }
      }

      friend std::ostream& operator<<(std::ostream& os, const Approx<T>& x)
      {
        os << "Approx(" << x._x << ")";
        return os;
      }
    
    private:

      const T _x;
      const double _eps;
  };
  
}