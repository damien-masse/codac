/** 
 *  \file codac2_Zonotope.h
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Maël Godard
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include <vector>
#include <codac2_Vector.h>
#include <codac2_Matrix.h>
#include <codac2_IntervalVector.h>

using namespace std;

namespace codac2
{
  /**
   * \class Zonotope
   * \brief Class representing a zonotope z + A*[-1,1]^m
   * 
   * This class represents a zonotope in n-dimensional space, defined by a center point `z` and a shape matrix `A`.
   * 
   * The Vector `z` and each column of the Matrix `A` must have the same dimension `n`, but the matrix `A` can have any number of columns `m`.
   */
  class Zonotope
  {
    public:

      /**
       * \brief Constructs a n-zonotope object with a given center and shape matrix
       * 
       * \param z Center of the zonotope (n-dimensional vector)
       * \param A Shape matrix of the zonotope (n x m matrix)
       */
      Zonotope(const Vector& z, const Matrix& A);

      /**
       * \brief Center of the zonotope
       */
      Vector z;

      /**
       * \brief Shape matrix of the zonotope
       */
      Matrix A;
      
  };
}
