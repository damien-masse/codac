/** 
 *  \file codac2_GaussJordan.h
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Luc Jaulin, Simon Rohou, Damien Massé
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include <ostream>
#include "codac2_Matrix.h"

namespace codac2
{
  // Gauss Jordan band diagonalization preconditioning
  /** \brief Gauss Jordan band diagonalization preconditioning :
   *  takes a Matrix A with A.cols >= A.rows and returns
   *  P such that P A is a ``band matrix''. The square left part of A
   *  must be full-rank.
   *  \param A matrix with A.cols >= A.rows
   *  \return P such that P A is ``band matrix''. Return Id if the
   *  left part of A is not well-conditioned.
   */
  Matrix gauss_jordan(const Matrix& A);

  // From https://www.ensta-bretagne.fr/jaulin/centered.html
}
