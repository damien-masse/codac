/** 
 *  \file codac2_IntLDLT.h
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Damien Massé
 *  \copyright  Copyright 2025 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include <ostream>
#include "codac2_Matrix.h"
#include "codac2_IntervalMatrix.h"
#include "codac2_IntervalVector.h"
#include "codac2_BoolInterval.h"

namespace codac2
{
  /** \brief Robust Cholesky decomposition for a symetric matrix of intervals,
   * with pivoting and without square root, based on Eigen decomposition.
   * The decomposition is of the form M = P^T [L][D][L]^T P where P is a
   * permutation matrix, [L] is unit-lower-triangular and [D] is diagonal.
   * Note that the eigen decomposition may fail (see info()).
   */
  class IntLDLT
  {
     public:
     /** \brief constructor from Matrix of double. Note that only
      *  the lower part of M is used.
      *  \param M the matrix for which the decomposition is computed
      */
     explicit IntLDLT(const Matrix &M);

     /** \brief constructor from Matrix of intervals. Eigen 
      *  decomposition is done on M.mid(). Once again, only the lower part
      *  of M is used.
      *
      *  \param M the matrix of intervals.
      */
     explicit IntLDLT(const IntervalMatrix &M);

     /** \brief returns if the eigen computation has been successful 
      *  on not. LDLT decomposition can fail with a 0 pivot (e.g.
      *  non-zero matrix with a 0 diagonal). The failure to extend
      *  the decomposition to matrix of intervals is not treated here:
      *  the result is the appearance of oo coefficients. 
      *
      *  \return true if the computation worked, false otherwise.
      */ 
     bool info() const;

     /** \brief check if the current matrix is negative (semidefinite)
      *  returns EMPTY if the decomposition failed.
      *  \return a BoolInterval
      */
     BoolInterval isNegative() const;
     /** \brief check if the initial matrix is positive (semidefinite)
      *  returns EMPTY if the decomposition failed.
      * 
      *  \return the BoolInterval
      */
     BoolInterval isPositive() const;
     /** \brief return an interval enclosing the determinant
      *  
      *  \return the product of the diagonal elements of [U]
      */
     Interval determinant() const;
     /** \brief equation solving M X = rhs 
      *  use successive operations 
      *  P^T X0 = rhs, L X1 = X0, D X2 = X1, L^T X3 = X2, P X = X3
      *  all matrices except D are invertible, 
      *  if D is (possibly) not invertible, compute a enclosing
      *  of the least-square solution of D X2 = X1 (does not mean
      *  that the overall solution is a least-square solution of M X = rhs).
      *  
      *  \param rhs right-hand side of the equation 
      *  \return a potential solution of the equation M X = rhs
      */
     IntervalMatrix solve(const IntervalMatrix &rhs) const;


     /** \brief rebuilding of the matrix, ie compute P^T [L][D][L]^T P
      *  can be used to evaluate the precision of the decomposition
      *
      *  \return the reconstructed matrix
      */
     IntervalMatrix ReconstructedMatrix() const;
     
     /** \brief the permutation P in the decomposition P^T [L][D][L]^T P.
      *  Note that Eigen returns it as a sequence of transpositions.
      *
      * \return the permutation P, as defined by Eigen
      */
     const Eigen::LDLT<Matrix>::TranspositionType
			&transpositionsP () const;
     /** \brief the Eigen decomposition of M.mid()
      *
      * \return the Eigen decomposition
      */
     const Eigen::LDLT<Matrix> &eigenLDLT() const;

     /** \brief returns the matrix storing [L] (stricly lower),
      * [D] (diagonal) and [L]^T (stricly upper)
      *
      * \return the interval matrix
      */
     const IntervalMatrix &matrixLDLT() const;

     /** \brief returns the vector of diagonals
      * note: Eigen returns a type Diagonal<const IntervalMatrix>.
      *
      * \return the vector
      */
     IntervalVector vectorD() const;
     /** \brief returns the matrix L
      * note: Eigen returns a view of the storage matrix.
      *
      * \return the matrix
      */
     IntervalMatrix matrixL() const;
     /** \brief returns the matrix L^T
      * note: Eigen returns a view of the tranpose of the storage matrix.
      *
      * \return the matrix
      */
     IntervalMatrix matrixU() const;

     private:
        int status; /* 0 : uninitialised 
		      -1 : eigen3 decomposition failed i
		      -2 : interval failed  1 : ok */
        Eigen::LDLT<Matrix> _LDLT;

        IntervalMatrix matrixLDLT_;   /* LDLT matrix */
     
        void computeMatrixLDLT(const IntervalMatrix &M);
        IntervalMatrix buildLDLTbounds
		(const IntervalMatrix &E);
  };

inline const IntervalMatrix &IntLDLT::matrixLDLT() const {
     assert_release(this->info());
     return this->matrixLDLT_;
}
inline IntervalMatrix IntLDLT::matrixL() const {
     assert_release(this->info());
     return IntervalMatrix
	(this->matrixLDLT_.triangularView<Eigen::UnitLower>());
}
inline IntervalMatrix IntLDLT::matrixU() const {
     assert_release(this->info());
     return this->matrixL().transpose();
}
inline IntervalVector IntLDLT::vectorD() const {
     assert_release(this->info());
     return IntervalVector
	(this->matrixLDLT_.diagonal());
}
inline const Eigen::LDLT<Matrix>::TranspositionType
			&IntLDLT::transpositionsP() const {
    assert_release(this->info());
    return this->_LDLT.transpositionsP();
}
inline const Eigen::LDLT<Matrix> &IntLDLT::eigenLDLT() const {
    return this->_LDLT;
}

}


