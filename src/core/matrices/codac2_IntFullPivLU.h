/** 
 *  \file codac2_IntFullPivLU.h
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Damien Massé
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include <ostream>
#include "codac2_Matrix.h"
#include "codac2_IntervalMatrix.h"
#include "codac2_BoolInterval.h"

namespace codac2
{
  // Full pivot LU decomposition applied by Eigen on a m*n matrix,
  // modified afterwards to handle interval of matrices

  class IntFullPivLU
  {
     public:
     /** constructor from Matrix of double 
      */
     explicit IntFullPivLU(const Matrix &M);

     /** constructor from Matrix of Intervals */
     explicit IntFullPivLU(const IntervalMatrix &M);

     /** determines if the initial matrix is injective */
     BoolInterval isInjective() const;
     /** determines if the initial matrix is invertible */
     BoolInterval isInvertible() const;
     /** determines if the initial matrix is surjective */
     BoolInterval isSurjective() const;
     /** overapproximation of the determinant. Only for square matrices */
     Interval determinant() const;
     /** rank **/
     Interval rank() const;
     /** approximation of the size of the kernel space.
      *  based on the result of rank(). As such, this is not
      *  the exact size of the kernel space as build by kernel() */
     Interval dimensionOfKernel() const;
     /** overapproximation of the kernel space. vectors
      *  outside the linear combinations of the columns vectors
      *  are guaranteed to be outside the kernel of the matrix */
     IntervalMatrix kernel() const;
     /** "underapproximation" of the column space of the matrix,
      *  ie return a set of independant columns of the original matrix
      *  which is possibly maximal. As with Eigen, you must 
      *  provide the original matrix used for the decomposition,
      *  hence two versions of the function are provided */
     IntervalMatrix image(const IntervalMatrix &M) const;
     /** "underapproximation" of the column space of the matrix,
      *  ie return a set of independant columns of the original matrix
      *  which is possibly maximal. As for Eigen, you must 
      *  provide the original matrix used for the decomposition,
      *  hence the use of a template */
     template <typename Derived> Derived
	image(const Eigen::MatrixBase<Derived> &M) const;
     /** equation solving M X = rhs 
      *  precisely look for solutions where the only non-zero
      *  values are those on non-zero pivots 
      *  if the matrix is full-rank and surjective (cols >= rows),
      *  it gives an overapproximation of the solutions (for each 
      *  possible values of rhs)
      *  if the matrix is full-rank and injective (rows >= cols),
      *  it returns empty if no solution is possible, and a
      *  possible overapproximation of the solutions otherwise
      *  (but there _may_ still be no solution)
      *  if the matrix is not full-rank, 
      *     empty means that no solution is possible with the
      *        initial precondition (non-zero values for non-zero pivots)
      *     non empty presents the possible solutions found */
     IntervalMatrix solve(const IntervalMatrix &rhs) const;


     /** rebuilding of the matrix for debug purposes
      *  (or estimation of the precision of the decomposition) */
     IntervalMatrix ReconstructedMatrix() const;
     /** maximum magnitude of the diagonal (maximum pivot) */
     double maxPivot() const;
     
     /** the permutation P in the decomposition P{-1}LUQ{-1} */
     const Eigen::FullPivLU<Matrix>::PermutationPType 
			&permutationP () const;
     /** the permutation Q in the decomposition P{-1}LUQ{-1} */
     const Eigen::FullPivLU<Matrix>::PermutationQType 
			&permutationQ () const;
     /** the eigen LU decomposition for M.mid() */
     const Eigen::FullPivLU<Matrix> &eigenLU() const;

     /** Enclosure of the adapted matrix, if possible */
     const IntervalMatrix &matrixLU() const;


     private:
        Eigen::FullPivLU<Matrix> _LU;

        IntervalMatrix matrixLU_;   /* LU matrix */
     
        void computeMatrixLU(const IntervalMatrix &M, double nonzero);
        static IntervalMatrix buildLUbounds(const IntervalMatrix &E);
  };

inline const IntervalMatrix &IntFullPivLU::matrixLU() const {
     return this->matrixLU_;
}
inline const Eigen::FullPivLU<Matrix>::PermutationPType 
			&IntFullPivLU::permutationP() const {
    return this->_LU.permutationP();
}
inline const Eigen::FullPivLU<Matrix>::PermutationQType 
			&IntFullPivLU::permutationQ() const {
    return this->_LU.permutationQ();
}
inline const Eigen::FullPivLU<Matrix> &IntFullPivLU::eigenLU() const {
    return this->_LU;
}

template<typename Derived>
inline  Derived IntFullPivLU::image
	(const Eigen::MatrixBase<Derived> &M) const
{
   int rk = this->rank().lb();
   if (rk==0) {
      return Derived::Zero(M.rows(),1); 
		/* NdDamien : où est-ce qu'on dit que Derived est une
		   matrice ??? je crois que je hais le C++ ;)  */
   }
   Derived ret = 
		Derived::Zero(M.rows(),rk);
   Index p = 0;
   Index dim = std::min(matrixLU_.rows(),matrixLU_.cols());
   auto Q = this->_LU.permutationQ();
   for (Index i = 0; i<dim; i++) {
      if (!matrixLU_(i,i).contains(0.0)) {
         ret.col(p) = M.col(Q.indices().coeff(i));
         p++;
      }
   }
   return ret;
}


}


