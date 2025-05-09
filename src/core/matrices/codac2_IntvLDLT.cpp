/** 
 *  \file codac2_LDLT.cpp
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Damien Massé
 *  \copyright  Copyright 2025 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <iostream>
#include <cstdlib>
#include "codac2_Matrix.h"
#include "codac2_Interval.h"
#include "codac2_IntervalMatrix.h"
#include "codac2_IntervalVector.h"
#include "codac2_BoolInterval.h"
#include "codac2_IntvLDLT.h"

namespace codac2 {

/* utility function : compute the bound using Gauss-Jordan successive
   algorithm. (note: quite different from buildLUbounds...)
   The residual matrix L^{-1} PMP^T L^T^{-1} is a parameter as E
   which is (initially) symmetric (only the lower part is used initially) */

IntervalMatrix IntvLDLT::buildLDLTbounds(const IntervalMatrix &E) {

    const Index dim = E.cols();
    IntervalMatrix IL = E;
    IntervalMatrix res = IntervalMatrix::Identity(dim,dim);
  
    for (Index k=0;k<dim;k++) {
       /* IL stores :
           (sum L^n)           | sum L^n E^T
          ----------------------------
           E sumL^T^n D_k^-1   |   E

          res stores L (stricly lower) and D (diagonal) */

       if (k>0) {
          /* storing new L */
          res.block(k,0,1,k) = IL.block(k,0,1,k);
          IL.block(k,0,1,k) = - res.block(k,0,1,k)*IL.topLeftCorner(k,k);
       }
       /* new diagonal coefficient */
       Interval d = IL(k,k) - (res.block(k,0,1,k)*IL.block(0,k,k,1)).value();
       if (d.contains(0.0) && k<dim-1) status=-2;
       
       if (k<dim-1) {
          IL.block(k,k+1,1,dim-k-1) =
			IL.block(k+1,k,dim-k-1,1).transpose()
			-res.block(k,0,1,k)*IL.topRightCorner(k,dim-k-1);
          IL.block(k+1,k,dim-k-1,1) = 
		IL.block(k,k+1,1,dim-k-1).transpose() /d;
       }
       res(k,k)=d;
    }
    return res;
}


/************************************************************/


/** constructor from Matrix of double 
 */
IntvLDLT::IntvLDLT(const Matrix &M) :
   status(0), _LDLT(M), 
   matrixLDLT_(M.rows(), M.rows())
{
   this->computeMatrixLDLT(IntervalMatrix(M));
}
  
/** constructor from Matrix of Intervals
 */
IntvLDLT::IntvLDLT(const IntervalMatrix &M) :
   status(0),_LDLT(M.mid()), 
   matrixLDLT_(M.rows(), M.rows())
{
   this->computeMatrixLDLT(M);
}

void IntvLDLT::computeMatrixLDLT(const IntervalMatrix &M) {
    if (_LDLT.info()!=Eigen::Success) {
       status=-1; return;
    }
    const Index dim = M.cols();
    IntervalMatrix build= 
       _LDLT.transpositionsP() * M * _LDLT.transpositionsP();
    std::cout << "build : " << build << "\n";

    /* construction de L{-1} */
    IntervalMatrix ImL=IntervalMatrix::Identity(dim,dim);
    Matrix L = _LDLT.matrixL();
//    L.triangularView<Eigen::StrictlyLower>() = 1*L;
    std::cout << "L : " << L << "\n";
    /* inversion of L */
    for (int c=0;c<dim;c++) {
      for (int r=c+1;r<dim;r++) {
         Interval s(L(r,c));
         for (int k=c+1;k<r;k++) {
            s += L(r,k)*ImL(k,c);
         }
         ImL(r,c)=-s;
      }
    }
    build = ImL * build * ImL.transpose();
   
    std::cout << " build : " << build << "\n";
    matrixLDLT_ = buildLDLTbounds(build);
    std::cout << " matrix LDLT_(0) : " << matrixLDLT_ << "\n";
    matrixLDLT_.triangularView<Eigen::StrictlyLower>()
        = IntervalMatrix(L) * (matrixLDLT_).triangularView<Eigen::UnitLower>();
    std::cout << " matrix LDLT_ : " << matrixLDLT_ << "\n";
    if (status==0) status=1;
}

bool IntvLDLT::info() const {
    return (status==1 || status==-2);
}

BoolInterval IntvLDLT::isNegative() const {
    if (status==0 || status==-1) return BoolInterval::EMPTY;
    BoolInterval ret = BoolInterval::TRUE;
    for (int i=0;i<matrixLDLT_.cols();i++) {
       if (matrixLDLT_(i,i).ub()<0.0) continue;
       if (matrixLDLT_(i,i).lb()>0.0) return BoolInterval::FALSE;
       ret = BoolInterval::UNKNOWN;
    }
    return ret;
}

BoolInterval IntvLDLT::isPositive() const {
    if (status==0 || status==-1) return BoolInterval::EMPTY;
    BoolInterval ret = BoolInterval::TRUE;
    for (int i=0;i<matrixLDLT_.cols();i++) {
       if (matrixLDLT_(i,i).lb()>0.0) continue;
       if (matrixLDLT_(i,i).ub()<0.0) return BoolInterval::FALSE;
       ret = BoolInterval::UNKNOWN;
    }
    return ret;
}

Interval IntvLDLT::determinant() const {
    if (status==0 || status==-1) return Interval::empty();
    Interval ret(1.0);
    for (int i=0;i<matrixLDLT_.cols();i++) {
       ret *= matrixLDLT_(i,i);
    }
    return ret;
}

IntervalMatrix IntvLDLT::solve(const IntervalMatrix &rhs) const {
    assert(status!=0);
    assert(rhs.rows() == matrixLDLT_.rows());
    if (status==-1) {
       return IntervalMatrix(matrixLDLT_.cols(),rhs.cols());
    }
    IntervalMatrix rhs2 = this->_LDLT.transpositionsP() * rhs;
    for (int i=0;i<matrixLDLT_.rows()-1;i++) {
      for (int j=i+1;j<matrixLDLT_.rows();j++) {
         rhs2.row(j) -= matrixLDLT_(j,i) * rhs2.row(i);
      }
    }
    for (int i=0;i<matrixLDLT_.rows();i++) {
       rhs2.row(i) = rhs2.row(i) / matrixLDLT_(i,i);
    }
    for (int i=matrixLDLT_.rows()-1;i>0;i--) {
      for (int j=0;j<i;j++) {
         rhs2.row(j) -= matrixLDLT_(i,j) * rhs2.row(i);
      }
    }
    return this->_LDLT.transpositionsP().inverse() * rhs2;
}

IntervalMatrix IntvLDLT::ReconstructedMatrix() const {
    assert(status==1 || status==-2);
    IntervalMatrix build = matrixLDLT_.diagonal().asDiagonal();
    build = this->matrixL() * build * this->matrixU();
    build = this->_LDLT.transpositionsP().inverse() * build *
		this->_LDLT.transpositionsP().inverse();
    return build;
}


}
