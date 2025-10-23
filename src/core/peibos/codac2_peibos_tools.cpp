/** 
 *  codac2_peibos_tools.cpp
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Maël Godard
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include "codac2_peibos_tools.h"
#include "codac2_inversion.h"
#include "codac2_IntvFullPivLU.h"

using namespace std;
using namespace codac2;

namespace codac2
{
  double split (const IntervalVector& X, double eps, vector<IntervalVector>& boxes)
  {
    if (X.max_diam()<=eps)
    {
      boxes.push_back(X);
      return X.max_diam();
    }
    else
    {
      auto p = X.bisect_largest(0.5);
      double diam1 = split(p.first,eps,boxes);
      double diam2 = split (p.second,eps,boxes);
      return std::max(diam1,diam2);
    }
  }

  double error(const IntervalVector& Y, const Vector& z, const IntervalMatrix& Jf, const Matrix& A, const IntervalVector& X)
  {
    Vector xc = X.mid();

    IntervalVector dX = X-xc;

    IntervalVector E = (Y - z) + (Jf - A)*dX;

    return E.norm().ub();
  }

  Matrix inflate_flat_parallelepiped(const Matrix& A, const Vector& e_vec, double rho)
  {
    Index m = A.cols();
    Index n = A.rows();

    IntervalMatrix A_int (A);

    IntvFullPivLU LUdec((IntervalMatrix) A_int.transpose());
    IntervalMatrix N = LUdec.kernel();

    vector<int> col_indices;

    for (int i = 0; i < m; i++)
      if (! (A.col(i).isZero()) )
        col_indices.push_back(i);

    // A_reduced is composed of the non empty vectors of A_int
    IntervalMatrix A_reduced (n, col_indices.size());
    for (int i = 0; i < (int) col_indices.size(); i++)
      A_reduced.col(i) = A_int.col(col_indices[i]);

    // Construct the new matrix A_tild = [A_reduced | N]
    IntervalMatrix A_tild (n,n);
    A_tild << A_reduced,N;

    IntervalMatrix Q = inverse_enclosure(A_tild.transpose() * A_tild);

    IntervalMatrix mult = IntervalMatrix::zero(n,n);
    for (int i = 0; i < n; i++)
      mult(i,i) = rho*sqrt(Q(i,i));

    for (int i = 0; i < m; i++)
      mult(i,i) += e_vec(i);

    // From here we have an IntervalMatrix A_inf, meaning that each generator is an interval vector
    IntervalMatrix A_inf = A_tild * mult;

    // The initial parallelepiped is <y> = Y*[-1,1]^n
    Matrix Y = A_inf.smig();

    // The following is similar to the previous operations. The difference is that here the Matrix Y is square and not interval

    // rho2 is the sum of the radiuses of the circle enclosing each interval generator
    Interval rho2 (0.);

    for (int i = 0; i < n; i++)
      rho2 += A_inf.col(i).rad().norm();

    IntervalMatrix Q2 = inverse_enclosure(Y.transpose() * Y);

    IntervalMatrix Y2 = IntervalMatrix::zero(n,n);

    for (int i = 0; i < n; i++)
      Y2.col(i) = Y.col(i)*(1+rho2*sqrt(Q2(i,i)));
    
    return Y2.smag();
  }
}