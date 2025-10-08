/** 
 *  codac2_peibos.cpp
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Maël Godard
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include "codac2_peibos.h"

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

  double error(const IntervalVector& z, const IntervalMatrix& JJg, const Matrix& B, const IntervalVector& X)
  {
    auto xc = X.mid();

    auto dX = X-xc;

    auto a =  z.mid();

    auto E = (z - a) + (JJg - B)*dX;

    return E.norm().ub();
  }

  Matrix inflate_flat_parallelepiped(const Matrix& A, const Vector& e_vec, double rho)
  {
    Index m = A.cols();
    Index n = A.rows();

    Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(A.transpose());
    Eigen::MatrixXd N = lu_decomp.kernel();

    Matrix A_tild (n,n);
    A_tild << A, N;

    Matrix Q = (A_tild.transpose() * A_tild).inverse();

    Matrix mult = Matrix::Zero(n,n);
    for (int i = 0; i < n; i++)
      mult(i,i) = rho*std::sqrt(Q(i,i));

    for (int i = 0; i < m; i++)
      mult(i,i) += e_vec(i);
    
    return A_tild*mult;
  }

  Parallelepiped parallelepiped_inclusion(const AnalyticFunction<VectorType>& f, const IntervalVector& X)
  {
    assert_release (f.input_size() < f.output_size() && "input size must be strictly lower than output size");
    assert_release (X.size() == f.input_size() && "dimension of X must match input size of f");

    auto z = f.eval(X.mid());

    Matrix B = f.diff(X.mid()).mid();

    // Maximum error computation
    double rho = error(z, f.diff(X), B, X);

    // Inflation of the parallelepiped
    auto A = inflate_flat_parallelepiped(B, X.rad(), rho);

    return Parallelepiped(z.mid(), A);
  }

  Parallelepiped parallelepiped_inclusion(const IntervalVector& z, const IntervalMatrix& JJf, const Matrix& Jf_tild, const AnalyticFunction<VectorType>& psi_0, const OctaSym& sigma, const IntervalVector& X)
  {
    // Computation of the Jacobian of g = f o sigma(psi_0)
    IntervalMatrix JJg = JJf * (sigma.permutation_matrix().template cast<Interval>()) * psi_0.diff(X);

    // B is an approximation of the Jacobian of g at the center of X
    Matrix B = (Jf_tild * sigma.permutation_matrix() * (psi_0.diff(X.mid()).mid()));

    // Maximum error computation
    double rho = error(z, JJg, B, X);

    // Inflation of the parallelepiped
    auto A = inflate_flat_parallelepiped(B, X.rad(), rho);

    return Parallelepiped(z.mid(), A);
  }

  vector<Parallelepiped> PEIBOS(const AnalyticFunction<VectorType>& f, const AnalyticFunction<VectorType>& psi_0, const vector<OctaSym>& symmetries, double epsilon, bool verbose)
  {
    return PEIBOS(f, psi_0, symmetries, epsilon, Vector::Zero(psi_0.output_size()), verbose);
  }

  vector<Parallelepiped> PEIBOS(const AnalyticFunction<VectorType>& f, const AnalyticFunction<VectorType>& psi_0, const vector<OctaSym>& symmetries, double epsilon, const Vector& offset, bool verbose)
  {
    Index m = psi_0.input_size();

    assert_release (f.input_size() == psi_0.output_size() && "output size of psi_0 must match input size of f");
    assert_release (offset.size() == psi_0.output_size() && "offset size must match output size of psi_0");
    assert_release (m < psi_0.output_size());
    assert_release (symmetries.size() > 0 && (int) symmetries[0].size() == psi_0.output_size() && "no generator given or wrong dimension of generator (must match output size of psi_0)");

    clock_t t_start = clock();

    vector<Parallelepiped> output;

    vector<IntervalVector> boxes;
    double true_eps = split(IntervalVector::constant(m,{-1,1}), epsilon, boxes);

    for (const auto& symmetry : symmetries)
    {
      VectorVar x(m);
      AnalyticFunction g_i ({x}, f(symmetry(psi_0(x))+offset));

      for (const auto& X : boxes)        
        output.push_back(parallelepiped_inclusion(g_i, X));
    }

    if (verbose)
    {
      printf("\nPEIBOS statistics:\n");
      printf("------------------\n");
      printf("Number of symmetries: %ld\n", symmetries.size());
      printf("Real epsilon: %.4f\n", true_eps);
      printf("Computation time: %.4fs\n\n", (double)(clock()-t_start)/CLOCKS_PER_SEC);
    }

    return output;
  }
}