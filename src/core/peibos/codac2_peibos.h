/** 
 *  \file codac2_peibos.h
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Maël Godard
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include "codac2_Parallelepiped.h"
#include "codac2_AnalyticFunction.h"
#include "codac2_OctaSym.h"
#include <codac2_OctaSym_operator.h>


namespace codac2
{
  double split (const IntervalVector& X, double eps, vector<IntervalVector>& boxes);

  double error(const IntervalVector& z, const IntervalMatrix& JJg, const Matrix& B, const IntervalVector& X);

  Matrix inflate_flat_parallelepiped (const Matrix& A, const Vector& e_vec, double rho);

  Parallelepiped parallelepiped_inclusion(const AnalyticFunction<VectorType>& f, const IntervalVector& X);
  Parallelepiped parallelepiped_inclusion(const IntervalVector& z, const IntervalMatrix& JJf, const Matrix& Jf_punc, const AnalyticFunction<VectorType>& psi_0, const OctaSym& symmetry, const IntervalVector& X);

  vector<Parallelepiped> PEIBOS(const AnalyticFunction<VectorType>& f, const AnalyticFunction<VectorType>& psi_0, const vector<OctaSym>& symmetries, double epsilon, bool verbose = false);
  vector<Parallelepiped> PEIBOS(const AnalyticFunction<VectorType>& f, const AnalyticFunction<VectorType>& psi_0, const vector<OctaSym>& symmetries, double epsilon, const Vector& offset, bool verbose = false);  
}