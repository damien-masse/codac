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
  /**
   * \brief Recursively split a box until the maximum diameter of each box is less than eps. Note that the resulting boxes are stored in place in the vector boxes.
   * 
   * \param X The box to split
   * \param eps The maximum diameter of each box
   * \param boxes The vector to store the resulting boxes
   * 
   * \return The maximum diameter of the resulting boxes (lower or equal to eps)
   */
  double split (const IntervalVector& X, double eps, vector<IntervalVector>& boxes);

  /**
   * \brief Compute the error term for the parallelepiped inclusion. The error is later used to inflate the flat parallelepiped \f$\mathbf{a}+\mathbf{B}\cdot(\left[\mathbf{x}\right]-\bar{\mathbf{x}})\f$
   * 
   * \param z The box enclosing \f$\mathbf{f}(\bar{\mathbf{x}})\f$. \f$\mathbf{a}\f$ is its center
   * \param JJg The interval Jacobian matrix containing \f$\frac{d\mathbf{f}}{d\mathbf{x}}([\mathbf{x}])\f$
   * \param B Any Matrix with the same dimension as JJg. Usually we pick an approximation of \f$\frac{d\mathbf{f}}{d\mathbf{x}}(\bar{\mathbf{x}})\f$
   * \param X The box \f$[\mathbf{x}]\f$
   */
  double error(const IntervalVector& z, const IntervalMatrix& JJg, const Matrix& B, const IntervalVector& X);

  /**
   * \brief Inflate the flat parallelepiped \f$\mathbf{A}\cdot e_{vec}\f$ by \f$\rho\f$.
   * 
   * \param A The shape matrix of the flat parallelepiped (n x m matrix with m < n)
   * \param e_vec The vector of scaling along each generator of the flat parallelepiped (m-dimensional vector)
   * \param rho The inflation factor
   * 
   * \return The shape matrix of the inflated parallelepiped (n x n matrix)
   */
  Matrix inflate_flat_parallelepiped (const Matrix& A, const Vector& e_vec, double rho);

  /**
   * \brief Compute a parallelepiped enclosing of \f$\mathbf{f}([\mathbf{x}])\f$ when \f$\mathbf{f}:\mathbb{R}^m\rightarrow\mathbb{R}^n\f$ is an analytic function with \f$m < n\f$.
   * 
   * \param f The analytic function \f$\mathbf{f}\f$
   * \param X The box \f$[\mathbf{x}]\f$
   * 
   * \return A Parallelepiped enclosing \f$\mathbf{f}([\mathbf{x}])\f$
   */
  Parallelepiped parallelepiped_inclusion(const AnalyticFunction<VectorType>& f, const IntervalVector& X);

  /**
   * \brief Used in PEIBOS. Compute a parallelepiped enclosing of \f$\mathbf{g}([\mathbf{x}])\f$ where \f$\mathbf{g} = \mathbf{f}\circ \sigma \circ \psi_0\f$.
   * 
   * \param z The box enclosing \f$\mathbf{f}(\sigma(\psi_0(\bar{\mathbf{x}})))\f$.
   * \param JJf The interval Jacobian matrix containing \f$\frac{d\mathbf{f}}{d\mathbf{y}}([\mathbf{y}])\f$ where \f$\mathbf{y} = \sigma(\psi_0(\left[\mathbf{x}\right]))\f$
   * \param Jf_tild An approximation of \f$\frac{d\mathbf{f}}{d\mathbf{y}}(\sigma(\psi_0(\bar{\mathbf{x}})))\f$
   * \param psi_0 The analytic function \f$\psi_0:\mathbb{R}^m\rightarrow\mathbb{R}^n\f$ to construct the atlas
   * \param sigma The symmetry operator to construct the atlas
   * \param X The box \f$[\mathbf{x}]\f$
   * 
   * \return A Parallelepiped enclosing \f$\mathbf{g}([\mathbf{x}])\f$
   */
  Parallelepiped parallelepiped_inclusion(const IntervalVector& z, const IntervalMatrix& JJf, const Matrix& Jf_tild, const AnalyticFunction<VectorType>& psi_0, const OctaSym& sigma, const IntervalVector& X);

  vector<Parallelepiped> PEIBOS(const AnalyticFunction<VectorType>& f, const AnalyticFunction<VectorType>& psi_0, const vector<OctaSym>& symmetries, double epsilon, bool verbose = false);
  vector<Parallelepiped> PEIBOS(const AnalyticFunction<VectorType>& f, const AnalyticFunction<VectorType>& psi_0, const vector<OctaSym>& symmetries, double epsilon, const Vector& offset, bool verbose = false);  
}