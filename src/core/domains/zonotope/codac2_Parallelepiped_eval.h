/** 
 *  \file codac2_Parallelepiped_eval.h
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Maël Godard
 *  \copyright  Copyright 2025 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include "codac2_AnalyticFunction.h"
#include "codac2_Wrapper.h"
#include "codac2_Parallelepiped.h"
#include "codac2_peibos.h"

namespace codac2
{
  template<typename T>
  template<typename... Args>
    requires ((!std::is_same_v<MatrixType,typename ExprType<Args>::Type>) && ...)
  inline Parallelepiped AnalyticFunction<T>::parallelepiped_eval(const Args&... x) const
    requires std::is_same_v<VectorType,T>
  {
    this->check_valid_inputs(x...);
    assert_release(this->input_size() < this->output_size() &&
                   "Parallelepiped evaluation requires more outputs than inputs.");
    assert_release(this->input_size() > 0 &&
                "Parallelepiped evaluation requires at least one input.");

    IntervalVector Y = this->eval(((typename Wrapper<Args>::Domain)(x)).mid()...);
    Vector z = Y.mid();

    Matrix A = this->diff(((typename Wrapper<Args>::Domain)(x)).mid()...).mid();

    // Maximum error computation
    double rho = error(Y, z, this->diff(x...), A, cart_prod(x...));

    // Inflation of the parallelepiped
    Matrix A_inf = inflate_flat_parallelepiped(A, (cart_prod(x...).template cast<Interval>()).rad(), rho);

    return Parallelepiped(z, A_inf);
  }
}