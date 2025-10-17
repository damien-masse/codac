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
#include "codac2_Parallelepiped.h"
#include "codac2_peibos.h"

namespace codac2
{
  template<typename T>
    requires std::is_base_of_v<AnalyticTypeBase,T>
  template<typename... Args>
  inline Parallelepiped AnalyticFunction<T>::parallelepiped_eval(const Args&... x) const
  {
    this->check_valid_inputs(x...);
    assert_release(this->input_size() < this->output_size());

    auto Y = this->eval(x.mid()...);
    auto z = Y.mid();

    Matrix A = this->diff(x.mid()...).mid();

    // Maximum error computation
    double rho = error(Y, z, this->diff(x...), A, cart_prod(x...));

    // Inflation of the parallelepiped
    auto A_inf = inflate_flat_parallelepiped(A, cart_prod(x...).rad(), rho);

    return Parallelepiped(z, A_inf);
  }
}