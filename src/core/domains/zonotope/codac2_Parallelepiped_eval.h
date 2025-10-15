/** 
 *  \file codac2_Parallelepiped_eval.h
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Maël Godard
 *  \copyright  Copyright 2025 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include "codac2_Parallelepiped.h"

namespace codac2
{
  template<typename T>
    requires std::is_base_of_v<AnalyticTypeBase,T>
  inline Parallelepiped AnalyticFunction<T>::parallelepiped_eval(const IntervalVector& x) const
  {
    this->check_valid_inputs(x);
    assert_release(this->input_size() < this->output_size());
    assert_release(this->args().size() == 1 && "f must have only one arg");

    // todo
    
    return Parallelepiped({{0.}},{{0.}});
  }
}