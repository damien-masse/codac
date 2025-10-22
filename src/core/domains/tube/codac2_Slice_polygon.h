/** 
 *  \file codac2_Slice_polygon.h
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include <iostream>
#include <iomanip>
#include "codac2_cart_prod.h"
#include "codac2_CtcDeriv.h"

namespace codac2
{
  template<typename T>
  ConvexPolygon Slice<T>::polygon_slice(const Slice<T>& v) const
    requires std::is_same_v<T,Interval>
  {
    const Interval& t = this->t0_tf();
    Interval input = this->input_gate(), output = this->output_gate();
    Interval proj_output =  input + t.diam() * v;
    Interval proj_input  = output - t.diam() * v;

    return CtcDeriv::polygon_slice(
      t, this->codomain(),
      input, proj_input,
      output, proj_output,
      v.codomain());
  }
}