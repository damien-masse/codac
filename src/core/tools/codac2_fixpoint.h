/** 
 *  \file codac2_fixpoint.h
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Simon Rohou
 *  \copyright  Copyright 2025 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

namespace codac2
{
  template<typename F, typename... X>
  void fixpoint(const F& contract, const X&... x)
  {
    return fixpoint(contract, 1e-2, x...);
  }

  template<typename F, typename... X>
  void fixpoint(const F& contract, double contraction_ratio, const X&... x)
  {
    double vol = -1., prev_vol;

    do
    {
      prev_vol = vol;
      contract();
      vol = 0.;
      ((vol += [&x]() {
        double w = x.volume(); return w == oo ? 0. : w;
      }()), ...);

    } while(prev_vol != vol);
    //} while((vol == 0. && prev_vol < 0.) || (prev_vol/vol) > contraction_ratio);
  }
}