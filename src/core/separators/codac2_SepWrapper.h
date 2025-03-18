/** 
 *  \file codac2_SepWrapper.h
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include "codac2_Sep.h"
#include "codac2_SepCtcPair.h"
#include "codac2_CtcWrapper.h"
#include "codac2_CtcUnion.h"

namespace codac2
{
  /// \cond SepWrapper<T>

  template<typename T>
  class SepWrapper : public Sep<SepWrapper<T>>
  {
    public:

      SepWrapper(const T& y)
        : Sep<SepWrapper<T>>(y.size()), _y(y)
      { }

    protected:

      const T _y;
  };

  /// \endcond

  template<>
  class SepWrapper<IntervalVector> : public SepCtcPair
  {
    public:

      SepWrapper(const IntervalVector& y)
        : SepCtcPair(complementary_union(y), CtcWrapper<IntervalVector>(y))
      { }

      BoxPair separate(const IntervalVector& x) const;

    protected:

      CtcUnion<IntervalVector> complementary_union(const IntervalVector& y) const
      {
        CtcUnion<IntervalVector> cu(y.size());
        for(const auto& complem_y : y.complementary())
          cu |= CtcWrapper<IntervalVector>(complem_y);
        return cu;
      }
  };
}