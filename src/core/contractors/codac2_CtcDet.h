/** 
 *  \file codac2_CtcDet.h
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Quentin Brateau
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include "codac2_Ctc.h"
#include "codac2_IntervalVector.h"

namespace codac2
{
  class CtcDet : public Ctc<CtcDet,IntervalVector>
  {
    public:
    
      //  Contractor on the constraint det([x] - [a] | [b] - [a]) = [y]
      CtcDet(const IntervalVector& a, const IntervalVector &b, const Interval &y)
        : Ctc<CtcDet,IntervalVector>(2), _a(a), _b(b), _y(y)
      {
        assert_release(a.size() == 2 && b.size() == 2 && "only 2d segments are supported");
      }

      void contract(IntervalVector& x) const;

    protected:

      const IntervalVector _a, _b;
      const Interval _y;
  };
}