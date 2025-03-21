/** 
 *  SepCross.cpp
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Quentin Brateau
 *  \copyright  Copyright 2025 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include "codac2_SepCross.h"

using namespace std;

namespace codac2
{
  BoxPair SepCross::separate(const IntervalVector& x) const
  {
    return SepCtcPair::separate(x);
  }
}