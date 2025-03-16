/** 
 *  codac2_SepWrapper.cpp
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include "codac2_SepWrapper.h"

using namespace std;
using namespace codac2;

BoxPair SepWrapper<IntervalVector>::separate(const IntervalVector& x) const
{
  return SepCtcPair::separate(x);
}

BoxPair SepWrapper_<PavingInOut>::separate(const IntervalVector& x) const
{
  assert(x.size() == P.size());

  BoxPair result = BoxPair(IntervalVector::empty(x.size()), IntervalVector::empty(x.size()));

  // Inner computation
  for (const auto &box: P.intersecting_boxes(x, PavingInOut::outer_complem))
  {
    result.inner |= (box & x);
  }

  // Outer computation
  for (const auto &box: P.intersecting_boxes(x, PavingInOut::inner))
  {
    result.outer |= (box & x);
  }

  // THe boundary belongs to both inner and outer
  for (const auto &box: P.intersecting_boxes(x, PavingInOut::bound))
  {
    result.inner |= (box & x);
    result.outer |= (box & x);
  }

  return result;
}