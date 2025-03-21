/** 
 *  codac2_CtcDet.cpp
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Quentin Brateau
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include "codac2_CtcDet.h"
#include "codac2_det.h"
#include "codac2_arith_sub.h"

using namespace std;
using namespace codac2;

void CtcDet::contract(IntervalVector& x) const
{
  assert_release(x.size() == 2 && "only 2d segments are supported");

  // Forward contraction
  IntervalVector a(_a), b(_b); // variables that may be contracted
  Interval y(_y);
  Interval i1 = x[0] - a[0];
  Interval i2 = Interval(b[1] - a[1]);
  Interval i3 = i1 * i2;
  Interval i4 = x[1] - a[1];
  Interval i5 = Interval(b[0] - a[0]);
  Interval i6 = i4 * i5;
  Interval i7 = i3 - i6;

  // Backward contraction
  i7 &= y;
  SubOp::bwd(i7, i3, i6);
  MulOp::bwd(i6, i4, i5);
  SubOp::bwd(i5, b[0], a[0]);
  SubOp::bwd(i4, x[1], a[1]);
  MulOp::bwd(i3, i1, i2);
  SubOp::bwd(i2, b[1], a[1]);
  SubOp::bwd(i1, x[0], a[0]);
}
