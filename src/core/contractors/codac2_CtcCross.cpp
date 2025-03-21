/** 
 *  codac2_CtcCross.cpp
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Quentin Brateau
 *  \copyright  Copyright 2025 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include "codac2_CtcCross.h"
#include "codac2_CtcDet.h"
#include "codac2_CtcUnion.h"
#include "codac2_CtcInter.h"

using namespace std;
using namespace codac2;

void CtcCross::contract(IntervalVector& x) const
{
  assert_release(x.size() == 2 && "only 2d segments are supported");
  IntervalVector a(_e[0]), b(_e[1]), r(_r); // variables that may be contracted

  CtcDet C1(a, b, Interval(-oo, 0));
  CtcDet C2(a-r, b-r, Interval(0, oo));
  CtcDet C3(a-r, a, Interval(-oo, 0));
  CtcDet C4(b-r, b, Interval(0, oo));

  CtcUnion C = C1 | C2 | C3 | C4;
  C.contract(x);
}

void CtcNoCross::contract(IntervalVector& x) const
{
  assert_release(x.size() == 2 && "only 2d segments are supported");
  IntervalVector a(_e[0]), b(_e[1]), r(_r); // variables that may be contracted

  CtcDet C1(a, b, Interval(0, oo));
  CtcDet C2(a-r, b-r, Interval(-oo, 0));
  CtcDet C3(a-r, a, Interval(0, oo));
  CtcDet C4(b-r, b, Interval(-oo, 0));

  CtcInter C = C1 & C2 & C3 & C4;
  C.contract(x);
}