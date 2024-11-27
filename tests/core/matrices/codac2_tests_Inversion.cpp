/** 
 *  Codac tests
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <catch2/catch_test_macros.hpp>
#include <codac2_Vector.h>
#include <codac2_Matrix.h>
#include <codac2_IntervalMatrix.h>
#include <codac2_Inversion.h>
#include <codac2_Approx.h>

using namespace std;
using namespace codac2;

TEST_CASE("Matrix")
{
  IntervalMatrix x({
  	{ {0.0,0.0}, {-0.1,-0.1}, {0.2,0.2}  },
  	{ {0.0,0.0}, {-0.2,-0.2}, {0.1,0.1} },
  	{ {0.1,0.1}, {-0.1,-0.1}, {0.1,0.1} }
  });

  double a=0.0;

  IntervalMatrix y = infinite_sum_enclosure(x,a);
  CHECK(Approx(Interval(a))==Interval(5.0)/Interval(2.0));

  IntervalMatrix z = inverse_enclosure(x);
  CHECK(z(1,2)==Interval(0));
  CHECK(Approx<IntervalMatrix>(z*x)==IntervalMatrix::Identity(3,3));

  Matrix u( {
     { 1,2,3 },
     { 1,3,5 },
     { 3,4,5 }
   });

   IntervalMatrix v=inverse_enclosure(u);
  CHECK((IntervalMatrix(u)*v).contains(Matrix::Identity(3,3)));

}
