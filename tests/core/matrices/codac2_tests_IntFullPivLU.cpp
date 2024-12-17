/** 
 *  Codac tests
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <catch2/catch_test_macros.hpp>
#include <codac2_IntFullPivLU.h>
#include <codac2_BoolInterval.h>

using namespace std;
using namespace codac2;

TEST_CASE("IntFullPivLU")
{
  Matrix M = Matrix::Random(14,18); /* we suppose that M is full rank... */
  IntFullPivLU LUdec(M);
  CHECK(LUdec.isInjective()==BoolInterval::FALSE);
  CHECK(LUdec.isSurjective()==BoolInterval::TRUE);
  CHECK(LUdec.rank()==Interval(14));
  IntervalMatrix K = LUdec.kernel();
  CHECK(K.cols()==4);
  CHECK((IntervalMatrix(M)*K).lpNorm<Eigen::Infinity>().ub()<=1e-10);

  IntervalMatrix I1 = LUdec.solve(IntervalMatrix::Identity(14,14));
  CHECK((IntervalMatrix(M)*I1-IntervalMatrix::Identity(14,14)).lpNorm<Eigen::Infinity>().ub()<=1e-10);
  

  Matrix M2 = Matrix::Random(18,14); /* we suppose that M2 is full rank too */
  IntFullPivLU LUdec2(M2);
  CHECK(LUdec2.isInjective()==BoolInterval::TRUE);
  CHECK(LUdec2.isSurjective()==BoolInterval::FALSE);
  CHECK(LUdec2.rank()==Interval(14));
  IntervalMatrix K2 = LUdec2.kernel();
  CHECK(K2.cols()==0);
 
  

  Matrix M3 
     { { 1,1,0,1,2,0,3,0 },
       { 1,1,0,2,2,0,4,0 },
       { 1,2,3,4,5,6,7,8 },
       { 2,2,0,3,4,0,7,0 } };
  IntFullPivLU LUdec3(M3);
  Matrix M4 
     { { 0,3,2 },
       { 0,4,2 },
       { 8,7,5 },
       { 0,7,4 } };
  CHECK(LUdec3.image(M3)==M4); /* completely "ad hoc" test */
}
