/** 
 *  Codac tests
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Damien Massé
 *  \copyright  Copyright 2025 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <catch2/catch_test_macros.hpp>
#include <codac2_IntLDLT.h>
#include <codac2_BoolInterval.h>
#include <iostream>

using namespace std;
using namespace codac2;

TEST_CASE("IntLDLT")
{
  /* ``classic'' invertible matrix */
  {
    IntervalMatrix M 
    { { 16, {1,2},  5, -1 },
      { {1,2},  14, -4, 6 },
      { 5,  -4,  15, 9 },
      { -1, 6,  9, 19 }};
    IntLDLT LDLTdec(M);
    std::cout << LDLTdec.matrixLDLT() << "\n";
    std::cout << LDLTdec.ReconstructedMatrix() << "\n";
    IntervalMatrix M2
    { { 3, 0, -2, 3},
      { 2, 6, 0, 3},
      { 0, 2, 1, 3},
      { 0,0, 9, 3 } };
    IntervalMatrix M3 = LDLTdec.solve(M2);
    std::cout << (M*M3) << "\n";
  
  }
  {
    Matrix M 
    { { 6, 1, -1},
      { 1, 12, 1},
      { -1,1, 0 } };
    IntLDLT LDLTdec(M);
    std::cout << LDLTdec.matrixLDLT() << "\n";
    std::cout << LDLTdec.ReconstructedMatrix() << "\n";
    IntervalMatrix M2
    { { 3, 0, 0, 3},
      { 0, 6, 0, 3},
      { 0,0, 9, 3 } };
    std::cout << (IntervalMatrix(M)*LDLTdec.solve(M2)) << "\n";
  }
}
