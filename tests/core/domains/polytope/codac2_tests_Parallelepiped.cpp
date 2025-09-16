/** 
 *  Codac tests
 *  
 *  Most of these tests come from the IBEX library (Gilles Chabert)
 *  See more: https://ibex-lib.readthedocs.io
 *  They have been revised to fit with Codac (v2)
 *
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Gilles Chabert, (Simon Rohou)
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <catch2/catch_test_macros.hpp>
#include <codac2_Parallelepiped.h>
#include <limits>

using namespace std;
using namespace codac2;

TEST_CASE("Parallelepiped")
{
  Parallelepiped p(Vector({0,2,4}), Matrix({{0.5,0.,0.},{0.,1.,0.},{0.,1.,1.}}));
  CHECK(p.bounding_box() == IntervalVector({{-0.5,0.5},{1.,3.},{2.,6.}}));

  Zonotope z = p.project({2,1,0});
  CHECK(z.z == Vector({4,2,0}));
  CHECK(z.A == Matrix({{0.,1.,1.},{0.,1.,0.},{0.5,0.,0.}}));
}

