/** 
 *  codac2_Zonotope.cpp
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Maël Godard
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include "codac2_Zonotope.h"

using namespace std;
using namespace codac2;

Zonotope::Zonotope(const Vector& c_, const Matrix& A_)
  : c(c_), A(A_)
{
  assert_release(c.size() == A.rows());
}

