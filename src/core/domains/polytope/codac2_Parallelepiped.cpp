/** 
 *  codac2_Parallelepiped.cpp
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Maël Godard
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include "codac2_Parallelepiped.h"

using namespace std;
using namespace codac2;

Parallelepiped::Parallelepiped(const Vector& c_, const Matrix& A_)
  : c(c_), A(A_)
{
  assert_release(c.size() == A.rows());
  assert_release(A.cols() <= c.size() && "too many vectors, you are describing a zonotope");
}

Zonotope Parallelepiped::project(const vector<int>& indices) const
{
  assert_release(*std::min_element(indices.begin(), indices.end()) >= 0 && "indices out of range");
  assert_release(*std::max_element(indices.begin(), indices.end()) <= c.size() && "indices out of range");

  Matrix A_cropped (indices.size(), A.cols());
  Vector c_cropped (indices.size());

  for (size_t i = 0; i < indices.size(); ++i)
  {
    A_cropped.row(i) = A.row(indices[i]);
    c_cropped[i] = c[indices[i]];
  }

  return Zonotope(c_cropped, A_cropped);
}

void generate_vertices(int i, int n, const Vector& c, const Matrix& A, vector<Vector>& L_v)
{
  if (i == n)
  {
    L_v.push_back(c);
  }
  else if (i<n)
  {
    generate_vertices(i+1, n, c + A.col(i), A, L_v);
    generate_vertices(i+1, n, c - A.col(i), A, L_v);
  }
}

vector<Vector> Parallelepiped::vertices() const
{
  vector<Vector> L_v;
  generate_vertices(0, c.size(),c,A,L_v);

  return L_v;
}

IntervalVector Parallelepiped::bounding_box() const
{
  IntervalVector box (c);
  for (auto vertice : vertices())
  {
    box |= IntervalVector(vertice);
  }
  return box;
}

