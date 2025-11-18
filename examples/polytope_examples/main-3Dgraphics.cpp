// Author : Damien Massé
// Construct different polytopes and draw them

#include <codac>
#include <codac-polytope.h>
#include <vector>
#include <memory>

using namespace std;
using namespace codac2;

void draw_polytope(Figure3D &fig, const Polytope &P, const StyleProperties &st) {
   std::vector<std::vector<Vector>> facets3D=P.compute_3Dfacets();
   Vector center = Vector::zero(3);
   Matrix transfo = Matrix::Identity(3,3);
   for (const std::vector<Vector> &vec : facets3D) {
      fig.draw_polygon(center,transfo,vec,st);
   }
}

int main(int argc, char *argv[])
{
   Figure3D fig("3D polytopes");

   /* box */
   IntervalVector box1 { { 1, 3 }, { 1, 2 }, { 2, 2.5 } };
   Polytope pbox1(box1);
   std::cout << "Box : " << pbox1 << "\n";
   draw_polytope(fig,pbox1,StyleProperties(Color::dark_green(),"box"));

   /* flat box */
   IntervalVector box2 { { 1, 3 }, { -1, -1 }, { 2, 2.5 } };
   Polytope pbox2(box2);
   std::cout << "Flat box : " << pbox2 << "\n";
   draw_polytope(fig,pbox2,StyleProperties(Color::dark_red(),"flat box"));

   /* parallelepiped */
   Polytope paral(Parallelepiped({-1,1,1}, {{ 0.5,0.6,0.0 },
					  { -0.3,0.0,0.4 },
					  { 0.0, -0.6, 0.8 }}));
   std::cout << "Parallelotop : " << paral << "\n";
   draw_polytope(fig,paral,StyleProperties(Color::dark_blue(),"parallelepiped"));

   /* zonotope */
   Polytope zono(Zonotope({-4,1,1}, {{ 0.5,0.6,0.0,0.4,0.0,1.1 },
					  { -0.3,0.0,0.4,-0.5,0.3,-0.3 },
					  { -0.2,0.0,0.3,-0.4,-0.2,-0.2 }}));
   std::cout << "Zonotope : " << zono << "\n";
   draw_polytope(fig,zono,StyleProperties(Color::dark_yellow(),"zonotope"));

   /* intersection */
   zono.meet_with_polytope(paral);

   std::cout << "Intersection : " << zono << "\n";
   draw_polytope(fig,zono,StyleProperties(Color::dark_gray(),"intersection"));
   zono.minimize_constraints();
   std::cout << "Intersection (minimized) : " << zono << "\n";

   return 0;
}
