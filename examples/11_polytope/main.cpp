// Author : Maël GODARD
// Adapted from CAPD examples

#include <codac>
#include <codac-polytope.h>
#include <vector>
#include <set>

using namespace std;
using namespace codac2;

int main()
{
   std::cout << std::scientific << std::setprecision(9);
  

   std::vector<Facet> fc = read_ineFile("cubocta.ine");
   std::cout << "fichier lu avec " << fc.size() << " facettes.\n";
   IntervalVector box = IntervalVector::constant(fc[0].row.size(),
		Interval(-10,10));
   std::vector<Facet> vide;
   DDbuildF2V build(box.size(),box,vide);
   Index i=0;
   for (int i=0;i<fc.size();i++) {
      build.add_facet(i+1,fc[i]);
   }
   std::cout << "built \n";
   std::vector<std::vector<Vector>> facets3D=build_3Dfacets(build);
   std::cout << "nb facets : " << facets3D.size() << "\n";

   Figure3D fig_poly("3D polyhedron");

   Vector center = Vector::zero(3);
   Matrix transfo = Matrix::Identity(3,3);
   for (const std::vector<Vector> &vec : facets3D) {
      fig_poly.draw_polygon(center,transfo,vec);
   }

   return 0;
}
