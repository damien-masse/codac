// Author : Damien Massé
// Read a description of 3D polyhedron in inequalities form (.ine file)
// and draw it

#include <codac>
#include <codac-polytope.h>
#include <vector>
#include <memory>

using namespace std;
using namespace codac2;

int main(int argc, char *argv[])
{
   std::cout << std::scientific << std::setprecision(9);
  

   const char *namefile = "cubocta.ine";
   if (argc>1) namefile= argv[1];
   std::shared_ptr<CollectFacets> fc = 
		std::make_shared<CollectFacets>(read_ineFile(namefile));
   std::cout << "fichier lu avec " << fc->nbfcts() << " facettes.\n";
   IntervalVector box = IntervalVector::constant((*fc)[0]->first.row.size(),
		Interval(-10,10));
   std::vector<Facet> vide;
   DDbuildF2V build(box.size(),box,fc,false);
   Index i=0;
   for (int i=0;i<fc->nbfcts();i++) {
      build.add_facet(i);
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
