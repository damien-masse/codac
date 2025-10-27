// Author : Maël GODARD
// Adapted from CAPD examples

#include <codac>
#include <codac-polytope.h>
#include <vector>
#include <set>

using namespace std;
using namespace codac2;


struct FacetComp {
bool operator()(const Facet &f1, const Facet &f2) const {
   for (int i=0;i<f1.row.size();i++) {
      if (f1.row[i]<f2.row[i]) return true;
      if (f2.row[i]<f1.row[i]) return false;
   }
   return f1.rhs<f2.rhs;
}
};

int main()
{
   std::cout << std::scientific << std::setprecision(9);
  

   std::vector<Facet> fc = read_ineFile("sampleh8.ine");
   std::cout << "fichier lu avec " << fc.size() << " facettes.\n";
   IntervalVector box = IntervalVector::constant(fc[0].row.size(),
		Interval(-10,10));
   std::vector<Facet> vide;
   DDbuildF2V build2(box.size(),box,vide);
   std::cout << "Initial build : " << build2 << std::endl;
   Index i=0;
/*
   for (const auto &f : fc) {
      std::cout << "add facet : " << f << std::endl;
      build2.add_facet(++i,f);
      std::cout << build2;
   }
*/
   for (int i=0;i<fc.size();i++) {
      std::cout << "add facet : " << fc[i] << std::endl;
      build2.add_facet(i+1,fc[i]);
      std::cout << build2;
   }
   return 0;
}
