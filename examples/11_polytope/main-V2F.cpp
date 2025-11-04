// Author : Damien Massé
// Transformation vertices -> facets from a .ext file

#include <codac>
#include <codac-polytope.h>
#include <vector>
#include <set>
#include <memory>

using namespace std;
using namespace codac2;


int main(int argc, char *argv[])
{
   std::cout << std::scientific << std::setprecision(9);
  

   const char *namefile = "irbox20-4.ext";
   if (argc>1) namefile= argv[1];
   std::vector<Vector> vx = read_extFile(namefile);
   std::cout << "File read with " << vx.size() << " vertices.\n";
   std::vector<Facet> vide;
   DDbuildV2F build2(1,vx[0]);
   std::cout << "Initial build : " << build2 << std::endl;
   for (int i=1;i<vx.size();i++) {
      std::cout << "add vertice " << (i+1) << " : " << vx[i] << std::endl;
      build2.add_vertex(i+1,vx[i]);
      std::cout << build2;
   }
      std::cout << build2;
   return 0;
}
