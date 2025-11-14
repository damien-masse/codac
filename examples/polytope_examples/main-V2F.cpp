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
   Polytope pol = Polytope::from_extFile(namefile);
   
   std::vector<IntervalVector> lst = pol.compute_vertices();
   std::cout << "nb vertices : " << lst.size() << std::endl;
   for (int i=0;i<lst.size();i++) {
     std::cout << lst[i] << std::endl;
   }
   return 0;
}
