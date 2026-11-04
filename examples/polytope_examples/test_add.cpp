// Author : Damien Massé
// Transformation facets to vertex from a .ine file

#include <codac>
#include <vector>
#include <set>
#include <memory>

using namespace std;
using namespace codac2;

int main(int argc, char *argv[])
{
   int64_t n;
   size_t a = SetBitset_unit::base_offset(-6,n);
   std::cout << poset_size_bitset << "\n";
   std::cout << ((-6)-((int64_t)poset_size_bitset-1))/(int64_t)poset_size_bitset << "\n";
   std::cout << (-6) << ":" << n << "," << a << " " << SetBitset_unit::value(n,a) << "\n"; 
   a = SetBitset_unit::base_offset(-33,n);
   std::cout << (-33) << ":" << n << "," << a << " " << SetBitset_unit::value(n,a) << "\n"; 
   a = SetBitset_unit::base_offset(0,n);
   std::cout << (0) << ":" << n << "," << a << " " << SetBitset_unit::value(n,a) << "\n"; 
   a = SetBitset_unit::base_offset(4,n);
   std::cout << (4) << ":" << n << "," << a << " " << SetBitset_unit::value(n,a) << "\n"; 
   a = SetBitset_unit::base_offset(97,n);
   std::cout << (97) << ":" << n << "," << a << " " << SetBitset_unit::value(n,a) << "\n"; 
   PosetElem p;
   std::cout << p << "\n";
   p.add(-6);
   std::cout << p << "\n";
   p.add(-33);
   std::cout << p << "\n";
   p.add(0);
   std::cout << p << "\n";
   p.add(4);
   std::cout << p << "\n";
   p.add(97);
   std::cout << p << "\n";
   PosetElem q;
   std::cout << q << "\n";
   q.add(-6);
   std::cout << q << "\n";
   q.add(-33);
   std::cout << q << "\n";
   q.add(0);
   std::cout << q << "\n";
   q.add(4);
   std::cout << q << "\n";
   q.add(97);
   std::cout << q << "\n";
   return 0;
}
