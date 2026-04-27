// Author : Damien Massé
// Graphical illustration of the polytope test

#include <codac>
#include <vector>
#include <memory>

using namespace std;
using namespace codac2;

int main(int argc, char *argv[])
{
//   std::cout << std::scientific << std::setprecision(20);

   // definition using vertices
   std::vector<Vector> vertices 
	{ {0,0,0}, {0.2,0.1,0}, {-0.1,0.1,0}, {0,0.15,0.2}, {-0.15,0,-0.15}, {0,0.3,0} };
   Polytope p1(vertices);
   // definition using facets
   std::vector<std::pair<Row,double>> facets 
       { { {1,1,0.75}, 0.3 },
         { {-1,-0.4,0.6} ,0.06 },
         { {-1,0.5,0.375} ,0.15 },
         { {-1,0.5,0} ,0.15 },
         { {0.5,-1,0.75} ,0 },
         { {-0.75,-1,0.75} ,0 },
         { {0.5,-1,-0.5} ,0 },
         { {1.0/3.0,1.0/3.0,-1} ,0.1 } };
   // the first argument is a bounding box, here the whole space
   Polytope p2(IntervalVector::Constant(3,Interval()), facets);

   // p1 and p2 are ``almost'' the same polytope

   std::cout << p1 << std::endl;

   std::vector<IntervalVector> vp2 = p2.vertices();
   for (auto &v : vp2) std::cout << v << std::endl;
//   fig.draw_polytope(p2,StyleProperties(Color::dark_blue(0.8),"p2"));

   Figure3D fig("manual_polytope");
   fig.draw_axes(1.0);
   fig.draw_polytope(p1,StyleProperties(Color::dark_red(0.8),"p1"));

   VectorVar x(3);
   AnalyticFunction f { {x},
    {
          (1.0-x[0])/sqrt(2.0+sqr(1.0-x[0])+sqr(x[1])+sqr(x[2])),
          x[1]/sqrt(2.0+sqr(1.0-x[0])+sqr(x[1])+sqr(x[2])),
          x[2]/sqrt(2.0+sqr(1.0-x[0])+sqr(x[1])+sqr(x[2]))
    }
   };

   Polytope q = f.polytope_eval(p1);
   fig.draw_polytope(q,StyleProperties(Color::dark_blue(0.8),"q"));
   std::cout << q << "\n";
   q.minimize_constraints();
   std::cout << q << "\n";


   SepPolytope sep_q(q);
   SepInverse sep_f(f,sep_q);
   auto pv1 = pave(IntervalVector::Constant(3,Interval(-3,3)),sep_f,0.02); 
   fig.draw_paving(pv1);


   
   return 0;
}
