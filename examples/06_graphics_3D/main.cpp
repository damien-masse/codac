// The generated .obj files can be visualized on https://3dviewer.net/

#include <codac>

using namespace std;
using namespace codac2;

int main()
{
  VectorVar x(3);

  // Example from the publication:
  // https://www.ensta-bretagne.fr/jaulin/paper_centeredActa.pdf

  AnalyticFunction f { {x},
    {
      -sqr(x[2])+2*x[2]*sin(x[2]*x[0])+cos(x[2]*x[1]),
      2*x[2]*cos(x[2]*x[0])-sin(x[2]*x[1])
    }
  };
  

  CtcInverse_ ctc(f, IntervalVector::zero(2));
  auto p_ctc = pave({{0,2},{2,4},{0,10}}, ctc, 0.02);
  Figure3D fig_ctc("Paving contractor");
  fig_ctc.draw_paving(p_ctc);


  SepInverse sep_ellipsoid(
      AnalyticFunction({x}, 2*sqr(x[0])+x[0]*x[1]+x[0]*x[2]+sqr(x[1])+sqr(x[2])),
      Interval(0.4,1));

  auto p_sep = pave({{-1.1,1.1},{-1.1,1.1},{-1.1,1.1}}, sep_ellipsoid&IntervalVector({{-1.5,1.5},{-0.1,0.1},{-1.5,1.5}}), 0.04);

  Figure3D fig_sep("Paving separator");
  fig_sep.draw_axes(0.4);
  fig_sep.draw_paving(p_sep);

}