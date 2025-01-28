from codac import *

x = VectorVar(3)

# Example from the publication:
# https://www.ensta-bretagne.fr/jaulin/paper_centeredActa.pdf
f = AnalyticFunction([x], [
  -sqr(x[2])+2*x[2]*sin(x[2]*x[0])+cos(x[2]*x[1]),
  2*x[2]*cos(x[2]*x[0])-sin(x[2]*x[1])
])

ctc = CtcInverse(f, [0,0])
p_ctc = pave([[0,2],[2,4],[0,10]], ctc, 0.02)
fig_ctc=Figure3D("Paving contractor")
fig_ctc.draw_paving(p_ctc)


sep_ellipsoid=SepInverse (AnalyticFunction([x], 2*sqr(x[0])+x[0]*x[1]+x[0]*x[2]+sqr(x[1])+sqr(x[2])),Interval(0.4,1))
p_sep = pave([[-1.1,1.1],[-1.1,1.1],[-1.1,1.1]], sep_ellipsoid&IntervalVector([[-1.5,1.5],[-0.1,0.1],[-1.5,1.5]]), 0.04)
fig_sep = Figure3D("Paving separator")
fig_sep.draw_axes(0.4)
fig_sep.draw_paving(p_sep)