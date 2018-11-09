#include "ibex.h"
#include "tubex.h"
#include "tubex-solver.h"
#include <iomanip>

using namespace std;
using namespace ibex;
using namespace tubex;

void contract(TubeVector& x)
{
  tubex::Function f("x", "-sin(x)");

  CtcPicard ctc_picard;
  ctc_picard.contract(f, x, FORWARD);

  CtcDeriv ctc_deriv;
  ctc_deriv.set_fast_mode(true);
  ctc_deriv.contract(x, f.eval_vector(x), FORWARD | BACKWARD);
}

int main()
{
  /* =========== PARAMETERS =========== */

    Tube::enable_syntheses(false);
    Vector epsilon(1, 0.05);
    Interval domain(0.,10.);
    TubeVector x(domain, 0.01, 1);
    TrajectoryVector truth(domain, tubex::Function("2.*atan(exp(-t)*tan(0.5))"));
    x.set(IntervalVector(truth(Interval(0.))), 0.); // initial condition
    // Note: use truth(Interval(0.)) instead of truth(0.) for a reliable evaluation

  /* =========== SOLVER =========== */

    tubex::Solver solver(epsilon);
    solver.set_refining_fxpt_ratio(0.9);
    solver.set_propa_fxpt_ratio(1.);
    solver.set_cid_fxpt_ratio(0.);
    solver.figure()->add_trajectoryvector(&truth, "truth");
    list<TubeVector> l_solutions = solver.solve(x, &contract);


  // Checking if this example still works:
  return (solver.solutions_contain(l_solutions, truth) != NO) ? EXIT_SUCCESS : EXIT_FAILURE;
}