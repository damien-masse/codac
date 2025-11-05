/**
 *  codac2_peibos_capd.cpp
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Maël Godard
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <cassert>
#include "codac2_peibos_capd.h"

using namespace std;

namespace codac2
{
  std::vector<Parallelepiped> PEIBOS(const capd::IMap& i_map, double tf, const AnalyticFunction<VectorType>& psi_0, const std::vector<OctaSym>& Sigma, double epsilon, bool verbose)
  {
    return PEIBOS(i_map, tf, psi_0, Sigma, epsilon, Vector::zero(psi_0.output_size()), verbose);
  }

  std::vector<Parallelepiped> PEIBOS(const capd::IMap& i_map, double tf, const AnalyticFunction<VectorType>& psi_0, const std::vector<OctaSym>& Sigma, double epsilon, const Vector& offset, bool verbose)
  {
    int m = psi_0.input_size();
    int n = psi_0.output_size();

    assert(offset.size() == n);
    assert (m < n);
    assert (Sigma.size() > 0 && (int) Sigma[0].size() ==  n);

    clock_t t_start = clock();

    std::vector<Parallelepiped> output;
    
    // CAPD solver setup
    capd::IMap g (i_map);
    capd::IOdeSolver solver(g, 30);
    
    solver.setAbsoluteTolerance(1e-20);
    solver.setRelativeTolerance(1e-20);

    capd::ITimeMap timeMap(solver);
    capd::ITimeMap timeMap_punct(solver);

    capd::interval initialTime(0.);
    capd::interval finalTime(tf);

    std::vector<IntervalVector> boxes;
    double true_eps = split(Interval(-1.,1.)*IntervalVector::Ones(m), epsilon, boxes);
    
    for (const auto& sigma : Sigma)
    {
      for (const auto& X : boxes)
      {

        // To get the flow function and its Jacobian (monodromy matrix) for [x]
        IntervalVector Y = sigma(psi_0.eval(X)) + offset;

        capd::IMatrix monodromyMatrix(n,n);
        capd::ITimeMap::SolutionCurve solution(initialTime); 
        capd::IVector c = to_capd(Y);

        capd::C1Rect2Set s(c);
        timeMap(finalTime, s, solution);
        capd::IVector result = timeMap(finalTime, s, monodromyMatrix);
        IntervalVector Jf = to_codac(monodromyMatrix);

        // To get the flow function and its Jacobian (monodromy matrix) for x_hat
        auto xc = X.mid();
        auto yc = (sigma(psi_0.eval(xc)) + offset).mid();

        capd::IMatrix monodromyMatrix_punct(n,n);
        capd::ITimeMap::SolutionCurve solution_punct(initialTime);
        capd::IVector c_punct = to_capd(IntervalVector(yc));

        capd::C1Rect2Set s_punct(c_punct);
        timeMap_punct(finalTime, s_punct, solution_punct);      
        capd::IVector result_punct = timeMap_punct(finalTime, s_punct, monodromyMatrix_punct);
        Vector Jf_tild = to_codac(monodromyMatrix_punct).mid();

        // Center of the parallelepiped
        IntervalVector z = to_codac(result);
        
        auto p = parallelepiped_inclusion(z, Jf, Jf_tild, psi_0, sigma, X);

        output.push_back(p);

      }
    }
    
    if (verbose)
    {
      printf("\nPEIBOS statistics:\n");
      printf("------------------\n");
      printf("Real epsilon: %.4f\n", true_eps);
      printf("Computation time: %.4fs\n\n", (double)(clock()-t_start)/CLOCKS_PER_SEC);
    }

    return output;
  }
}
