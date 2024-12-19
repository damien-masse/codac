/** 
 *  codac2_Trajectory_operator.cpp
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou, Luc Jaulin
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include "codac2_Trajectory_operator.h"

using namespace std;
using namespace codac2;


IntervalVector TrajectoryOp::fwd(const SampledTrajectory<Vector>& x1, const Interval& x2)
{
  return x1(x2);
}

VectorOpValue TrajectoryOp::fwd(const SampledTrajectory<Vector>& x1, const ScalarOpValue& x2)
{
  IntervalMatrix d(x1.size(),x2.da.cols());
  
  return {
    fwd(x1,x2.m),
    fwd(x1,x2.a),
    d,
    x2.def_domain
  };
}

void TrajectoryOp::bwd(
  [[maybe_unused]] const SampledTrajectory<Vector>& x1,
  [[maybe_unused]] const IntervalVector& y,
  [[maybe_unused]] Interval& x2)
{
  // todo
}