/** 
 *  Codac tests
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Damien Massé
 *  \copyright  Copyright 2025 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <catch2/catch_test_macros.hpp>
#include <codac2_Vector.h>
#include <codac2_Row.h>
#include <codac2_Matrix.h>
#include <codac2_clp.h>

using namespace std;
using namespace codac2;

TEST_CASE("Polytope-CLP-nobox")
{
	  /* non empty, bounded polytope : deformed octahedron not with // constraints */
      {
	  Matrix x({
	  	{ 1.5, 1.5, 1.4 },
	  	{ 0.5, 1.4, -0.5 },
	  	{ 1.4, -0.5, 0.5 },
	  	{ -0.5, 0.5, 1.5 },
	  	{ 0.5, -0.5, -1.5 },
	  	{ -0.5, -1.5, 0.5 },
	  	{ -1.5, 0.5, -0.5 },
	  	{ -1.5, -1.5, -1.5 }
	  });
	  Vector r({3.0,3.0,-1.0,3.0,-1.0,-1.0,3.0,-1.0});
	
	  LPclp clp(x,r);
	  LPclp::lp_result_stat res = clp.solve(true);
	  CHECK(res[LPclp::NOTEMPTY]);
	  IntervalVector diff = x*clp.getFeasiblePoint()-r;
	  CHECK(diff.ub().maxCoeff()<0.0);
	  CHECK(res[LPclp::BOUNDED]);
	  IntervalRow dual = clp.getBoundedRow().head(x.rows());
	  CHECK(dual.lb().maxCoeff()>=0.0);
	  CHECK((dual*x).mig().maxCoeff()<=1e-10);
	
	  Row obj({1.0,1.0,1.0});
	  clp.setObjective(obj);
	  res = clp.solve(false);
	  CHECK(res[LPclp::NOTEMPTY]);
	  diff = x*clp.getFeasiblePoint()-r;
	  CHECK(diff.ub().maxCoeff()<1e-8); /* can be >0 near the basis */
	  CHECK(res[LPclp::BOUNDED]);
	  dual = clp.getBoundedRow().head(x.rows());
	  CHECK(dual.lb().maxCoeff()>=0.0);
	  CHECK((dual*x-obj).mig().maxCoeff()<=1e-10);
  }

  {
	  /* non empty, bounded polytope : deformed octahedron with // constraints  */
	  Matrix x({
	  	{ 1.5, 1.5, 1.5 },
	  	{ 0.5, 1.5, -0.5 },
	  	{ 1.5, -0.5, 0.5 },
	  	{ -0.5, 0.5, 1.5 },
	  	{ 0.5, -0.5, -1.5 },
	  	{ -0.5, -1.5, 0.5 },
	  	{ -1.5, 0.5, -0.5 },
	  	{ -1.5, -1.5, -1.5 }
	  });
	  Vector r({3.0,3.0,-1.0,3.0,-1.0,-1.0,3.0,-1.0});
	
	  LPclp clp(x,r);
	  LPclp::lp_result_stat res = clp.solve(true);
	  CHECK(res[LPclp::NOTEMPTY]);
	  IntervalVector diff = x*clp.getFeasiblePoint()-r;
	  CHECK(diff.ub().maxCoeff()<0.0);
	  CHECK(res[LPclp::BOUNDED_APPROX]);
	  IntervalRow dual = clp.getBoundedRow().head(x.rows());
	  CHECK(dual.lb().maxCoeff()>=0.0);
	  CHECK((dual*x).mig().maxCoeff()<=1e-10);
	
	  Row obj({1.0,1.0,1.0});
	  clp.setObjective(obj);
	  res = clp.solve(false);
	  CHECK(res[LPclp::NOTEMPTY]);
	  diff = x*clp.getFeasiblePoint()-r;
	  CHECK(diff.ub().maxCoeff()<0.0);
	  CHECK(res[LPclp::BOUNDED]);
	  dual = clp.getBoundedRow().head(x.rows());
	  CHECK(dual.lb().maxCoeff()>=0.0);
	  CHECK((dual*x-obj).mig().maxCoeff()<=1e-10);
  }

  {
	  /* empty polytope : deformed octahedron with non// constraints */
	  Matrix x({
	  	{ 1.5, 1.5, 1.4 },
	  	{ 0.5, 1.4, -0.5 },
	  	{ 1.4, -0.5, 0.5 },
	  	{ -0.5, 0.5, 1.5 },
	  	{ 0.5, -0.5, -1.5 },
	  	{ -0.5, -1.5, 0.5 },
	  	{ -1.5, 0.5, -0.5 },
	  	{ -1.5, -1.5, -1.5 }
	  });
	  Vector r({3.0,3.0,-1.0,0.9,-1.0,-1.0,-2.0,-1.0});
	
	  LPclp clp(x,r);
	  LPclp::lp_result_stat res = clp.solve(true);
	  CHECK(res[LPclp::EMPTY]); 
          IntervalRow dual = clp.getEmptyRow().head(x.rows());
          CHECK(dual.lb().minCoeff()>=0.0);
          CHECK((dual*x).mig().maxCoeff()<=1e-10);
          CHECK((dual.dot(r)).ub()<0.0);
  }

  {
	  /* empty polytope : deformed octahedron with // constraints  */
	  Matrix x({
	  	{ 1.5, 1.5, 1.5 },
	  	{ 0.5, 1.5, -0.5 },
	  	{ 1.5, -0.5, 0.5 },
	  	{ -0.5, 0.5, 1.5 },
	  	{ 0.5, -0.5, -1.5 },
	  	{ -0.5, -1.5, 0.5 },
	  	{ -1.5, 0.5, -0.5 },
	  	{ -1.5, -1.5, -1.5 }
	  });
	  Vector r({3.0,3.0,-1.0,0.9,-1.0,-1.0,-2.0,-1.0});
	
	  LPclp clp(x,r);
	  LPclp::lp_result_stat res = clp.solve(true);
	  CHECK(res[LPclp::EMPTY_APPROX]); /* TO CHANGE */
          IntervalRow dual = clp.getEmptyRow().head(x.rows());
//          CHECK(dual.lb().minCoeff()>=0.0);
          CHECK((dual*x).mig().maxCoeff()<=1e-10);
          CHECK((dual.dot(r)).ub()<0.0);
  }

  {
	  /* unbounded polytope : deformed octahedron  with non // constraints */
	  Matrix x({
	  	{ 1.5, 1.5, 1.6 },
	  	{ 0.5, 1.6, -0.5 },
	  	{ 1.6, -0.5, 0.5 },
	  	{ -0.5, 0.5, 1.5 },
	  	{ 0.5, -0.5, -1.6 },
	  	{ -1.6, 0.5, -0.6 },
	  });
	  Vector r({3.0,3.0,-1.0,3.0,-1.0,3.0});
	
	  LPclp clp(x,r);
	  LPclp::lp_result_stat res = clp.solve(true);
	  CHECK(res[LPclp::NOTEMPTY]);
	  IntervalVector diff = x*clp.getFeasiblePoint()-r;
	  CHECK(diff.ub().maxCoeff()<0.0);
	  CHECK(res[LPclp::UNBOUNDED]);
	  IntervalVector vect = clp.getUnboundedVect();
	  CHECK((x*vect).ub().maxCoeff()<=0.0);
 
	  Row obj({0.0,-1.0,-1.0});
	  clp.setObjective(obj);
	  res = clp.solve(false);
	  CHECK(res[LPclp::NOTEMPTY]);
	  diff = x*clp.getFeasiblePoint()-r;
	  CHECK(diff.ub().maxCoeff()<1e-10);
	  CHECK(res[LPclp::UNBOUNDED]); 
	  diff = x*clp.getUnboundedVect();
          CHECK(diff.ub().maxCoeff()<=1e-10);
  }

  {
	  /* unbounded polytope : deformed octahedron  with // constraints */
	  Matrix x({
	  	{ 1.5, 1.5, 1.5 },
	  	{ 0.5, 1.5, -0.5 },
	  	{ 1.5, -0.5, 0.5 },
	  	{ -0.5, 0.5, 1.5 },
	  	{ 0.5, -0.5, -1.5 },
	  	{ -1.5, 0.5, -0.5 },
	  });
	  Vector r({3.0,3.0,-1.0,3.0,-1.0,3.0});
	
	  LPclp clp(x,r);
	  LPclp::lp_result_stat res = clp.solve(true);
	  CHECK(res[LPclp::NOTEMPTY]);
	  IntervalVector diff = x*clp.getFeasiblePoint()-r;
	  CHECK(diff.ub().maxCoeff()<0.0);
	  CHECK(res[LPclp::BOUNDED_APPROX]);
			/* note : polytope is unbounded, but
			   the emptiness solving did not gave
			   an infinite direction because of // of some
			   constraints */
	  IntervalRow dual = clp.getBoundedRow().head(x.rows());
	  CHECK(dual.lb().maxCoeff()>=0.0);
	  CHECK((dual*x).mig().maxCoeff()<=1e-10);
 
	  Row obj({-1.0,-1.0,-1.0});
	  clp.setObjective(obj);
	  res = clp.solve(false);
	  CHECK(res[LPclp::NOTEMPTY]);
	  diff = x*clp.getFeasiblePoint()-r;
	  CHECK(diff.ub().maxCoeff()<1e-10);
	  CHECK(res[LPclp::UNBOUNDED_APPROX]); 
		/* can't prove that it is unbounded, once again due to //
		   of some constraints */
	  diff = x*clp.getUnboundedVect();
          CHECK(diff.ub().maxCoeff()<1e-10);
  }
}

TEST_CASE("Polytope-CLP-box")
{
	  /* non empty, bounded polytope */
      {
	  Matrix x({
	  	{ 1.5, 1.5, 1.5 },
	  	{ 0.5, 1.5, -0.5 },
	  	{ 1.5, -0.5, 0.5 },
	  	{ -0.5, 0.5, 1.5 }
	  });
	  Vector r({3.0,3.0,-1.0,3.0});
          IntervalVector b { { -5,5 } , {-5,5}, {-5,5 } };
	
	  LPclp clp(x,r);
          clp.set_bbox(b);
	  LPclp::lp_result_stat res = clp.solve(true);
	  CHECK(res[LPclp::NOTEMPTY]);
	  IntervalVector diff = x*clp.getFeasiblePoint()-r;
	  CHECK(diff.ub().maxCoeff()<0.0);
	  CHECK(res[LPclp::BOUNDED]);
	
	  Row obj({1.0,-1.0,1.0});
	  clp.setObjective(obj);
	  res = clp.solve(false);
	  CHECK(res[LPclp::NOTEMPTY]);
	  diff = x*clp.getFeasiblePoint()-r;
	  CHECK(diff.ub().maxCoeff()<1e-8); /* can be >0 near the basis */
	  CHECK(res[LPclp::BOUNDED]);
  }

	  /* empty polytope */
      {
	  Matrix x({
	  	{ 1.5, 1.5, 1.4 },
	  	{ 0.5, 1.4, -0.5 },
	  	{ 1.4, -0.5, 0.5 },
	  	{ -0.5, 0.5, 1.5 }
	  });
	  Vector r({3.0,3.0,-1.0,3.0});
          IntervalVector b { { 3,6 } , {0,6}, {3,6 } };
	
	  LPclp clp(x,r);
          clp.set_bbox(b);
	  LPclp::lp_result_stat res = clp.solve(true);
	  CHECK(res[LPclp::EMPTY]);
  }


}
