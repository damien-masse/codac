/**
 * \file codac2_polytope_clp.h
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Damien Massé
 *  \copyright  Copyright 2025
 *  \license    This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 */

#pragma once

#include <vector>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <ctime>
#include <cmath>
#include <bitset>
#include <memory>

#include "codac2_Index.h"
#include "codac2_Matrix.h"
#include "codac2_Vector.h"
#include "codac2_Row.h"
#include "codac2_IntervalRow.h"
#include "codac2_IntervalVector.h"
#include "codac2_IntervalMatrix.h"

class ClpSimplex; /* or #include <ClpSimplex.hpp> */

namespace codac2 {

/** \class LPclp
 *  \brief LPclp class, to modelise different LP problems
 *  
 */
class LPclp {
    public:
      /** \brief Constructor of the whole problem
       *   
       *  LP problem : max objvect X s.t. mat X <= rhsvect  
       *
       *  \param mat Matrix of the constraints (constraints in rows)
       *  \param objvect Vector of objectives (length = mat.cols())
       *  \param rhsvect RHS of the contraints  (length = mat.rows())
       */
      LPclp (const Matrix &mat, const Row &objvect,
		const Vector &rhsvect);

      /** \brief Constructor of the constraints, objective = 0
       *   
       *  LP constraint : mat X <= rhsvect  
       *  Can be used either to enter the objective afterwards
       *  or the check non-emptiness, or minimize the set of constraints
       *
       *  \param mat Matrix of the constraints (constraints in rows)
       *  \param rhsvect RHS of the contraints  (length = mat.rows())
       */
      LPclp (const Matrix &mat, const Vector &rhsvect);

      /** \brief Empty constructors, with only the size of the problem
       *  \param dim number of variables
       *  \param nbcsts number of constraints
       */
      LPclp(Index dim, Index nbcsts);

      /** \brief destructor */
      ~LPclp();

      /** \brief Add a constraint
       *  \param vst the (row) vector 
       *  \param rhs its RHS 
       *  \return the row number of the constraint */
      Index addConstraint(const Row &vst, double rhs);
      /** \brief Set the objvective
       *  \param objvect the new objective
       *  \param rhs its RHS */
      void setObjective(const Row &objvect);
      /** \brief Set the ''activity'' of a constraint 
       *  inactive constraints <=> rhs put to oo
       *  \param cst row number of the constraint
       *  \param state new state (default true)
       */
      void setActive(Index cst, bool state=true);
      /** \brief Get the ''activity'' of a constraint 
       *  inactive constraints <=> rhs put to oo
       *  \param cst row number of the constraint
       *  \return the ''activity'' of the constraint
       */
      bool getActive(Index cst) const;
      /** \brief Get the ''redundant'' status of a constraint 
       *  can be used after a call to minimize()
       *  \param cst row number of the constraint
       *  \return the ''activity'' of the constraint
       */
      bool isRedundant(Index cst) const;

      /** \brief Get the matrix of constraint 
         \return the matrix 
      */
      const Matrix &getMat() const;
      /** \brief Get the objective row
       *  \return the objective 
       */
      const Row &getObjvect() const;
      /** \brief Get the RHS vector 
       *  \return the rhs vector
       */
      const Vector &getRhsvect() const;
      /** \brief Get a single constraint, as an IntervalRow
       *  \return the rhs vector
       */
      Row getCst(Index i) const;

      /** \brief Status a constraint (internal)
       *  \arg \c REMOVED Constraint removed or never added
       *  \arg \c INACTIVE Inactive constraint 
       *  \arg \c REDUNDANT Detected as redundant by \c minimize()
       *  \arg \c LOCKED Forced in (not used)
       *  \arg \c SETIN Equality (not used)
       */
      enum CSTSTATUS {
          REMOVED,   /* constraint has been removed (or never added) */
 	  INACTIVE,  /* inactive, can be made active */
          REDUNDANT, /* shown redundant wrt other constraints */
          LOCKED,    /* must be used for this objective */
          SETIN,     /* equality forced */
          SIZE_CSTSTATUS
      };
      using lp_cststatus = std::bitset<SIZE_CSTSTATUS>;

      /** \brief Results of a solving process
       *  Different status are defined, with respect to 
       *  primal results (existence of points in the polytope) 
       *  and dual results (existence of constaints matching the 
       *  objective). Also, results about the basis are computed.
       *  \c LPRESULT is used in a bitset to describes the different
       *  results.
       *  \arg \c CHANGED No result because not computed or modified
       *  \arg \c EMPTY Polytope is empty 
       *  \arg \c EMPTY_APPROX LP-coin says empty, but the verification
       *          was not conclusive
       *  \arg \c EMPTY_BBOX Emptiness was proved thanks to the bbox
       *  \arg \c NOTEMPTY The analysis gave a guaranteed point
       *  \arg \c NOTEMPTY_APPROX LP-coin says non-empty, but the 
       *          verification of its result is not clear
       *  \arg \c UNBOUNDED Polytope is unbounded
       *  \arg \c UNBOUNDED_APPROX LP-coin says unbounded, but the
       *          verification was not conclusive
       *  \arg \c BOUNDED A bound is found for the objective
       *  \arg \c BOUNDED_APPROX LP-coin gave a bound, but the
       *          verification was not conclusive
       *  \arg \c BOUNDED_BBOX Bound was proved thanks to the bbox
       *  \arg \c ERROR_LPCOIN LP-coin returned an error
       *  \arg \c ERROR_PRIMAL_CHECK The checking of LP-coin primal (empty or not) result failed.
       *  \arg \c ERROR_DUAL_CHECK The checking of LP-coin dual (unbounded or not) result failed.
       */
      enum LPRESULT {
            CHANGED,
	    EMPTY,
            EMPTY_APPROX,
	    EMPTY_BBOX,
	    NOTEMPTY,
            NOTEMPTY_APPROX,
	    UNBOUNDED,
	    UNBOUNDED_APPROX,
	    BOUNDED,
	    BOUNDED_APPROX,
	    BOUNDED_BBOX,
	    ERROR_LPCOIN,
	    ERROR_PRIMAL_CHECK,
	    ERROR_DUAL_CHECK,
	    SIZE_LPRESULT
      };
      /** \code{.cpp} 
          using lp_result_stat = std::bitset<SIZE_LPRESULT>;
          \endcode */
      using lp_result_stat = std::bitset<SIZE_LPRESULT>;

      /** \brief Solve the problem, i.e. build it if needed, call lp-coin
       *  and check the result.
       *  \param checkempty check for emptiness: add a column in the model
       *       and maximize on this column
       *  \param option for now, the only values are 0 (normal) or 4 (restart
       *  from previous solutions if possible). */
      lp_result_stat solve(bool checkempty=false, int option=0);

      /** \brief detect redundant constraints, with a tolerance.
       * return the number of remaining constraints, -1 if the polytope
       * is empty.
       * if checkempty=true, we check emptiness first will all constraints.
       * tolerance works as follows : 
       *   max is the interval computed by maximizing without the constraint
       *   rhs the rhs of the constraint
       *       max_l <= max_u <= rhs => redundant
       *       max_l <= rhs <= max_u => redundant if rhs + tol_up >= max_u
       *       rhs <= max_l <= max_u => redundant if 
       *		rhs+tol_up >= max_u AND rhs+tol_down >= max_l
       */
      int minimize_polytope(const Interval &tolerance, bool checkempty=true);

      /** \brief Returns the value of the objective after
       *  solving. The value is guaranteed as an overapproximation.
       *  E.g. empty means that the polytope is empty, but +oo means only
       *  that the polytope may be unbounded */
      const Interval &getValobj() const;
      /** \brief Returns a solution point (as a ``small'' box)
       *  (with \c NOTEMPTY and \c NOTEMPTY_APPROX result). 
       * If a specific emptiness checking was done, then 
       * the point is in the ``middle'' of the polytope */
      const IntervalVector &getFeasiblePoint() const;
      /** \brief Returns a dual solution (as a ``small'' box)
       *  (with \c BOUNDED and \c BOUNDED_APPROX result). 
       * For a emptiness checking, the vector may prove the
       * global boundedness of the polytope */
      const IntervalRow &getBoundedRow() const;
      /** \brief Returns a emptiness vector
       *  (with \c EMPTY and \c EMPTY_APPROX result).  */
      const IntervalRow &getEmptyRow() const;
      /** \brief Returns a unboundness vector
       *  (with \c UNBOUNDED and \c UNBOUNDED_APPROX result).  */
      const IntervalVector &getUnboundedVect() const;

   private:
      /* initial problem */
      Index nbRows, nbCols;  /* dimension of the matrix */
      Matrix Amat;     /* the constraint matrix.*/
      Row objvect;  /* objective vector */
      Vector rhsvect;  /* rhs values */
      std::vector<lp_cststatus> cststat;

      lp_result_stat status;
      std::vector<bool> rowBasis;  /* limited to the rows */
      bool built;
      bool built_emptytest;
      ClpSimplex *model=nullptr;
      Interval Valobj;
      IntervalVector primalSol;
      IntervalRow dualSol;
      IntervalVector primalRay;
      IntervalRow dualRay;
      
      IntervalVector bbox;
      double timeout=4.0;
      double maxIteration=200;
      double tolerance=1e-9;

      void buildModel(bool emptytest); 
      void setModel(bool emptytest); 

      void correctBasisInverse(IntervalMatrix &basisInverse,
           const IntervalMatrix &basis) const;
      void reset_solution();
      lp_result_stat testGeneralApproach(int stat, bool checkempty);

      static IntervalMatrix infinite_sum_enclosure (const IntervalMatrix& M);
};


std::ostream& print_lp_stat(std::ostream& os, const LPclp::lp_result_stat x);
// std::ostream& operator<<(std::ostream& os, const LPclp &lp);

inline const Matrix &LPclp::getMat() const { return Amat; }
inline const Row &LPclp::getObjvect() const { return objvect; } 
inline const Vector &LPclp::getRhsvect() const { return rhsvect; }
inline Row LPclp::getCst(Index i) const { return Amat.row(i); }
inline const Interval &LPclp::getValobj() const { return Valobj; }
inline bool LPclp::isRedundant(Index cst) const { return cststat[cst][REDUNDANT]; }
inline const IntervalVector &LPclp::getFeasiblePoint() const { return primalSol; }
inline const IntervalRow &LPclp::getBoundedRow() const { return dualSol; }
inline const IntervalRow &LPclp::getEmptyRow() const { return dualRay; }
inline const IntervalVector &LPclp::getUnboundedVect() const { return primalRay; }

}

