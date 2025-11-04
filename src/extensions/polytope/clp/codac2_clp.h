/**
 * \file codac2_clp.h
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
#include "codac2_BoolInterval.h"
#include "codac2_Facet.h"

class ClpSimplex;  /* instead of including ClpSimplex.hpp */
class CoinBuild;   

namespace codac2 {


/** \class LPclp
 *  \brief LPclp class, to modelise different LP problems
 *  
 */
class LPclp {
    public:
      /** \brief Constructor with the collection of facets 
       *
       * LP problem : max objvect X s.t. mat X <= rhsvect  
       *
       * \param dim dimension of the space
       * \param facets the facets 
       * \param objvect the objective */
      LPclp (Index dim, std::shared_ptr<CollectFacets> &facets,
	const Row &objvect);

      /** \brief Constructor with the collection of facets, no objective 
       *
       * LP problem : max objvect X s.t. mat X <= rhsvect  
       *
       * \param dim dimension of the space
       * \param facets the facets */
      LPclp (Index dim, std::shared_ptr<CollectFacets> &facets);

      /** \brief Constructor with the collection of facets, bbox, no objective 
       *
       * LP problem : max objvect X s.t. mat X <= rhsvect  
       *
       * \param dim dimension of the space
       * \param facets the facets
       * \param box the bounding box for Neumaier's approx */
      LPclp (Index dim, std::shared_ptr<CollectFacets> &facets,
		const IntervalVector &box);

      /** \brief Constructor from a matrix : create a new CollectFacets
       *   
       *  LP problem : max objvect X s.t. mat X <= rhsvect  
       *
       *  \param mat Matrix of the constraints (constraints in rows)
       *  \param objvect Vector of objectives (length = mat.cols())
       *  \param rhsvect RHS of the contraints  (length = mat.rows())
       *  \param eqSet rows which are equalities
       */
      LPclp (const Matrix &mat, const Row &objvect,
		const Vector &rhsvect, 
		const std::vector<Index> &eqSet=std::vector<Index>());

      /** \brief Constructor of the constraints, objective = 0
       *   
       *  LP constraint : mat X <= rhsvect  
       *  Can be used either to enter the objective afterwards
       *  or the check non-emptiness, or minimize the set of constraints
       *
       *  \param mat Matrix of the constraints (constraints in rows)
       *  \param rhsvect RHS of the contraints  (length = mat.rows())
       *  \param eqSet rows which are equalities
       */
      LPclp (const Matrix &mat, const Vector &rhsvect, 
		const std::vector<Index> &eqSet=std::vector<Index>());

      /** \brief Copy constructor
       *  Do not build the model, even if the original model is built
       *  also do not copy all results related to the last computation :
       *  status, row_basis, valObj, primalSol, primalRay, dualSol, dualRay
       */
      LPclp (const LPclp& P);

      /** \brief Empty constructor
       *  \param dim number of variables
       */
      LPclp(Index dim);

      /** \brief destructor */
      ~LPclp();

      /** \brief set bbox
       *  \param bbox the new bbox
       */
      void set_bbox(const IntervalVector &box);

      /** \brief Add a constraint (warning : only for "local" CollectFacets)
       *  \param facet the facet
       *  \return the row number of the constraint, or -1 if the
       *  insertion failed (an constraint with the same base exists) */
      Index addConstraint(const Facet &facet);
      /** \brief Add a constraint (warning : only for "local" CollectFacets)
       *  \param vst the (row) vector 
       *  \param rhs its RHS 
       *  \param isEq true if it is an equality
       *  \return the row number of the constraint, or the number of the
       *  existing constraint if a constraint with the same base was modified */
      Index addConstraint(const Row &vst, double rhs, bool isEq=false);
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

      /** \brief Get the set of facets 
         \return the facets
      */
      const CollectFacets &getFacets() const;
      /** \brief Get the collections of facets as pointer
         \return the facets
      */
      std::shared_ptr<CollectFacets> getFacetsPtr() ;
      /** \brief Get the objective row
       *  \return the objective 
       */
      const Row &getObjvect() const;
      /** \brief Build the matrix
       *  \return the matrix
       */
      Matrix getMat() const;
      /** \brief Build the RHS vector 
       *  \return the rhs vector
       */
      Vector getRhsvect() const;
      /** \brief Get a single constraint, as an facet
       *  \return the rhs vector
       */
      CollectFacets::mapCIterator getCst(Index i) const;

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

      /** \brief detect redundant equality constraints
       * returns the number of remaining equality constraints, -1 if empty
       * tolerance works as follows : 
       *   each new equality is added to a (pseudo)-triangular matrix
       *    of IntervalVector
       *   when a row with all coefficients containing 0 appears,
       *      a) if 0 is in the rhs, declare the row as redundant
       *      b) if not, and all coefficients are exactly 0, returns empty
       *      c) if bounding box, use it to try to prove emptiness, if not
       *         possible keep the row (and hope for emptiness later?)
       */
      int minimize_eqpolytope();

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
      std::shared_ptr<CollectFacets> Afacets; /* the facets */
      Row objvect;  /* objective vector */
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
 
      void addFacetToCoinBuild(CoinBuild &buildObject,
			const Facet &facet, lp_cststatus status,
			int *row2Index, double *row2Vals) const;

      void buildModel(bool emptytest); 
      void setModel(bool emptytest); 

      void correctBasisInverse(IntervalMatrix &basisInverse,
           const IntervalMatrix &basis) const;
      void reset_solution();
      Interval dotprodrhs(const IntervalRow &d) const;

      BoolInterval checkDualForm(IntervalRow &dualVect,
                IntervalRow &initSum, const IntervalMatrix &basisInverse,
                const std::vector<Index> &wholeBasis, Index nbRowsInBasis,
                bool checkempty) const;
      lp_result_stat testGeneralApproach(int stat, bool checkempty);


      static IntervalMatrix infinite_sum_enclosure (const IntervalMatrix& M);

};


std::ostream& print_lp_stat(std::ostream& os, const LPclp::lp_result_stat x);
// std::ostream& operator<<(std::ostream& os, const LPclp &lp);

// inline const Matrix &LPclp::getMat() const { return Amat; } /* TODO */
inline const Row &LPclp::getObjvect() const { return objvect; } 
// inline const Vector &LPclp::getRhsvect() const { return rhsvect; } /* TODO */
inline CollectFacets::mapCIterator 
	LPclp::getCst(Index i) const { return (*Afacets)[i]; }
inline const Interval &LPclp::getValobj() const { return Valobj; }
inline bool LPclp::isRedundant(Index cst) const { return cststat[cst][REDUNDANT]; }
 inline const IntervalVector &LPclp::getFeasiblePoint() const { return this->primalSol; }
inline const IntervalRow &LPclp::getBoundedRow() const { return this->dualSol; }
inline const IntervalRow &LPclp::getEmptyRow() const { return this->dualRay; }
inline const IntervalVector &LPclp::getUnboundedVect() const { return this->primalRay; }
inline void LPclp::set_bbox(const IntervalVector &box) { this->bbox=box; }

}

