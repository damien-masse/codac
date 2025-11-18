/**
 * \file codac2_Polytope.h
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
#include "codac2_clp.h"
#include "codac2_Polytope_util.h"
#include "codac2_Facet.h"
#include "codac2_Polytope_dd.h"


namespace codac2 {
  /**
   * \class Polytope
   * \brief Represents a bounded convex polytope as a set of constraints
   * and a bounding box (the bounding box is part of the constraints)  */
  class Polytope 
  {
    public:
        /**
         * \brief basic constuctor
         */ 
        Polytope();

        /**
         * \brief basic constuctor with a dimension, full
         */ 
        Polytope(Index dim);

        /**
         * \brief basic constuctor with a dimension, empty or full
         */ 
        Polytope(Index dim, bool empty);

        /**
         * \brief Constructs a bounded polytope representing a box
         * from a bounded interval vector
         */
        Polytope(const IntervalVector &box);

        /**
         * \brief Copy (no minimization)
         */
        Polytope(const Polytope &P);

        /** 
         * \brief Constructs a polytope from a set of vertices
         * for now, just take the bounding box of the vertices
         */
        Polytope(const std::vector<Vector> &vertices);

        /** 
         * \brief Constructs a polytope from a set of vertices (intervalVector)
         * and a set of linear forms, described as a CollectFacets
         * \param vertices the vertices
         * \param facetforms the facets description (the CollectFacets)
         */
        Polytope(const std::vector<IntervalVector> &vertices, 
		const CollectFacets &facetsform);

        /**
         * \brief from a parallelepiped
         * \param par the parallelepiped
         */
        Polytope(const Parallelepiped &par);

        /**
         * \brief from a zonotope
         * \param zon the zonotope
         */
        Polytope(const Zonotope &zon);

        /**
         * \brief Constructs a polytope from a bounding box and 
         * a set of constraints 
         */
        Polytope(const IntervalVector &box, 
                 const std::vector<std::pair<Row,double>> &facets, 
	         bool minimize=false);

        /** 
         * \brief build a polytope from an ine File (cddlib format)
         * \param filename the file */
	static Polytope from_ineFile(const char *filename);
        /** 
         * \brief build a polytope from an ext File (cddlib format)
         * \param filename the file */
	static Polytope from_extFile(const char *filename);

        /**
         * \brief Destructor
         */
        ~Polytope();

        /***************** ACCESS **************************/
        /**
         * \brief Get the dimension of the space 
         * 
         * \return the dimension of the space
         */
        Index dim() const;
        Index size() const;

        /**
         * \brief Get the dimension of the polytope (ie the dimension
         * of the space minus the equalities). -1 means that the polytope
         * is empty.
         * 
         * \return the dimension of the polytope.
         */
//        Index dim_polytope() const;

        /**
         * \brief Get the number of constraints
         * 
         * \return the number of facets (including equalities)
         */
        Index nbFacets() const;

        /**
         * \brief Get the number of equality constraints
         * 
         * \return the number of equalities
         */
        Index nbEqFacets() const;

        /**
         * \brief is empty
         */
        bool is_empty(bool check=true) const;

        /**
         * \brief has a flat dimension
         */
        bool is_flat() const;

        /** \brief ``component'' : interval of the bounding box 
         * (only const !!!). Update the bbox if needed.
         * \return enclosing of the component
         */
        const Interval& operator[](Index i) const;

        /**
         * ``middle'' (vector inside the polytope, if possible)
         */
        Vector mid() const;

        /** 
         * \brief bound a constraint (fast), either by a "close" one (+box)
         * or by the vertices if they are computed.
         * no update or LP are done 
         * return an interval I : row X is guaranteed to be less than I.ub()
         * and I.lb gives an indication of the precision of the computation
         * (i.e. if diam(I) is low, a constraint close to row was used
         * \param fbase the constraint, expressed as a FacetBase
         * \return the interval I */
        Interval fast_bound(const FacetBase &base) const;

        /** 
         * \brief bound a constraint, using either clp or DD
         * \param r the row
         * \return the max
         */
        double bound_row(const Row &r) const;
  
#if 0
        /** \brief ``distance'' from the satisfaction of a constraint 
         *  ie the inflation of the rhs needed to satisfy the constraint
         *  \param fc the constraint
         *  \return an interval
         */
        double distance_cst (const Facet &fc) const;

        /**
         * \brief relationships with a box (fast check)
         * \param p the box
         * \return polytope_inclrel checking inclusion and intersection */
        Facet_::polytope_inclrel relation_Box(const IntervalVector& p) const;
#endif

        /**
         * \brief Checks whether the polytope contains a given point, or
         *  includes a box
         *
         * \param p The point or box to check, enclosed in an IntervalVector.
         * \return BoolInterval indicating possible containment.
         */
        BoolInterval contains(const IntervalVector& p) const;

        /**
         * \brief Checks inclusion in another polytope
         * \param P the polytope
         * \param checkF2V checks using DD
         * \param checkCLP checks using CLP
         * \return true if inclusion is guaranteed
         */
        bool is_included(const Polytope &P, 
			bool checkF2V=true, bool checkCLP=true) const;

        /** \brief intersects a box 
         * \param x the box (IntervalVector) 
         * \return if the polytope intersects the box
         */
        BoolInterval intersects(const IntervalVector& x) const;
        
        /** \brief intersects a polytope
         * \param x the polytope
         * \return if the polytope intersects the box
         */
        BoolInterval intersects(const Polytope &p) const;

        /** minimize the constraints, removing (possibly) redundant
         *  constraints.
         */
        void minimize_constraints() const;

        /************* Box access *************************/
 
        /**
         * \brief get the current bounding box of the polytope
         * (note: it may not be tightest bounding box, it is the 
         *  bounding box used for approximations results) 
         *
         * \return The current bounding box
         */
        const IntervalVector &box(bool tight=true) const;

        /**
         * \brief test if bounding box is included
         * \return true if the polytope is included in x
         */
        bool box_is_included(const IntervalVector& x) const;

        /** 
         * \brief tests if the polytope is bisectable 
	 * (i.e. its bounding box is bisectable).
         * do not update the box
         * \return true if the bbox is bisectable
         */
        bool is_bisectable() const;


        /************* Modification **********************/
        /** modification creates a new polytope         **/

        /** \brief set to empty */
        void set_empty();

        /** \brief set to (singleton) 0^dim */
        void clear();

        /** \brief add a inequality (pair row X <= rhs)
         *  do basic checks, but do not minimize the system
         *  \param facet the constraint
         *  \return false if the (basic) checks showed the constraint to
         *  be redundant (or inside tolerance), true if added */
        bool add_constraint(const std::pair<Row,double>& facet,
		double tolerance=0.0);

        /** \brief add a inequality with intervalVector (cst X <= rhs)
         *  using the bounding box
         *  do basic checks, but do not minimize the system
         *  \param cst the row constraint
         *  \param x the rhs
         *  \return false if the (basic) checks showed the constraint to
         *  be redundant (or inside tolerance), true if added */
        bool add_constraint(const IntervalRow &cst, double rhs,
		double tolerance=0.0);

        /** \brief two inequalities with intervalVector (cst X in rhs)
         *  using the bounding box
         *  do basic checks, but do not minimize the system
         *  \param cst the row constraint
         *  \param x the rhs
         *  \return pair of booleans (one for each constraints possibly added */
        std::pair<bool,bool> add_constraint_band(const IntervalRow &cst,
		 const Interval &rhs, double tolerance=0.0);

        /** \brief intersect with a box.
         *  \param b the box
         *  \return 0 if nothing done, 1 changed made, -1 results is empty */
        int meet_with_box(const IntervalVector &b);

        /** \brief intersect with a box.
         *  \param b the box
         *  \return *this */
        Polytope &operator&=(const IntervalVector &b);

        /** \brief intersect with a polytope.
         *  \param P the polytope
         *  \return 0 if nothing done, 1 changed made, -1 results is empty */
        int meet_with_polytope(const Polytope &P);

        /** \brief intersect with a polytope.
         *  \param P the polytope
         *  \return *this */
        Polytope &operator&=(const Polytope &P);

        /** \brief inflation by a cube
         *  this <- this + [-rad,rad]^d
         *  \param rad radius of the box
         *  \return the
         */
        Polytope inflate(double rad);

        /** \brief inflation by a ball while keeping the set of facets
         *  this <- this + B_d(rad)
         *  \param rad radius of the ball
         *  \return *this
         */
        Polytope inflate_ball(double rad);

        /** \brief expansion of a dimension (flat or not) 
         *  \param dm index of the dimension
         *  \param rad radius 
         *  \return *this */
        Polytope unflat(Index dm, double rad);

        /** \brief centered homothety
         *  x <- [c] + delta*([x]-[c]) ( (1-delta)[c] + delta*[x] )
         *  \param c ``center''
         *  \param delta expansion
         *  \return this
         */
        Polytope homothety(IntervalVector c, double delta);


        /*********** Printing and other ``global access'' *********/
        /**
         * \brief output the bbox and the set of facets 
         * \param os the output stream
         * \param P the polytope
         * \return the stream
         */
        friend std::ostream& operator<<(std::ostream& os,
			const Polytope &P);

        /**
         * \brief Computes a set of vertices enclosing the polytope
         * 
         * \return set of vertices, as IntervalVectors
         */
        std::vector<IntervalVector> compute_vertices() const;

        /**
         * \brief Computes the set of 3D facets
         * 
         * \return set of set of vertices, as Vectors
         */
        std::vector<std::vector<Vector>> compute_3Dfacets() const;


        /** \brief state of the polytope */
        enum POLSTATE {
            EMPTY,       /* is empty */
            NOTEMPTY,    /* is NOT empty */
            MINIMIZED,   /* minimized : redundant constraints removed */
            BOXUPDATED,  /* box is minimal */
            CLPFORM,     /* has an up-to-date clp form */
            F2VFORM,     /* has an up-to-date F2V form */
            V2FFORM,     /* has an up-to-date V2F form */
            INVALID,     /* not initialised */
            SIZE_POLSTATE
        };
        using pol_state = std::bitset<SIZE_POLSTATE>;
        static constexpr pol_state pol_state_empty = 
			1<<EMPTY | 1<<MINIMIZED | 1<<BOXUPDATED;
        static constexpr pol_state pol_state_init = 
			1<<NOTEMPTY | 1<<MINIMIZED | 1<<BOXUPDATED;
        
             

   private:
      Index _dim;                      /* dimension */
      mutable IntervalVector _box;     /* bounding box */
      mutable std::shared_ptr<CollectFacets> _facets; 
	/* "native" formulation , may be shared by other formulations */
      mutable std::unique_ptr<LPclp> _clpForm; /* LPclp formulation, if used */
      mutable std::unique_ptr<DDbuildF2V> _DDbuildF2V; 
		/* DDbuildF2V formulation, if used */
      mutable std::unique_ptr<DDbuildV2F> _DDbuildV2F; 
		/* DDbuildV2F formulation, generally not used */
      mutable pol_state state = pol_state_init;
      void minimize_constraints_clp(const Interval &tolerance=Interval(0.0)) const;
      void minimize_constraints_F2V() const;
      void update_box_clp() const;
      void update_box_F2V() const;
      void update_box() const;
      bool check_empty_clp() const;
      bool check_empty_F2V() const;
      bool check_empty() const;
      double bound_row_clp(const Row &r) const;
      double bound_row_F2V(const Row &r) const;
      void build_DDbuildF2V() const;
      void build_clpForm() const;
      void set_empty_private() const;
};


inline const IntervalVector &Polytope::box(bool tight) const { 
      if (tight) this->update_box();
      return _box; 
}
inline Index Polytope::dim() const { return _dim; }
inline Index Polytope::size() const { return _dim; }
inline bool Polytope::is_empty(bool check) const { 
	if (!check) return state[EMPTY]; 
	return this->check_empty(); }
inline bool Polytope::is_flat() const { return state[EMPTY] ||
		_box.is_degenerated() || (_facets && _facets->nbeqfcts()>0); }
inline Index Polytope::nbFacets() const { 
 	if (state[EMPTY]) return -1;
        if (!_facets) return (2*_dim);
        return (2*_dim + _facets->nbfcts());
     }

inline void Polytope::set_empty_private() const {
   state = pol_state_empty;
   _box.set_empty();
   _clpForm=nullptr;
   _DDbuildF2V=nullptr;
   _DDbuildV2F=nullptr;
   _facets=nullptr;
}


inline void Polytope::set_empty() {
   this->set_empty_private();
}


}
