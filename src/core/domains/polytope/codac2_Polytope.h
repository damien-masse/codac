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
#include <initializer_list>

#include "codac2_Index.h"
#include "codac2_Matrix.h"
#include "codac2_Vector.h"
#include "codac2_Row.h"
#include "codac2_IntervalRow.h"
#include "codac2_IntervalVector.h"
#include "codac2_IntervalMatrix.h"
#include "codac2_BoolInterval.h"
#include "codac2_Polytope_util.h"
#include "codac2_Facet.h"
#include "codac2_Polytope_dd.h"
#include "codac2_Parallelepiped.h"
#include "codac2_Zonotope.h"
#ifdef WITH_CLP
#include "../../../extensions/clp/codac2_clp.h"
#endif


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
         * \param vertices the vertices
         */
        Polytope(const std::vector<Vector> &vertices);

        /** 
         * \brief Constructs a polytope from a set of vertices as intervalboxes
         * \param vertices the vertices
         */
        Polytope(const std::vector<IntervalVector> &vertices);

        /** 
         * \brief Constructs a polytope from a set of vertices (intervalVector)
         * and a set of linear forms, described as a CollectFacets
         * \param vertices the vertices
         * \param facetsform the facets description (the CollectFacets)
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
         * \brief Constructs a polytope from a box and a CollectFacets which
	 * will be moved in the polytope. Perform a "box extraction" and a 
         * renumbering on the collection.
         * \param box the box
         * \param facets the collection of facets 
         */
        Polytope(IntervalVector &&box, CollectFacets &&facets);

        /**
         * \brief Constructs a polytope from a CollectFacets which
	 * will be moved in the polytope. Perform a "box extraction" and a 
         * renumbering on the collection.
         * \param facets the collection of facets 
         */
        Polytope(CollectFacets &&facets);

        /**
         * \brief Constructs a polytope from a bounding box and 
         * a set of constraints (inequalities)
 	 * 
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

        /** \brief convex union of polytopes
         *  use vertices to join the polytopes, as such the result may be
         *  too large
         *  \param lst a non-empty initializer_list of valid polytopes of 
	 *  same dimension
         *  \return the polytope */
        static Polytope union_of_polytopes(std::initializer_list<Polytope> lst);

        /** \brief convex union of polytopes
         *  use vertices to join the polytopes, as such the result may be
         *  too large
         *  \param lst a non-empty initializer_list of valid polytopes of 
	 *  same dimension
         *  \return the polytope */
        static Polytope union_of_polytopes(const std::vector<Polytope> &lst);

        /**
         * \brief Destructor
         */
        ~Polytope();

        /**
         * \brief copy assignment operator 
         * \param P copy
         */
        Polytope &operator=(const Polytope &P);
     


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
         * \param base the constraint, expressed as a FacetBase
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
        bool is_subset(const Polytope &P, 
			bool checkF2V=true, bool checkCLP=true) const;

        /** \brief intersects a box 
         * \param x the box (IntervalVector) 
         * \return if the polytope intersects the box
         */
        BoolInterval intersects(const IntervalVector& x) const;
        
        /** \brief intersects a polytope
         * \param p the polytope
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
        bool box_is_subset(const IntervalVector& x) const;

        /** 
         * \brief tests if the polytope is bisectable 
	 * (i.e. its bounding box is bisectable).
         * do not update the box
         * \return true if the bbox is bisectable
         */
        bool is_bisectable() const;


        /************* Modification **********************/
        /* keeping the current polytope                  */

        /** \brief set to empty */
        void set_empty();

        /** \brief set to (singleton) 0^dim */
        void clear();

        /** \brief add a inequality (pair row X <= rhs)
         *  do basic checks, but do not minimize the system
         *  \param facet the constraint
         *  \param tolerance tolerance for redundancy checking (CLP only)
         *  \return false if the (basic) checks showed the constraint to
         *  be redundant (or inside tolerance), true if added */
        bool add_constraint(const std::pair<Row,double>& facet,
		double tolerance=0.0);

        /** \brief add a inequality with intervalVector (cst X <= rhs)
         *  using the bounding box
         *  do basic checks, but do not minimize the system
         *  \param cst the row constraint
         *  \param rhs the rhs
         *  \param tolerance tolerance for redundancy checking (CLP only)
         *  \return false if the (basic) checks showed the constraint to
         *  be redundant (or inside tolerance), true if added */
        bool add_constraint(const IntervalRow &cst, double rhs,
		double tolerance=0.0);

        /** \brief two inequalities with intervalVector (cst X in rhs)
         *  using the bounding box
         *  do basic checks, but do not minimize the system. Also,
         *  should not be used for equalities (cst and rhs punctual)
         *  \param cst the row constraint
         *  \param rhs the rhs
         *  \param tolerance tolerance for redundancy checking (CLP only)
         *  \return pair of booleans (one for each constraints possibly added */
        std::pair<bool,bool> add_constraint_band(const IntervalRow &cst,
		 const Interval &rhs, double tolerance=0.0);

        /** \brief add a equality (pair row X = rhs)
         *  do basic checks, but do not minimize the system
         *  \param facet the constraint
         *  \return -1 if empty, 0 no change (not probable), 1 change */
        int add_equality(const std::pair<Row,double>& facet);

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

        /** \brief union with a polytope
         *  \param P the polytope
         *  \return 0 if nothing done, 1 changed made, -1 results is empty */
        int join_with_polytope(const Polytope &P);

        /** \brief union with a polytope
         *  \param P the polytope
         *  \return *this */
        Polytope &operator|=(const Polytope &P);

        /** \brief inflation by a cube, keeping the shape (not optimal)
         *  this <- this + [-rad,rad]^d
         *  \param rad radius of the box
         *  \return *this (keep minimized, except rounding)
         */
        Polytope &inflate(double rad);

        /** \brief inflation with a box, keeping the shape (not optimal)
         *  this <- this + box
         *  \param box the box
         *  \return *this (keep minimized, except rounding)
         */
        Polytope &inflate(const IntervalVector& box);

        /** \brief inflation by a ball, keeping the shape (not optimal)
         *  this <- this + B_d(rad) (note: vector norm is computed on double)
         *  \param rad radius of the ball
         *  \return *this (keep minimized, except rounding)
         */
        Polytope &inflate_ball(double rad);

        /** \brief expansion of a dimension (not optimal)
         *  equivalent to inflation by a box 0,...,[-rad,rad],...0
         *  \param dm index of the dimension
         *  \param rad radius 
         *  \return *this (keep minimized, except rounding) */
        Polytope &unflat(Index dm, double rad);

        /** \brief centered homothety (optimal if c is punctual)
         *  x <- [c] + delta*(x-[c]) ( or (1-delta)[c] + delta*x )
         *  note : delta must be > 0
         *  \param c ``center''
         *  \param delta expansion
         *  \return *this
         */
        Polytope &homothety(const IntervalVector &c, double delta);

        /* creating a new polytope                  */

        /** \brief return a polytope intersected with an hyperplan
	 * in a specific coordinate (rebuild the constraints without
         * any dependence on this coordinate in a new polytope)
         *  \param dm index of the coordinate
         *  \param x the value of the coordinate
         *  \return the polytope */
        Polytope meet_with_hyperplane(Index dm, double x) const;

        /** \brief reverse affine transformation 
	 *  compute { x' in Box s.t. ([M] x' + [P]) in the initial polytope }
         *  The bounding box is needed to handle interval in [M].
         *  \param M non-empty matrix
         *  \param P non-empty vector 
         *  \param bbox bouding-box of the new polytope
         *  \return a new polytope */
        Polytope reverse_affine_transform(const IntervalMatrix &M,
		const IntervalVector &P, const IntervalVector &bbox) const;

        /** \brief linear bijective affine transformation 
	 *  compute { [M] x + [P] s.t. x in the initial polytope }
         *  with MInv encloses the inverse of M 
         *  M is used only to approximate the bounding box of the result,
         *  Minv is used to compute the new constraints
         *  \param M matrix
         *  \param Minv its inverse
         *  \param P non-empty vector 
         *  \return a new polytope */
        Polytope bijective_affine_transform(const IntervalMatrix &M,
			const IntervalMatrix &Minv, 
			const IntervalVector &P) const;

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
#ifdef WITH_CLP
            CLPFORM,     /* has an up-to-date clp form */
#endif
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
#ifdef WITH_CLP
      mutable std::unique_ptr<LPclp> _clpForm; /* LPclp formulation, if used */
#endif
      mutable std::unique_ptr<DDbuildF2V> _DDbuildF2V; 
		/* DDbuildF2V formulation, if used */
      mutable std::unique_ptr<DDbuildV2F> _DDbuildV2F; 
		/* DDbuildV2F formulation, generally not used */
      mutable pol_state state = pol_state_init;
      void minimize_constraints_F2V() const;
      void update_box_F2V() const;
      void update_box() const;
      bool check_empty_F2V() const;
      bool check_empty() const;
      double bound_row_F2V(const Row &r) const;
      void build_DDbuildF2V() const;
      void set_empty_private() const;
#ifdef WITH_CLP
      void minimize_constraints_clp(const Interval &tolerance=Interval(0.0)) const;
      void update_box_clp() const;
      bool check_empty_clp() const;
      double bound_row_clp(const Row &r) const;
      void build_clpForm() const;
#endif
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
#ifdef WITH_CLP
   _clpForm=nullptr;
#endif
   _DDbuildF2V=nullptr;
   _DDbuildV2F=nullptr;
   _facets=nullptr;
}


inline void Polytope::set_empty() {
   this->set_empty_private();
}

std::ostream& operator<<(std::ostream& os, const Polytope &P);

}
