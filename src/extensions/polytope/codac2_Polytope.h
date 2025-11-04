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

#if 0

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
         * \brief Constructs a polytope from a set of vertices
         * and a set of linear forms, i.e. the minimal polytope 
         * built on the linear forms that includes the vertices.
         * \param vertices the vertices
         * \param facetforms the facets description (must NOT be empty)
         */
        Polytope(const std::vector<Vector> &vertices, 
		const std::vector<Row> &facetforms);

        /**
         * \brief from a slanted box (M.V with reverse known)
	 * \param M matrix
         * \param rM inverse matrix
         * \param V vector
         */
        Polytope(const IntervalMatrix& M, const IntervalMatrix &rM,
		const IntervalVector &V);

        /**
         * \brief from a slanted box (M.V, reverse computed)
         * \param M matrix
         * \param V vector
         */
        Polytope(const IntervalMatrix& M, const IntervalVector &V);

        /**
         * \brief Constructs a polytope from a bounding box and 
         * a set of constraints 
         */
        Polytope(const IntervalVector &box, 
                 const std::vector<Facet> &facets, 
	         bool minimize=false);

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
        bool is_empty() const;

        /**
         * \brief has a flat dimension
         */
        bool is_flat() const;

        /** \brief ``component'' : interval of the bounding box 
           * (only const !!!) */
        const Interval& operator[](Index i) const;

        /**
         * ``middle'' (vector inside the polytope, if possible)
         */
        Vector mid() const;

        /** \brief ``distance'' from the satisfaction of a constraint 
         *  ie the inflation of the rhs needed to satisfy the constraint
         *  \param fc the constraint
         *  \return the distance
         */
        double distance_cst (const Facet &fc) const;


        /**
         * \brief relationships with a box (fast check)
         * \param p the box
         * \return polytope_inclrel checking inclusion and intersection */
        polytope_inclrel relation_Box(const IntervalVector& p) const;

        /**
         * \brief Checks whether the polytope contains a given point, or
         *  includes a box
         *
         * \param p The point or box to check, enclosed in an IntervalVector.
         * \return BoolInterval indicating possible containment.
         */
        BoolInterval contains(const IntervalVector& p) const;

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

        /************* Box access *************************/
 
        /**
         * \brief get the current bounding box of the polytope
         * (note: it may not be tightest bounding box, it is the 
         *  bounding box used for approximations results) 
         *
         * \return The current bounding box
         */
        const IntervalVector &box() const;

        /**
         * \brief compute and update the bounding box
         * (using LPs on the polytope) 
         *
         * \return the updated hull box
         */
        const IntervalVector &update_box();

        /**
         * \brief test if bounding box is included
         * (mainly useful if current bounding box is tight)
         */
        bool box_is_included(const IntervalVector& x) const;


        /************* Modification **********************/
        /** modification creates a new polytope         **/

        /** \brief set to empty */
        void set_empty();

        /** \brief set to (singleton) 0 */
        void clear();

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


        /**
         * \brief Computes the sets of vertices of the polytope
         * 
         * WARNING: the results is unsafe: no guarantee is given
         * that the convex hull of the result fully encloses the polytope
         * works only for bounded polyhedron, suppose _box exists and is bounded
         *
         * \return set of vertices, as Vectors
         */
        std::vector<Vector> vertices() const;

   private:
      Index _dim;                      /* dimension */
      Index _nbEqcsts;                 /* nb of equality constraints */
      Index _nbineqCsts;               /* nb of inequality constraints */
      bool  _empty; 		       /* is empty */
      IntervalVector _box;             /* bounding box */
      mutable std::unique_ptr<LPclp> _clpForm; /* LPclp formulation */
      mutable std::unique_ptr<DDbuildF2V> _DDbuildF2V; 
					/* DDbuildF2V formulation */
      mutable std::unique_ptr<DDbuildV2F> _DDbuildV2F; 
					/* DDbuildV2F formulation */
      bool _box_updated;
      bool _clpFrom_updated;
      bool _buildF2V_updated;
      bool _buildV2F_updated;
      void build_clpForm_from_box();
      void minimize_constraints();
};

inline const IntervalVector &Polytope::box() const { return _box; }
inline Index Polytope::dim() const { return _dim; }
inline Index Polytope::size() const { return _dim; }
inline bool Polytope::is_empty() const { return _empty; }
inline bool Polytope::is_flat() const { return _nbEqcsts!=0; }
inline Index Polytope::nbFacets() const { return _nbEqcsts!=0; }



}

#endif
