/**
 * \file codac2_Polytope_clp.h
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
#include "codac2_Polytope.h"
#include "codac2_Parallelepiped.h"
#include "codac2_Zonotope.h"

namespace codac2 {
  /**
   * \class Polytope_clp extends polytope with a clp Form
   *  Represents a bounded convex polytope as a set of constraints
   * and a bounding box (the bounding box is part of the constraints)  */
  class Polytope_clp : public Polytope {

   public:
      Polytope_clp() : Polytope();
      explicit Polytope_clp(Index dim) : Polytope(dim) {};
      explicit Polytope_clp(Index dim, bool empty) : Polytope(dim,empty) {};
      explicit Polytope_clp(const IntervalVector &box) : Polytope(box) {};
      Polytope_clp(const Polytope_clp &P) : Polytope(P) {}; /* FIXME */
      Polytope_clp(Polytope_clp &&P) = default;
      explicit Polytope_clp(const std::vector<Vector> &vertices) : Polytope(vertices) {};
      explicit Polytope_clp(const std::vector<IntervalVector> &vertices) : Polytope(vertices) {};
      Polytope_clp(const std::vector<IntervalVector> &vertices,
                const CollectFacets &facetsform) : Polytope(vertices, facetsform) {};
      explicit Polytope_clp(const Parallelepiped &par) : Polytope(par) {};
      explicit Polytope_clp(const Zonotope &zon) : Polytope(zon) {};
      Polytope_clp(const IntervalVector &box,
                 const std::vector<std::pair<Row,double>> &facets,
                 bool minimize=false) : Polytope(box,facets,minimize) {};

      /**
       *  copy assignment operator 
       * \param P copy
       */
      Polytope_clp &operator=(const Polytope_clp &P);


      void minimize_constraints() const override;
     
   private:
      mutable std::unique_ptr<LPclp> _clpForm = nullptr; /* LPclp formulation */

      mutable bool clpUptodate = nullptr;

      void update_box() const override;
      bool check_empty() const override;
      void set_empty_private() const override;

      void minimize_constraints_clp(const Interval &tolerance=Interval(0.0))
			 const;
      void update_box_clp() const;
      bool check_empty_clp() const;
      double bound_row_clp(const Row &r) const;
      void build_clpForm() const;

      void update_box_F2V() const;
      bool check_empty_F2V() const;
      double bound_row_F2V(const Row &r) const;
      void build_DDbuildF2V() const;
  }

