/**
 * \file codac2_clp_Polytope.cpp (Polytope methods which use CLP)
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Damien Massé
 *  \copyright  Copyright 2025
 *  \license    This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 */

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
#include "codac2_inversion.h"
#include "codac2_Parallelepiped.h"
#include "codac2_Zonotope.h"
#include "codac2_Polytope.h"
#include "codac2_Polytope_clp.h"
#include "codac2_Facet.h"
#include "codac2_clp.h"

using namespace codac2;


namespace codac2 {

Polytope_clp &Polytope_clp::operator=(const Polytope_clp &P) {
   assert(!P.state[INVALID]);
   this->_dim = P._dim;
   this->_box = P._box;
   this->_facets = P.state[EMPTY] ? std::make_shared<CollectFacets>(_dim) :
                           std::make_shared<CollectFacets>(*(P._facets));
   this->_clpForm=nullptr; /* FIXME */
   this->_DDbuildF2V=nullptr;
   this->_DDbuildV2F=nullptr;
   this->state = P.state & (pol_state_empty | 
                pol_state(1<<NOTEMPTY | 1<<INVALID));
   return *this;
}

Polytope_clp &Polytope_clp::operator=(const Polytope &P) {
   *this = std::move(Polytope_clp(P));
   return *this;
}
Polytope_clp &Polytope_clp::operator=(Polytope &&P) {
   *this = std::move(Polytope_clp(P));
   return *this;
}

double Polytope_clp::bound_row(const Row &r) const {
   assert(r.size()==_dim);
   if (state[EMPTY]) return -oo;
//   if (!state[F2VFORM]) { 
        return this->bound_row_clp(r); 
//   }
//   return Polytope::bound_row_F2V(r); 
}

double Polytope_clp::bound_row_clp(const Row &r) const {
   this->build_clpForm();
   _clpForm->setObjective(r);
   LPclp::lp_result_stat ret = _clpForm->solve();
   if (ret[EMPTY]) return -oo;
   return _clpForm->getValobj().ub();
}

void Polytope_clp::build_clpForm() const {
   if (clpUptodate) return;
   if (state[EMPTY]) { 
        _clpForm=nullptr;
        return;
   }
   _clpForm = std::make_unique<LPclp>(_dim,_facets,_box);
   clpUptodate=true;
}

void Polytope_clp::minimize_constraints_clp(const Interval &tolerance) const { 
   assert_release(!state[MINIMIZED]);
   this->build_clpForm();
   int ret = _clpForm->minimize_polytope(tolerance, false, !state[BOXUPDATED]);
   if (ret==-1) {
      this->set_empty_private();
      return;
   }
   state[NOTEMPTY]=true;
   if (!state[BOXUPDATED]) {
     _box &= _clpForm->get_bbox();
     state[BOXUPDATED]=true;
   }
   bool changed=false;
   for (Index i=0;i<_facets->nbfcts();i++) {
      if (_clpForm->isRedundant(i)) {
	_facets->removeFacetById(i+1);
        changed=true;
      }
   } 
   if (changed) {
      std::vector<Index> corresp = _facets->renumber();
/*      if (state[F2VFORM]) {
         _DDbuildF2V->update_renumber(corresp);
      } */
      state[F2VFORM]=false;
      clpUptodate=false; state[V2FFORM]=false;
      _clpForm=nullptr;
   }
   state[MINIMIZED]=true;
}

void Polytope_clp::minimize_constraints() const {
   if (state[MINIMIZED]) return;
   if (state[EMPTY]) return;
   /* default : F2V */
//   if (!state[F2VFORM] && state[CLPFORM]) { 
        this->minimize_constraints_clp(); 
//        return;
//   }
//   this->minimize_constraints_F2V(); 
}


void Polytope_clp::update_box_clp() const {
   assert_release(!state[BOXUPDATED]);
   this->build_clpForm();
   int ret = _clpForm->minimize_box();
   if (ret==-1) {
      this->set_empty_private();
      return;
   }
   _box &= _clpForm->get_bbox();
   state[NOTEMPTY]=true;
   state[BOXUPDATED]=true;
}

void Polytope_clp::update_box() const {
   if (state[BOXUPDATED]) return;
   if (state[EMPTY]) return;
//   if (!state[F2VFORM] && state[CLPFORM]) { 
        this->update_box_clp(); 
//        return;
//   }
//   this->update_box_F2V(); 
}


bool Polytope_clp::check_empty_clp() const {
   assert_release(!state[NOTEMPTY] && !state[EMPTY]);
   this->build_clpForm();
   int ret = _clpForm->check_emptiness();
   if (ret==-1) {
      this->set_empty_private();
      return true;
   }
   state[NOTEMPTY]=true;
   return false;
}

bool Polytope_clp::check_empty() const {
   if (state[NOTEMPTY]) return false;
   if (state[EMPTY]) return true;
//   if (!state[F2VFORM] && state[CLPFORM]) { 
        return this->check_empty_clp(); 
//   }
//   return this->check_empty_F2V(); 
}



}

