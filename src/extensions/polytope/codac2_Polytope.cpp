/**
 * \file codac2_Polytope.cpp
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
#include "codac2_clp.h"
#include "codac2_Polytope.h"
#include "codac2_Facet.h"

using namespace codac2;


namespace codac2 {

Polytope::Polytope() : _dim(-1), _empty(true)
{
}

Polytope::Polytope(Index dim) : Polytope(dim,false)
{
}

Polytope::Polytope(Index dim, bool empty) : _dim(dim), _empty(empty),
		_box(IntervalVector::constant(_dim, 
				empty ? Interval::empty() : Interval())),
		_facets(std::make_shared<CollectFacets>(dim))
{
}

Polytope::Polytope(const IntervalVector &box) : _dim(box.size()),
		_empty(box.is_empty()), _box(box), 
		_facets(std::make_shared<CollectFacets>(_dim)) 
{
}

Polytope::Polytope(const Polytope &P) : _dim(P._dim), _empty(P._empty),
        _box(P._box),
	_facets(P._empty ? std::make_shared<CollectFacets>(_dim) :
		           std::make_shared<CollectFacets>(*(P._facets))),
	_minimized(P._minimized) 
{
}

Polytope::Polytope(const std::vector<Vector> &vertices) : 
  Polytope()
{
   if (vertices.empty()) return;
   _dim=vertices[0].size();
   _empty= false;
   /* build the V2F form */
   _DDbuildV2F = std::make_unique<DDbuildV2F>(1,vertices[0]);
   for (Index i=1;i<(Index)vertices.size();i++) {
       _DDbuildV2F->add_vertex(i+1,vertices[i]);
   }
   /* get the CollectFacets */
   _facets=_DDbuildV2F->getFacets();
   _box = _facets->extractBox();
   _facets->encompass_vertices(vertices,_box,true);
   _minimized=true; /* considered */
   _DDbuildV2F.release(); /* encompass may have added constraints */
}   

Polytope::Polytope(const std::vector<IntervalVector> &vertices,
		const CollectFacets &facetsform) :  Polytope()
{
  _dim = facetsform.getDim();
  if (_dim==-1) return;
  _facets=std::make_shared<CollectFacets>(facetsform);
  if (vertices.empty()) {
     _box = IntervalVector::constant(_dim,Interval::empty());
     _empty=true;
     return;
  }
  _facets=_DDbuildV2F->getFacets();
  _box = _facets->extractBox();
  _facets->encompass_vertices(vertices,_box,true);
  _minimized=false;
}


Polytope::Polytope(const IntervalVector &box, 
		const std::vector<std::pair<Row,double>> &facets,
		bool minimize) : Polytope(box) {
   for (const std::pair<Row,double> &cst : facets) {
      this->add_constraint(cst);
   }
   if (minimize) this->minimize_constraints_clp();
   else _minimized=false;
}



Polytope::Polytope(const Parallelepiped &par) :
	Polytope(par.box()) {
   /* we can use the vertices, or inverse the shape,
      we consider the inversion */
   IntervalMatrix u = inverse_enclosure(par.A);
   for (Index i=0;i<u.rows();i++) {
      this->add_constraint_band(u.row(i),
			u.row(i).dot(par.z)+Interval(-1.0,1.0),0.0);
   }   
   _minimized=true;
}

Polytope::Polytope(const Zonotope &zon) :
	Polytope() {
   _dim=zon.z.size();
   _empty= false;
   /* we use the vertices (note : can do better) */
   IntervalVector v = IntervalVector::constant(zon.A.cols(),-1.0);
   this->_DDbuildV2F = std::make_unique<DDbuildV2F>(1,(zon.z+zon.A*v).mid());
   Index idV = 2;
   Index a = 0;
   while (a<v.size()) {
     while (a<v.size() && v[a]==1.0) {
       v[a]=-1.0;
       a++;
     } 
     if (a<v.size()) {
       v[a]=1.0; a=0;
       this->_DDbuildV2F->add_vertex(idV++,IntervalVector(zon.z+zon.A*v));
     }
   }
   /* get the CollectFacets */
   _facets=_DDbuildV2F->getFacets();
   IntervalVector range = IntervalVector::constant(zon.A.cols(),
					Interval(-1.0,1.0));
   _box = zon.z+zon.A*range; /* FIXME : use box */
   _facets->encompass_zonotope(IntervalVector(zon.z), 
		IntervalMatrix(zon.A), range, true);
   _DDbuildV2F.release(); /* encompass may have added constraints */
   _minimized=true;
}

Polytope Polytope::from_ineFile(const char *filename) {
   std::shared_ptr<CollectFacets> fcts
	= std::make_shared<CollectFacets>(read_ineFile(filename));
   if (fcts->getDim()==-1) return Polytope();
   Polytope p(fcts->getDim(),false);
   p._facets = fcts;
   p._minimized=false;
   return p;
}

Polytope Polytope::from_extFile(const char *filename) {
   std::vector<Vector> vts = read_extFile(filename);
   return Polytope(vts);
}

Polytope::~Polytope() {
}


Interval Polytope::fast_bound(const FacetBase &base) const {
   if (base.isNull()) return Interval::zero();
   Index gdim = base.gtDim();
   Interval res;
   if (base.row[gdim]>0.0) {
      res= Interval(_box[gdim].ub())*base.row[gdim];
   } else {
      res= Interval(-_box[gdim].lb())*base.row[gdim];
   }
   if (!base.isCoord()) {
      Row bcopy = base.row;
      bcopy[gdim]=0.0;
      res += bcopy.dot(_box);
      if (_facets->get_map().empty()) return res;
      std::pair<CollectFacets::mapCIterator,CollectFacets::mapCIterator>
         pIt = _facets->get_map().equal_range(base);
      if (pIt.first!=pIt.second) { /* we found an exact key */
         return Interval(pIt.first->second.rhs);
      } else {
         if (pIt.first!=_facets->get_map().begin()) --pIt.first;
         Interval a1 = Facet_::bound_linear_form(*(pIt.first), base.row, _box);
         if (pIt.second!=_facets->endFacet()) {
             Interval a2 = Facet_::bound_linear_form(*(pIt.second), base.row, _box);
             if (a2.ub()<a1.ub()) a1=a2;
         }
         if (a1.ub()<res.ub()) res=a1;
      }
      /* check the equalities ? */
      for (Index i=0;i<_facets->nbeqfcts();i++) {
         CollectFacets::mapIterator it = _facets->get_eqFacet(i);
         Interval a = Facet_::bound_linear_form(*it, base.row, _box);
         if (a.ub()<res.ub()) res=a;
      }
   }
   /* check the vertices if _buildF2V is up-to-date */
   if (_DDbuildF2V && _buildF2V_updated) {
       const std::forward_list<DDvertex> &vertices=_DDbuildF2V->get_vertices();
       for (const DDvertex &vt : vertices) {
            IntervalVector vtB = _DDbuildF2V->compute_vertex(vt.vertex);
            Interval a = base.row.dot(vtB);
            if (a.ub()<res.ub()) res=a;
       }
   }
   return res;
}

void Polytope::build_DDbuildF2V() const {
   if (_buildF2V_updated) return;
   if (_empty) { 
        _DDbuildF2V=nullptr;
	_buildF2V_updated=true; 
        return;
   }
   /* TODO : insert box progressively ? */
   _DDbuildF2V = std::make_unique<DDbuildF2V>(_dim,_box,_facets,true);
   for (CollectFacets::mapCIterator it = _facets->get_map().begin();
		it!=_facets->get_map().end(); ++it) {
      if ((*it).second.eqcst) continue;
//      std::cout << " add facet " << *it << "\n";
      _DDbuildF2V->add_facet(it);
//      std::cout << " result " << *_DDbuildF2V << "\n";
      if (_DDbuildF2V->is_empty()) { 
         _empty=true;
         _box.set_empty();
         _DDbuildF2V=nullptr;
         break;
      }
   }
   _buildF2V_updated=true;
}

void Polytope::build_clpForm() const {
   if (_clpForm_updated) return;
   if (_empty) { 
        _clpForm=nullptr;
	_clpForm_updated=true; 
        return;
   }
   _clpForm = std::make_unique<LPclp>(_dim,_facets,_box);
   _clpForm_updated=true;
}

void Polytope::minimize_constraints_clp(const Interval &tolerance) const { 
   if (_empty) return;
   this->build_clpForm();
   int ret = _clpForm->minimize_polytope(tolerance, true);
   if (ret==-1) {
      _empty=true;
      _box.set_empty();
      _minimized=true;
      return;
   }
   _box &= _clpForm->get_bbox();
   bool changed=false;
   for (Index i=0;i<_facets->nbfcts();i++) {
      if (_clpForm->isRedundant(i)) {
	_facets->removeFacetById(i+1);
        changed=true;
      }
   } 
   if (changed) {
      _facets->renumber();
      _clpForm_updated=_buildF2V_updated=_buildV2F_updated=false;
      _clpForm=nullptr;
   }
   _minimized=true;
}


#if 0
const IntervalVector& Polytope::update_box() {
   if (_dim==-1 || _empty) return this->_box;
   if (!_clpForm) return this->_box;
   Row obj = Row::zero(_dim);
   for (Index i=0;i<_dim;i++) {
       obj[i]=1.0;
       _clpForm->setObjective(obj);
       LPclp::lp_result_stat ret =  _clpForm->solve(false,4);
       if (ret[LPclp::EMPTY]) { _empty=true; _box.set_empty(); return _box; }
       if (ret[LPclp::BOUNDED]) _box[i] = min(_box[i],_clpForm->getValobj());
       obj[i]=-1.0;
       _clpForm->setObjective(obj);
       ret =  _clpForm->solve(false,4);
       if (ret[LPclp::EMPTY]) { _empty=true; _box.set_empty(); return _box; }
       if (ret[LPclp::BOUNDED]) _box[i] = max(_box[i],-_clpForm->getValobj());
   }
   _clpForm->set_bbox(_box);
   return this->_box;
}

polytope_inclrel Polytope::relation_Box(const IntervalVector& p) const {
   assert_release(_dim>=0);
   if (p.is_empty()) return (1<<INCLUDES) | (1<<DISJOINT);
   if (_empty) return (1<<NOTINCLUDE) | (1<<DISJOINT);
   if (!_clpForm) return (1<<INCLUDES) | (1<<INTERSECTS);
   /* to get INCLUDES, must have INCLUDES for all 
      to get NOTINCLUDE, must have NOTINCLUDE once
      to get INTERSECTS, must have INCLUDES for all but one, and INTERSECTS
                         for this one 
      to get DISJOINT, must have DISJOINT once */
   bool notinclude=false;
   bool allinclude=true;
   bool intersects=true;
   const std::vector<Facet>& facets = _clpForm->getFacets();
   for (const Facet &facet : facets) {
       polytope_inclrel res2 = facet.relation_Box(p);
       if (res2[DISJOINT]) return (1<<NOTINCLUDE) | (1<<DISJOINT);
       if (res2[INCLUDES]) continue;
       if (allinclude && res2[INTERSECTS]) { allinclude=false; }
       else { allinclude=false; intersects=false; }
       if (res2[NOTINCLUDE]) notinclude=true;
   }
   return (allinclude ? (1<<INCLUDES | 1<<INTERSECTS)
	: ((intersects ? 1<<INTERSECTS : 1<<MAYINTERSECT) |
           (notinclude ? 1<<NOTINCLUDE : 1<<MAYINCLUDE)));
}

BoolInterval Polytope::contains(const IntervalVector& p) const {
   polytope_inclrel r = this->relation_Box(p);
   return (r[INCLUDES] ? BoolInterval::TRUE : 
		(r[NOTINCLUDE] ? BoolInterval::FALSE : BoolInterval::UNKNOWN));
}

BoolInterval Polytope::intersects(const IntervalVector& p) const {
   polytope_inclrel r = this->relation_Box(p);
   return (r[INTERSECTS] ? BoolInterval::TRUE : 
		(r[DISJOINT] ? BoolInterval::FALSE : BoolInterval::UNKNOWN));
}

const Interval& Polytope::operator[](Index i) const {
   assert_release(i<=_dim);
   return _box[i];
}

Vector Polytope::mid() const {
   assert_release(_dim>=0);
   if (_empty) return 
	Vector::Constant(_dim,std::numeric_limits<double>::quiet_NaN());
   if (!_clpForm) return _box.mid();
   LPclp::lp_result_stat stat = _clpForm->solve(true,0);
   if (stat[LPclp::NOTEMPTY] || stat[LPclp::NOTEMPTY_APPROX])
       return _clpForm->getFeasiblePoint().mid();
   return Vector::Constant(_dim,std::numeric_limits<double>::quiet_NaN());
}

double Polytope::distance_cst(const Facet &fc) const {
   assert_release(_dim==fc.row.size());
   if (_empty) return -oo;
   if (!_clpForm) return (fc.row.dot(_box)-fc.rhs).ub();
   _clpForm->setObjective(fc.row);
   LPclp::lp_result_stat stat = _clpForm->solve(false,0);
   if (stat[LPclp::EMPTY]) return -oo;
   return (_clpForm->getValobj()-fc.rhs).ub();
}

bool Polytope::box_is_included(const IntervalVector& x) const {
   return _box.is_subset(x);
}


void Polytope::clear() {
   assert_release(_dim>=0);
   _empty=false;
   _box = IntervalVector::Zero(_dim);
   this->build_clpForm_from_box();
}

#endif

bool Polytope::add_constraint(const std::pair<Row,double>& facet, 
	double tolerance) {
   if (_empty) return false;
   FacetBase base = FacetBase(facet.first);
   if (base.isNull()) {
      if (facet.second>=tolerance) return false;
      this->set_empty(); return true;
   }
   Interval act = this->fast_bound(base);
   if (act.ub()<=facet.second+tolerance) return false;
   if (base.isCoord()) {
      Index gdim = base.gtDim();
      Interval val = Interval(facet.second)/facet.first[gdim];
      if (facet.first[gdim]>0.0) {
         _box[gdim] &= Interval(-oo,val.ub());
         if (_box[gdim].is_empty()) { 
	   this->set_empty(); 
           return true;
	 }
      } else {
         _box[gdim] &= Interval(-val.ub(),+oo);
         if (_box[gdim].is_empty()) { 
	   this->set_empty(); 
           return true;
	 }
      }
      if (_clpForm_updated) {
         this->_clpForm->set_bbox(_box);
      }
      _buildF2V_updated = _buildV2F_updated = false;
			/* TODO : update _DDbuildF2V */
      _minimized=false;
      return true; 
   }
   auto res= _facets->insert_facet(facet.first, facet.second, 
		false, CollectFacets::MIN_RHS);
   if (res.first==0) return false;
   if (_clpForm_updated) {
      this->_clpForm->updateConstraint(res.first);
   }
   if (_buildF2V_updated) {
      this->_DDbuildF2V->add_facet((*_facets)[res.first-1]);
   }
   _buildV2F_updated = false;
   _minimized=false;
   return true;
}

bool Polytope::add_constraint(const IntervalRow &cst, double rhs,
	double tolerance) {
   if (_empty) return false;
   Row cstmid = cst.mid();
   IntervalRow rem = cst-cstmid;
   Interval d = rhs+rem.dot(_box);
   if (d.ub()==+oo) return false;
   return  this->add_constraint(std::pair(cstmid, d.ub()), tolerance);
}

std::pair<bool,bool> Polytope::add_constraint_band(const IntervalRow &cst,
	 const Interval &rhs, double tolerance) {
   if (_empty) return { false, false };
   Row cstmid = cst.mid();
   IntervalRow rem = cst-cstmid;
   Interval d = rhs+rem.dot(_box);
   bool rub, rlb;
   if (d.ub()==+oo) rub=false;
   else rub = this->add_constraint(std::pair(cstmid, d.ub()), tolerance);
   if (d.lb()==-oo) rlb=false;
   else rlb = this->add_constraint(std::pair(-cstmid, -d.lb()), tolerance);
   return { rlb, rub };
}
     
void Polytope::set_empty() const {
   _empty=true;
   _box.set_empty();
}

void Polytope::set_empty() {
   _empty=true;
   _box.set_empty();
}

std::vector<IntervalVector> Polytope::compute_vertices() const {
   this->build_DDbuildF2V();
   if (_empty) return std::vector<IntervalVector>();
   std::vector<IntervalVector> ret;
   for (const DDvertex &vtx : _DDbuildF2V->get_vertices()) {
       ret.push_back(_DDbuildF2V->compute_vertex(vtx.vertex));
   }
   return ret;
}

std::vector<std::vector<Vector>> Polytope::compute_3Dfacets() const {
   this->build_DDbuildF2V();
   if (_empty) return std::vector<std::vector<Vector>>();
   return build_3Dfacets(*_DDbuildF2V);
} 

}

