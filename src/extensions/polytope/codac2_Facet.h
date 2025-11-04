/**
 * \file codac2_Facet.h  
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Damien Massé
 *  \copyright  Copyright 2025
 *  \license    This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 */

#pragma once

#include <utility> /* uses swap */
#include <bitset>
#include <cmath>
#include <set>

#include "codac2_Index.h"
#include "codac2_Vector.h"
#include "codac2_Row.h"
#include "codac2_Interval.h"
#include "codac2_IntervalVector.h"
#include "codac2_IntervalRow.h"
#include "codac2_BoolInterval.h"
#include "codac2_operators.h"

namespace codac2 {

/* base of comparison for Facet (note: only a preorder) */
struct FacetBase {
   Index bdim; /* if A is the index of the greatest abs value and B
                  the second greatest (=A if 0 outside)
		    v[A]>0 v[B]>0     A*(2*dim)+A-B (+2*dim if B>A)
                    v[A]>0 v[B]<0     A*(2*dim)+A-(B+dim)+2*dim
                    v[A]<0 v[B]>0     (A+dim)*(2*dim)+A+dim-B
                    v[A]<0 v[B]<0     (A+dim)*(2*dim)+(A+dim)-(B+dim)
							 (+2*dim if B>A)
                  noting Adim = A if v[A]>0 and A+dim if v[A]<0
                  and Bdim = B if v[B]>0 and B+dim if v[B]<0
                  bdim = Adim*(2*dim+1) - Bdim (+2*dim if Adim>Bdim)
                  if null vector, value = -1 */
   double vdim; /* |v[B]|/|v[A]| (0 if no value for B) */
   Row row;

   FacetBase(const Row &row) : bdim(-1), vdim(0.0), row(row) {
      this->compute_key();
   }

   FacetBase(Row &&row) : bdim(-1), vdim(0.0), row(row) {
      this->compute_key();
   }
   FacetBase(const FacetBase &fc) : bdim(fc.bdim), vdim(fc.vdim), row(fc.row) {
   }
   FacetBase(FacetBase &&fc) : bdim(fc.bdim), vdim(fc.vdim), row(fc.row) {
   }

   inline void swap(FacetBase &fb) {
      this->row.swap(fb.row);
      std::swap(this->bdim,fb.bdim);
      std::swap(this->bdim,fb.bdim);
   }
   
   void change_row(Row &&row) {
      this->row=row; this->bdim=-1; this->vdim=0.0;
      this->compute_key();
   }
   void change_row(const Row &row) {
      this->row=row; this->bdim=-1; this->vdim=0.0;
      this->compute_key();
   }

   void negate_row() {
      this->row=-row;
      Index u = 2*row.size()*row.size();
      if (this->bdim>=u) this->bdim-=u; else this->bdim+=u;
   }

   private:

   void compute_key() {
      if (row.size()==0) return; 
      double val=row[0];
      double valabs=std::abs(val);
      if (val<0.0) bdim=row.size(); else if (val>0.0) bdim=0;
      Index b2dim = -1;
      double val2abs=0.0;
      for (Index i=1;i<row.size();++i) {
          double val_i=std::abs(row[i]);
          if (val_i>valabs) {
	     val2abs=valabs;
	     valabs=val_i;
             val=row[i];
 	     b2dim=bdim;
             if (val<0.0) bdim=row.size()+i; else bdim=i;
          } else if (val_i>val2abs) {
             val2abs=val_i;
             if (val_i<0.0) b2dim=row.size()+i; else b2dim=i;
          }
      }
      if (valabs==0.0) return; /* bdim=-1 */
      if (b2dim==-1) b2dim=bdim;
      if (b2dim>bdim) 
          bdim = bdim * (2*row.size()+1) + 2*row.size() - b2dim;
      else
          bdim = bdim * (2*row.size()+1) - b2dim;
      vdim = valabs/val2abs;
   }
};

struct FacetRhs {
   double rhs;
   bool eqcst;
   Index Id;

   FacetRhs(double rhs, bool eqcst, Index Id) : rhs(rhs), eqcst(eqcst), Id(Id)
	 {};
   FacetRhs(const FacetRhs rhs, Index Id) : rhs(rhs.rhs), eqcst(rhs.eqcst),
	 Id(Id)
	 {};

   void swap(FacetRhs &r) {
     std::swap(rhs,r.rhs);
     std::swap(eqcst,r.eqcst);
   }
};

using Facet = std::pair<FacetBase,FacetRhs>;

/* I want specific functions for Facet, but I do not want to create
   a new class inheriting Facet (as the Map will have the pair)
   => using namespace */
namespace Facet_ {

  inline Facet make(const Row &row, double rhs, bool eqcst, Index Id) {
     return std::make_pair(FacetBase(row),FacetRhs(rhs,eqcst,Id));
  }
  inline Facet make(Row &&row, double rhs, bool eqcst, Index Id) {
     return std::make_pair(FacetBase(row),FacetRhs(rhs,eqcst,Id));
  }
  inline Facet make(const Row &row, double rhs, bool eqcst) {
     return std::make_pair(FacetBase(row),FacetRhs(rhs,eqcst,0));
  }
  inline Facet make(Row &&row, double rhs, bool eqcst) {
     return std::make_pair(FacetBase(row),FacetRhs(rhs,eqcst,0));
  }


  inline void swap_Facet(Facet &f1,Facet &f2) {
     f1.first.swap(f2.first);
     f1.second.swap(f2.second);
  }

   
  /** \brief Potential inclusion relation between sets 
   * INCLUDES : A is included in B
   * MAYINCLUDE : A may be included in B
   * NOTINCLUDE : A\\B is not empty
   * INTERSECTS : A inter B is non empty
   * MAYINTERSECT : A inter B may be non empty
   * DISJOINT : A inter B is empty
   */
  enum INCLREL {
    INCLUDES,
    MAYINCLUDE,
    NOTINCLUDE,
    INTERSECTS,
    MAYINTERSECT,
    DISJOINT,
    SIZEINCLREL
  };
  using polytope_inclrel = std::bitset<SIZEINCLREL>;
  const polytope_inclrel inclrel_includes = 1<<INCLUDES;
  const polytope_inclrel inclrel_mayinclude = 1<<MAYINCLUDE;
  const polytope_inclrel inclrel_notinclude = 1<<NOTINCLUDE;
  const polytope_inclrel inclrel_intersects = 1<<INTERSECTS;
  const polytope_inclrel inclrel_mayintersect = 1<<MAYINTERSECT;
  const polytope_inclrel inclrel_disjoint = 1<<DISJOINT;

  /* test only includes and intersects */
  polytope_inclrel relation_Box(const Facet &f, const IntervalVector &b);

   /* give the bound for a Row and a box, based on the largest index (bdim).
      ie compute R(bdim)/r(bdim) rhs + (R - R(bdim)/r(bdim) r) B,
      assuming R(bdim)/r(bdim) is positive, or the constraint is an equality */
  Interval bound_linear_form(const Facet &f, 
		const Row &row2, const IntervalVector &b);

  void contract_Box(const Facet &f, IntervalVector &b);
  void contract_out_Box(const Facet &f, IntervalVector &b);
}

std::ostream& operator<<(std::ostream& os, const Facet& f); 

struct FacetBaseComp {
   inline bool operator()(const FacetBase& lhs, const FacetBase& rhs) const
    {
        if (lhs.bdim<rhs.bdim) return true;
        else if (lhs.bdim>rhs.bdim) return false;
        if (lhs.vdim<rhs.vdim) return true;
        else if (lhs.vdim>rhs.vdim) return false;
        for (int i=0;i<lhs.row.size();i++) 
           if (lhs.row[i]<rhs.row[i]) return true;
           else if (lhs.row[i]>rhs.row[i]) return false;
        return false;
    }
};

/** encapsulate a collection of facets with double index :
    by the base vector, and by index, with storing of eqcsts */
class CollectFacets {
   public:
     using mapType = std::map<FacetBase,FacetRhs,FacetBaseComp>;
     using mapIterator = mapType::iterator;
     using mapCIterator = mapType::const_iterator;

     /** generate an empty collection */
     CollectFacets() {};
     /** generate the set of facets from a matrix and a vector, i.e.
      *  mat X <= rhs (eqSet states which rows of mat are equalities)
      *  @param mat the matrix
      *  @rhsvect the vector of rhs
      *  @eqSet the set of equalities */
     CollectFacets(const Matrix &mat, const Vector &rhsvect, 
			const std::vector<Index> &eqSet);
  
     /** return the number of facets, as the size of _allfacets */
     Index nbfcts() const;
     /** return the number of equality facets, as the size of _eqfacets */
     Index nbeqfcts() const;
     /** return the map */
     const mapType &get_map() const;

     enum ACTION_DUPLICATE { /* operation to do when a modification
				create a duplicate facet */
       KEEP_RHS, /* change nothing (default behaviour) */
       CHANGE_RHS, /* change the value in all cases */
       MAX_RHS, /* maximize the RHS */
       MIN_RHS  /* minimize the RHS */
     };
       
   
     /** insert a new facet row X <= rhs
      * @param row Row
      * @param rhs right-hand side 
      * @param eqcst equality constraint
      * @param act action if an rhs with the same row exists
      * @return the Id of the constraint, and true if a new facet
      * has been inserted, false if old facet */
     std::pair<Index,bool> insert_facet(const Row &row, 
			double rhs, bool eqcst, ACTION_DUPLICATE act=KEEP_RHS);
     /** insert a new facet row X <= rhs
      * @param row Row
      * @param rhs right-hand side 
      * @param eqcst equality constraint
      * @param act action if an rhs with the same row exists
      * @return the Id of the new constraint, and true if a new facet
      * has been inserted, false if old facet. If old facet and no change,
      * returns 0 */
     std::pair<Index,bool> insert_facet(Row &&row, double rhs, bool eqcst,
			ACTION_DUPLICATE act=KEEP_RHS);
     /** insert a new facet row X <= rhs, using the Facet constuction
      * @param facet the facet
      * @param act action if an rhs with the same row exists
      * @return the Id of the constraint, and true if a new facet
      * has been inserted, false if old facet. If old facet and no change
      * returns 0 */
     std::pair<Index,bool> insert_facet(const Facet &facet,
		ACTION_DUPLICATE act=KEEP_RHS);

     /** get a facet by Id, starting from 0 (i.e. facet of index Id+1)
      * @param id the Id of the facet, minus 1
      * @return an iterator on the facet */
     mapIterator operator[](Index id);
     /** get a equality facet by the equality Id,
      * starting from 0 (i.e. facet of index eqIndex(Id+1))
      * @param eqId the equality Id of the facet
      * @return an iterator on the facet */
     mapIterator get_eqFacet(Index eqId);
     /** get the end iterator of the facet
      * @return an iterator on the facet */
     mapIterator endFacet(); 

     /** change a eqfacet, keeping its Id (unless the new row creates
      * a duplicate)
      *  @param eqId (the eqId of the facet)
      *  @param nrow (the new row of the facet)
      *  @param nrhs (the new right-hand side)
      *  @param act action if duplicate
      *  @return an iterator on the new facet, true if the modified facet
      *  has been inserted, false if there was a duplicate. In this case
      *  if the existing facet is not changed, return endFacet() */
     std::pair<mapIterator,bool> change_eqFacet(Index eqId, Row&& nrow,
			 double nrhs, ACTION_DUPLICATE act=KEEP_RHS);
     /** change a facet, keeping its Id (unless the new row creates
      * a duplicate)
      *  @param Id (the Id of the facet starting from 1)
      *  @param nrow (the new row of the facet)
      *  @param nrhs (the new right-hand side)
      *  @param act action if duplicate
      *  @return an iterator on the new facet, true if the modified facet
      *  has been inserted, false if there was a duplicate. In this
      *  cas, without change, return endFacet() */
     std::pair<mapIterator,bool> change_ineqFacet(Index id, Row&& nrow,
			 double nrhs, ACTION_DUPLICATE act=KEEP_RHS);

     /** change the rhs of a facet
      *  @param it the iterator of the facet
      *  @param nrhs the new right-hand side
      *  @return the iterator of the facet (=it) */
     mapIterator change_ineqFacet_rhs(mapIterator it, double nrhs);

     /** transform an equality facet to an unequality facet, and add
        a new facet with a new bound nbound
        if nbound<rhs, keep row x <= rhs and add -row x <= -nbound
        if nbound>rhs, put row x <= nbound and add -row x <= -rhs
        if nbound = -oo or +oo, do not modify/add the new facet
        @param eqId the equality Id of the equality facet
        @param nbound the new bound
        @return depends on the construction
            nbound=-oo => current facet iterator
            -oo<nbound<rhs => new facet iterator (for nbound), or endFacet if
                insertion "failed"
            rhs<nbound<+oo => new facet iterator (for rhs), or endFacet if 
                insertion "failed"
            nbound=+oo => new facet iterator (for rhs), or endFacet if 
	        insertion "failed" */
     mapIterator dissociate_eqFacet(Index eqId , double nbound,
		ACTION_DUPLICATE act=KEEP_RHS);

   private:
      mapType _map;
      std::vector<mapIterator> _allFacets;
      std::vector<Index> _eqFacets;   /* index in _allFacets */

      /* look for and remove a value in eqFacet */
      void remove_in_eqFacets(Index id);
      /* finish an insertion process */
      std::pair<Index,bool> result_insertion(Index id,
		std::pair<mapIterator,bool> res,
		 double rhs, bool eqcst, ACTION_DUPLICATE act);
};

inline const CollectFacets::mapType &CollectFacets::get_map() const
   { return _map; }

inline CollectFacets::mapIterator CollectFacets::operator[](Index id)  
   { return this->_allFacets[id]; }
inline CollectFacets::mapIterator CollectFacets::get_eqFacet(Index id) 
   { return this->_allFacets[this->_eqFacets[id]]; }

inline Index CollectFacets::nbfcts() const {
   return this->_allFacets.size(); 
}

inline Index CollectFacets::nbeqfcts() const {
   return this->_eqFacets.size(); 
}

inline CollectFacets::mapIterator CollectFacets::endFacet() { 
   return this->_map.end();
}

}
