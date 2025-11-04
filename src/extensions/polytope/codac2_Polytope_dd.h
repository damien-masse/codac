/**
 * \file codac2_Polytope_dd.h implementation of the dd algorithm
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
#include <map>
#include <forward_list>
#include <set>
#include <unordered_set>
#include <stdexcept>

#include "codac2_Index.h"
#include "codac2_Vector.h"
#include "codac2_Row.h"
#include "codac2_Facet.h"


namespace codac2 {


#if 0
struct DDlink; /* link between two vertices */
#endif

struct DDvertex {
   Index Id;
   IntervalVector vertex;
   Interval lambda;
   std::vector<Index> fcts; 
	/* note : double redirection here, a ``real'' facet may have several
           Ids */
   
   enum ddvertexstatus {
     VERTEXNEW, /* created */
     VERTEXLT, /* < (GT,LT) => new vertex */
     VERTEXON, /* <,=,> (GT,LE) => for now, equivalent to LT */
     VERTEXGE, /* = and possible small >
		  keep the vertex. no change with (GT,GE) */
     VERTEXGT, /* vertex GT, untreated and not in stack */
     VERTEXSTACK, /* vertex in stack for later treatment */
     VERTEXREM /* vertex GT treated, to be removed */
   };
   ddvertexstatus status;

   std::vector<std::forward_list<DDvertex>::iterator> links;
 
   DDvertex(const IntervalVector &v, Index Id) : Id(Id), 
	 vertex(v), status(VERTEXNEW) {
	reduce_vector();
   }

   DDvertex(IntervalVector &&v, Index Id, 
		std::vector<std::forward_list<DDvertex>::iterator> &nlinks) :
	 Id(Id), 
	 vertex(v), status(VERTEXNEW) {
         this->links.swap(nlinks);
	 reduce_vector();
   }

   bool addLnk(const std::forward_list<DDvertex>::iterator &a,
		 bool check=false) {
      if (check) {
         for (auto b : links) {
            if (a==b) return false;
         } 
      }
      this->links.push_back(a);
      return true;
   }

//   void addLnkVertex(const std::forward_list<DDvertex>::iterator &V1);

   void removeLnk(const std::forward_list<DDvertex>::iterator &a) {
      for (Index i=0;i<(Index)links.size();i++) {
         if (links[i]==a) { 
	    links[i]=links[links.size()-1]; 
	    links.pop_back();
            return;
         }
      }
   }

   void changeLnk(const std::forward_list<DDvertex>::iterator &a,
		const std::forward_list<DDvertex>::iterator &b) {
      for (Index i=0;i<(Index)links.size();i++) {
         if (links[i]==a) { links[i]=b; return; }
      }
      /* should not happen */
      links.push_back(b);
   }

   bool addFct(Index a) {
      auto rev_it = this->fcts.rbegin();
      while (rev_it != this->fcts.rend() && (*rev_it)>a) rev_it++;
      if (rev_it!=this->fcts.rend() && (*rev_it)==a) return false;
      this->fcts.insert(rev_it.base(),a);
      return true;
   }
 
   inline void reduce_vector() {
      int fx;
      frexp(vertex[0].mid(),&fx);
      for (Index i=1;i<vertex.size();i++) {
         int fx2;
         frexp(vertex[i].mid(),&fx2);
         if (fx2>fx) fx=fx2;
      }
      if (std::abs(fx)<10) return;
      double mult=exp2(-fx);
      for (Index i=0;i<vertex.size();i++) {
         vertex[i] *= mult;
      }
   }
};

#if 0
struct DDlink {
   std::forward_list<DDvertex>::iterator V1;
   std::forward_list<DDvertex>::iterator V2;
//   std::unordered_set<std::vector<Index>> fcts;
  
   DDlink(const std::forward_list<DDvertex>::iterator &V1,
	  const std::forward_list<DDvertex>::iterator &V2
		/*, const std::vector<Index> &fcts*/) : V1(V1), V2(V2) {
/*       this->fcts.insert(fcts); */
   }

   void changeVertex(const std::forward_list<DDvertex>::iterator &Vold,
                const std::forward_list<DDvertex>::iterator &Vnew) {
       if (this->V1==Vold) { this->V1=Vnew; return; }
       this->V2=Vnew;
   }
   
   const std::forward_list<DDvertex>::iterator &otherVertex
		(const std::forward_list<DDvertex>::iterator &V) const {
       if (V==V1) return V2;
       return V1;
   }
  
#if 0
   void addFct(const std::vector<Index> &fcts) {
      this->fcts.insert(fcts);
   }
#endif
};
#endif

class DDbuildF2V {
  public:
     /* create the build with a collection of facet, but
        use only the equality facets */
     DDbuildF2V(Index dim, const IntervalVector &box,
 		std::shared_ptr<CollectFacets> facetsptr,
		bool include_box=true);
     
//     void add_constraint_box(const IntervalVector &box);
   
     /* add an inequality facet in the build */
     int add_facet(Index id);
     int add_facet(CollectFacets::mapCIterator itFacet);

     IntervalVector compute_vertex(const IntervalVector &vect) const;
     Index get_dim() const;
     Index get_fdim() const;

     bool is_empty() const;
     const IntervalMatrix &get_M_EQ() const;

     const std::vector<IntervalVector> &get_lines() const;
     const std::forward_list<DDvertex> &get_vertices() const;
     const std::vector<Index> &get_reffacets() const;

     friend std::ostream& operator<<(std::ostream& os, const DDbuildF2V& build);

  private:
     bool addDDlink(const std::forward_list<DDvertex>::iterator& V1,
   		      const std::forward_list<DDvertex>::iterator& V2);

     std::forward_list<DDvertex>::iterator addDDvertex(const IntervalVector &v);
     std::forward_list<DDvertex>::iterator addDDvertex(IntervalVector &&v);
/*     std::map<Index,DDvertex>::iterator addVertexSon(const DDvertex &p1,
		 DDvertex &p2, double rhs, Index idFacet);  */
/*     void removeVertex(Index Id, bool removeLnks); */

#if 0
     /* x_dim = eqcst x + rhs */
     struct EqFacet {
        Index dim;
        Row eqcst;
        Interval valcst;
        
        EqFacet (Index dim, const Row &eqcst, const Interval &valcst) :
		dim(dim), eqcst(eqcst), valcst(valcst) { };

        /* computation of the value */
        void computeVal(IntervalVector &x) const {
            x[dim] = this->eqcst.dot(x) + valcst;
        }

        /* adapt box */
        void adapt_box(IntervalVector &bbox) const;
       
        /* adaptation of an equality constraint  cst x = rhs */
        /* with alpha = cst[dim], 
           we have (cst' + alpha*eqcst).mid x' = rhs - alpha*valcst +
				[-1,1] * (cst' + alpha*eqcst).rad * box */
        void adapt_eqconstraint(const IntervalVector &bbox,
			        Row &cst, Interval &rhs) const;

        /* adaptation of an inequality constraint  cst x <= rhs */
        /* with alpha = cst[dim], 
           we have (cst' + alpha*eqcst) x' <== (rhs - alpha*valcst).ub */
        void adapt_ineqconstraint(IntervalRow &cst, double &rhs) const;

        /* construction of a new eqfacet from an equality constraint 
                  cst x = rhs */
        static EqFacet build_EqFacet(const IntervalVector &bbox,
			const Row &cst, Interval &rhs);
     };
#endif

     Index dim;
     Index fdim; /* dimension - nb independant equalities constraints */
     IntervalVector bbox;

     std::shared_ptr<CollectFacets> facetsptr;

     /* if fdim<dim :
        constraint a X <= b becomes (-b 0) + (a M_EQ) 
        vertex : (c Y) becomes X = M_EQ (1 Y/c) 
        ray : (0 Y) becomes a ray X = M_EQ (0 Y) ?
        lines : (0 Y) becomes a line X = M_EQ (0 Y) ?
     */
     IntervalMatrix M_EQ;

     std::vector<IntervalVector> lines; /* lines (for unbounded polyhedron) */
               /* rays can be handled as classical vertices */

     std::forward_list<DDvertex> vertices;
#if 0
     std::forward_list<DDlink> total_links;
#endif
     std::vector<Index> reffacets; /* indirection for facets */

     bool empty, flat;
     Index nbIn;
     double tolerance = 1e-9;
};

inline Index DDbuildF2V::get_dim() const { return this->dim; }
inline Index DDbuildF2V::get_fdim() const { return this->fdim; }
inline bool DDbuildF2V::is_empty() const { return this->empty; }
inline const IntervalMatrix &DDbuildF2V::get_M_EQ() const { return this->M_EQ; }
inline const std::vector<IntervalVector> &DDbuildF2V::get_lines() const 
	{ return this->lines; }
inline const std::forward_list<DDvertex> &DDbuildF2V::get_vertices() const 
	{ return this->vertices; }
inline const std::vector<Index> &DDbuildF2V::get_reffacets() const 
	{ return this->reffacets; }


struct DDfacet {
   CollectFacets::mapIterator facetIt;
   double lambda;
   
   enum ddfacetstatus {
     FACETNEW,
     FACETIN,
     FACETON,
     FACETOUT,
     FACETREM
   };
   ddfacetstatus status;

   std::vector<std::forward_list<DDfacet>::iterator> links;
   std::vector<Index> vtx;
 
   DDfacet(CollectFacets::mapIterator facetIt,
		std::vector<std::forward_list<DDfacet>::iterator> &nlinks,
		std::vector<Index> &nvtx) : facetIt(facetIt), status(FACETNEW) {
         this->links.swap(nlinks);
         this->vtx.swap(nvtx);
   }

   bool addLnk(std::forward_list<DDfacet>::iterator a, bool check=false) {
      if (check) {
         for (auto b : links) {
            if (a==b) return false;
         } 
      }
      this->links.push_back(a);
      return true;
   }

   void removeLnk(const std::forward_list<DDfacet>::iterator &a) {
      for (Index i=0;i<(Index)links.size();i++) {
         if (links[i]==a) { 
	    links[i]=links[links.size()-1]; 
	    links.pop_back();
            return;
         }
      }
   }

   void changeLnk(const std::forward_list<DDfacet>::iterator &a,
		const std::forward_list<DDfacet>::iterator &b) {
      for (Index i=0;i<(Index)links.size();i++) {
         if (links[i]==a) { links[i]=b; return; }
      }
      /* should not happen */
      links.push_back(b);
   }

   void addVtx(Index a) {
      auto it = this->vtx.begin();
      while (it != this->vtx.end() && (*it)<a) it++;
      if (it!=this->vtx.end() && (*it)==a) return;
      this->vtx.insert(it,a);
   }

};

class DDbuildV2F {
  public:
     DDbuildV2F(Index idVertex, const Vector &vertex);
   
     int add_vertex(Index idVertex, const Vector& vertex);
     int add_vertex(Index idVertex, const IntervalVector& vertex);

     friend std::ostream& operator<<(std::ostream& os, const DDbuildV2F& build);

  private:
     /* can NOT return facets.end() */
     std::forward_list<DDfacet>::iterator addDDfacet(
		const CollectFacets::mapIterator &it,
		std::vector<std::forward_list<DDfacet>::iterator> &links,
		std::vector<Index> &vtx);
     /* can return facets.end() if insertion failed */
     std::forward_list<DDfacet>::iterator addDDfacet(
		std::pair<Row,double> &&row_rhs,
		std::vector<std::forward_list<DDfacet>::iterator> &links,
		std::vector<Index> &vtx);
     /* can return facets.end() if insertion failed */
     std::forward_list<DDfacet>::iterator addDDfacet(
 		const std::pair<Row,double> &row_rhs,
		std::vector<std::forward_list<DDfacet>::iterator> &links,
		std::vector<Index> &vtx);
     /* can return facets.end() if insertion failed */
     std::forward_list<DDfacet>::iterator
	 addFacetSon(const std::forward_list<DDfacet>::iterator &It1, 
		std::forward_list<DDfacet>::iterator &It2, Index idVertex); 
     void removeFacet(std::forward_list<DDfacet>::iterator It,
			 bool removeLnks);

     inline std::pair<Row,double> reduce_facet(const Row &row, double rhs) {
        int fx;
        frexp(row[0],&fx);
        for (Index i=1;i<row.size();i++) {
            int fx2;
            frexp(row[i],&fx2);
            if (fx2>fx) fx=fx2;
         }
         if (std::abs(fx)<10) return std::pair(row,rhs);
         double mult=exp2(-fx);
         return std::pair(row*mult,rhs*mult);
     }
     
     std::shared_ptr<CollectFacets> facetsptr;
     std::forward_list<DDfacet> facets;
     std::vector<Index> idVertices; /* list of vertices entered, useful as
                long as there are eqfacets */
     Index nbIn=0;
     double tolerance =1e-9;
};

}

