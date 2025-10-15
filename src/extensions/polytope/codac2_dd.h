/**
 * \file codac2_dd.h implementation of the dd algorithm
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
#include <set>
#include <unordered_set>
#include <stdexcept>

#include "codac2_Index.h"
#include "codac2_Vector.h"
#include "codac2_Row.h"
#include "codac2_Facet.h"


namespace codac2 {

struct DDlink; /* link between two vertices */

struct DDvertex {
   Index Id;
   IntervalVector vertex;
   Interval lambda;
   
   enum ddvertexstatus {
     VERTEXNEW,
     VERTEXLT, /* < (GT,LT) => new vertex */
     VERTEXON, /* <,=,> (GT,LE) => for now, equivalent to LT */
     VERTEXGE, /* = and possible small >
		  keep the vertex. no change with (GT,GE) */
     VERTEXGT, /* remove vertex, add new vertices for links */
   };
   ddvertexstatus status;

   std::vector<std::map<Index,DDlink>::iterator> links;
 
   DDvertex(const IntervalVector &v, Index Id) : Id(Id), 
	 vertex(v), status(VERTEXNEW) {
   }

   DDvertex(IntervalVector &&v, Index Id, 
		std::vector<std::map<Index,DDlink>::iterator> &nlinks) :
	 Id(Id), 
	 vertex(v), status(VERTEXNEW) {
         this->links.swap(nlinks);
   }

   void addLnk(const std::map<Index,DDlink>::iterator &a) {
      this->links.push_back(a);
   }

   void addLnkVertex(const std::map<Index,DDvertex>::iterator &V1,
		 const std::vector<Index> &fcts);

   void removeLnk(const std::map<Index,DDlink>::iterator &a);

#if 0
   void addFct(Index a) {
      auto rev_it = this->fcts.rbegin();
      while (rev_it != this->fcts.rend() && (*rev_it)>a) it++;
      if (rev_it!=this->fcts.rend() && (*rev_it)==a) return;
      this->fcts.insert(rev_it.base(),a);
   }
#endif
};

struct DDlink {
   std::map<Index,DDvertex>::iterator V1;
   std::map<Index,DDvertex>::iterator V2;
   std::unordered_set<std::vector<Index>> fcts;
  
   DDlink(Index dim, const std::map<Index,DDvertex>::iterator &V1,
		const std::map<Index,DDvertex>::iterator &V2,
		const std::vector<Index> &fcts) : V1(V1), V2(V2) {
       this->fcts.insert(fcts);
   }

   void changeVertex(const std::map<Index,DDvertex>::iterator &Vold,
                const std::map<Index,DDvertex>::iterator &Vnew) {
       if (this->V1==Vold) { this->V1=Vnew; return; }
       this->V2=Vnew;
   }
   
   const std::map<Index,DDvertex>::iterator &otherVertex
		(const std::map<Index,DDvertex>::iterator &V) const {
       if (V==V1) return V2;
       return V1;
   }
  
   void addFct(const std::vector<Index> &fcts) {
      this->fcts.insert(fcts);
   }
};

class DDbuildF2V {
  public:
     DDbuildF2V(Index dim, const IntervalVector &box,
		const std::vector<Facet> &eqconstraints);
   
     int add_facet(Index idFacet, const Facet& fct);

     friend std::ostream& operator<<(std::ostream& os, const DDbuildF2V& build);

  private:
     std::map<Index,DDvertex>::iterator initDDvertex(Index id,
		 const IntervalVector &v, std::vector<Index> &fcts);
     Index initDDlinks(const std::map<Index,DDvertex>::iterator& V1,
   		      const std::map<Index,DDvertex>::iterator& V2,
   		      std::vector<Index> &fct);

     std::map<Index,DDvertex>::iterator addDDvertex(const IntervalVector &v,
		std::vector<Index> &fcts);
     std::map<Index,DDvertex>::iterator addDDvertex(IntervalVector &&v,
	        std::vector<Index> &fcts);
     std::map<Index,DDvertex>::iterator addVertexSon(const DDvertex &p1,
		 DDvertex &p2, double rhs, Index idFacet); 
     void removeVertex(Index Id, bool removeLnks);

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

     IntervalVector recompute_vertex(const IntervalVector &vect) const; 

     Index dim;
     Index fdim;
     IntervalMatrix M_EQ1, M_EQ2;
     std::vector<EqFacet> eqfacets;
     std::vector<Index> dim_facets;

     std::map<Index,DDvertex> vertices;
     std::map<Index,DDlink> total_links;
     std::map<Index,std::pair<IntervalRow,double>> ineqfacets;

     bool empty;
     Index nbIdVertex;
     Index nbIdLinks;
     double tolerance = 1e-9;
};

struct DDfacet {
   Index Id;
   Facet facet;
   double lambda;
   
   enum ddfacetstatus {
     FACETNEW,
     FACETIN,
     FACETON,
     FACETOUT
   };
   ddfacetstatus status;

   std::set<Index> links;
   std::vector<Index> vtx;
 
   DDfacet(const Facet &v, Index Id, 
		std::set<Index> &nlinks,
		std::vector<Index> &nvtx) : Id(Id), facet(v), status(FACETNEW) {
         this->links.swap(nlinks);
         this->vtx.swap(nvtx);
   }
   DDfacet(Facet &&v, Index Id, 
		std::set<Index> &nlinks,
		std::vector<Index> &nvtx) : Id(Id), facet(v), status(FACETNEW) {
         this->links.swap(nlinks);
         this->vtx.swap(nvtx);
   }

   void replaceLnks(Index old, Index nw) {
      this->links.erase(old);
      this->links.insert(nw);
   }

   void addLnk(Index a) {
      this->links.insert(a);
   }

   void addVtx(Index a) {
      auto it = this->vtx.begin();
      while (it != this->vtx.end() && (*it)<a) it++;
      if (it!=this->vtx.end() && (*it)==a) return;
      this->vtx.insert(it,a);
   }

   void removeLnk(Index a) {
      this->links.erase(a);
   }
};

class DDbuildV2F {
  public:
     DDbuildV2F(Index idVertex, const Vector &vertex);
   
     int add_vertex(Index idVertex, const Vector& vertex);
     int add_vertex(Index idVertex, const IntervalVector& vertex);

     friend std::ostream& operator<<(std::ostream& os, const DDbuildV2F& build);

  private:
     Index addDDfacet(const Facet &f, std::set<Index> &links,
		std::vector<Index> &vtx);
     Index addDDfacet(Facet &&f, std::set<Index> &links,
		std::vector<Index> &vtx);
     Index addFacetSon(const DDfacet &p1, DDfacet &p2, Index idVertex); 
     void removeFacet(Index Id, bool removeLnks);
     
     std::map<Index,DDfacet> facets;
     std::vector<Facet> eqfacets; /* no need to 'link' eqfacets */
     std::vector<Index> idVertices; /* list of vertices entered, useful as
                long as there are eqfacets */
     Index nbIn=0;
     double tolerance =1e-9;
};

}

