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

#include "codac2_Index.h"
#include "codac2_Vector.h"
#include "codac2_Row.h"
#include "codac2_Facet.h"


namespace codac2 {

struct DDvertex {
   Index Id;
   IntervalVector vertex;
   Interval lambda;
   
   enum ddvertexstatus {
     VERTEXNEW,
     VERTEXLT,
     VERTEXGT,
     VERTEXEQ,
   };
   ddvertexstatus status;

   std::set<Index> links;
   std::vector<Index> fcts;
 
   DDvertex(const IntervalVector &v, Index Id, 
		std::set<Index> &nlinks,
		std::vector<Index> &nfcts) : Id(Id), 
	 vertex(v), status(VERTEXNEW) {
         this->links.swap(nlinks);
         this->fcts.swap(nfcts);
   }
   DDvertex(IntervalVector &&v, Index Id, 
		std::set<Index> &nlinks,
		std::vector<Index> &nfcts) : Id(Id), 
	 vertex(v), status(VERTEXNEW) {
         this->links.swap(nlinks);
         this->fcts.swap(nfcts);
   }

   void replaceLnks(Index old, Index nw) {
      this->links.erase(old);
      this->links.insert(nw);
   }

   void addLnk(Index a) {
      this->links.insert(a);
   }

   void addFct(Index a) {
      auto it = this->fcts.begin();
      while (it != this->fcts.end() && (*it)<a) it++;
      if (it!=this->fcts.end() && (*it)==a) return;
      this->fcts.insert(it,a);
   }

   void removeLnk(Index a) {
      this->links.erase(a);
   }
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

class DDbuildF2V {
  public:
     DDbuildF2V(const IntervalVector &box);
   
     int add_facet(Index idFacet, const Facet& fct);

     friend std::ostream& operator<<(std::ostream& os, const DDbuildF2V& build);

  private:
     Index addDDvertex(const IntervalVector &v, std::set<Index> &links,
		std::vector<Index> &fcts);
     Index addDDvertex(IntervalVector &&v, std::set<Index> &links,
		std::vector<Index> &fcts);
     Index addVertexSon(const DDvertex &p1, DDvertex &p2, 
				Index idFacet, bool isEq); 
     void removeVertex(Index Id, bool removeLnks);
     
     std::map<Index,DDvertex> vertices;
     Index nbIn;
     double tolerance = 1e-9;
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

