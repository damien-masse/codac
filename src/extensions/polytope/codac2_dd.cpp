/**
 * \file codac2_dd.h implementation of the dd algorithm
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
#include <list>

#include "codac2_Index.h"
#include "codac2_Vector.h"
#include "codac2_Row.h"
#include "codac2_Facet.h" 
#include "codac2_dd.h"

#undef DEBUG_CODAC2_DD

namespace codac2 {

DDbuildF2V::DDbuildF2V(const IntervalVector &box) : firstIn(-1) {
   assert_release(!box.is_unbounded());
   /* we use "negative" number for the box constraints :
      one for each negative bounds (-s to -1), one for all positive bounds 
	(-s-1) */
   std::set<Index> links;
   std::vector<Index> fcts;
   for (Index i=-box.size();i<0;i++) {
      fcts.push_back(i);
   }
   for (Index i=1;i<=box.size();i++) {
      links.insert(i);
   }
   IntervalVector l = box.lb();
   this->addDDvertex(l,links,fcts);
   Interval sum = box.sum();
   for (Index i=0;i<box.size();i++) {
      l=box.lb();
      l[i]=l[i]+sum.ub()-sum.lb();
      for (Index j=-box.size()-1;j<0;j++) {
         if (j==-i-1) continue;
         fcts.push_back(j);
      }
      for (Index j=0;j<=box.size();j++) {
         if (j==i+1) continue;
         links.insert(j);
      }
      this->addDDvertex(l,links,fcts);
   }
}

Index DDbuildF2V::addDDvertex(const IntervalVector &v,
		 std::set<Index> &links, std::vector<Index> &fcts) {
    nbIn++;
    vertices.emplace(nbIn,
         DDvertex(v,nbIn,links,vtx));
    return nbIn;
}
Index DDbuildF2V::addDDvertex(IntervalVector &&v,
		 std::set<Index> &links, std::vector<Index> &fcts) {
    nbIn++;
    vertices.emplace(nbIn,
         DDvertex(v,nbIn,links,vtx));
    return nbIn;
}

/* p1 has lambda>0, p2 has lambda<0 ;
   p1 is outside the new constraint, p2 is inside.
   For inequalities, we need to get "as much" p1 as possible, 
   that is maximize |lambda2| and minimize |lambda1| */
Index DDbuildF2V::addVertexSon(const DDvertex &p1, DDvertex &p2,
        Index idFacet, bool isEq) {
    Interval ratio2=1.0/(1.0-p1.lambda/p2.lambda); /* maximize ratio */
    IntervalVector newV = ratio2.ub()*p1.vertex + 
		(Interval(1.0)-ratio2.ub())*p2.vertex;
    /* computation of facets */
    std::vector<Index> nfcts;
    Index a1=0; Index a2=0;
    bool idFacetAdded=false;
    while (a1<(Index)p1.fcts.size() && a2<(Index)p2.fcts.size()) {
       if (p1.fcts[a1]<p2.fcts[a2]) {
          a1++; if (a1==(Index)p1.fcts.size()) break;
       } else if (p2.fcts[a2]<p1.fcts[a1]) {
          a2++; if (a2==(Index)p2.fcts.size()) break;
       } else {
         if (p1.fcts[a1]>idFacet && !idFacetAdded) {
            nfcts.push_back(idFacet);
            idFacetAdded=true;
         }
         nfcts.push_back(p1.fcts[a1]);
         a1++; a2++;
       }
    }
    if (!idFacetAdded) nfcts.push_back(idFacet);
    std::vector<Index> lnks;
    if (!isEq) lnks.push_back(p2.Id);
    Index newId=this->addDDvertex(newV,nfcts,lnks);
    if (!isEq) p2.replaceLnks(p1.Id,newId);
    return newId;
}


struct PairMaxSets {
   Index a,b;
   std::vector<Index> elems;

   PairMaxSets() {}

   PairMaxSets(Index a, Index b, const std::vector<Index> e1, const std::vector<Index> e2) : a(a), b(b) {
       Index i=0,j=0;
       while(i<(Index)e1.size() && j<(Index)e2.size()) {
         if (e1[i]==e2[j]) { 
           elems.push_back(e1[i]);
           i++; j++;
         } else if (e1[i]<e2[j]) i++;
         else j++;
       }
   }

   void swap(PairMaxSets &p2) {
      std::swap(this->a,p2.a);
      std::swap(this->b,p2.b);
      this->elems.swap(p2.elems);
   }

   int comp(const PairMaxSets &pm) const { 
		/* -1 if elems included in pm.elems; 0 if ==,
		   +1 elems includes pm.elems, +2 if not comparable */
       const std::vector<Index> &e1=elems;
       const std::vector<Index> &e2=pm.elems;
       std::vector<long int>::size_type i=0,j=0;
       bool leq=(e1.size()<=e2.size()), geq=(e1.size()>=e2.size()); 
			/* leq=F,geq=T => +1
                           leq=F,geq=F => 2
                           leq=T,geq=F => -1
                           leq=T,geq=T => 0 */
       while(i<e1.size() && j<e2.size()) {
         if (e1[i]<e2[j]) { 
            if (!geq) return 2;
            leq=false; i++;
         } else if (e1[i]>e2[j]) {
            if (!leq) return 2;
            geq=false; j++;
         } else { i++; j++; }
       }
       if (i<e1.size()) {
          if (!geq) return 2;
          leq=false;
       } 
       if (j<e2.size() || !geq) {
          if (!leq) return 2;
          return -1;
       } 
       return (leq ? 0 : 1);
   }
};


int DDbuildF2V::add_facet(Index idFacet, const Facet& fct) {
   IntervalRow rowI = fct.row;
   int cnt=0;
   /* the first step is the modification of the status */
   for (DDvertex &d : vertices) {
      if (d.status==DDvertex::VERTEXREMOVED) continue;
      Interval a = rowI*d.vertex - fct.rhs;
      d.lambda=a.lb();
      if (a.ub()<0.0) { d.status=DDvertex::VERTEXLT; }
      else if (a.lb()<=0.0) { d.status=DDvertex::VERTEXEQ; }
      else { d.status=DDvertex::VERTEXGT;  }
   }
   Index i = firstIn;
   Index fNew=-1;
   while (i!=fNew && i!=-1) {
     DDvertex &d = vertices[i];
     i = d.nextId;
     if (d.status==DDvertex::VERTEXEQ) {
        d.addFct(idFacet);
     }
     if (d.status==DDvertex::VERTEXEQ || d.status==DDvertex::VERTEXLT) continue;
       /* in case fct.eqcst is true, we remove VERTEXLT later */
     for (Index idL : d.links) {
        DDvertex &d2 = vertices[idL];
        if (d2.status==DDvertex::VERTEXEQ) {
           d2.removeLnk(d.Id);
           continue;
        }
        if (d2.status!=DDvertex::VERTEXLT) { 
	   continue; /* may be GT or REMOVED, for each do nothing */
        }
        Index newId = this->addVertexSon(d,d2,idFacet,fct.eqcst);
        cnt++;
        if (fNew==-1) fNew=newId;
     }
     this->removeVertex(d.Id,false);
   }
   /* elimination of old vertices */
   i = firstIn;
   while (i!=fNew) {
      DDvertex &d = vertices[i];
      i = d.nextId;
      if (fct.eqcst && d.status==DDvertex::VERTEXLT) {
          this->removeVertex(d.Id,true);
      }
   } 
   /*construction of new links for the eq vertices */
   i = firstIn;
   std::vector<PairMaxSets> maxset;
   while (i!=-1) {
      DDvertex &d = vertices[i];
      i = d.nextId;
      if (d.status==DDvertex::VERTEXEQ || d.status==DDvertex::VERTEXNEW) {
         Index j=i;
         while (j!=-1) {
            DDvertex &d2 = vertices[j];
            j=d2.nextId;
            if (d2.status==DDvertex::VERTEXEQ || d2.status==DDvertex::VERTEXNEW) {
               PairMaxSets p(d.Id,d2.Id,d.fcts,d2.fcts);
               PairMaxSets &actuel = p;
               bool placed=false;
               unsigned long int k=0, l=0;
               while (k<maxset.size()) {
                  int u=actuel.comp(maxset[k]);
                  if (u==2) { 
		    if (k!=l) maxset[k].swap(maxset[l]);
	    	    k++; l++;
		    continue;
		  }
                  if (u==-1 || u==0) break; 
			/* note : if u==-1 or 0, 
			   placed should be false so k==l */
                  if (u==1) {
                     if (!placed) { 
			maxset[k]=p; 
			placed=true;
			l++;
                     }
                     k++;
                  } 
               }
               if (l<k) maxset.resize(l);
            }
         }
      }
   }
   for (const PairMaxSets &pm : maxset) {
      vertices[pm.a].addLnk(pm.b);
      vertices[pm.b].addLnk(pm.a);
   }
   return cnt;
}

void DDbuildF2V::removeVertex(Index Id, bool removeLnks) {
   DDvertex &d = vertices[Id];
   if (d.prevId!=-1) vertices[d.prevId].nextId=d.nextId;
   if (d.nextId!=-1) vertices[d.nextId].prevId=d.prevId;
   if (firstFree!=-1) vertices[firstFree].prevId=Id;
   d.status=DDvertex::VERTEXREMOVED;
   d.nextId=firstFree;
   firstFree=Id;
   if (removeLnks) {
      for (Index l : d.links) {
         vertices[l].removeLnk(Id);
      }
   }
}

std::ostream& operator<<(std::ostream& os, const DDbuildF2V& build) {
   os << "DDbuildF2V (" << build.vertices.size() << " vtx) : " << std::endl;
   for (const DDvertex &vt : build.vertices) {
      if (vt.status==DDvertex::VERTEXREMOVED) continue;
      os << vt.vertex << std::endl;
      os << " fct: (";
      for (Index a : vt.fcts) {
         os << a << ",";
      }
      os << ")" << std::endl;
   }
   os << "endDDbuildF2V" << std::endl;
   return os;
}

DDbuildV2F::DDbuildV2F(Index idVertex, const Vector &vertex) {
   Row r = Row::Zero(vertex.size());
   for (Index i=0;i<vertex.size();i++) {
      r[i]=1.0;
      eqfacets.emplace_back(r,vertex[i],true);
      r[i]=0.0;
   }
   idVertices.push_back(idVertex);
}

Index DDbuildV2F::addFacetSon(const DDfacet &p1, DDfacet &p2,
        Index idVertex) {
#ifdef CODAC2_DD
    std::cout << " p1.facet.row " <<  p1.facet.row << std::endl;
    std::cout << " p2.facet.row " <<  p2.facet.row << std::endl;
#endif
    Row newV = - p2.lambda*p1.facet.row + p1.lambda*p2.facet.row;
    double newRhs = - p2.lambda*p1.facet.rhs + p1.lambda*p2.facet.rhs;
#ifdef CODAC2_DD
    std::cout << " newV " << newV << std::endl;
    std::cout << " newRhs " << newRhs << std::endl;
#endif
    /* computation of vertices */
    std::vector<Index> nvtx;
    Index a1=0; Index a2=0;
    bool idVertexAdded=false;
    while (a1<(Index)p1.vtx.size() && a2<(Index)p2.vtx.size()) {
       if (p1.vtx[a1]<p2.vtx[a2]) {
          a1++; if (a1==(Index)p1.vtx.size()) break;
       } else if (p2.vtx[a2]<p1.vtx[a1]) {
          a2++; if (a2==(Index)p2.vtx.size()) break;
       } else {
         if (p1.vtx[a1]>idVertex && !idVertexAdded) {
            nvtx.push_back(idVertex);
            idVertexAdded=true;
         }
         nvtx.push_back(p1.vtx[a1]);
         a1++; a2++;
       }
    }
    if (!idVertexAdded) nvtx.push_back(idVertex);
    std::set<Index> lnks;
    lnks.insert(p2.Id);
    Index newId=this->addDDfacet(Facet(newV,newRhs,false),lnks,nvtx); 
    p2.replaceLnks(p1.Id,newId);
    return newId;
}

int DDbuildV2F::add_vertex(Index idVertex, const Vector& vertex) {
    return this->add_vertex(idVertex,IntervalVector(vertex));
}

int DDbuildV2F::add_vertex(Index idVertex, const IntervalVector& vertex) {
   int cnt=0;
   /* first check equalities */
   if (!eqfacets.empty()) {
      Index eqViolated=-1;
      Interval delta;
      for (Index i=eqfacets.size()-1;i>=0;i--) {
          Interval calc = eqfacets[i].row.dot(vertex)-eqfacets[i].rhs;
          if (calc.mig()<this->tolerance) continue;
          if (eqViolated==-1) { eqViolated=i; delta=calc; continue; }
          else {
             double correction = -(calc/delta).mid();
             eqfacets[i].row += correction*eqfacets[eqViolated].row;
             eqfacets[i].rhs += correction*eqfacets[eqViolated].rhs;
          }
      }
      if (eqViolated!=-1) {
         if (facets.empty()) {
            this->nbIn=0;
            std::set<Index> lnk = { 2 };
            std::vector<Index> vtx = { idVertex };
            if (delta.lb()>0.0) {
               this->addDDfacet(Facet(eqfacets[eqViolated].row,
			eqfacets[eqViolated].row.dot(vertex).ub(),false),
		     lnk,vtx);
               vtx = idVertices;
               lnk = { 1 };
               this->addDDfacet(Facet(-eqfacets[eqViolated].row,
			-eqfacets[eqViolated].rhs,false),
		     lnk,vtx);
            } else {
               this->addDDfacet(Facet(-eqfacets[eqViolated].row,
			(-eqfacets[eqViolated].row.dot(vertex)).ub(),false),
		     lnk,vtx);
               vtx = idVertices;
               lnk = { 1 };
               this->addDDfacet(Facet(eqfacets[eqViolated].row,
			eqfacets[eqViolated].rhs,false),
		     lnk,vtx);
           }
           cnt=2;
         } else {
           /* add the remaining inequality, linked with all existing facets */
           std::set<Index> lnk;
           for (std::pair<const Index, DDfacet> &f : facets) {
                 f.second.status=DDfacet::FACETOUT; /* to remove FACETNEW */
		 lnk.insert(f.first);
           }
           Index u;
           std::vector<Index> vtx = idVertices;
           if (delta.lb()>0.0) {
              u = this->addDDfacet(Facet(-eqfacets[eqViolated].row,
				-eqfacets[eqViolated].rhs,false),lnk,vtx);
           } else {
              u = this->addDDfacet(Facet(eqfacets[eqViolated].row,
			eqfacets[eqViolated].rhs,false),lnk,vtx);
           }
           /* change all existing facets to take the new vertex into account */
    
           for (std::pair<const Index, DDfacet> &f : facets) {
               if (f.second.status!=DDfacet::FACETOUT) continue;
               Interval q = f.second.facet.rhs-
				f.second.facet.row.dot(vertex);
	       double correction = (q/delta).mid();
               f.second.facet.row += correction*eqfacets[eqViolated].row;
               f.second.facet.rhs += correction*eqfacets[eqViolated].rhs;
               f.second.addLnk(u);
               f.second.addVtx(idVertex);
           }
           cnt=1;
         }
         if (eqViolated!=(Index)eqfacets.size()-1) {
            eqfacets[eqViolated]=eqfacets[eqfacets.size()-1];
         } 
         eqfacets.pop_back();
         if (!eqfacets.empty()) idVertices.push_back(idVertex);
         return cnt;
      } else idVertices.push_back(idVertex); /* and follows */
   }
   /* no equalities removed */
   for (std::pair<const Index, DDfacet> &f : facets) {
       DDfacet &d = f.second;
       Interval calc = d.facet.row.dot(vertex)-d.facet.rhs;
       d.lambda = calc.lb();
       if (calc.mig()<this->tolerance) { d.status=DDfacet::FACETON; }
       else if (calc.lb()>0.0) { d.status=DDfacet::FACETOUT; }
       else { d.status=DDfacet::FACETIN; }
   } 
   Index endF = nbIn;
   for (std::map<Index,DDfacet>::iterator it = facets.begin();
                  it!=facets.end() && it->first<=endF;) {
     DDfacet &d = it->second;
     ++it;
     if (d.status==DDfacet::FACETON) {
        d.addVtx(idVertex);
        continue;
     }
     if (d.status==DDfacet::FACETIN || d.status==DDfacet::FACETNEW) continue;
     if (d.links.empty()) { /* not supposed to happen */
        /* very specific case of a 1-dimension polyhedra */
        d.vtx.clear();
        d.vtx.push_back(idVertex);
        d.facet.rhs = d.facet.row.dot(vertex).ub();
        continue;
     }
     for (Index idL : d.links) {
        auto it2 = facets.find(idL); 
        if (it2==facets.end()) continue;
        DDfacet &d2 = it2->second;
        if (d2.status==DDfacet::FACETON) {
           d2.removeLnk(d.Id);
           continue;
        }
        if (d2.status!=DDfacet::FACETIN) { 
	   continue; /* do nothing for out */
        }
#ifdef CODAC2_DD
        std::cout << "    addFacet " << d.Id << " " << d2.Id << " " << d.status << " " << d2.status << std::endl;
        std::cout << "    (" << d.lambda << " " << d2.lambda << ")    " << std::endl;
#endif
#ifdef CODAC2_DD
        Index newid = this->addFacetSon(d,d2,idVertex);
        std::cout << " ... done " << newid <<std::endl;
#else
        this->addFacetSon(d,d2,idVertex);
#endif
        cnt++;
     }
     this->removeFacet(d.Id,false);
     
   }
   /*construction of new links for the eq facets */
   std::vector<PairMaxSets> maxset;
   for (std::map<Index,DDfacet>::iterator it = facets.begin();
                  it!=facets.end();++it) {
      DDfacet &d = it->second;
      if (d.status==DDfacet::FACETON || d.status==DDfacet::FACETNEW) {
         for (std::map<Index,DDfacet>::iterator it2 = next(it);
		it2!=facets.end(); ++it2) {
            DDfacet &d2 = it2->second;
            if (d2.status==DDfacet::FACETON || d2.status==DDfacet::FACETNEW) {
               PairMaxSets actuel(d.Id,d2.Id,d.vtx,d2.vtx);
               bool placed=false;
               unsigned long int k=0, l=0;
               while (k<maxset.size()) {
                  int u=actuel.comp(maxset[k]);
                  if (u==2) { 
		    if (k!=l) maxset[k].swap(maxset[l]);
	    	    k++; l++;
		    continue;
		  }
                  if (u==-1 || u==0) { placed=true; break; }
			/* note : if u==-1 or 0, 
			   placed should be false so k==l */
                  if (u==1) {
                     if (!placed) { 
			maxset[k]=actuel;
			placed=true;
			l++;
                     } 
                     k++;
                  } 
               }
               if (l<k) maxset.resize(l);
               else if (!placed) {
                   maxset.push_back(actuel);
               }
            }
         }
      }
   }
   for (const PairMaxSets &pm : maxset) {
//      std::cout << " maxset : " << pm.a << " " << pm.b << std::endl;
      facets.at(pm.a).addLnk(pm.b);
      facets.at(pm.b).addLnk(pm.a);
   }
   return cnt;
}

Index DDbuildV2F::addDDfacet(const Facet &f, std::set<Index> &links,
		std::vector<Index> &vtx) {
    
    nbIn++;
    facets.emplace(nbIn, 
	 DDfacet(f,nbIn,links,vtx));
    return nbIn;
}

Index DDbuildV2F::addDDfacet(Facet &&f, std::set<Index> &links,
		std::vector<Index> &vtx) {
    
    nbIn++;
    facets.emplace(nbIn, 
	 DDfacet(f,nbIn,links,vtx));
    return nbIn;
}

void DDbuildV2F::removeFacet(Index Id, bool removeLnks) {
   DDfacet &d = facets.at(Id);
   if (removeLnks) {
      for (Index l : d.links) {
         auto it = facets.find(l);
         if (it==facets.end()) continue;
         it->second.removeLnk(Id);
      }
   }
   this->facets.erase(Id);
}

std::ostream& operator<<(std::ostream& os, const DDbuildV2F& build) {
   os << "DDbuildV2F (" << build.eqfacets.size() << " eq, " <<
		build.facets.size() << " ineq) : " << std::endl;
   for (const Facet &eqf : build.eqfacets) {
      os << eqf.row << " = " << eqf.rhs << std::endl;
   }
   for (const auto &ifacet : build.facets) {
      auto &ieqf = ifacet.second;
      os << ieqf.Id << " : " << ieqf.facet.row << " <= " << ieqf.facet.rhs << std::endl;
      os << " vtx: (";
      for (Index a : ieqf.vtx) {
         os << a << ",";
      }
      os << ")" << std::endl;
      os << " lnks: (";
      for (Index a : ieqf.links) {
         os << a << ",";
      }
      os << ")" << std::endl;
   }
   os << "endDDbuildV2F" << std::endl;
   return os;
}

}

