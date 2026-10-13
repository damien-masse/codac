/**
 * \file codac2_PosetMax.h   Sets of (mostly contiguous) integers
 * ----------------------------------------------------------------------------
 *  \date       2026
 *  \author     Damien Massé
 *  \copyright  Copyright 2026
 *  \license    This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 */

#pragma once

#include <utility>
#include <bitset>
#include <cmath>
#include <set>
#include <cstdint>
#include <map>
#include <unordered_map>
#include <vector>
#include <memory>

#include "codac2_Index.h"

//class DDvertex;

namespace codac2 {
   /* the poset is an (ordered) vector of Index */
   /* used in .h because of its templated version */

/*** elements ****************************************************/
/*** we may store more that set of elems to speed up computation */

/* each set is stored as a vector of pairs (unsigned integer / bitset)
 * each bitset is used to describe 64 elements */
const std::size_t poset_size_bitset = 
		std::numeric_limits<unsigned long long>::digits;
using poset_bitset = std::bitset<poset_size_bitset>;
struct SetBitset_unit : std::pair<int64_t,poset_bitset> {

   SetBitset_unit(Index b,poset_bitset s) :
         std::pair<int64_t,poset_bitset>(b,s) { }

   static inline size_t base_offset(Index n, int64_t &cbase) {
     if (n>=0) { cbase=n/poset_size_bitset; return n%poset_size_bitset; }
     cbase = (n-(poset_size_bitset-1))/poset_size_bitset;
     return n-(cbase/poset_size_bitset);
   }
   static inline Index value(int64_t cbase, size_t offset) {
     return cbase*poset_size_bitset+offset;
   }
};


/** \brief Pair of an element of type T and an (ordered) vector of
 *  Index, along with the partial (subset) order on the vector */
/*  A total order is also defined to compare poset elems and create */
/*  a hashmap structure */
struct PosetElem {
   std::vector<SetBitset_unit> elems;
   Index nb_elems;

   /** \brief empty element */
   PosetElem() : nb_elems(0) {}

   /** \brief build an element (e1 inter e2)
    *  \param e1 ordered set
    *  \param e2 ordered set */
   PosetElem(const PosetElem &p1, const PosetElem &p2);

   /** \brief build an element (e1 inter e2)
    *  \param e1 ordered set
    *  \param e2 ordered set */
   friend PosetElem operator&(const PosetElem &p1, const PosetElem &p2);

   void clear() {
     elems.clear();
     nb_elems=0;
   }

   /** \brief has an element 
    *  \param e the element
    *  \return true if e is in the set
    */  
   bool has(Index e) const;

   /** \brief add an element 
    *  \param e the element
    *  \return true if added, false if already exists 
    */  
   bool add(Index e); 

   /** \brief removed an element 
    *  \param e the element
    *  \return true if removed, false if absent
    */  
   bool remove(Index e); 

   /** \brief swap two PosetElem
    *  \param p2 the element to be swapped with (*this)
    */
   void swap(PosetElem &p2) {
      std::swap(this->nb_elems,p2.nb_elems);
      this->elems.swap(p2.elems);
   }
 
   /** \brief comparison between the set of indices
    *  \param pm the element to be compared
    *  \return -1 if (*this).elems is included in pm.elems
    *           0 if both sets are equal
    *           1 if (*this).elems is a superset of pm.elems
    *           2 if both sets are not comparable
    */
   int comp(const PosetElem &pm) const;

   /** \build a new PosetElem by renumbering an existing one :
    *  negative values are kept. non-negative values associated (by vect)
    *  to negative are removed. other are replaced by the vect value 
    *  \param vect the renumbering vector
    */
   PosetElem renumber(const std::vector<Index> vect) const;

   /** \brief ``reverse'' relation, i.e. add id to mapping[s] for all
    * s in the posetElem
    *  \param Id the Index to add
    *  \param map the mapI (Index -> PosetElem)
    */
   void reverse_relation(Index id, std::unordered_map<Index,PosetElem> &mapping) const;

   /** \brief display (for debug) */
   friend std::ostream& operator<<(std::ostream &out, const PosetElem &p);
};

/* comparisons of PosetElem for map. Note: we give a priority order, i.e.
   bigger sets are first */
struct PosetElemCmp
{
   bool operator()(const PosetElem &e1, const PosetElem &e2) const {
      if (e1.nb_elems>e2.nb_elems) return true;
      if (e1.nb_elems<e2.nb_elems) return false;
      std::vector<long int>::size_type i=0,j=0;
      while(i<e1.elems.size() && j<e2.elems.size()) {
        if (e1.elems[i].first<e2.elems[j].first)
          return true;
        if (e1.elems[i].first>e2.elems[j].first)
          return false;
        if (e1.elems[i].second.to_ullong()>e2.elems[j].second.to_ullong()) 
          return true;
        if (e1.elems[i].second.to_ullong()<e2.elems[j].second.to_ullong()) 
          return false;
        i++; j++;
      }
      /* normally, the only possible remaining case is equality */
      return false;
   }
}; 

template <typename T>
using PosetElemMap = std::multimap<PosetElem,T,
	PosetElemCmp>;

/** \brief stucture to handle sets of maximal elements */
template <typename T>
struct PosetElemMaxSet : std::vector<std::pair<T,PosetElem>> {

 
   bool addElement(std::pair<const T,PosetElem> &&p) {
      bool ismax=true;
      PosetElem &elem = p.second;
      unsigned long int k=0,l=this->size();
      while (k<l) {
          int u=elem.comp((*this)[k].second);
          if (u==2) { 
             k++;
             continue;
          }
          if (u==-1 || u==0) { ismax=false; break;  }
          if (u==1) {
             l--;
             if (l==k) break;
             std::swap((*this)[k],(*this)[l]);
          }
      }
      if (ismax) {
         this->resize(l+1);
         (*this)[l].first=p.first;
         (*this)[l].second=std::move(p.second);
      } else if (l!=this->size()) { /* should not happen */
         this->resize(l);
      }
      return ismax;
   }
    
   /** \brief move in a new element (if maximal)
    *  \param p the pair (T, poset) (note: the poset may be swapped)
    *  \return true if the element is added, false otherwise
    */
   bool addElement(std::pair<const T,PosetElem> &p) {
      bool ismax=true;
      PosetElem &elem = p.second;
      unsigned long int k=0,l=this->size();
      while (k<l) {
          int u=elem.comp((*this)[k].second);
          if (u==2) { 
             k++;
             continue;
          }
          if (u==-1 || u==0) { ismax=false; break;  }
          if (u==1) {
             l--;
             if (l==k) break;
             std::swap((*this)[k],(*this)[l]);
          }
      }
      if (ismax) {
         this->resize(l+1);
         (*this)[l].first=p.first;
         std::swap((*this)[l].second,elem);
      } else if (l!=this->size()) { /* should not happen */
         this->resize(l);
      }
      return ismax;
   }
};


}
