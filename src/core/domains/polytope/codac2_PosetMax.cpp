/**
 * \file codac2_PosetMax.cpp   Sets of integers for polytope constraints
 * ----------------------------------------------------------------------------
 *  \date       2026
 *  \author     Damien Massé
 *  \copyright  Copyright 2026
 *  \license    This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 */

#include <utility>
#include <bitset>
#include <cmath>
#include <set>
#include <cstdint>

#include "codac2_Index.h"
#include "codac2_PosetMax.h"

using namespace codac2;

namespace codac2 {

PosetElem::PosetElem(const PosetElem &p1, const PosetElem &p2) : nb_elems(0) {
    std::vector<SetBitset_unit>::size_type i=0,j=0;
    while(i<p1.elems.size() && j<p2.elems.size()) {
      if (p1.elems[i].first==p2.elems[j].first) {
         poset_bitset inter = p1.elems[i].second & p2.elems[j].second;
         if (inter.any()) {
            this->nb_elems += inter.count();
            this->elems.emplace_back(p1.elems[i].first,inter);
         }
         i++; j++;
      } else if (p1.elems[i].first<p2.elems[j].first) {
          i++;
      } else {
          j++;
      }
   }

}
PosetElem operator&(const PosetElem &p1, const PosetElem &p2) {
    return PosetElem(p1,p2);
}

bool PosetElem::has(Index e) const {
     int64_t n;
     size_t a = SetBitset_unit::base_offset(e,n);
     auto it = elems.rbegin(); /* note : reverse iteration */
     while (it!=elems.rend() && it->first>n) {
         ++it;
     }
     if (it==elems.rend() || it->first<n) return false;
     return (it->second[a]);
}

bool PosetElem::add(Index e) {
     int64_t n;
     size_t a = SetBitset_unit::base_offset(e,n);
     auto it = elems.rbegin(); /* note : reverse iteration */
     while (it!=elems.rend() && it->first>n) {
         ++it;
     }
     if (it==elems.rend() || it->first<n) {
         elems.emplace(it.base(),n,(1<<a));
         this->nb_elems++;
         return true;
     } 
     if (it->second[a]) return false;
     it->second[a]=true; 
     this->nb_elems++;
     return true;
}

bool PosetElem::remove(Index e) {
     int64_t n;
     size_t a = SetBitset_unit::base_offset(e,n);
     auto it = elems.rbegin();
     while (it!=elems.rend() && it->first>n) {
         ++it;
     }
     if (it==elems.rend() || it->first<n) return false;
     if (!it->second[a]) return false;
     it->second[a]=false;
     this->nb_elems--;
     if (it->second.none()) {
        elems.erase(it.base()-1);
     }
     return true;
}

void PosetElem::reverse_relation(Index id, std::unordered_map<Index,PosetElem> &mapping) const {
    for (const auto &s : this->elems) {
       for (size_t i=0;i<poset_size_bitset;i++) {
         if (s.second[i]) 
             mapping[SetBitset_unit::value(s.first,i)].add(id);
       }
    }
}

PosetElem PosetElem::renumber(const std::vector<Index> vect) const {
    PosetElem ret;
    for (const auto &s : this->elems) {
       if (s.first<0) continue;
       for (size_t i=0;i<poset_size_bitset;i++) {
         if (s.second[i]) {
             Index v = SetBitset_unit::value(s.first,i);
             if (vect[v]<0) continue;
             ret.add(v);
         }
       }
    }
    return ret;
}


int PosetElem::comp(const PosetElem &pm) const {
              /* -1 if elems included in pm.elems; 0 if ==,
                 +1 elems includes pm.elems, +2 if not comparable */
     const std::vector<SetBitset_unit> &e1=elems;
     const std::vector<SetBitset_unit> &e2=pm.elems;
     std::vector<SetBitset_unit>::size_type i=0,j=0;
     bool leq=(this->nb_elems<=pm.nb_elems), geq=(this->nb_elems>=pm.nb_elems);
     bool same_length=(this->nb_elems==pm.nb_elems);
                      /* leq=F,geq=T => +1
                         leq=F,geq=F => 2
                         leq=T,geq=F => -1
                         leq=T,geq=T => 0 */
     while(i<e1.size() && j<e2.size()) {
       if (e1[i].first==e2[j].first) {
          poset_bitset inter = e1[i].second & (~e2[j].second);
          if (inter.any()) {
             if (!geq || same_length) return 2;
             leq=false;
          }
          inter = (~e1[i].second) & e1[j].second;
          if (inter.any()) {
             if (!leq || same_length) return 2;
             geq=false;
          }
          i++; j++;
       } else if (e1[i].first<e2[j].first) {
          if (!geq || same_length) return 2;
          leq=false;
          i++;
       } else {
          if (!leq || same_length) return 2;
          geq=false;
          j++;
       }
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


std::ostream& operator<<(std::ostream &out, const PosetElem &p) {
   out << p.nb_elems << "(" ;
   for (const auto &s : p.elems) {
      for (size_t i=0;i<poset_size_bitset;i++) {
        if (s.second[i]) 
            out << SetBitset_unit::value(s.first,i) << ",";
      }
   }
   out << ")";
   return out;
}

}
