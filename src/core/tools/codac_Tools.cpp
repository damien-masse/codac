/** 
 *  Tools
 * ----------------------------------------------------------------------------
 *  \date       2020
 *  \author     Simon Rohou
 *  \copyright  Copyright 2021 Codac Team
 *  \license    This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 */

#include <iostream>
#include <sstream>
#include <algorithm>
#include <functional>
#include "codac_Tools.h"

using namespace std;
using namespace ibex;

namespace codac
{
  void Tools::replace_all(string& input, const string& search, const string& format)
  {
    size_t index = 0;

    while(index != string::npos)
    {
      index = input.find(search, index);
      if(index != string::npos)
      {
        input.replace(index, search.length(), format);
        index += format.length();
      }
    }
  }

  void Tools::trim(string& str)
  {
    str.erase(str.begin(), find_if(str.begin(), str.end(), not1(ptr_fun<int, int>(isspace))));
    str.erase(find_if(str.rbegin(), str.rend(), not1(ptr_fun<int, int>(isspace))).base(), str.end());
  }

  double Tools::rand_in_bounds(const Interval& itv)
  {
    // The seed of the pseudo-random number generator is voluntarily initialized
    // outside this function, on demand.
    return max(itv.lb(),min(itv.ub(),rand()/double(RAND_MAX)*itv.diam()+itv.lb()));
  }
}