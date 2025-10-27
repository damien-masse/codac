/**
 * \file codac2_Polytope_util.h utilities for polytope 
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
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <ctime>
#include <cmath>
#include <bitset>
#include <memory>
#include <map>
#include <forward_list>
#include <set>
#include <stdexcept>

#include "codac2_Index.h"
#include "codac2_Vector.h"
#include "codac2_Row.h"
#include "codac2_Facet.h"
#include "codac2_dd.h"


namespace codac2 {

/** read a .ine file and return a list of facets
 *  this is a very crude function, which assumes the format of the 
 *  file is correct 
 *  @param filename name of the file */
std::vector<Facet> read_ineFile(const char *filename);

/** construct the list of facets of a 3D-polyhedron from the buildF2V 
 * structure 
 * use the mid of each vertex. For rays/line use an arbitrary bound 
 * @param build the buildF2V structure */
std::vector<std::vector<Vector>> build_3Dfacets(const DDbuildF2V &build, double bound=50.0);
   
}

