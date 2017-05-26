
//
// This implements half-space intersection
//

#pragma once


#include "util/mathtool/Point.h"
#include "util/mathtool/Vector.h"
#include "util/CD3d_vsa.h" //needed for VAS_plane
#include <list>

//
//compute the intersection of the half-spaces defined by the give planes and 
//the point o that is inside the intersection
//
//implemented using qhull
//
cd_m * Hspace_intersection(std::list<VAS_plane>& planes, const mathtool:: Point3d & o);