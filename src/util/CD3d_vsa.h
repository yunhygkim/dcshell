#pragma once

//implment Variational Shape Approximation

#include "util/mathtool/Point.h"
#include "util/mathtool/Vector.h"
#include <list>


struct VAS_plane
{
	mathtool::Vector3d n;
	mathtool::Point3d p;
};

//
// compute Variational Shape Approximation
//
class cd_m;
std::list<VAS_plane> VSA(cd_m * m, int k);
