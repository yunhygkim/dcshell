#pragma once

//implment Variational Shape Approximation

#include "util/mathtool/Point.h"
#include "util/mathtool/Vector.h"
#include <list>

#include "measure/Cd3d_hull.h"

//
// compute progress hull
//

//
// max_fsize: maximum number of faces left in the hull
// max_vol_inc : max vol increase allowed, if max_vol_inc is violated, the function will return the hull even though it has more faces than "max_fsize"
// enclosing_constrains: a list of planes that enforce the new hull the be constrained inside
//
// !! Note: this function assumes that the input_hull is inside enclosing_constrains
//
bool ProgressiveHull(const cd_hull& input_hull, cd_hull& simplified_hull, int max_fsize, float max_vol_inc, const list<hull_plane>& enclosing_constrains);
bool ProgressiveHullQP(const cd_hull& input_hull, cd_hull& remeshed_hull, float max_vol_inc, const list<hull_plane>& enclosing_constrains);
