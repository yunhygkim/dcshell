#pragma once
#include <CD3d_hull.h>

//compute the volume of intersection between hull and p
//without computing the intersection explicitly
double compute_volume(const cd_hull& hull, hull_plane& hp, bool print=false);

//compute the volume of derivative between hull and p
//without computing the intersection explicitly
void compute_cone_volume_derivative_coeeficient
(cd_v * v, hull_plane& hp, double& c2, double& c1, double& c0);

// test only
double compute_volume_explicit(const cd_hull& hull, hull_plane& hp);
