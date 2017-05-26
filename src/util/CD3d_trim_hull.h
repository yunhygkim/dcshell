#pragma once

//implment Variational Shape Approximation

#include "util/mathtool/Point.h"
#include "util/mathtool/Vector.h"
#include <list>

#include "measure/Cd3d_hull.h"

//
// trim overlapping hulls
// given two convex hulls of A and B called hullA an hullB, return the plane "trim_plane"
// so that the volume of A that is outside "hullA\trim_plane" and the volume of B that is outside "hullB\trim_plane"
// is minimized!
//
// if the hulls do not overlap return false
//

//using svm
bool trim_overlapping_hulls_svm(const cd_hull& hullA, const cd_hull& hullB, hull_plane& trim_plane, double svm_C);
bool trim_overlapping_hulls_svm_opt(const cd_hull& hullA, const cd_hull& hullB, hull_plane& trim_plane, double svm_C);

//using simple heuristic
bool trim_overlapping_hulls_heuristic(const cd_hull& hullA, const cd_hull& hullB, hull_plane& trim_plane);


//
// given a hull and a list of trim planes, return a trimed hull
//
//
bool trim_hull_by_planes(const cd_hull& hull, const list<hull_plane>& trim_planes, cd_hull& trimed_hull, const Point3d& o);
bool trim_hull_by_planes(const cd_hull& hull, const list<hull_plane>& trim_planes, cd_hull& trimed_hull, double vertex_collapse_tolerance);

//
// Expand hull by extra points
//
// max_vol_inc: max volume difference allowed between hull and exp_hull
// new_pts_sets: input points, in sets
//
bool expand_hull_by_points(const cd_hull& hull, const list<hull_plane>& trim_planes, vector< vector<Point3d> >& new_pts_sets, cd_hull& exp_hull, double max_vol_inc);
