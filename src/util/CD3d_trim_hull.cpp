#include "util/CD3d_trim_hull.h"
#include "util/CD3d_util.h"
#include "util/CD3d_hspace_intersection.h"
#include "util/CD.h"
#include "intersection.h"
#include "hull_vol.h"

#include "CD3d_model.h"
#include "svm.h"
#include<map>
#include<set>
#include<algorithm>

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>

using namespace mathtool;
using namespace std;
//
//
//void sample_from_triangle(cd_f * f, double d, vector<Point3d>& samples)
//{
//	const Point3d& p1 = f->v[0]->pos;
//	const Point3d& p2 = f->v[1]->pos;
//	const Point3d& p3 = f->v[2]->pos;
//
//	//compute triangle area
//	Vector3d v1 = p2 - p1;
//	Vector3d v2 = p3 - p1;
//	Vector3d v3 = p3 - p2;
//	Vector3d v4 = p1 - p3;
//	float area = (v1%v2).norm() / 2;
//
//
//	//compute number of points needed
//	int size = (int)ceil(area / (d*d));
//
//	double r[2];
//
//	//randomly generate these points
//	while (size>0)
//	{
//		//halton(r);
//		//hammersley_sequence(1,r);
//
//		r[0] = drand48();
//		r[1] = drand48();
//		Point3d pos = p1 + v1*r[0] + v2*r[1];
//
//		//check if inside the triangle
//		if ((v1 % (pos - p1))*f->n <= 0) continue;
//		if ((v3 % (pos - p2))*f->n <= 0) continue;
//		if ((v4 % (pos - p3))*f->n <= 0) continue;
//		//in tri
//
//		samples.push_back(pos);
//		size--;
//	}//j
//}
//
////
////make k samples uniformlly inside the hull
////
////this is NOT uniform so we will have to improve this later
////
//void sample_in_hull(const cd_hull& hull, cd_m * hull_m, vector<Point3d>& samples, int k)
//{
//	//
//	float bbox[6];
//	computeBoundingBox(hull_m->v(), bbox);
//
//	//sample inside the hull
//	for (int i = 0; i < k; i++)
//	{
//		Point3d sample;
//		//draw a sample from box
//		for (int j = 0; j < 3; j++)
//			sample[j] = (bbox[j * 2 + 1] - bbox[j * 2])*drand48() + bbox[j * 2];
//		//check if the point is inside the hull
//		if (hull.isInsideHull(sample))
//			samples.push_back(sample);
//	}
//
//}
//
//void sample_from_hull(const cd_hull& hull, vector<Point3d>& samples, int k)
//{
//	cd_m * m = hull.to_m();
//
//	//sample the points from inside the hull
//	//sample_in_hull(hull, m, samples, k);
//
//	//sample the points from surface of the hull
//	for (auto& f : m->f())
//	{
//		sample_from_triangle(f, 1, samples);
//	}
//
//	//now add the vertice of the hull into the sample as well.
//	for (auto& p : m->v()) samples.push_back(p->pos);
//
//	m->destroy();
//	delete m;
//}

//Red/Blue points
struct RBPoint
{
  RBPoint(const Point3d& _p, bool _blue):pos(_p), isblue(_blue), v(NULL){ dist=0; inIntersection=false;}
  RBPoint(cd_v * _v, bool _blue):pos(_v->pos), isblue(_blue), v(_v){ dist=0; inIntersection=false;}

  cd_v * v;
  Point3d pos;
  bool isblue;
  bool inIntersection;
  float dist;
  bool operator<(const RBPoint& p) const { return dist<p.dist; }
};

vector<Point3d> all_sur_pts;
vector<Point3d> all_sv_A;
vector<Point3d> all_sv_B;

//make sure that the point p and q must have different lables
Point3d bisection(svm_model * model, const Point3d& p, const Point3d& q, double tau_sqr);

double least_square_fit(const vector<Point3d>& pts, Point3d& o, Vector3d& n); //implemented below

Point3d getPointInIntersection(const cd_hull& hull, const list<hull_plane>& trim_planes);

//given twp facets, compute their interesction as line segment
void tri_tri_intersection(cd_f * f1, cd_f* f2, vector<Point3d>& intersections, float resolution=1e10);

//given two overlapping convex object, find the best position o so that
//(the volume of B on positive side of the plane)+(the volume of A on on the negtive side of the plane)
//is miminized.
void optimize_position_exact(const cd_hull& hullA, const cd_hull& hullB, Point3d& o, const Vector3d& n);

//given red/blue point sets, find the best position o so that
//# of blue points on the positive side of the plane and # of red points on the negtive side of the plane
//is miminized.
void optimize_position(vector<RBPoint>& pts, Point3d& o, const Vector3d& n);

//given red/blue point sets, find the best orientation n so that
//the sum of # of blue points on the positive side of the plane and # of red
//points on the negtive side of the plane
//is miminized.
void optimize_orientation(vector<RBPoint>& pts, const Point3d& o,  Vector3d& n);
void optimize_orientation_breute_force(vector<RBPoint>& pts, const Point3d& o,  Vector3d& n);

inline void build_red_blue(const cd_hull& hullA, const cd_hull& hullB, vector<RBPoint>& pts)
{
  pts.reserve(hullA.getSampledPoints().size()+hullB.getSampledPoints().size());

  //red points
	for(auto& pt:hullA.getSampledPoints()){
    RBPoint pt2(pt,false);
    pt2.inIntersection=hullB.isInsideHull(pt);
    pts.push_back(pt2);
  }

  //blue points
	for(auto& pt:hullB.getSampledPoints()){
    RBPoint pt2(pt,true);
    pt2.inIntersection=hullA.isInsideHull(pt);
    pts.push_back(pt2);
  }

}

inline void build_red_blue(cd_m * A, cd_m * B, vector<RBPoint>& pts)
{
  pts.reserve(A->v().size()+B->v().size());

  //red points
	for(cd_v * v:A->v()){
    RBPoint pt2(v,false);
    pts.push_back(pt2);
  }

  //blue points
  for(cd_v * v:B->v()){
    RBPoint pt2(v,true);
    pts.push_back(pt2);
  }
}

inline void getIntersectionPoints(const cd_hull& hullA, const cd_hull& hullB, vector<Point3d>& intersections, double res=1e10)
{
  //find all pairs of intersecting triangles between A & B
  cd_m * hullA_m = hullA.to_m();
  cd_m * hullB_m = hullB.to_m();
  list< pair<int,int> > contacts;
  {
    CD_RAPID cd;
    int Acd = cd.buildModel(*hullA_m);
    int Bcd = cd.buildModel(*hullB_m);
    cd.isInCollision(Acd,Bcd,contacts);
  }

  //find the actual intersection between all pairs of triangles
  {
    vector<cd_f *> Af(hullA_m->f().begin(),hullA_m->f().end());
    vector<cd_f *> Bf(hullB_m->f().begin(),hullB_m->f().end());

    for(auto& contact: contacts)
    {
      auto & ta = Af[contact.first];
      auto & tb = Bf[contact.second];

      vector<Point3d> myintersections;
      //find intersection between ta & tb
      tri_tri_intersection(ta, tb, myintersections,res);
      for(auto& pt:myintersections)
      {
        bool duplicated=false;
        for(auto& pt2:intersections)
        {
          if(AlmostEqual3(pt.get(),pt2.get()))
          {
            duplicated=true;
            break;
          }
        }//end for pt2

        if(duplicated==false)
          intersections.push_back(pt);
      }//end for pt
    }
  }

  hullA_m->destroy();
  hullB_m->destroy();
}

bool trim_overlapping_hulls_svm(vector<RBPoint>& pts, hull_plane& trim_plane, double svm_C)
{
	//use libsvm to find the separating plane
	struct svm_parameter param;

	// default values
	param.svm_type = C_SVC; //C_SVC, NU_SVC
	param.kernel_type = LINEAR;
	param.cache_size = 100;
	param.C = svm_C;
	param.eps = 1e-3;
	param.shrinking = 1;
	param.probability = 0;
	param.nr_weight = 0;
	param.weight_label = NULL;
	param.weight = NULL;
	param.gamma=0.6;
	param.degree = 3;

	//prepare the svm problem
	svm_problem prob;
	prob.l = pts.size();
	svm_node *x_space = new svm_node[4 * prob.l];
	assert(x_space);
	prob.x = new svm_node *[prob.l];
	prob.y = new double[prob.l];

	int i = 0;
	for (auto& pt : pts)
	{
		x_space[4 * i].index = 1;
		x_space[4 * i].value = pt.pos[0];
		x_space[4 * i + 1].index = 2;
		x_space[4 * i + 1].value = pt.pos[1];
		x_space[4 * i + 2].index = 3;
		x_space[4 * i + 2].value = pt.pos[2];
		x_space[4 * i + 3].index = -1;
		prob.x[i] = &x_space[4 * i];
		prob.y[i] = (pt.isblue)?-1:1;
		i++;
	}

	// for (auto& pt : hullB_pts)
	// {
	// 	x_space[4 * i].index = 1;
	// 	x_space[4 * i].value = pt[0];
	// 	x_space[4 * i + 1].index = 2;
	// 	x_space[4 * i + 1].value = pt.pos[1];
	// 	x_space[4 * i + 2].index = 3;
	// 	x_space[4 * i + 2].value = pt.pos[2];
	// 	x_space[4 * i + 3].index = -1;
	// 	prob.x[i] = &x_space[4 * i];
	// 	prob.y[i] = (pt.isblue)?1:-1;
	// 	i++;
	// }

	const char * error_msg = svm_check_parameter(&prob, &param);

	if (error_msg)
	{
		cerr << "! libsvm ERROR: " << error_msg << endl;
		exit(1);
	}

	//train
	svm_model *model = svm_train(&prob, &param);

	//retreive the separating plane.
	vector<Point3d> surface_points;
	Point3d pt;

	//
	int sv_size = svm_get_nr_sv(model);
	int * sv_indices = new int[sv_size]; //Each sv_indices[i] is in the range of [1, ..., num_traning_data].
	svm_get_sv_indices(model, sv_indices);

	vector<Point3d> hullA_sv, hullB_sv;
	//int hullA_pt_size = hullA_pts.size();

	//cout << "sv_size=" << sv_size << endl;
	//cout << "hullA_pt_size=" << hullA_pt_size << endl;
	//cout << "hullB_pt_size=" << hullB_pts.size() << endl;

	//retrive th support vectors
	svm_node x[4];
	x[0].index = 1;
	x[1].index = 2;
	x[2].index = 3;
	x[3].index = -1;
	for (int i = 0; i < sv_size; i++)
	{
		int id = sv_indices[i] - 1;
		auto & s = pts[id]; //(id < hullA_pt_size) ? hullA_pts[id] : hullB_pts[id - hullA_pt_size];
		x[0].value = s.pos[0];
		x[1].value = s.pos[1];
		x[2].value = s.pos[2];

		if (svm_predict(model, x)>0) hullA_sv.push_back(s.pos); //this sv is on the positive side
		else hullB_sv.push_back(s.pos); //otherwise
	}


	//cout << "hullA_sv size=" << hullA_sv.size() << endl;
	//cout << "hullB_sv size=" << hullB_sv.size() << endl;

  if(hullB_sv.empty() || hullA_sv.empty()) return false;

	all_sv_A.insert(all_sv_A.end(), hullA_sv.begin(), hullA_sv.end());
	all_sv_B.insert(all_sv_B.end(), hullB_sv.begin(), hullB_sv.end());;

	//A has lable=+1 and B has label=-1
	//find the closest pairs and then do bisection
	//get k-closest pairs
	const int max_sv_size=100;

	if(hullA_sv.size() > max_sv_size) //reduce hullA size if it is too big
	{
		//get a subset of hullA_sv
		random_shuffle(hullA_sv.begin(),hullA_sv.end());
		vector<Point3d> tmp(hullA_sv.begin(),hullA_sv.begin()+max_sv_size);
		hullA_sv.swap(tmp);
	}

	if(hullB_sv.size() > max_sv_size) //reduce hullB size if it is too big
	{
		//get a subset of hullA_sv
		random_shuffle(hullB_sv.begin(),hullB_sv.end());
		vector<Point3d> tmp(hullB_sv.begin(),hullB_sv.begin()+max_sv_size);
		hullB_sv.swap(tmp);
	}

	//find the cut between all pairs
	for (auto& pa : hullA_sv)
	{
		for (auto& pb : hullB_sv)
		{
			Point3d decision_pt = bisection(model, pa, pb, 1e-10);
			surface_points.push_back(decision_pt);
		}//for pb
	}//end pa


	//cout << "surface_points size=" << surface_points.size() << endl;
	all_sur_pts.insert(all_sur_pts.end(), surface_points.begin(), surface_points.end());

	least_square_fit(surface_points, trim_plane.o, trim_plane.n);

	if ((hullA_sv.front() - trim_plane.o)*trim_plane.n < 0)
	{
		trim_plane.n = -trim_plane.n;
	}

	//this is more for drawing purpose, can be delete later
	Point3d proj_com;
	for (auto& pt : surface_points)
	{
		auto vec = (pt - trim_plane.o);
		auto proj = trim_plane.o + (vec - trim_plane.n*(vec*trim_plane.n));
		for(int j=0;j<3;j++) proj_com[j] += proj[j];
	}
	for (int j = 0; j < 3; j++) proj_com[j] /= surface_points.size();

	double raidus = 0;
	for (auto& pt : surface_points)
	{
		auto vec = (pt - trim_plane.o);
		auto proj = trim_plane.o + (vec - trim_plane.n*(vec*trim_plane.n));
		auto d = (proj - proj_com).normsqr();
		if (d>raidus) raidus = d;
	}
	raidus = sqrt(raidus);

	trim_plane.o = proj_com;
	trim_plane.error = raidus;

	//

	delete[] sv_indices;
	delete[] prob.x;
	delete[] prob.y;
	delete[] x_space;

	return true;
}



bool trim_overlapping_hulls_svm_opt(const cd_hull& hullA, const cd_hull& hullB, hull_plane& trim_plane, double svm_C)
{
   cout<<"- calling svm-opt"<<endl;
   vector<Point3d> more_A;

   vector<RBPoint> pts;
   build_red_blue(hullA, hullB, pts);

  //  getIntersectionPoints(hullA, hullB, more_A);
  //  for(auto& pt: more_A)
  //  {
  //    RBPoint pt2(pt,false);
  //    pt2.inIntersection=true;
  //    pts.push_back(pt2);
  //  }

   bool r=trim_overlapping_hulls_svm(pts,trim_plane,svm_C);
   if(!r) //if failed, we try again without the new points
   {
     cerr<<"! Error: trim_overlapping_hulls_svm failed"<<endl;
     return false;
    //  pts.clear();
    //  build_red_blue(hullA, hullB, pts);
    //  r=trim_overlapping_hulls_svm(pts,trim_plane,svm_C);
    //  if(!r){
    //    cerr<<"! Error: trim_overlapping_hulls_svm failed"<<endl;
    //    return false;
    //  }
     //
    //  for(auto& pt: more_A) //do this for optimization.
    //  {
    //    RBPoint pt2(pt,false);
    //    pt2.inIntersection=true;
    //    pts.push_back(pt2);
    //  }
   }

   cout<<"- optimize orientation"<<endl;
   //optimize_orientation(pts, trim_plane.o, trim_plane.n);

   cout<<"- optimize position"<<endl;
   auto o_bkup=trim_plane.o;
	 //optimize_position(pts, trim_plane.o, trim_plane.n);
   optimize_position_exact(hullA, hullB, trim_plane.o, trim_plane.n);

   //reproject o_bkup to the new plane, for visualization purpose
   auto v=(o_bkup-trim_plane.o);
   trim_plane.o=trim_plane.o+(v-trim_plane.n*(trim_plane.n*v));

	 return r;
}

//
// trim overlapping hulls using svm
// if the hulls do not overlap return false
//
bool trim_overlapping_hulls_svm(const cd_hull& hullA, const cd_hull& hullB, hull_plane& trim_plane, double svm_C)
{
	//Samlpe hullA and hullB based on their volumes

	const vector<Point3d> & hullA_pts = hullA.getSampledPoints();
	const vector<Point3d> & hullB_pts = hullB.getSampledPoints();

	cout << "hullA_pts size=" << hullA_pts.size() << endl;
	cout << "hullB_pts size=" << hullB_pts.size() << endl;

	//use libsvm to find the separating plane

	struct svm_parameter param;

	// default values
	param.svm_type = C_SVC; //C_SVC, NU_SVC
	param.kernel_type = LINEAR;
	param.cache_size = 100;
	param.C = svm_C;
	param.eps = 1e-3;
	param.shrinking = 1;
	param.probability = 0;
	param.nr_weight = 0;
	param.weight_label = NULL;
	param.weight = NULL;
	param.gamma=0.6;
	param.degree = 3;

	//prepare the svm problem
	svm_problem prob;
	prob.l = hullA_pts.size() + hullB_pts.size();
	svm_node *x_space = new svm_node[4 * prob.l];
	assert(x_space);
	prob.x = new svm_node *[prob.l];
	prob.y = new double[prob.l];

	int i = 0;
	for (auto& pt : hullA_pts)
	{
		x_space[4 * i].index = 1;
		x_space[4 * i].value = pt[0];
		x_space[4 * i + 1].index = 2;
		x_space[4 * i + 1].value = pt[1];
		x_space[4 * i + 2].index = 3;
		x_space[4 * i + 2].value = pt[2];
		x_space[4 * i + 3].index = -1;
		prob.x[i] = &x_space[4 * i];
		prob.y[i] = +1;
		i++;
	}

	for (auto& pt : hullB_pts)
	{
		x_space[4 * i].index = 1;
		x_space[4 * i].value = pt[0];
		x_space[4 * i + 1].index = 2;
		x_space[4 * i + 1].value = pt[1];
		x_space[4 * i + 2].index = 3;
		x_space[4 * i + 2].value = pt[2];
		x_space[4 * i + 3].index = -1;
		prob.x[i] = &x_space[4 * i];
		prob.y[i] = -1;
		i++;
	}

	const char * error_msg = svm_check_parameter(&prob, &param);

	if (error_msg)
	{
		cerr << "! libsvm ERROR: " << error_msg << endl;
		exit(1);
	}

	//train
	svm_model *model = svm_train(&prob, &param);

	//retreive the separating plane.
	vector<Point3d> surface_points;
	Point3d pt;

	//
	int sv_size = svm_get_nr_sv(model);
	int * sv_indices = new int[sv_size]; //Each sv_indices[i] is in the range of [1, ..., num_traning_data].
	svm_get_sv_indices(model, sv_indices);

	vector<Point3d> hullA_sv, hullB_sv;
	int hullA_pt_size = hullA_pts.size();

	cout << "sv_size=" << sv_size << endl;
	cout << "hullA_pt_size=" << hullA_pt_size << endl;
	cout << "hullB_pt_size=" << hullB_pts.size() << endl;

	//retrive th support vectors
	svm_node x[4];
	x[0].index = 1;
	x[1].index = 2;
	x[2].index = 3;
	x[3].index = -1;
	for (int i = 0; i < sv_size; i++)
	{
		int id = sv_indices[i] - 1;
		auto & s = (id < hullA_pt_size) ? hullA_pts[id] : hullB_pts[id - hullA_pt_size];
		x[0].value = s[0];
		x[1].value = s[1];
		x[2].value = s[2];

		if (svm_predict(model, x)>0) hullA_sv.push_back(s); //this sv is on the positive side
		else hullB_sv.push_back(s); //otherwise
	}


	cout << "hullA_sv size=" << hullA_sv.size() << endl;
	cout << "hullB_sv size=" << hullB_sv.size() << endl;

  if(hullB_sv.empty() || hullA_sv.empty()) return false;

	all_sv_A.insert(all_sv_A.end(), hullA_sv.begin(), hullA_sv.end());
	all_sv_B.insert(all_sv_B.end(), hullB_sv.begin(), hullB_sv.end());;

	//A has lable=+1 and B has label=-1
	//find the closest pairs and then do bisection
	//get k-closest pairs
	const int max_sv_size=100;

	if(hullA_sv.size() > max_sv_size) //reduce hullA size if it is too big
	{
		//get a subset of hullA_sv
		random_shuffle(hullA_sv.begin(),hullA_sv.end());
		vector<Point3d> tmp(hullA_sv.begin(),hullA_sv.begin()+max_sv_size);
		hullA_sv.swap(tmp);
	}

	if(hullB_sv.size() > max_sv_size) //reduce hullB size if it is too big
	{
		//get a subset of hullA_sv
		random_shuffle(hullB_sv.begin(),hullB_sv.end());
		vector<Point3d> tmp(hullB_sv.begin(),hullB_sv.begin()+max_sv_size);
		hullB_sv.swap(tmp);
	}

	//find the cut between all pairs
	for (auto& pa : hullA_sv)
	{
		for (auto& pb : hullB_sv)
		{
			Point3d decision_pt = bisection(model, pa, pb, 1e-10);
			surface_points.push_back(decision_pt);
		}//for pb
	}//end pa


	//cout << "surface_points size=" << surface_points.size() << endl;
	all_sur_pts.insert(all_sur_pts.end(), surface_points.begin(), surface_points.end());

	least_square_fit(surface_points, trim_plane.o, trim_plane.n);

	if ((hullA_sv.front() - trim_plane.o)*trim_plane.n < 0)
	{
		trim_plane.n = -trim_plane.n;
	}

	//this is more for drawing purpose, can be delete later
	Point3d proj_com;
	for (auto& pt : surface_points)
	{
		auto vec = (pt - trim_plane.o);
		auto proj = trim_plane.o + (vec - trim_plane.n*(vec*trim_plane.n));
		for(int j=0;j<3;j++) proj_com[j] += proj[j];
	}
	for (int j = 0; j < 3; j++) proj_com[j] /= surface_points.size();

	double raidus = 0;
	for (auto& pt : surface_points)
	{
		auto vec = (pt - trim_plane.o);
		auto proj = trim_plane.o + (vec - trim_plane.n*(vec*trim_plane.n));
		auto d = (proj - proj_com).normsqr();
		if (d>raidus) raidus = d;
	}
	raidus = sqrt(raidus);

	trim_plane.o = proj_com;
	trim_plane.error = raidus;

	//

	delete[] sv_indices;
	delete[] prob.x;
	delete[] prob.y;
	delete[] x_space;

	return true;
}

//
// trim overlapping hulls using simple heuristic
// if the hulls do not overlap return false
//
bool trim_overlapping_hulls_heuristic(const cd_hull& hullA, const cd_hull& hullB, hull_plane& trim_plane)
{
	// //find all pairs of intersecting triangles between A & B
  vector<Point3d> intersections;
  getIntersectionPoints(hullA, hullB, intersections, 0.1);

  cd_m * hullA_m = hullA.to_m();
  cd_m * hullB_m = hullB.to_m();

	// list< pair<int,int> > contacts;
	// {
	//   CD_RAPID cd;
	// 	int Acd = cd.buildModel(*hullA_m);
	//   int Bcd = cd.buildModel(*hullB_m);
	// 	cd.isInCollision(Acd,Bcd,contacts);
  // }
  //
	// //find the actual intersection between all pairs of triangles
	// vector<Point3d> intersections;
	// {
	// 	vector<cd_f *> Af(hullA_m->f().begin(),hullA_m->f().end());
	// 	vector<cd_f *> Bf(hullB_m->f().begin(),hullB_m->f().end());
  //
	// 	for(auto& contact: contacts)
	// 	{
	// 		auto & ta = Af[contact.first];
	// 		auto & tb = Bf[contact.second];
  //
	// 		//find intersection between ta & tb
	// 		tri_tri_intersection(ta, tb, intersections);
	// 	}
  // }

	least_square_fit(intersections, trim_plane.o, trim_plane.n); //implemented below

	all_sur_pts.insert(all_sur_pts.end(), intersections.begin(), intersections.end());

	//make sure that n is pointing outside B
  double max_A_d=-FLT_MAX;
	Point3d max_A_pt;

	for(auto& v : hullA_m->v())
	{
		if(hullB.isInsideHull(v->pos)==false) //found the furthese vertex of A outside B
		{
			double d = fabs((v->pos - trim_plane.o)*trim_plane.n);
			if(d>max_A_d){ max_A_d=d; max_A_pt=v->pos;}
    }
	}

  //---------------------------------------

  double max_B_d=-FLT_MAX;
	Point3d max_B_pt;
  for(auto& v : hullB_m->v())
  {
    if(hullA.isInsideHull(v->pos)==false) //found the furthese vertex of B outside A
    {
      double d = fabs((v->pos - trim_plane.o)*trim_plane.n);
      if(d>max_B_d){ max_B_d=d; max_B_pt=v->pos;}
    }
  }

	if(max_A_d==-FLT_MAX && max_B_d==-FLT_MAX)
	{
		delete hullA_m;
    delete hullB_m;

	  cerr<<"! Error: Cannot find any point of A outside B or point of B outside A!!!"<<endl;
		return false;
	}

  if(max_A_d!=-FLT_MAX && max_B_d!=-FLT_MAX)
  {
    auto da=(max_A_pt - trim_plane.o)*trim_plane.n;
    auto db=(max_B_pt - trim_plane.o)*trim_plane.n;
    if(da*db<0) //different signs, easy case
    {
      if (da < 0)
      {
        trim_plane.n = -trim_plane.n;
      }
    }
    else //same sign??!!
    {
      if(da<0 && max_A_d>=max_B_d)
      {
        trim_plane.n = -trim_plane.n;
      }
      else if(db>0 && max_B_d>=max_A_d)
      {
        trim_plane.n = -trim_plane.n;
      }
    }
    ///------------------------------------------
  }
  else if (max_A_d!=-FLT_MAX)
  {
    auto da=(max_A_pt - trim_plane.o)*trim_plane.n;
    if (da < 0)
    {
      trim_plane.n = -trim_plane.n;
    }
  }
  else if (max_B_d!=-FLT_MAX)
  {
    auto d=(max_B_pt - trim_plane.o)*trim_plane.n;
    if (d > 0)
    {
      trim_plane.n = -trim_plane.n;
    }
  }

	//this is more for drawing purpose, can be delete later
	Point3d proj_com;
	for (auto& pt : intersections)
	{
		auto vec = (pt - trim_plane.o);
		auto proj = trim_plane.o + (vec - trim_plane.n*(vec*trim_plane.n));
		for(int j=0;j<3;j++) proj_com[j] += proj[j];
	}
	for (int j = 0; j < 3; j++) proj_com[j] /= intersections.size();

	double raidus = 0;
	for (auto& pt : intersections)
	{
		auto vec = (pt - trim_plane.o);
		auto proj = trim_plane.o + (vec - trim_plane.n*(vec*trim_plane.n));
		auto d = (proj - proj_com).normsqr();
		if (d>raidus) raidus = d;
	}
	raidus = sqrt(raidus);

	trim_plane.o = proj_com;
	trim_plane.error = raidus;


	delete hullA_m;
  delete hullB_m;

	return true;

}

//repair model
void repair_trimed_hull(cd_m * convex, cd_hull& hull, double vertex_collapse_tolerance);

//
// given a hull and a list of trim planes, return a trimed hull
//
//
bool trim_hull_by_planes(const cd_hull& hull, const list<hull_plane>& trim_planes, cd_hull& trimed_hull, const Point3d& o)
{
  //convert all facet of the hull into half-space.
	std::list<VAS_plane> planes;
	VAS_plane vsp;
	for (auto& plane : trim_planes)
	{
		vsp.n = plane.n;
		vsp.p = plane.o;

		if ( !isfinite(vsp.n[0]) ||  !isfinite(vsp.n[1]) ||  !isfinite(vsp.n[2])
		     || !isfinite(vsp.p[0]) ||  !isfinite(vsp.p[1]) ||  !isfinite(vsp.p[2]) )
		{
			cout << "vsp.n=(" << vsp.n << ") vsp.p=(" << vsp.p << ")" << endl;
			cout << "something is not right, ignore this plane"<<endl;;
			continue;
		}

		planes.push_back(vsp);
	}

	//cout<<"Found "<<planes.size()<<" cutting planes"<<endl;

	cd_m * m=hull.to_m();
	for (auto& f : m->f())
	{
    if(f->area()<1e-10) continue; //too small....
		vsp.n = -f->n;
		vsp.p = f->v[0]->pos;
		planes.push_back(vsp);
	}

	//cout << "======================== Hspace_intersection start ==========================" << endl;
	//compute the intersection
	cd_m * trimed_hull_m = Hspace_intersection(planes, o);
	//cout << "======================== Hspace_intersection end ==========================" << endl;

  cd_hull newhull;
  newhull.buildhull(trimed_hull_m);
  trimed_hull=newhull;

  trimed_hull_m->destroy();
  delete trimed_hull_m;

  m->destroy();
  delete m;

  return true;
}

bool trim_hull_by_planes(const cd_hull& hull, const list<hull_plane>& trim_planes, cd_hull& trimed_hull, double vertex_collapse_tolerance)
{

  //cout<<"trim_hull_by_planes start"<<endl;

	//get the point inside the hull
	auto o=getPointInIntersection(hull, trim_planes);
  if(o[0]==FLT_MAX)
  {
    cerr<<"! Error: trim_hull_by_planes cannot find point inside hull."<<endl;
    return false;
  }

	//convert all facet of the hull into half-space.
	std::list<VAS_plane> planes;
	VAS_plane vsp;
	for (auto& plane : trim_planes)
	{
		vsp.n = plane.n;
		vsp.p = plane.o;

		if ( !isfinite(vsp.n[0]) ||  !isfinite(vsp.n[1]) ||  !isfinite(vsp.n[2])
		     || !isfinite(vsp.p[0]) ||  !isfinite(vsp.p[1]) ||  !isfinite(vsp.p[2]) )
		{
			cout << "vsp.n=(" << vsp.n << ") vsp.p=(" << vsp.p << ")" << endl;

			cout << "something is not right, ignore this plane"<<endl;;
			continue;
		}

		planes.push_back(vsp);
	}

	//cout<<"Found "<<planes.size()<<" cutting planes"<<endl;

	cd_m * m=hull.to_m();
	for (auto& f : m->f())
	{
    if(f->area()<1e-10) continue; //too small....
		vsp.n = -f->n;
		vsp.p = f->v[0]->pos;
		planes.push_back(vsp);
	}

	//cout << "======================== Hspace_intersection start ==========================" << endl;
	//compute the intersection
	cd_m * trimed_hull_m = Hspace_intersection(planes, o);
	//cout << "======================== Hspace_intersection end ==========================" << endl;

  if(vertex_collapse_tolerance==0)
  {
    cd_hull newhull;
    newhull.buildhull(trimed_hull_m);
    trimed_hull=newhull;
  }
  else
  {
  	//repair model (remove 0 area triangles)
  	repair_trimed_hull(trimed_hull_m, trimed_hull, vertex_collapse_tolerance);

  	while (trimed_hull_m->f().size() - trimed_hull.getFaceSize() != 0 && trimed_hull.getFaceSize()!=0)
  	{
  		//destroy the model
  		trimed_hull_m->destroy();
  		delete trimed_hull_m;

  		//
  		trimed_hull_m = trimed_hull.to_m();
  		repair_trimed_hull(trimed_hull_m, trimed_hull, vertex_collapse_tolerance);
  	}
  }//end if

	//destroy the model
	trimed_hull_m->destroy();
	delete trimed_hull_m;

	m->destroy();
	delete m;

	//cout<<"trim_hull_by_planes done"<<endl;

	return true;
}

double least_square_fit(const vector<Point3d>& pts, Point3d& o, Vector3d& n)
{
	//compute X2, Y2, XY, XZ, YZ, X, Y, Z from pts
	double X2 = 0, Y2 = 0, XY = 0, XZ = 0, YZ = 0, X = 0, Y = 0, Z = 0;

	for (auto& pt : pts)
	{
    //cout<<"pt="<<pt<<endl;
		X2 += (pt[0] * pt[0]);
		Y2 += (pt[1] * pt[1]);
		XY += (pt[0] * pt[1]);
		XZ += (pt[0] * pt[2]);
		YZ += (pt[1] * pt[2]);
		X += (pt[0]);
		Y += (pt[1]);
		Z += (pt[2]);
	}

	Eigen::Matrix3d m;
	Eigen::Vector3d b;
	m << X2, XY, X, XY, Y2, Y, X, Y, pts.size();
	b << XZ, YZ, Z;

	Eigen::Vector3d x = m.fullPivHouseholderQr().solve(b); //m.colPivHouseholderQr().solve(b);


	//cout << "solution = \n"<< x << endl;
	//cout << "residual = \n" << (m*x - b) << endl;
	n.set(x(0),x(1),-1);
	n = n.normalize();
	o.set(0,0,x(2));

	//compute the error
	double max_e = -FLT_MAX;
	for (auto& pt : pts)
	{
		double e = x(0)*pt[0] + x(1)*pt[1] + x(2) - pt[2];
		e=e*e;
    if(e>max_e) max_e=e;
	}

	//cout << "max error=" << max_e << endl;
	return  max_e;
}


inline bool isInside(const list<hull_plane>& trim_planes, const Point3d& pt)
{
	for (auto& plane : trim_planes)
	{
		double  dot = (pt - plane.o)*plane.n;
		if (dot <= 0) return false;
	}
	return true;
}

Point3d getPointInIntersection(const cd_hull& hull, const list<hull_plane>& trim_planes)
{
	vector<Point3d> hull_pts;
	hull.to_pvector(hull_pts);
	auto com=computeCenter(hull_pts);

	if (isInside(trim_planes, com)) return com;

  const int max_iterations = 1000000;
	int iter=0;
	while (iter++ < max_iterations)
	{
		//draw samples
		Point3d sample;
		double sum_w = 0;
		for (auto& p : hull_pts)
		{
			auto w = drand48();//+1e-10;
			sum_w += w;
			for (int j = 0; j<3; j++) sample[j] += w*p[j];
		}

		//cout << "sum_w=" << sum_w << endl;

		for (int j = 0; j<3; j++) sample[j] /= sum_w;
		if (isInside(trim_planes, sample)) return sample;
	}

//  for(auto& pt : hull_pts)
//		if (isInside(trim_planes, pt)) return pt;

  return Point3d(FLT_MAX,FLT_MAX,FLT_MAX);
}

//check if c is between a and b
inline bool inbetween(const Point3d& a, const Point3d& b, const Point3d& c)
{
	auto v = (b - a);
	auto u = (c - a);
	auto dot = u*v;
	if (dot<0) return false;
	if (dot>v.normsqr()) return false;
	return true;
}

//repair model
void repair_trimed_hull(cd_m * convex, cd_hull& hull, double vertex_collapse_tolerance)
{
	//init
	const float small = vertex_collapse_tolerance;
	map<cd_v*, bool> v2dead;
	for (auto v : convex->v()) v2dead[v] = false;

	//check all edges
	for (auto e : convex->e())
	{
		if (v2dead[e->v[0]]) continue;
		if (v2dead[e->v[1]]) continue;
		if (e->length < small){ v2dead[e->v[0]] = true;}
	}

	//check all triangles
	for (auto f : convex->f())
	{
		auto & a = f->v[0];
		auto & b = f->v[1];
		auto & c = f->v[2];

		if (v2dead[a]) continue;
		if (v2dead[b]) continue;
		if (v2dead[c]) continue;

		if (f->area() > small) continue;

		if (inbetween(a->pos, b->pos, c->pos)) v2dead[c] = true;
		else if (inbetween(b->pos, c->pos, a->pos)) v2dead[a] = true;
		else if (inbetween(c->pos, a->pos, b->pos)) v2dead[b] = true;
	}

	//collect points that are not deleted
	vector<Point3d> pts;
	for (auto v : v2dead)
	{
		if (v.second == false) pts.push_back(v.first->pos);
	}

	//check if a new model should be build
	//if (pts.size() == convex->v().size()) return; //nothing is changed....

	//cout << "pts.size()=" << pts.size() << " convex->v().size()=" << convex->v().size() << endl;

	//build the new model using pts
	cd_hull newhull;

	if (pts.size()>=4) newhull.buildhull(pts);

	//cout << "old hull face size=" << convex->f().size() << endl;

	//done
	hull=newhull;

	//
	//cout << "new hull face size=" << hull.getFaceSize() << endl;
}


///
/// Create segments between two intersecting triangles
///

void tri_tri_intersection(cd_f* f1, cd_f* f2, vector<Point3d>& intersections, float resolution)
{
    if(f1->area()<1e-10 || f2->area()<1e-10) return; //nothing to see here

    //init
    Point3d u[3]; //vertices of f1
    Point3d v[3]; //vertices of f2

    for(int j=0;j<3;j++) u[j]=f1->v[j]->pos;
    for(int j=0;j<3;j++) v[j]=f2->v[j]->pos;

    double ndot=f1->n*f2->n;
    if(fabs(fabs(ndot)-1)<1e-20){ //parallel

       if(ndot<0) //facing opposite direction
			 {
         cerr<<"! Warning: build_segment degenerate case: f1 and f2 are parallel and facing opposite directions"<<endl;
				 for(int j=0;j<3;j++) intersections.push_back(u[j]);
		     for(int j=0;j<3;j++) intersections.push_back(v[j]);
         return;
			 }
			 else{
			   //ignore this case for now...
				 cerr<<"! Warning: build_segment degenerate case: f1 and f2 are parallel and facing the same direction"<<endl;
			 }

    }

    //not parallel
    //compute intersection
    Point3d x3d[2];
    bool   ex[2][6]; //record which edge contributes to which vertex
    char r=tri_tri_intersection(u,f1->n,v,f2->n,x3d,ex[0],ex[1]);
    if(r=='0'){ //odd, no intersection found
        return;
    }

    //interpolate...
    //if( fabs((x3d[0]-u[0])*f1->n)<1e-10 && fabs((x3d[0]-v[0])*f2->n)<1e-10 )
    intersections.push_back(x3d[0]);
    //if( fabs((x3d[1]-u[0])*f1->n)<1e-10 && fabs((x3d[1]-v[0])*f2->n)<1e-10 )
    intersections.push_back(x3d[1]);

    auto vec=(x3d[1]-x3d[0]);
    double vec_norm=vec.norm();
    if(resolution<vec_norm)
    {
      int steps=(int)ceil(vec_norm/resolution);
      vec=vec/steps;
      for(int i=1;i<steps;i++) intersections.push_back(x3d[0]+vec*i);
    }
}

Point3d bisection(svm_model * model, const Point3d& p, const Point3d& q, double tau_sqr)
{
	Point3d s=p;
	Point3d t=q;

	svm_node x[4];
	x[0].index = 1;
	x[1].index = 2;
	x[2].index = 3;
	x[3].index = -1;

	double d_sqr = (s - t).normsqr();
	while (d_sqr > tau_sqr)
	{
		Point3d mid((s[0] + t[0]) / 2, (s[1] + t[1]) / 2, (s[2] + t[2]) / 2);
		x[0].value = mid[0];
		x[1].value = mid[1];
		x[2].value = mid[2];
		double d = svm_predict(model,x);

		if (d > 0) s = mid;
		else t = mid;

		d_sqr = (s - t).normsqr();
	}

	//now s and t are close to each other, return the mid point
	return Point3d((s[0] + t[0]) / 2, (s[1] + t[1]) / 2, (s[2] + t[2]) / 2);
}

pair<double,double> optimize_position_exact(const vector<RBPoint>& pts, hull_plane& hp)
{
  double red_c2=0, red_c1=0,red_c0=0;
  double blue_c2=0, blue_c1=0, blue_c0=0;

  //
  double lowest_d  = pts.back().dist;

  //cout<<"hp.o="<<hp.o<<" hp.n="<<hp.n<<endl;

  //accumulate coefficients
  for(auto& pt:pts)
  {
    if(pt.isblue)
    {
      compute_cone_volume_derivative_coeeficient(pt.v, hp, blue_c2, blue_c1, blue_c0);
    }
    else
    {
      compute_cone_volume_derivative_coeeficient(pt.v, hp, red_c2, red_c1, red_c0);
    }
  }

  //solve quadratic equation...
  double a=(blue_c2-red_c2);
  double b=(blue_c1-red_c1);
  double c=(blue_c0-red_c0);

  //cout<<"a="<<a<<" b="<<b<<" c="<<c<<endl;

  double d=0;
  if(!isfinite(a) || !isfinite(b) || !isfinite(c)) return make_pair(FLT_MAX,FLT_MAX);
  if(a==0 && b==0) return make_pair(FLT_MAX,FLT_MAX);
  if(a==0) //LINEAR
  {
    d=-c/b;
  }
  else{
    double tmp=b*b-4*a*c;
    if(tmp<1e-10) return make_pair(FLT_MAX,FLT_MAX); //1e-10 for stability...
    tmp=sqrt(tmp);
    double x1=(-b+tmp)/(2*a);
    double x2=(-b-tmp)/(2*a);
    //cout<<"x1="<<x1<<" x2="<<x2<<endl;

    //smaller x means more in the direction of n
    //d=min(x1,x2);
    return make_pair(x1,x2);
  }

  // if(d>=lowest_d)
  //   return make_pair(FLT_MAX,FLT_MAX);
  //
  // hp.o=Point3d(0,0,0)+hp.n*(-d);
  //
  // return d;
    cout<<"yh-remove later!!!"<<endl;
    return make_pair(FLT_MAX,FLT_MAX);    //yh-error
}

//given two overlapping convex object, find the best position o so that
//(the volume of B on positive side of the plane)+(the volume of A on on the negtive side of the plane)
//is miminized.
void optimize_position_exact(const cd_hull& hullA, const cd_hull& hullB, Point3d& o, const Vector3d& n)
{
  // if(debug)
  // {
  //   cd_hull tmp = hullA;
  //   cout<<"A vol="<<tmp.computeHullVol()<<endl;
  //   cd_hull tmp2 = hullB;
  //   cout<<"B vol="<<tmp2.computeHullVol()<<endl;
  // }

  //Point3d my_o(21.264556840247553993, -22.637143479825031278, 28.004212727551966111);
  //Vector3d my_n(-0.8862423356744910663, -0.21464119924175334431, 0.41049199510621975362);
  //bool debug = ((my_o-o).norm()<1e-10) && ((my_n-n).norm()<1e-10);
  //cout.precision(20);
  //cout<<"o="<<o<<" n="<<n<<endl;

  cd_m * A=hullA.to_m();
  cd_m * B=hullB.to_m();
  vector<RBPoint> pts;
  build_red_blue(A,B,pts);

	for(auto& pt: pts)
	{
    pt.dist=-(pt.pos[0]*n[0]+pt.pos[1]*n[1]+pt.pos[2]*n[2]);
    //pt.dist=-((pt.pos-o)*n);
  }

  //most red are at the front and blue are at the end
  sort(pts.begin(),pts.end());
  double lowest_d=pts.back().dist, highest_d=pts.front().dist;

  // {
  //   RBPoint last_blue=pts.back();
  //   while(pts.back().isblue){
  //     last_blue=pts.back();
  //     pts.pop_back();
  //   }
  //   if(last_blue.v != pts.back().v) pts.push_back(last_blue);
  // }


  //
  // {
  //   //find first blue (fb) and last red (lr)
  //   int fb_id=0, lr_id=0;
  //   uint size=pts.size();
  //   for(lr_id=size-1;lr_id>=0;lr_id--) if(!pts[lr_id].isblue) break;
  //   for(;fb_id<size;fb_id++) if(pts[fb_id].isblue) break;
  //   if(fb_id>0) fb_id--;
  //   if(lr_id<size-1) lr_id+=1;
  //   lowest_d=pts[lr_id].dist;
  //   highest_d=pts[fb_id].dist;
  //
  //   for(;size>lr_id+1; size--) pts.pop_back();
  //   //vector<RBPoint> tmp(pts.begin()+fb_id, pts.begin()+lr_id);
  //   //tmp.swap(pts);
  // }

  double min_vol=FLT_MAX;

  //cout<<"------"<<endl;
  while(pts.empty()==false)
  {
    //if(pts.back().dist<highest_d) break;

    //cout<<"~~~~~~"<<endl;
    //cout<<"lowest d="<<lowest_d<<" highest d="<<highest_d<<endl;
    hull_plane hp;
    //hp.o=o; //BUG: this should not be o...
    hp.n=n;

    auto d_pair=optimize_position_exact(pts, hp);
    //if(d!=FLT_MAX && d>=highest_d)
    double d[2]={d_pair.first,d_pair.second};

//exit(1);

    for(int i=0;i<2;i++)
    {
      if(d[i]!=FLT_MAX && d[i]<=lowest_d && d[i]>=highest_d)
      {
        hp.o=Point3d(0,0,0)+hp.n*(-d[i]);
        double vol_b=compute_volume(hullB, hp);
        double vol_a=compute_volume(hullA, hp);



        double vol=vol_b-vol_a;
        // if(debug)
        // {
        //   cout<<"d="<<d[i]<<" vol="<<vol<<" vol_a="<<vol_a<<" vol_b="<<vol_b<<endl;
        // }

        if(vol_b<0 || vol_a<0) //this should never happend, ignore it.
        {
          cout<<"hp.o="<<hp.o<<" hp.n="<<hp.n<<endl;
          cout<<"! error: something wrong"<<endl;
          //compute_volume(hullB, hp, true);
          continue;
        }

        if(vol<min_vol)
        {
          //cout<<"!!! new min"<<endl;
          //compute_volume(hullA, hp, true);
          min_vol=vol;
          o=hp.o;
        }
      }
    }//end for i


    pts.pop_back();
  }

  A->destroy();
  B->destroy();
  delete A;
  delete B;

  //if(debug) exit(1);
}

//given red/blue point sets, find the best position o so that
//the sum of # of blue points on the positive side of the plane and # of red points on the negtive side of the plane
//is miminized.
void optimize_position(vector<RBPoint>& pts, Point3d& o, const Vector3d& n)
{
  uint min_count=0;
	for(auto& pt:pts)
	{
		pt.dist=(pt.pos-o)*n;
    if(pt.isblue) min_count++;
	}

	sort(pts.begin(),pts.end());

	// start to count...
	if(!pts.front().isblue){
		reverse(pts.begin(), pts.end());
	}

	uint best_id=UINT_MAX;
  uint blue_red_diff=min_count;
  int blue_count=min_count;
  int red_count=0;

	uint count=min_count;
	uint psize=pts.size();

	for(uint i=0;i<psize;i++)
	{
		auto& pt=pts[i];
		if(pt.isblue){
      count--;
      blue_count--;
    }
		else{
      count++;
      red_count++;
    }

    if(pt.inIntersection) //to update the position, pt must be in interesction
    {
  		if(count<min_count)
  		{
  			min_count=count;
  			best_id=i;
        blue_red_diff=abs(red_count-blue_count);
  		}
      else if(count==min_count)
      {
        uint diff=abs(red_count-blue_count);
        if(diff<blue_red_diff)
        {
          blue_red_diff=diff;
          best_id=i;
        }
      }
    }//(pt.isIntersection)

	}//end for i

//cout<<"min_count="<<min_count<<endl;

  if(best_id==psize-1)
	{
		 o=pts[best_id].pos;
	}
	else if(best_id<psize)
	{
    if(pts[best_id+1].inIntersection)
    {
  		const Point3d & s=pts[best_id].pos;
  		const Point3d & t=pts[best_id+1].pos;
  		o.set( (s[0]+t[0])/2, (s[1]+t[1])/2, (s[2]+t[2])/2);
    }
    else
    {
      o=pts[best_id].pos;
    }
	}
}

//given red/blue point sets, and an extra point to constrain the orientation,
//find the best orientation n so that
//the sum of # of blue points on the positive side of the plane and # of red
//points on the negtive side of the plane
//is miminized.

//return the sum
inline uint optimize_orientation
(vector<RBPoint>& pts,
 const Point3d& o,  Vector3d& n,
 const Point3d& p)
{
  Vector3d v=(p-o);
  {
    float vnorm=v.norm();
    if(vnorm==0) return UINT_MAX;
    v=v/vnorm;
  }

  //update n so that p is included on the plane
  const Vector3d u=n%v;
  Vector3d new_n=v%u;
  Vector3d old_n=n;
  if(new_n*n<0) new_n=-new_n;

  //now, for every point q, compute the angle that this plane needs to rotate
  //around axis v to touch point q
  uint blue_count=0; //# of blue points on the positive side
  uint red_count=0;  //# of red points on the negative side
  uint pts_size=pts.size();

  double max_y=-FLT_MAX;
  bool is_max_Y_red=false;

  //note: the furthest point should be red...
  for(auto& pt: pts)
  {

    //
    double x=(pt.pos-o)*u;
    double y=(pt.pos-o)*new_n;
    double theta=atan2(y,x);
    pt.dist=(theta>0)?theta:(theta+PI);
    //
    if(pt.isblue && y>0) blue_count++;
    else if(!pt.isblue && y<0) red_count++;

    if(y>max_y)
    {
      max_y=y;
      is_max_Y_red=(!pt.isblue);
    }
  }//end for pt

  if(!is_max_Y_red) return UINT_MAX; //not valid

  sort(pts.begin(),pts.end());

  uint min_sum=blue_count+red_count;
  uint cur_sum=min_sum;

  //
  for(auto& pt: pts)
  {
    //cout<<"pt.dist="<<pt.dist<<endl;

    if(AlmostEqual3(pt.pos.get(),p.get())) continue;

    double y=(pt.pos-o)*new_n;
    if(pt.isblue) cur_sum+=(y>0)?-1:1;
    else cur_sum+=(y>0)?1:-1;

    double tmp=min(cur_sum, pts_size-cur_sum);

    if(tmp<min_sum)
    {
      min_sum=tmp;
      //compute the new n
      n=(v%(pt.pos-o)).normalize();
      if(n*old_n<0) n=-n;
    }
  }//end for pt

  return min_sum;
}

inline uint count_errors(vector<RBPoint>& pts, const Point3d& o,  Vector3d& n)
{
  uint count=0;
  for(auto& p: pts)
  {
    double d=(p.pos-o)*n;
    if(p.isblue && d>0) count++;
    else if(!p.isblue && d<0) count++;
  }
  return count;
}

//return the sum
void optimize_orientation_breute_force
(vector<RBPoint>& pts,
 const Point3d& o,  Vector3d& n)
{
  uint min_count=count_errors(pts,o,n);
  Vector3d best_n=n;

  uint psize=pts.size();
  for(uint i=0;i<psize;i++)
  {
    auto& p1=pts[i];
    if(p1.inIntersection==false) continue; //only consider points in intersection for potential second point
    auto v1=p1.pos-o;

    for(uint j=i+1;j<psize;j++)
    {
      auto& p2=pts[j];
      if(p2.inIntersection==false) continue; //only consider points in intersection for potential second point
      if(AlmostEqual3(p1.pos.get(),p2.pos.get())) continue;

      //compute the new n
      auto v2=p2.pos-o;
      Vector3d nn=(v1%v2).normalize();
      if(nn*n<0) nn=-nn;
      uint count=count_errors(pts,o,nn);
      if(count<min_count)
      {
        min_count=count;
        best_n=nn;
      }
    }

  }//end for p

  n=best_n;
}



//given red/blue point sets, find the best orientation n so that
//the sum of # of blue points on the positive side of the plane and # of red
//points on the negtive side of the plane
//is miminized.
void optimize_orientation(vector<RBPoint>& pts, const Point3d& o,  Vector3d& n)
{
  uint min_count=0;
  Vector3d best_n=n;
  for(auto& p: pts)
  {
    double d=(p.pos-o)*n;
    if(p.isblue && d>0) min_count++;
    else if(!p.isblue && d<0) min_count++;
  }

  for(auto& p: pts)
  {
    if(p.inIntersection==false) continue; //only consider points in intersection for potential second point

    Vector3d my_n=n;
    uint count=optimize_orientation(pts, o, my_n, p.pos);
    if(count<min_count)
    {
      min_count=count;
      best_n=my_n;
    }
  }//end for p

  n=best_n;
}

bool expand_hull_by_points(const cd_hull& hull,
                           const list<hull_plane>& trim_planes,
                           vector< vector<Point3d> >& new_pts_sets,
                           cd_hull& exp_hull, double max_vol_inc)
{
  exp_hull=hull;

  //
  //collect points from the hull
  vector<Point3d> hull_pts;
  double init_vol = exp_hull.computeHullVol();
  {
    cd_m * m = exp_hull.to_m();
    for(auto & v : m->v())
    {
      hull_pts.push_back(v->pos);
    }
    m->destroy();
    delete m;
  }

  //
  //this is a subset sum problem
  //select the largest setset of new_pts_sets so that the total is maxmized button
  //smaller than max_vol_inc
  //

  //filter points that are in the trim_planes
  //and also inside a subset of hull faces that can enclose a finite volumn...
  vector< pair<double, uint> > sorted_by_vol_inc;
  uint new_pts_size=new_pts_sets.size();
  for(uint i=0;i<new_pts_size;i++)
  {
    auto& pts=new_pts_sets[i];
    vector<Point3d> tmp;
    for(auto& pt: pts)
    {
      //make sure pt is in the contains
      if( isInside(trim_planes, pt) ) tmp.push_back(pt);
    }//end pt
    pts.swap(tmp);
  }//end i

  for(uint i=0;i<new_pts_size;i++)
  {
    auto& pts=new_pts_sets[i];

    //build the new model using pts
    vector<Point3d> allpts=hull_pts;
    allpts.insert(allpts.end(), pts.begin(), pts.end());

    cd_hull newhull;
    if (allpts.size()>=4) newhull.buildhull(allpts);

    //
    double vol_inc=(newhull.computeHullVol()-init_vol);
    sorted_by_vol_inc.push_back(make_pair(vol_inc, i));
  }
  //end for i

  //now we can process the point sets from small vol_inc to large vol_inc
  sort(sorted_by_vol_inc.begin(),sorted_by_vol_inc.end());

  for(auto& pts_id : sorted_by_vol_inc)
  {
    auto& pts=new_pts_sets[pts_id.second];
    //
    vector<Point3d> allpts=hull_pts;
    allpts.insert(allpts.end(), pts.begin(), pts.end());
    cd_hull newhull;
    if (allpts.size()>=4) newhull.buildhull(allpts);
    //
    double vol_inc=(newhull.computeHullVol()-init_vol);

    if(vol_inc>max_vol_inc) //too big
    {
      //cout<<"vol_inc="<<vol_inc<<endl;
      return true;
    }

    //cout<<"! vol_inc="<<vol_inc<<" pts size="<<pts.size()<<endl;

    //absorb all the points
    hull_pts.swap(allpts);
    exp_hull=newhull;
  }

  double vol_inc=(exp_hull.computeHullVol()-init_vol);
  cout<<"! Wow: took all vol_inc="<<vol_inc<<endl;

  return true;
}
