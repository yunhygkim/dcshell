#ifndef _CD_HULL_H_
#define _CD_HULL_H_

///////////////////////////////////////////////////////////////////////////////
// QHull Headers
extern "C"{
#include <stdio.h>
//#include "qhull.h"
#include "libqhull.h"
#include "poly.h"
#include "qset.h"
}

///////////////////////////////////////////////////////////////////////////////
#include "util/model/IPolygonal.h"
#include <cassert>
#include <list>
#include <map>
using namespace std;

///////////////////////////////////////////////////////////////////////////////
struct cd_v;
struct cd_e;
struct cd_f;
class cd_m;
class cd_hull;

///////////////////////////////////////////////////////////////////////////////
struct hull_plane
{
	hull_plane(){ error = DBL_MAX; }
	Point3d o;
	Vector3d n;
	double error;
};

///////////////////////////////////////////////////////////////////////////////
//
//  Following classes are private classes which have only private constructors
//
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
class cd_hull{
public:

    ///////////////////////////////////////////////////////////////////////////
    // Construction
    cd_hull();
    ~cd_hull(){ free(); }

    ///////////////////////////////////////////////////////////////////////////
    // Core Functions

    ///build the convex hull from the model
    void buildhull(cd_m * m);

    ///build the convex hull from tips and features
    void buildhull(const vector<cd_v*>& vs);

	  ///build the convex hull from a set of points
  	void buildhull(const vector<Point3d>& pts);

    //get the id of the i-th vertex of facet fid
    int vid(int fid, int i) const { return hull_vid[fid*3+i]; }

    // check whether a vertex is inside the convex hull
    bool isInsideHull(cd_v * vertex);
  	bool isInsideHull(const Point3d& p);
  	bool isInsideHull(const Point3d& p) const;

    //get normal direction of the vertex on the hull
    //note that vertex must be on the convex hull
    const Vector3d& getHullNormal(cd_v * vertex);

    //operator =
		cd_hull& operator=(const cd_hull& other)
		{
			//Data for the facet of hull
			this->fsize = other.fsize;
			this->hull_vid = new int[fsize * 3];   			    //number of vertices will be fsize*3
			assert(this->hull_vid);

			memcpy(this->hull_vid, other.hull_vid, sizeof(int)*fsize * 3);

			this->hull_vertices = other.hull_vertices;	//number of vertices will be fsize*3
			this->facet_normals = other.facet_normals;	//number of normals will be fsize
			this->input_points = other.input_points;	//input vertices to compute the convex hull
			this->v2normal = other.v2normal;
			this->volume = other.volume;

			return *this;
		}

		//
		// convert this convex hull to a model object
		// this model does not reuse the cd_v from input and therefore is
		// brand new; you should destroy the model when it is not needed.
		//
		// offset is the expansion applied to the new model
		//
		cd_m * to_m(float offset=0) const;


		//
		// convert this convex hull to a vector of points on the convex hull
		//
		void to_pvector(vector<Point3d>& pts) const;

		// simplify this convex hull, returns a brand new cd_m
		//
		// k is the number of facets left
		// offset is the expansion applied to the simplified convex shape
		//
		cd_m * simplify(int k, float offset) const;  //VSA-based method, faster but no volume control

		//progress hull, LP-based method
		//enclosing_constrains ensure that the new hull are inside the given hull planes
		void simplify(int k, float max_vol_inc, const list<hull_plane>& enclosing_constrains, cd_hull & simplified_hull) const;

    // remesh this convex hull, returns a brand new cd_m
    //
    // k is the number of facets left
    // offset is the expansion applied to the remeshed convex shape
    //
    cd_m * remesh(int k, float offset) const;  //VSA-based method, faster but no volume control

    //remesh the hull, QP-based method
    //enclosing_constrains ensure that the new hull are inside the given hull planes
    void remesh(float max_vol_inc, const list<hull_plane>& enclosing_constrains, cd_hull & remeshed_hull) const;


		//compute the volume of the hull
		double computeHullVol();

		//sample points inside and on the hull
		void samplePoints(int k, double d); //d is surface coverage
		void clearSamplePoints() { sampled_points.clear(); }
		const vector<Point3d> & getSampledPoints() const { return sampled_points; }
		void addSampledPoints(const vector<Point3d>& samples)
		{
			sampled_points.insert(sampled_points.end(), samples.begin(),samples.end());
		}

    ///////////////////////////////////////////////////////////////////////////////
    // access functions

    //note that the id is NOT vertex id but fid*3+vid
    int getFaceSize() const { return fsize; }
    const Point3d& getVertex(int id) const { return hull_vertices[id]; }
    const Vector3d& getTriNormal(int fid) const { return facet_normals[fid]; }
    int * getHullVID() { return hull_vid; }
    bool isHullBuilt() { return hull_built; }

///////////////////////////////////////////////////////////////////////////////
// Proected Methods
protected:

    ///////////////////////////////////////////////////////////////////////////
    // General Create / Delete
    int callQhull(char* options, coordT *pt, int vsize);
    void copyPoints(coordT * pt, const list<cd_v*>& vs);
    void copyPoints(coordT * pt, const vector<cd_v*>& vs);
	  void copyPoints(coordT * pt, const vector<Point3d>& vs);
    void storeHull();
    void free();

    ///////////////////////////////////////////////////////////////////////////
    // Other stuff..
    bool isVisible(const Point3d& p, int fid) const;

///////////////////////////////////////////////////////////////////////////////
// Private Data
private:

    //Data for the facet of hull
    int fsize;
    int * hull_vid;   			    //number of vertices will be fsize*3
    vector<Point3d> hull_vertices;	//number of vertices will be fsize*3
    vector<Vector3d> facet_normals;	//number of normals will be fsize
	  vector<Point3d> input_points;	//input vertices to compute the convex hull

		vector<Point3d> sampled_points;	//points sampled in and on the hull

		double volume; //the volume of this convex hull

    //mapping vid to a normal direction of the normal
    //note that this id is point id, i.e. the order that points are added to the convex hull
    map<int,Vector3d> v2normal;

    // used for qhull
    int curlong, totlong;

    //for visibility test
    int last_fid; //face id

    //whether qhull return an error
    bool hull_built;
};

#endif //_CD_HULL_H_
