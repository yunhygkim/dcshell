#include "util/CD3d_hspace_intersection.h"
#include "CD3d_model.h"
#include "CD3d_hull.h"

using namespace mathtool;
using namespace std;

// QHull Headers
extern "C" {
#include <stdio.h>
#include "libqhull.h"
#include "poly.h"
#include "qset.h"
}

inline void dual(Point3d& pt, const facetT* f, const Point3d& o)
{
	pt[0] = o[0] - float(f->normal[0] / f->offset);
	pt[1] = o[1] - float(f->normal[1] / f->offset);
	pt[2] = o[2] - float(f->normal[2] / f->offset);
}

inline void dual(coordT * pt, const VAS_plane& p, const Point3d& o)
{
	//float offset=(v.pos-o)*v.n;
	float offset = (p.p - o)*p.n;

	pt[0] = p.n[0] / offset;
	pt[1] = p.n[1] / offset;
	pt[2] = p.n[2] / offset;
}

//
//compute the intersection of the half-spaces defined by the give planes and
//the point o that is inside the intersection
//
//implemented using qhull
//
cd_m * Hspace_intersection(list<VAS_plane>& planes, const Point3d & o)
{
	int psize = planes.size();
	coordT * pt = new coordT[psize * 3];
	assert(pt);

	//convert planes to dual
	coordT * ptr = pt;
	for (list<VAS_plane>::const_iterator i = planes.begin(); i != planes.end(); i++)
	{
		dual(ptr, *i, o);
		//cout<<ptr[0]<<","<<ptr[1]<<","<<ptr[2]<<endl;
		ptr += 3;
	}

	//build convex hull from dual
	//using qhull
	static char * options = "qhull Qt Pp";
	int curlong, totlong;
	qh_init_A(stdin, stdout, stderr, 0, NULL);
	qh_initflags(options);
	qh_init_B(pt, psize, 3, false);
	qh_qhull();
	qh_check_output();

	//convert the dual back, these are the vertices
	//of the intersection
	PtVector pv;

	//global variable for Qhull
	facetT *facet;

	int fid = 0;
	FORALLfacets
	{
		Point3d pt;
		dual(pt, facet, o);
		pv.push_back(pt);
		facet->id = fid++;
	}

	//free data
	qh_freeqhull(!qh_ALL);
	qh_memfreeshort(&curlong, &totlong);
	delete[] pt;

	cd_hull hull;
	hull.buildhull(pv);
	return hull.to_m();

	// use these points to build the convex hull again...
	/*
	psize = pv.size();
	pt = new coordT[psize * 3];
	assert(pt);

	//fill up pt
	ptr = pt;
	int vid = 0;
	FOREACH_Pt(pv)
		ptr[vid * 3]   = p[0];
		ptr[vid * 3+1] = p[1];
		ptr[vid * 3+2] = p[2];
		vid++;
	FOREACH_END

	qh_init_A(stdin, stdout, stderr, 0, NULL);
	static char * options2 = "qhull Qt Pp";
	qh_initflags(options2);
	qh_init_B(pt, psize, 3, false);
	qh_qhull();
	qh_check_output();

	//get facets of the intersection
	vertexT *vertex;
	vertexT **vertexp;
	setT *vertices;
	vid = 0;
	TriVector tris; //triangles
	FORALLfacets
	{
		vertices = qh_facet3vertex(facet);
		vector<int> ids;
		FOREACHvertex_(vertices)
		{
			ids.push_back( qh_pointid(vertex->point) );
		}

		if(ids.size()!=3) cerr<<"! Error: ids.size()="<<ids.size()<<endl;
		assert(ids.size()==3);

		reverse(ids.begin(), ids.end());
		//dummy triangulation
		for (unsigned int i = 0; i < ids.size() - 2; i++)
		{
			Tri tri;
			tri[0] = ids[0];
			tri[1] = ids[i + 1];
			tri[2] = ids[i + 2];
			tris.push_back(tri);
		}
		//done

		qh_settempfree(&vertices);
	}

	//free data
	qh_freeqhull(!qh_ALL);
	qh_memfreeshort(&curlong, &totlong);
	delete[] pt;

	//build model from pt and ids
	ML ml;
	PartVector partv;
	partv.push_back(Part(0, tris.size() - 1));
	buildModels(partv, pv, tris, ml);
	assert(ml.size()==1);

	return ml.front();
	*/
}
