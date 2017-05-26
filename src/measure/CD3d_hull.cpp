#include "measure/CD3d_hull.h"
#include "CD3d_model.h"


///////////////////////////////////////////////////////////////////////////////
//extern "C"{
//char qh_version[] = "cd hull";
//}

///////////////////////////////////////////////////////////////////////////////
#include <set>
using namespace std;

///////////////////////////////////////////////////////////////////////////////
// Hull F Functions
cd_hull::cd_hull()
{
    hull_vid=NULL;
    last_fid=0;
    fsize=0;
    curlong=totlong=0;
    hull_built=false;
	volume = 0;
}

void cd_hull::buildhull(cd_m * m)
{
    int vsize=m->v().size();
	if (vsize < 4) return; //nothing to do...

	//this->input_vertices = vector<cd_v*>(m->v().begin(), m->v().end());

    //allocate memory
    coordT * pt=new coordT[vsize*3]; //each point will have three coord
    assert(pt);

    //copy points
    copyPoints(pt,m->v());

    //using qhull
    static char * options=(char *)"qhull Qt Pp";
    callQhull(options, pt, vsize);

    //store results
    storeHull();


    ////add by Guilin 01/06/2014
    ////---------------------------------------//
    //vector<cd_v*> vecs(m->v().begin(), m->v().end());
    //int * hvid = this->getHullVID();
    //int hvid_size = this->getFaceSize()*3;
    //for(int i = 0; i < hvid_size; i++)
    //{
    //	vecs[hvid[i]]->isHullV = true;
    //}

    //not used
    delete [] pt; pt=NULL;
}

void cd_hull::buildhull(const vector<cd_v*>& vs)
{
	if(hull_vid)
	{
		this->free();
	}

	int vsize=vs.size();

	//allocate memory
	coordT * pt=new coordT[vsize*3]; //each point will have three coord
	assert(pt);

	//copy points
	copyPoints(pt, vs);

	//using Qhull
	static char * options=(char *)"qhull Qt Pp";
	if(!callQhull(options, pt, vsize))
	{
		//store Qhull
		storeHull();

		this->hull_built = true;
	}

	//not used
	delete [] pt; pt=NULL;
}


///build the convex hull from a set of points
void cd_hull::buildhull(const vector<Point3d>& pts)
{
	if (hull_vid)
	{
		this->free();
	}

	int vsize = pts.size();

	//allocate memory
	coordT * pt = new coordT[vsize * 3]; //each point will have three coord
	assert(pt);

	//copy points
	copyPoints(pt, pts);

	//using Qhull
	static char * options = (char *)"qhull QJ Pp";
	if (!callQhull(options, pt, vsize))
	{
		//store Qhull
		storeHull();

		this->hull_built = true;
	}

	//not used
	delete[] pt; pt = NULL;
}

int cd_hull::callQhull(char* options, coordT *pt, int vsize)
{
//	qh_init_A(stdin, stdout, stderr, 0, NULL);
//	qh_initflags(options);
//	qh_init_B (pt, vsize, 3, false);
//	qh_qhull();

	int exitcode = qh_new_qhull (3, vsize, pt, false, options, NULL, stderr);

	return exitcode;

	//qh_check_output();
}


bool cd_hull::isInsideHull(cd_v * vertex)
{
	return isInsideHull(vertex->pos);
}

bool cd_hull::isInsideHull(const Point3d& p)
{
	if(this->isVisible(p, last_fid))
		return false;

	for(int i=0;i<this->fsize;i++)
	{
		if(!this->isVisible(p, i))
			continue;

		// visible, not in the convex hull, store the facet id, return false
		last_fid = i;
		return false;
	}

	//all invisible, in the convex hull
	return true;
}

bool cd_hull::isInsideHull(const Point3d& p) const
{


	for (int i = 0; i<this->fsize; i++)
	{
		if (!this->isVisible(p, i))
			continue;

		// visible, not in the convex hull, store the facet id, return false
		return false;
	}

	//all invisible, in the convex hull
	return true;
}

///////////////////////////////////////////////////////////////////////////////
//
//
//Protected
//
//
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// General Create / Delete
///////////////////////////////////////////////////////////////////////////////

void cd_hull::copyPoints(coordT * pt, const VL& vs)
{
	this->input_points.clear();
    FOREACH_CV(vs)
        const Point3d& p=v->pos;
        for( int i=0;i<3;i++ ) pt[v->id*3+i]=p[i];
		this->input_points.push_back(p);
    FOREACH_END
}

void cd_hull::copyPoints(coordT * pt, const vector<cd_v*>& vs)
{
	int vid = 0;
	this->input_points.clear();
	typedef vector<cd_v*>::const_iterator CVVIT;
    for(CVVIT it = vs.begin(); it != vs.end(); ++it)
    {
        const Point3d& p=(*it)->pos;
        for( int i=0;i<3;i++ ) pt[vid++]=p[i];
		this->input_points.push_back(p);
    }
}

void cd_hull::copyPoints(coordT * pt, const vector<Point3d>& vs)
{
	int vid = 0;
	this->input_points.clear();
	for (const auto & p : vs)
	{
		for (int i = 0; i<3; i++) pt[vid++] = p[i];
		this->input_points.push_back(p);
	}
}

//Allocate space for storing and store normal of facets
void cd_hull::storeHull()
{
    //global variable for Qhull
    facetT *facet;
    vertexT *vertex;
    vertexT **vertexp;
    setT *vertices;

    this->fsize=qh num_facets;

    //allocate
    hull_vid=new int[fsize*3];
    assert(hull_vid!=NULL);

    this->hull_vertices.reserve(fsize*3);
    this->facet_normals.reserve(fsize);

    int fid=0;
    int vid=0;

    //mapping vid to normal
    v2normal.clear();

    FORALLfacets {
        facet->id=fid++;
        facet_normals.push_back(Vector3d(facet->normal).normalize());
        vertices= qh_facet3vertex (facet);
        FOREACHvertex_(vertices) {
        	hull_vertices.push_back(Point3d(vertex->point));
            hull_vid[vid]=qh_pointid(vertex->point);

            if(v2normal.find(hull_vid[vid])==v2normal.end())
            {
                v2normal[hull_vid[vid]]=facet_normals.back();
            }
            else
            {
                Vector3d& vec=v2normal[ hull_vid[vid] ];
                vec = vec +  facet_normals.back();
            }
        	vid++;
        }

        qh_settempfree(&vertices);
    }

    //normalize the normal directions
    for(map<int,Vector3d>::iterator i=v2normal.begin();i!=v2normal.end();i++)
    {
        i->second=i->second.normalize();
    }//end for i

}

void cd_hull::free()
{
	if(hull_vid!=NULL)
	{
		delete [] hull_vid;
		hull_vid = NULL;

		this->facet_normals.clear();
		this->hull_vertices.clear();
		this->input_points.clear();

	  //free mem
    qh_freeqhull(!qh_ALL);
    qh_memfreeshort (&curlong, &totlong);
	}
}


bool cd_hull::isVisible(const Point3d& p, int fid) const
{
    //for(int i=0;i<3;i++)
    {
        const double tiny_number=1e-8;

        const Point3d& pt_on_hull = this->hull_vertices[fid*3];
        Vector3d v = (p - pt_on_hull); //.normalize();

        //TODO: or we can simply don't normalize...
        double v_norm=v.norm();
        if(v_norm<tiny_number) return false;

        v=v/v_norm;

        float dot=(v*this->facet_normals[fid]);
        if(dot>tiny_number)
        {
            return true;
        }
    }

    return false;
}

const Vector3d& cd_hull::getHullNormal(cd_v * vertex)
{
    if(vertex!=NULL)
    {
        if(v2normal.find(vertex->id)!=v2normal.end())
            return v2normal[vertex->id];
    }

    throw cd_msg("cd_hull::getHullNormal Error");
}

//
// convert this convex hull to a model object
// this model does not reuse the cd_v from input and therefore is
// brand new; you should destroy the model when it is not needed.
//
// offset is the expansion applied to the new model
//
#include "util/CD3d_util.h"
cd_m * cd_hull::to_m(float offset) const
{
	int tsize = getFaceSize();
	TriVector tris; //triangles

	set<int> vidset;
	for (int i = 0; i<tsize; i++)
		for (int d = 2; d >= 0; d--) vidset.insert(hull_vid[3 * i + d]);

	PtVector pv;
	map<int, int> id2id;
	int vid = 0;
	for (set<int>::iterator i = vidset.begin(); i != vidset.end(); i++, vid++)
	{
		id2id[*i] = vid;
		pv.push_back(this->input_points[*i]);
	}

	for (int i = 0; i<tsize; i++)
	{
		Tri tri;
		for (int d = 2; d >= 0; d--) tri[d] = id2id[ hull_vid[3 * i + d] ];
		swap(tri[1], tri[2]);
		tris.push_back(tri);
	}

	//compute offset
	PtVector pv_offset;
	if (offset != 0)
	{
		Point3d com = computeCenter(pv);
		FOREACH_CPt(pv)
			Point3d newp = p + (p - com).normalize()*offset;
		pv_offset.push_back(newp);
		FOREACH_END
	}

	//build model from pt and ids
	ML ml;
	PartVector partv;
	partv.push_back(Part(0, tris.size() - 1));
	if (offset != 0)
		buildModels(partv, pv_offset, tris, ml);
	else
		buildModels(partv, pv, tris, ml);

  if(ml.size()!=1)
  {
    cout<<"vidset size="<<vidset.size()<<" tsize="<<tsize<<endl;
    cout<<"ml.size()="<<ml.size()<<endl;
    assert(ml.size() == 1);
  }

	return ml.front();
}


void cd_hull::to_pvector(vector<Point3d>& pts) const
{
	int tsize = getFaceSize();
	TriVector tris; //triangles

	set<int> vidset;
	for (int i = 0; i < tsize; i++)
	for (int d = 2; d >= 0; d--) vidset.insert(hull_vid[3 * i + d]);

	pts.reserve(vidset.size());
	for (set<int>::iterator i = vidset.begin(); i != vidset.end(); i++)
	{
		pts.push_back(this->input_points[*i]);
	}//
}


// simplify this convex hull, returns a brand new cd_m
// k is the number of facets left
// offset is the expansion applied to the simplified convex shape
#include "util/CD3d_vsa.h"
#include "util/CD3d_hspace_intersection.h"
#include "util/CD3d_progressive_hull.h"

cd_m * cd_hull::simplify(int k, float offset) const
{
	cout << "- VSA-based methods is called" << endl;

	//convex hull to cd_m first
	cd_m * hull_m = to_m(offset);

	Point3d com = computeCenter(hull_m->v());

	//call VSA
	list<VAS_plane> planes=VSA(hull_m, k);

	//compute variational-shape approximation of hull_m
	cd_m * simplified_hull = Hspace_intersection(planes, com);

	//delete hull_m
	hull_m->destroy();
	delete hull_m;

	return simplified_hull;
}

void cd_hull::simplify(int k, float max_vol_inc, const list<hull_plane>& enclosing_constrains, cd_hull & simplified_hull) const
{
    cout << "- ProgressiveHull is called" << endl;
    ProgressiveHull(*this, simplified_hull, k, max_vol_inc, enclosing_constrains);
}

void cd_hull::remesh(float max_vol_inc, const list<hull_plane>& enclosing_constrains, cd_hull & remeshed_hull) const
{
    // remesh the hull using qp
    ProgressiveHullQP(*this, remeshed_hull, max_vol_inc, enclosing_constrains);
}

inline double tetra_vol(const Point3d& p1, const Point3d& p2, const Point3d& p3, const Point3d& p4)
{
	auto N = ((p2 - p1) % (p3 - p1));
	return -((p4 - p1)*N) / 6;
}

//compute the volume of the hull
double cd_hull::computeHullVol()
{
	if (this->volume == 0)
	{

		double vol = 0;
		int tsize = getFaceSize();

		const Point3d orig(0, 0, 0);
		int id = 0;
		for (int i = 0; i < tsize; i++)
		{
			Tri tri;
			const auto & p3 = this->input_points[hull_vid[id++]];
			const auto & p2 = this->input_points[hull_vid[id++]];
			const auto & p1 = this->input_points[hull_vid[id++]];

			vol += tetra_vol(p1, p2, p3, orig);
		}

		this->volume = vol;
	}

	return this->volume;
}


//sample points inside and on the hull
void sample_in_hull(const cd_hull& hull, cd_m * hull_m, vector<Point3d>& samples, int k);
void sample_from_triangle(cd_f * f, double d, vector<Point3d>& samples);
void cd_hull::samplePoints(int k, double d)
{
  //cout<<"samplePoints k="<<k<<endl;
  srand48(k);
	cd_m * m = to_m();

	//sample the points from inside the hull
 	sample_in_hull(*this, m, sampled_points, k);
//
// 	for (auto& f : m->f())
// 	{
// 		sample_from_triangle(f, d, sampled_points);
// 	}
//

	//now add the vertice of the hull into the sample as well.
	//for (auto& p : m->v()) this->sampled_points.push_back(p->pos);

	m->destroy();

	delete m;
}


void sample_from_triangle(cd_f * f, double d, vector<Point3d>& samples)
{
	const Point3d& p1 = f->v[0]->pos;
	const Point3d& p2 = f->v[1]->pos;
	const Point3d& p3 = f->v[2]->pos;

	//compute triangle area
	Vector3d v1 = p2 - p1;
	Vector3d v2 = p3 - p1;
	Vector3d v3 = p3 - p2;
	Vector3d v4 = p1 - p3;
	float area = (v1%v2).norm() / 2;
  if(area==0) return;

	//compute number of points needed
	int size = (int)ceil(area / (d*d));

	double r[2];

	//randomly generate these points
	while (size>0)
	{
		//halton(r);
		//hammersley_sequence(1,r);

		r[0] = drand48();
		r[1] = drand48();
		Point3d pos = p1 + v1*r[0] + v2*r[1];

		//check if inside the triangle
		if ((v1 % (pos - p1))*f->n <= 0) continue;
		if ((v3 % (pos - p2))*f->n <= 0) continue;
		if ((v4 % (pos - p3))*f->n <= 0) continue;
		//in tri

		samples.push_back(pos);
		size--;
	}//j
}

//
//make k samples uniformlly inside the hull
//
void sample_in_hull(const cd_hull& hull, cd_m * hull_m, vector<Point3d>& samples, int k)
{
	//
	float bbox[6];
	computeBoundingBox(hull_m->v(), bbox);

	//sample inside the hull
	//for (int i = 0; i < k; i++)
  while(k>0)
	{
		Point3d sample;
		//draw a sample from box
		for (int j = 0; j < 3; j++)
			sample[j] = (bbox[j * 2 + 1] - bbox[j * 2])*drand48() + bbox[j * 2];
		//check if the point is inside the hull
		if (hull.isInsideHull(sample))
    {
			samples.push_back(sample);
      k--;
    }
	}//end while
}
