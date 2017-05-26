#include "util/CD3d_vsa.h"
#include "CD3d_model.h"
#include "util/CD3d_eigen.h"

#include<map>
#include<set>
#include <cfloat>

using namespace mathtool;
using namespace std;

//
// compute Variational Shape Approximation
//

struct VAS_proxy
{
	VAS_proxy()
	{
		seed = NULL;
		flag = 0;
	}

	float error(cd_f * f)
	{
		return (plane.n - f->n).normsqr();
	}

	float error()
	{
		float err = 0;
		FOREACH_F(fl) err += error(f)*f->area(); FOREACH_END
		return err;
	}

	//recenter the proxy
	void update()
	{
		if (seed == NULL) return;

		//compute the area-weighted normal of all facets
		Vector3d avg_n(0,0,0);
		FOREACH_F(fl)
			if (f->flag != this->flag) continue; //no longer belong to this proxy
			avg_n = avg_n + f->n*f->area();
		FOREACH_END
		plane.n = avg_n.normalize();

		//find the closest facet (re-seed)
		cd_f * best_f = NULL;
		float min_err = FLT_MAX;
		FOREACH_F(fl)
			if (f->flag != this->flag) continue; //no longer belong to this proxy
			float err = error(f);
			if (err < min_err)
			{
				min_err = err;
				best_f = f;
			}
		FOREACH_END

		//update
		seed = best_f;
		fl.clear();

		if (seed != NULL)
		{
			fl.push_back(seed);
			seed->flag = flag = generateID(); //get a new flag...
		}
	}

	cd_f * seed;
	FL fl;
	VAS_plane plane;
	flagT flag;
};

typedef vector<VAS_proxy>    PxL;
typedef PxL::iterator        PxIT;
#define FOREACH_Px(pxl)      {for(PxIT ipx=pxl.begin();ipx!=pxl.end();ipx++)  { VAS_proxy& px=*ipx;

inline Point3d find_extreme_pt(const VL& vl, const Vector3d& vec)
{
	float max_dist = -FLT_MAX;
	Point3d extreme_pt;
	const Point3d O(0, 0, 0);

	FOREACH_CV(vl)
		//if (f->flag != this->flag) continue; //no longer belong to this proxy
		//for (int i = 0; i < 3; i++)
	{
		const Point3d& pt = v->pos;
		float d = (pt - O)*vec;
		if (d>max_dist)
		{
			max_dist = d;
			extreme_pt = pt;
		}
	}//end for i
	FOREACH_END

		return extreme_pt;
}

inline PxL VSA_seeds(cd_m * m, int k)
{
	vector<cd_f*> fv(m->f().begin(),m->f().end());
	set<cd_f*> seeds;
	int fsize = fv.size();

	while (seeds.size() != k)
		seeds.insert(fv[fsize*drand48()]);

	//create proxy for each seed facet
	PxL pxl;
	for (set<cd_f*>::iterator i = seeds.begin(); i != seeds.end(); i++)
	{
		VAS_proxy px;
		px.seed = *i;
		px.plane.n = (*i)->n;
		px.fl.push_back(*i);
		(*i)->flag = px.flag = generateID();
		pxl.push_back(px);
	}

	return pxl;
}

inline void VSA_propogate(PxL& proxies)
{
	map<int, PxIT> flag2proxy;
	list<cd_f *> open;

	FOREACH_Px(proxies)
		if (px.seed == NULL) continue;
		flag2proxy[px.flag] = ipx;
		open.push_back(px.seed);
	FOREACH_END

	while (open.empty() == false) //loop when open is not empty
	{
		cd_f * f=open.front();
		open.pop_front();
		
		//check neighbors
		for (int i = 0; i < 3; i++)
		{
			cd_f * nei = f->nei(i);
			if (nei == NULL) continue;

			bool grab_nei = false;

			if (flag2proxy.find(nei->flag) == flag2proxy.end()) //nei is not assigned to ANY proxy
			{
				grab_nei = true;
			}
			else if (nei->flag!=f->flag) //nei is assigned to a different proxy, check if we find a better proxy here
			{
				float err_old = flag2proxy[nei->flag]->error(nei);
				float err_new = flag2proxy[f->flag]->error(nei);
				if (err_new < err_old)
				{
					grab_nei = true;
				}
			}

			if (grab_nei)
			{
				nei->flag = f->flag;
				flag2proxy[f->flag]->fl.push_back(nei);
				open.push_back(nei);
			}

		}//end for

	}//end while

}

inline float VSA_updateproxy(PxL& proxies)
{
	float err = 0;
	FOREACH_Px(proxies)
		if (px.seed == NULL) continue;
		err += px.error();
		px.update();
	FOREACH_END

	return err;
}

list<VAS_plane> VSA(cd_m * m, int k)
{
	list<VAS_plane> planes;
	k = k - 6; //reserve 6 for oriented bounding box
	if (k > 0) //start VAS
	{

		PxL pxl = VSA_seeds(m, k);

		const int max_iter = 50;
		const float err_threshold = 1e-4;
		float last_err = FLT_MAX;
		int iter = 0;

		while (iter++ < max_iter)
		{
			VSA_propogate(pxl);

			float err = VSA_updateproxy(pxl);
			if (fabs(err - last_err) < err_threshold) break;
			last_err = err;
		}//end while

		FOREACH_Px(pxl)
			if (px.fl.empty()) continue; //invalid proxy...
		px.plane.n = -px.plane.n;
		px.plane.p = find_extreme_pt(m->v(), px.plane.n);
		planes.push_back(px.plane);
		FOREACH_END
	}

	//add 6 PA planes
	PV pts;
	FOREACH_V(m->v()) pts.push_back(v->pos); FOREACH_END
	Vector3d PAs[3];
	Point3d com=EigenVectors(pts, PAs[0], PAs[1], PAs[2]);
	for (int i = 0; i < 6; i++)
	{
		VAS_plane plane;
		plane.n = PAs[i/2];
		if (i % 2 == 0) plane.n = -plane.n;
		plane.p = find_extreme_pt(m->v(), plane.n);
		float len = (plane.p - com).norm()-0.01; //shrink a very tiny bit
		plane.p = com + (plane.p - com).normalize()*len;
		planes.push_back(plane);
	}
	//done with PA planes

	return planes;
}
