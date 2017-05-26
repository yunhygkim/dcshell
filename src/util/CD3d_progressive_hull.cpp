#include "util/CD3d_progressive_hull.h"
#include "CD3d_model.h"

#include<map>
#include<set>
//#include<math.h>

#include <glpk.h> //linear programming solver
#include <CGAL/basic.h>     //quadratic programming solver
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>

#ifdef CGAL_USE_GMP
#include <CGAL/Gmpzf.h>
typedef CGAL::Gmpzf ET;
#else
#include <CGAL/MP_Float.h>
typedef CGAL::MP_Float ET;
#endif

typedef CGAL::Quadratic_program<double> Program;
typedef CGAL::Quadratic_program_solution<ET> Solution;

using namespace mathtool;
using namespace std;

//find the best position to replace the collaping edge e
bool SolveLP(list<cd_f*>& faces, cd_e * e, const list<hull_plane>& enclosing_constrains, Point3d& solution);

//find the best position to move v so that it is close to the average of v's neighbor
bool SolveLP(list<cd_f*>& faces, cd_v * v, const list<hull_plane>& enclosing_constrains, Point3d& solution);

double least_square_fit(const vector<Point3d>& pts, Point3d& o, Vector3d& n); //implemented in CD_trim_hull.cpp

//find the optimal pos to move v so that it is close to the average of v's neighbor
bool SolveQP(list<cd_f*>& faces, cd_v * v, const list<hull_plane>& enclosing_constrains, Point3d& solution);

struct collapsed_edge
{
	collapsed_edge(){ e = NULL;  vol_inc = DBL_MAX; }
	bool operator<(const collapsed_edge& e) const { return vol_inc < e.vol_inc; }

	cd_e * e;
	//list<cd_f*> visible_faces;
	Point3d p;
	double vol_inc; //increased volume
};

struct centered_vertex
{
    centered_vertex(){ v = NULL;  vol_inc = DBL_MAX; }
    bool operator<(const centered_vertex& v) const { return vol_inc < v.vol_inc; }

    cd_v * v;   //ori position
    Point3d p;  //new position
    double vol_inc; //increased volume
};

struct moved_vertex
{
	moved_vertex(){ v = NULL;  vol_inc = DBL_MAX; }
	bool operator<(const moved_vertex& v) const { return vol_inc < v.vol_inc; }

	cd_v * v;
	Point3d p; //new position
	double vol_inc; //increased volume
};


list<cd_f*> get_incident_faces(cd_e * e)
{
	set<cd_f*> faces;

	if (e->f[0] != NULL) faces.insert(e->f[0]);
	if (e->f[1] != NULL) faces.insert(e->f[1]);

	for (cd_e * ve : e->v[0]->edges)
	{
		if (ve == e) continue;
		if (ve->f[0] != NULL) faces.insert(ve->f[0]);
		if (ve->f[1] != NULL) faces.insert(ve->f[1]);
	}

	for (cd_e * ve : e->v[1]->edges)
	{
		if (ve == e) continue;
		if (ve->f[0] != NULL) faces.insert(ve->f[0]);
		if (ve->f[1] != NULL) faces.insert(ve->f[1]);
	}

	return list<cd_f*>(faces.begin(), faces.end());
}

list<cd_f*> get_incident_faces(cd_v * v)
{
	set<cd_f*> faces;

	for (cd_e * ve : v->edges)
	{
		if (ve->f[0] != NULL) faces.insert(ve->f[0]);
		if (ve->f[1] != NULL) faces.insert(ve->f[1]);
	}

	return list<cd_f*>(faces.begin(), faces.end());
}

list<cd_e*> get_incident_edges(list<cd_f *>& faces)
{
	set<cd_e*> edges;
	for (cd_f * f : faces)
	{
		if (f == NULL) continue;
		edges.insert(f->e[0]);
		edges.insert(f->e[1]);
		edges.insert(f->e[2]);
	}
	return list<cd_e*>(edges.begin(),edges.end());
}

double volume(cd_f * f, const Point3d& p)
{
	return (f->area()/3)*((p - f->v[1]->pos)*f->n);
}

void laplacian_smoothing(cd_m * m, cd_v * v, moved_vertex& mv, const list<hull_plane>& enclosing_constrains)
{
	mv.v = v;
	list<cd_f*> faces = get_incident_faces(v);

	//
	if (SolveLP(faces, v, enclosing_constrains, mv.p) == false)
	{
		mv.vol_inc = FLT_MAX;
		return;
	}

	//compute increase volume
	mv.vol_inc = 0;
	for (cd_f * f : faces)
	{
		double vol = volume(f, mv.p);
		if (vol > 0)
		{
			mv.vol_inc += vol;
		}
		else
		{
			mv.vol_inc = FLT_MAX;
			break;
		}
	}//end of f
}

inline void getvisiblefaces(const Point3d& pt, list<cd_f*>& faces)
{
	list<cd_f*> open=faces;
	flagT flag=generateID();
	for(cd_f* f:faces) f->flag=flag;

	while(!open.empty())
	{
		cd_f* f=open.front();
		open.pop_front();
		for(short i=0;i<3;i++)
		{
			cd_f * of=f->e[i]->otherf(f);
			if(of->flag==flag) continue; //visited
			of->flag=flag;
			if( (pt - of->v[1]->pos)*of->n >0 ) //visible from this face
			{
				faces.push_back(of);
				open.push_back(of);
			}
		}//end for i
	}//end while
}



//collapse an edge e to a point p
//return the increase volume
void collapse(cd_m * m, cd_e * e, collapsed_edge& ce, const list<hull_plane>& enclosing_constrains)
{

	ce.e = e; // yh-cmt- store the collapse_edge using e
	list<cd_f*> visible_faces;
	visible_faces = get_incident_faces(e); //yh-cmt set the faces which are neighbors of the edge e and vertices of the edge e


	if (SolveLP(visible_faces, e, enclosing_constrains, ce.p) == false) //yh-cmt- uses glpk library for linear programming // ce.p is a output which we want to find by solving linear problems // we will replace the edge e with the ce.p // we decrease edges to simplify the mesh model using this // yh?-cmt- enclosing_constraints are able to make the resulting simplified model is still a convex.
	{
		ce.vol_inc = FLT_MAX;
		return;
	}

	//validate the result, make sure if all faces are visible by ce.p
	for (cd_f * f : visible_faces)
	{
		if( (ce.p - f->v[1]->pos)*f->n <=0 ) //invisible from this face
		{
			ce.vol_inc = FLT_MAX;
			return;
		}
	}

    //yh-cmt if the solver finds ce.p
	//compute increased volume

	//get all faces that are visible by ce.p
	getvisiblefaces(ce.p, visible_faces);
	ce.vol_inc = 0; //yh-cmt initialize


	for (cd_f * f : visible_faces)
	{
		double vol = volume(f, ce.p);
		if (vol >= 0)
		{
			ce.vol_inc += vol;
		}
		else
		{ //this should not happen
			cerr << "! Error: invalid solution! vol=" << vol << " f area="<<f->area()
			     <<" dot="<<(ce.p - f->v[1]->pos)*f->n <<endl;
			exit(1);
		}
	}//end of f
}

//
// using QP, find an optimal vertex, which is close to a center point
//
void moveToCenter(cd_m * m, cd_v * v, centered_vertex& cv, const list<hull_plane>& enclosing_constrains)
{
    cv.v = v;
    list<cd_f*> faces = get_incident_faces(v); //neighbor faces

    //find an optimal vertex
    if(SolveQP(faces, v, enclosing_constrains, cv.p)==false)
    {
        cv.vol_inc = FLT_MAX; //QP fails
        return;
    }

		//JML: I think this is too slow...
    //store new v(p)
    // flagT flag = generateID();
    // cv.v->flag = flag;
		//
    // vector<Point3d> new_pts;
    // new_pts.push_back(cv.p);
		//
    // //collect the remaining points
    // for(cd_v * v: m->v())
    // {
    //     if(v->flag == flag) continue;
    //     new_pts.push_back(v->pos);
    // }
		//
    // //build a temp convex
    // cd_hull new_hull;
    // new_hull.buildhull(new_pts);
		//
    // //check vol_inc
    // float new_vol = new_hull.computeHullVol();
    // float vol_inc=new_vol-prev_vol;
		//
    // if(vol_inc > 0 ) cv.vol_inc=vol_inc;
    // else    cv.vol_inc=FLT_MAX;


    //compute increase volume
    cv.vol_inc =0;
    getvisiblefaces(cv.p,faces);

    for(auto f: faces)
    {
        double vol=volume(f, cv.p);

        if(vol>0){
					  cv.vol_inc += vol;  //vaild
				}

        else
        {
            cv.vol_inc=FLT_MAX;
            break;
        }
    }
}

//
// max_fsize: maximum number of faces left in the hull
// max_vol_inc : max vol increase allowed, if max_vol_inc is violated, the function will return the hull even though it has more faces than "max_fsize"
// enclosing_constrains: a list of planes that enforce the new hull the be constrained inside
//
// !! Note: this function assumes that the input_hull is inside enclosing_constrains
//
bool ProgressiveHull(const cd_hull& input_hull, cd_hull& simplified_hull, int max_fsize, float max_vol_inc, const list<hull_plane>& enclosing_constrains)
{
	cout << "input_hull.getFaceSize()=" << input_hull.getFaceSize() << endl;
	simplified_hull = input_hull;

	if (input_hull.getFaceSize() < max_fsize) //done
	{
		return true;
	}

	//convex hull to cd_m first
	cd_m * hull_m = input_hull.to_m(0);

	cout << "------------------------------------" << endl;
	float total_vol_inc = 0;
	float init_vol = simplified_hull.computeHullVol(); // Yun: check the current volume to keep tracking

		//Point3d com(0, 0, 0);
		//for (auto& v : hull_m->v())
		//{
		//	com[0] += v->pos[0];
		//	com[1] += v->pos[1];
		//	com[2] += v->pos[2];
		//}
		//com[0] /= hull_m->v().size();
		//com[1] /= hull_m->v().size();
		//com[2] /= hull_m->v().size();

		//for (auto& f : hull_m->f())
		//{
		//	double dot = f->n*(com - f->v[0]->pos);
		//	cout << "dot=" << dot << endl;
		//}

		//


#if 0
	while(true)
	{
		vector<Point3d> new_pts;
		collapsed_edge best_ce;

		for (auto & e : hull_m->e()) // yun: look through each edge of the convex shape
		{
			collapsed_edge ce;  // yun: track how much the volume increase
      collapse(hull_m, e, ce, enclosing_constrains); // yun: collapse an edge e to a point p and return the increase volume
			if (ce.vol_inc < best_ce.vol_inc) // yun: find out what is the bast edge to collapse (the best edge - the smallest volume increases)
			{
				best_ce = ce;
			}
		}//

		total_vol_inc += best_ce.vol_inc;
		if (total_vol_inc > max_vol_inc)
		{
			break;
		}

		flagT flag = generateID();
		best_ce.e->flag = flag;
		best_ce.e->v[0]->flag = flag;
		best_ce.e->v[1]->flag = flag;
		new_pts.push_back(best_ce.p);

		//collect the remaining points
		for (cd_v * v : hull_m->v())
		{
			if (v->flag == flag) continue;
			new_pts.push_back(v->pos);
		}

		//build a new convex shape out of new_pts
		cd_hull new_hull;
		new_hull.buildhull(new_pts);

		float new_vol = new_hull.computeHullVol();
		total_vol_inc=new_vol - init_vol;
		if (total_vol_inc > max_vol_inc)   //yun: check the maximum volume increases
		{
			break;
		}

		//cout << "new_hull.getFaceSize()=" << new_hull.getFaceSize() << "max_fsize=" << max_fsize<< endl;
		if (new_hull.getFaceSize() < max_fsize)
		{
			// if (simplified_hull.getFaceSize() == 0)
			// {
			// 	simplified_hull = new_hull;
			// }

			hull_m->destroy();
			delete hull_m;

			//cout << "Good new_hull.getFaceSize()=" << new_hull.getFaceSize() << endl;
			return true;
		}

		else
		{
			hull_m->destroy();
			delete hull_m;
			hull_m = new_hull.to_m(0);
			simplified_hull = new_hull;
		}

	}//end while(true)

#else

	int init_face_size=simplified_hull.getFaceSize();

	while(true)
	{
		vector<collapsed_edge> collapsed_edges;
		vector<Point3d> new_pts;

		for (auto & e : hull_m->e()) // yun: look through each edge of the convex shape
		{
			collapsed_edge ce;  // yun: track how much the volume increase
			collapse(hull_m, e, ce, enclosing_constrains); // yun: collapse an edge e to a point p and return the increase volume
			if(ce.vol_inc < max_vol_inc) collapsed_edges.push_back(ce);
		}//

		sort(collapsed_edges.begin(),collapsed_edges.end());

		flagT flag = generateID();
		uint count=0;
		for(auto& ce : collapsed_edges)
		{
			total_vol_inc += ce.vol_inc;
			if(total_vol_inc > max_vol_inc) break;
			ce.e->flag = flag;
			ce.e->v[0]->flag = flag;
			ce.e->v[1]->flag = flag;
			new_pts.push_back(ce.p);
			count++;
	  }

		if(count==0) break;

		//collect the remaining points
		for (cd_v * v : hull_m->v())
		{
			if (v->flag == flag) continue; // yun: the points which are not flagged and delete the flagged point which is the best collapsed edge
			new_pts.push_back(v->pos);
		}

		//build a new convex shape out of new_pts
		cd_hull new_hull;
		new_hull.buildhull(new_pts);
		simplified_hull = new_hull;

		if(init_face_size==new_hull.getFaceSize())
			break; //no improvement.
		//cout<<"new_hull.getFaceSize()="<<new_hull.getFaceSize()<<endl;
	}


#endif

	hull_m->destroy();
	delete hull_m;

	return false;
}

//
// Find an optimal vertex, which is close to the center of the neighboring vertices
// enclosing_constrains: a list of planes that enforce the new hull to be constrained inside
//
bool ProgressiveHullQP(const cd_hull& input_hull, cd_hull& remeshed_hull, float max_vol_inc, const list<hull_plane>& enclosing_constrains)
{
    remeshed_hull = input_hull;

    //convex hull to cd_m first
    cd_m * hull_m = input_hull.to_m(0);

    float total_vol_inc = 0;
    float init_vol = remeshed_hull.computeHullVol();

    Vector<double, 3> init_fatness;
    Vector<double, 3> curr_fatness;
    
    //while(true) // do work repeatedly before satisfying the termination condition
    {
        vector<centered_vertex> cvs;
        centered_vertex best_cv; //best centered v

        // find the best v
        for(auto& v: hull_m->v())
        {
            centered_vertex cv;
            moveToCenter(hull_m, v, cv, enclosing_constrains); //find the best new vertex
            //if(cv.vol_inc < max_vol_inc)    cvs.push_back(cv);
            cvs.push_back(cv);
            if(cv.vol_inc<best_cv.vol_inc)  best_cv=cv;
        }

        //sort vertices from small vol inc to larger vol in
		sort(cvs.begin(), cvs.end());

		vector<Point3d> new_pts;
		for(auto& cv : cvs)
		{
			if(cv.vol_inc > max_vol_inc) //the vol_inc is either too big or invalid
			{
				new_pts.push_back(cv.v->pos);
				continue;
            }

            //increase total vol
            total_vol_inc+=cv.vol_inc;

            if(total_vol_inc<max_vol_inc) //valid, still small enough
            {
                new_pts.push_back(cv.p);
            }
            else //too big
            {
                new_pts.push_back(cv.v->pos);
            }
        }//end cvs

        // //check max_vol_inc
        // total_vol_inc += best_cv.vol_inc;
        // if(total_vol_inc>max_vol_inc)
        // {
        //     cout<<"(n-3) stop the process: total_vol_inc > max_vol_inc"<<endl;
        //     //cout<<"* total_vol_inc: "<<total_vol_inc<<" / max_vol_inc: "<<max_vol_inc<<endl;
        //     break; //done
        // }
				//
        // //if total_vol_inc <= max_vol_inc, keep going
        // //store new v(p)
        // flagT flag = generateID();
        // best_cv.v->flag = flag;
				//
				//
        // new_pts.push_back(best_cv.p);
				//
				// cout<<"move "<<best_cv.v->pos<<" to "<<best_cv.p<<" vol_inc="<<best_cv.vol_inc<<endl;
				//
        // //colect the remaining points
        // for(cd_v * v: hull_m->v())
        // {
        //     if(v->flag == flag) continue;
        //     new_pts.push_back(v->pos);
        // }

        //build a new convex shape out of new_pts
        cd_hull new_hull;
        new_hull.buildhull(new_pts);

        
        
		// float new_vol = new_hull.computeHullVol();
        // if(new_vol - init_vol > max_vol_inc)
        // {
        //     cout<<"(n-4) stop the process: (new - init_vol) > max_vol_inc"<<endl;
        //     cout<<"* vol_inc: "<<new_vol-init_vol<<" / max_vol_inc: "<<max_vol_inc<<endl;
				//
				// 		hull_m->destroy();
				// 		delete hull_m;
        //     return false;
        // }

        hull_m->destroy();
		delete hull_m;
		hull_m=new_hull.to_m(0);
		remeshed_hull=new_hull;
        
        return true;
    }
    return false;
}

//find the best position to replace the collaping edge e
bool SolveLP(list<cd_f*>& faces, cd_e * e, const list<hull_plane>& enclosing_constrains, Point3d& solution)
{
	int constraint_size = faces.size() + enclosing_constrains.size();
	int variable_size = 3;

	//cout << "solveLP variable_size=" << variable_size << " constraint_size=" << constraint_size << endl;

	glp_prob * lp = glp_create_prob();
	assert(lp);

	//start to use lp
	glp_set_prob_name(lp, "lp");
	glp_set_obj_dir(lp, GLP_MIN);
	glp_add_rows(lp, constraint_size);
	glp_add_cols(lp, variable_size);

	//init rows
	int iaja_size = 0;
	int row_id = 1;

	//the location of the new point
	double obj_coeff[3] = { 0, 0, 0 };

	//constraints: the new point needs to see all faces to remain enclosing the original convex shape
	char tmp[64];
	for (auto & f : faces)
	{
		//sprintf(tmp, "r%08d", row_id);
		double v = f->n[0] * f->v[0]->pos[0] + f->n[1] * f->v[0]->pos[1] + f->n[2] * f->v[0]->pos[2]; // e) the unknown vertex x - one of the vertices of the face = unknown vector
		glp_set_row_name(lp, row_id, "");
		glp_set_row_bnds(lp, row_id, GLP_LO, v, 0.0);

		iaja_size += 3;
		row_id++;

		//cout << "c[" << tmp << "] " << f->n[0] << "x+" << f->n[1] << "y+" << f->n[2] << "z  >= " << v << "=" << f->n[0] <<"*"<< f->v[0]->pos[0]<<" + "<<f->n[1]<<"*"<< f->v[0]->pos[1]<<" + "<<f->n[2]<<"*"<<f->v[0]->pos[2] << endl;

		double f_area = f->area();
		for(short d=0;d<3;d++) obj_coeff[d] += (f_area*f->n[d]);
	}

	//additional enclosing constrains
	for (auto& p : enclosing_constrains)
	{
		double v = p.n[0] * p.o[0] + p.n[1] * p.o[1] + p.n[2] * p.o[2];
		glp_set_row_name(lp, row_id, "");
		glp_set_row_bnds(lp, row_id, GLP_LO, v, 0.0);
		iaja_size += 3;
		row_id++;
	}

	//cout << "iaja_size=" << iaja_size << endl;

	//objective function: minimize the increase volume
	//init cols
	for (int i = 1; i <= variable_size; i++)
	{
		//char tmp[64];
		//sprintf(tmp, "s%08d", i);
		glp_set_col_name(lp, i, "");
		//glp_set_col_kind(lp, i, GLP_CV);
		glp_set_col_bnds(lp, i, GLP_FR, 0, 0);
		{
			glp_set_obj_coef(lp, i, obj_coeff[i-1]);
			//cout << "x[" << i << "] > " << obj_coeff[i - 1] << endl;
		}

	}//end i

	//init ia, ja, and ar
	int * ia = new int[1 + iaja_size];
	int * ja = new int[1 + iaja_size];
	double * ar = new double[1 + iaja_size];
	assert(ia && ja && ar);

	int ia_id = 1;
	row_id = 1;
	for (auto & f : faces)
	{
		//for (auto j : t1c.second.first)
		//for each row
		for (short d = 0; d < 3;d++)
		{
			ia[ia_id] = row_id;
			ja[ia_id] = d+1;
			ar[ia_id] = f->n[d]; //coefficient of the constrains from: x*f->n > f->v[0]*f->n , where x is unknown
			ia_id++;
		}//end for j

		row_id++;
	}//end for i

	//additional enclosing constrains
	for (auto& p : enclosing_constrains)
	{
		for (short d = 0; d < 3; d++)
		{
			ia[ia_id] = row_id;
			ja[ia_id] = d + 1;
			ar[ia_id] = p.n[d]; //coefficient of the constrains from: x*f->n > f->v[0]*f->n , where x is unknown
			ia_id++;
		}//end for j

		row_id++;
	}

	glp_load_matrix(lp, iaja_size, ia, ja, ar);

	glp_smcp parm;
	glp_init_smcp(&parm);
	parm.msg_lev = GLP_MSG_OFF;
	parm.tm_lim = 100;
	//parm.presolve=GLP_ON;
	int err = glp_simplex(lp, &parm);
	//cout << "err=" << err << endl;

	double z = glp_get_obj_val(lp);
	//cout << "objective value=" << z << endl;

	//get mip status
	int glp_prim_stat = glp_get_status(lp);

	//switch (glp_prim_stat)
	//{
	//	//case GLP_OPT: cout << "solution is optimal;" << endl; break;
	//	//case GLP_FEAS: cout << "solution is feasible;" << endl; break;
	//	case GLP_INFEAS: cout << "solution is infeasible;" << endl;   break;
	//	case GLP_NOFEAS: cout << "problem has no feasible solution;" << endl;  break;
	//	case GLP_UNBND: cout << "problem has unbounded solution;" << endl;  break;
	//	case GLP_UNDEF: cout << "solution is undefined." << endl;  break;
	//}

	bool solution_found = glp_prim_stat == GLP_OPT || glp_prim_stat == GLP_FEAS;

	if (solution_found)
	{
		solution[0] = glp_get_col_prim(lp, 1);		solution[1] = glp_get_col_prim(lp, 2);
		solution[2] = glp_get_col_prim(lp, 3);
		if (isfinite(solution[0]) == false && isfinite(solution[1]) == false && isfinite(solution[2]) == false)
		{
			solution_found = false;
		}
	}

	glp_delete_prob(lp);
	delete[] ia;
	delete[] ja;
	delete[] ar;

	return solution_found;
}



//find the best position to move v so that it is close to the average of v's neighbor
bool SolveLP(list<cd_f*>& faces, cd_v * v, const list<hull_plane>& enclosing_constrains, Point3d& solution)
{
	//find the plane for all neighboring verteices
	//the plane will be used for objective function
	vector<Point3d> nei;
	Point3d nei_center;
	for (auto e : v->edges)
	{
		const auto& pos = e->otherv(v)->pos;
		nei.push_back(pos);
		for (short i = 0; i < 3; i++) nei_center[i] += pos[i];
	}
	for (short i = 0; i < 3; i++) nei_center[i] /= v->edges.size();

	Point3d o;
	Vector3d n;
	least_square_fit(nei, o, n); //implemented in CD_trim_hull.cpp

	//
	int constraint_size = faces.size() + enclosing_constrains.size() + 1; //the last one is the volume constraint
	int variable_size = 3;

	//cout << "solveLP variable_size=" << variable_size << " constraint_size=" << constraint_size << endl;

	glp_prob * lp = glp_create_prob();
	assert(lp);

	//start to use lp
	glp_set_prob_name(lp, "lp");
	glp_set_obj_dir(lp, GLP_MAX);
	glp_add_rows(lp, constraint_size);
	glp_add_cols(lp, variable_size);

	//init rows
	int iaja_size = 0;
	int row_id = 1;

	//the location of the new point
	double obj_coeff[3] = { 0, 0, 0 };

	//constraints: the new point needs to see all faces to remain enclosing the original convex shape
	//char tmp[64];
	for (auto & f : faces)
	{
		//sprintf(tmp, "r%08d", row_id);
		double v = f->n[0] * f->v[0]->pos[0] + f->n[1] * f->v[0]->pos[1] + f->n[2] * f->v[0]->pos[2];
		glp_set_row_name(lp, row_id, "");
		glp_set_row_bnds(lp, row_id, GLP_LO, v, 0.0);
		iaja_size += 3;
		row_id++;

		//cout << "c[" << tmp << "] " << f->n[0] << "x+" << f->n[1] << "y+" << f->n[2] << "z  >= " << v << "=" << f->n[0] <<"*"<< f->v[0]->pos[0]<<" + "<<f->n[1]<<"*"<< f->v[0]->pos[1]<<" + "<<f->n[2]<<"*"<<f->v[0]->pos[2] << endl;

		double f_area = f->area();
		for (short d = 0; d<3; d++) obj_coeff[d] += (f_area*f->n[d]);
	}

	//additional enclosing constrains
	for (auto& p : enclosing_constrains)
	{
		double v = p.n[0] * p.o[0] + p.n[1] * p.o[1] + p.n[2] * p.o[2];
		glp_set_row_name(lp, row_id, "");
		glp_set_row_bnds(lp, row_id, GLP_LO, v, 0.0);
		iaja_size += 3;
		row_id++;
	}

	//set up volume constraint here

	//cout << "iaja_size=" << iaja_size << endl;

	//objective function: maximize the alignment obtween n and (x-o)
	//init cols
	for (int i = 1; i <= variable_size; i++)
	{
		//char tmp[64];
		//sprintf(tmp, "s%08d", i);
		glp_set_col_name(lp, i, "");
		//glp_set_col_kind(lp, i, GLP_CV);
		glp_set_col_bnds(lp, i, GLP_FR, 0, 0);
		{
			glp_set_obj_coef(lp, i, obj_coeff[i - 1]);
			//cout << "x[" << i << "] > " << obj_coeff[i - 1] << endl;
		}

	}//end i

	//init ia, ja, and ar
	int * ia = new int[1 + iaja_size];
	int * ja = new int[1 + iaja_size];
	double * ar = new double[1 + iaja_size];
	assert(ia && ja && ar);

	int ia_id = 1;
	row_id = 1;
	for (auto & f : faces)
	{
		//for (auto j : t1c.second.first)
		//for each row
		for (short d = 0; d < 3; d++)
		{
			ia[ia_id] = row_id;
			ja[ia_id] = d + 1;
			ar[ia_id] = f->n[d]; //coefficient of the constrains from: x*f->n > f->v[0]*f->n , where x is unknown
			ia_id++;
		}//end for j

		row_id++;
	}//end for i

	//additional enclosing constrains
	for (auto& p : enclosing_constrains)
	{
		for (short d = 0; d < 3; d++)
		{
			ia[ia_id] = row_id;
			ja[ia_id] = d + 1;
			ar[ia_id] = p.n[d]; //coefficient of the constrains from: x*f->n > f->v[0]*f->n , where x is unknown
			ia_id++;
		}//end for j

		row_id++;
	}

	glp_load_matrix(lp, iaja_size, ia, ja, ar);

	glp_smcp parm;
	glp_init_smcp(&parm);
	parm.msg_lev = GLP_MSG_OFF;
	parm.tm_lim = 1000;
	int err = glp_simplex(lp, &parm);
	//cout << "err=" << err << endl;

	double z = glp_get_obj_val(lp);
	//cout << "objective value=" << z << endl;

	//get mip status
	int glp_prim_stat = glp_get_status(lp);

	//switch (glp_prim_stat)
	//{
	//	//case GLP_OPT: cout << "solution is optimal;" << endl; break;
	//	//case GLP_FEAS: cout << "solution is feasible;" << endl; break;
	//	case GLP_INFEAS: cout << "solution is infeasible;" << endl;   break;
	//	case GLP_NOFEAS: cout << "problem has no feasible solution;" << endl;  break;
	//	case GLP_UNBND: cout << "problem has unbounded solution;" << endl;  break;
	//	case GLP_UNDEF: cout << "solution is undefined." << endl;  break;
	//}

	bool solution_found = glp_prim_stat == GLP_OPT || glp_prim_stat == GLP_FEAS;

	if (solution_found)
	{
		solution[0] = glp_get_col_prim(lp, 1);		solution[1] = glp_get_col_prim(lp, 2);
		solution[2] = glp_get_col_prim(lp, 3);
		if (isfinite(solution[0]) == false && isfinite(solution[1]) == false && isfinite(solution[2]) == false)
		{
			solution_found = false;
		}
	}

	glp_delete_prob(lp);
	delete[] ia;
	delete[] ja;
	delete[] ar;

	return solution_found;
}

//
// find the best position to move v so that it is close to the center of v's neighbor using QP
//
bool SolveQP(list<cd_f*>& faces, cd_v * v, const list<hull_plane>& enclosing_constrains, Point3d& solution)
{
    vector<Point3d> nei_vs; //neighbor vertices
    Point3d c_pos(0,0,0);   //center pos

    for (auto e: v->edges)
    {
        const auto& pos = e->otherv(v)->pos;
        nei_vs.push_back(pos);

        for(int i = 0; i<3; i++)    c_pos[i] += pos[i]; //sum v pos
    }

    for(int i = 0; i<3; i++)    c_pos[i] /= v->edges.size(); // centered v

    //constructs a quadratic program with no variables and no constraints, ready for data to be added
    Program qp (CGAL::SMALLER, false, 0, false, 0);

    //set non-default entries
    const int X = 0;
    const int Y = 1;
    const int Z = 2;

    // objective function:
    // minimize the distance between c_pos and new vertex(unknown location), V
    // ignore square root
    // (V[x] - c_pos[x])^2 + (V[y] - c_pos[y])^2 + (V[z] -c_pos[z]) ^2
    {
        //v_x^2 + v_y^2 + v_z^2 !! specify 3D !!
        qp.set_d(X, X, 3);  qp.set_d(Y, Y, 3);  qp.set_d(Z,Z,3);

        double c[3] = {0, 0, 0};
        double c0 =0;
        for(int i=0; i<3; i++)
        {
            c[i]= -2*c_pos[i];
            c0 += c_pos[i]*c_pos[i];
        }

        //v_x^2 + v_y^2 + v_z^2 !! specify 3D !!
        qp.set_d(X, X, 2);  qp.set_d(Y, Y, 2);  qp.set_d(Z,Z,2);

        //(-2*c_pos[x])*X + (-2*c_pos[y])*Y +(-2*c_pos[z])*Z
        qp.set_c(X, c[0]);  qp.set_c(Y, c[1]);  qp.set_c(Z, c[2]);

        // c_pos[x]*c_pos[x] + c_pos[y]*c_pos[y] + c_pos[z]*c_pos[z]
        qp.set_c0(c0);
    }

    {
				// constraint 1
		    // : (V - v)* F.normal >= 0 (V:unknown vertex, v: current vertex, F:neighbor faces)
        int row_id=0;
        for(auto& f : faces)
        {
            double b = f->n[0]*v->pos[0] + f->n[1]*v->pos[1] + f->n[2]*v->pos[2];

            //F.norm[x]*X + F.norm[y]*Y + F.norm[z]*Z >= b
            qp.set_a(X, row_id, f->n[0]);    qp.set_a(Y, row_id, f->n[1]);    qp.set_a(Z, row_id, f->n[2]);    qp.set_b(row_id, b);
            qp.set_r(row_id, CGAL::LARGER);
            row_id++;
        }
				//constraint 2
		    //additional enclosing constrains
		    for (auto& p : enclosing_constrains)
		    {
		        double b = p.n[0] * p.o[0] + p.n[1] * p.o[1] + p.n[2] * p.o[2];
						qp.set_a(X, row_id, p.n[0]);    qp.set_a(Y, row_id,  p.n[1]);    qp.set_a(Z, row_id,  p.n[2]);    qp.set_b(row_id, b);
						qp.set_r(row_id, CGAL::LARGER);
		        row_id++;
		    }
    }

    //solve the program, using ET as the exact Type
    Solution s = CGAL::solve_quadratic_program(qp, ET());
    assert (s.solves_quadratic_program(qp));

    // output solution

    if(s.is_optimal() || !(s.is_infeasible()))
    {
        if(false) //print obj, var values
        {
            //std::cout << s <<endl; //show solutions.
            std::cout<<"obj val: "<< to_double(s.objective_value()) << endl;

            for(auto var = s.variable_values_begin(); var <s.variable_values_end(); ++var)
                cout<<to_double(*var)<<" ";

            cout<<endl;
        }

        int i=0;
        for(auto var=s.variable_values_begin(); var < s.variable_values_end(); ++var, ++i)
            solution[i] = to_double(*var);

        if(isfinite(solution[0])==false || isfinite(solution[1]) == false || isfinite(solution[2]) == false){
            cout<<"**solution is not finite"<<endl;
            return false;
        }

				// float d1=(c_pos-solution).normsqr();
				// float d2=(c_pos-v->pos).normsqr();
				// cout<<"dist to center ="<<d1<<" before d="<<d2<<((d1>d2)?" further":" closer")<<endl;

        return true; // done
    }
    else
		{
			  //cout<<"! Error: QP failed to find a solution"<<endl;
		    return false;
		}

}
