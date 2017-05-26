#include "GL/shadowMap_glsl.h"
#include "CD3d_draw_GL.h"
#include "measure/CD3d_hull.h"
#include "CD3d_stat.h"
#include "util/CD3d_util.h"
#include "GL/gliLight.h"
#include "util/CD.h"
#include "util/CD3d_trim_hull.h"
#include "io/CD3dexport.h"
#include <GL/gliFont.h>
#include <math.h>
#include <sstream>

///////////////////////////////////////////////////////////////////////////////

extern cd_state state;

// ==== temp ===== //
extern list<Point3d> cubes;
extern double cubeX;
extern double cubeY;
extern double cubeZ;

void drawCubes(list<Point3d>& cubes, double xsize, double ysize, double zsize);

//extern list<list<pair<Point3d,Point3d> > > g_minCutEdges;
// ================ //

///////////////////////////////////////////////////////////////////////////////

list<glmodel> g_glmodels;
CD_RAPID g_rapid;
Point3d g_model_center; //this is closest center of a g_model to state.model_com

typedef Vector<double, 6> FATNESS_STAT;
FATNESS_STAT getFatness(cd_m * m, FATNESS_STAT fatness);
void show_fatness(int i, FATNESS_STAT init_fatness, FATNESS_STAT curr_fatness);
void save_fatness(int i, FATNESS_STAT init_fatness, FATNESS_STAT curr_fatness);

///////////////////////////////////////////////////////////////////////////////
void drawPath( const VL& path, GLenum type )
{
    //draw path
    glBegin(type);
    FOREACH_CV(path)
        const Point3d& p=v->pos;
        glVertex3d(p[0],p[1],p[2]);
    FOREACH_END
    glEnd();
}

///////////////////////////////////////////////////////////////////////////////

inline void drawPolyText()
{
    //char value[128];
	const float xoffset = 3;
}

inline void drawTextInfo()
{
    if(!state.show_txt) return;
    glPushAttrib(GL_CURRENT_BIT);

    //draw reference axis
    glMatrixMode(GL_PROJECTION); //change to Ortho view
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0,20,0,20);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glDisable(GL_LIGHTING);

    //Draw text
    glTranslated(0,20,0);
    //drawPolyText();

    glPopMatrix();

    //pop GL_PROJECTION
    glMatrixMode(GL_PROJECTION); //change to Pers view
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopAttrib();
}

///////////////////////////////////////////////////////////////////////////////

inline Point3d compCenter(cd_f * f){
    Point3d c;
    for( int j=0;j<3;j++ )
        c[j]=(f->v[0]->pos[j]+f->v[1]->pos[j]+f->v[2]->pos[j])/3;
    return c;
}

///////////////////////////////////////////////////////////////////////////////
//
//
//  glmodel
//
//
///////////////////////////////////////////////////////////////////////////////

glmodel::glmodel(cd_m * m)
{
    solid_flat=solid_smooth=wire=weight=openings=hull_flat
		= face_number_gid = vertex_number_gid = hull_wire =
		hull_sample_gid = -1;

    this->m=m;
    com=computeCenter(m->v());
    radius=computeRadius(m->v(),com);

  	hull_model=NULL;
  	hull_completely_covered = false;
    trim_hull_failed=false;

	  if (m->v().size() < 4) hull_completely_covered = true; //not a 3d model...
}

void glmodel::draw_solid_flat()
{
	if( solid_flat==-1 )
        build();

    if (state.show_b_hull) glColor3f(.95, .25, .25);
    else if(state.b_grey ) glColor3f(.75,.75,.75);
    else glColor4f(color[0],color[1],color[2],state.opacity);

	if(state.b_lighting) glEnable(GL_LIGHTING);

	glCallList(solid_flat);
}

void glmodel::draw_solid_smooth()
{
	if(solid_smooth==-1) build_solid_smooth();
  if (state.show_b_hull) glColor3f(.95, .25, .25);
  else if(state.b_grey ) glColor3f(.75,.75,.75);
  else glColor4f(color[0],color[1],color[2],state.opacity);
	if(state.b_lighting) glEnable(GL_LIGHTING);
    glCallList(solid_smooth);
}

void glmodel::draw_wire()
{
	if( wire==-1 )
		build_wire();
	if(state.b_grey) glColor3f(0,0,0);
	else glColor4f(color[0]/2,color[1]/2,color[2]/2,state.opacity);
	glDisable(GL_LIGHTING);
	glCallList(wire);
}

void glmodel::draw_weight()
{

}


void glmodel::draw_samples()
{
	if (this->hull.getSampledPoints().empty()) return;

	if (hull_sample_gid == -1)
	{
		hull_sample_gid = glGenLists(1);
		glNewList(hull_sample_gid, GL_COMPILE);

		glBegin(GL_POINTS);
		for (auto& pt : this->hull.getSampledPoints())
		{
			glVertex3dv(pt.get());
		}
		glEnd();

		glEndList();
	}

	if (state.b_grey) glColor3f(.75, .75, .75);
	else glColor4f(color[0], color[1], color[2], state.opacity);

	glPointSize(3);
	glDisable(GL_LIGHTING);
	glCallList(hull_sample_gid);
}

void glmodel::draw_openings()
{
	if( openings==-1 ) build_openings();
    //glColor3f(color[0]/5,color[1]/5,color[2]/5);
    if(state.b_grey) glColor3f(0,0,0);
    else glColor3f(0.3f,0.3f,0.5f);
    glDisable(GL_LIGHTING);
    glCallList(openings);
}

void glmodel::draw_hull()
{
	if( hull_flat==-1 ) build_hull();

  if(state.b_grey){
		if (state.show_smooth || state.show_solid)  glColor4f(.75, .75, .75,1);// , 0.25);
		else glColor4f(.75, .75, .75, 1); // state.opacity);
	}
    else{
		if (state.show_smooth || state.show_solid) glColor4f(color[0], color[1], color[2],1);// , 0.5); //show both hull and model
		else glColor4f(color[0], color[1], color[2], 1);// state.opacity);
	}
	if(state.b_lighting) glEnable(GL_LIGHTING);

  glCallList(hull_flat);
}

void glmodel::draw_hull_wire()
{
	if( hull_wire==-1 )
        build_hull_wire();
    if(state.b_grey) glColor3f(0,0,0);
    else glColor4f(color[0]/2,color[1]/2,color[2]/2,state.opacity);
    glDisable(GL_LIGHTING);
    glCallList(hull_wire);
}

void glmodel::draw()
{

}

void glmodel::build_hull()
{
  hull_flat=glGenLists(1);
  glNewList(hull_flat,GL_COMPILE );
  glEnable( GL_POLYGON_OFFSET_FILL );
  glPolygonOffset( 1.0, 1.0 );
	::draw(hull);
  glDisable( GL_POLYGON_OFFSET_FILL );
  glEndList();
}

void glmodel::build_hull_wire()
{
    if( m->v().size()<4 )
    {
        build_wire();
        hull_wire=wire;
        return;
    }

    hull_wire=glGenLists(1);
    glNewList(hull_wire,GL_COMPILE );
    glDisable(GL_LIGHTING);
    glLineWidth(1.0);
	  int tsize = hull.getFaceSize();
    for(int i=0;i<tsize;i++)
    {
        glBegin(GL_LINE_LOOP);
		for (int d = 2; d >= 0; d--) glVertex3dv(hull.getVertex(3 * i + d).get());
        glEnd();
    }
    glEndList();
}

void glmodel::create_simplified_hull()
{

	if (this->hull.getFaceSize() == 0)
	{
		this->hull.buildhull(m);
		if (this->hull.getFaceSize() == 0) return; //nothing get created...
		this->hull_model = this->hull.to_m(-0.01); //thrink the model a bit
		hull_cd_id = g_rapid.buildModel(*this->hull_model);
	}

  if (this->hull.getFaceSize() <= state.hull_simplification_min_face_count) return;

  //if (this->hull.getFaceSize() <= state.hull_simplification_min_face_count) return; //do not simplify any more! already pretty small //yh-original

	int fsize = state.hull_simplification_min_face_count; //max(state.hull_simplification_min_face_count, this->hull.getFaceSize()*state.hull_simplification_face_reduction_ratio);
	auto vol_inc = hull.computeHullVol()*state.hull_simplification_max_vol_inc_ratio;

	cout << "vol_inc=" << vol_inc << "param.hull_simplification_max_vol_inc_ratio=" << state.hull_simplification_max_vol_inc_ratio << endl;

	cd_hull simplified_hull;
	this->hull.simplify(fsize, vol_inc, this->hull_trim_planes, simplified_hull);

	if (this->hull_model != NULL)
	{
		this->hull_model->destroy();
		delete this->hull_model;
		this->hull_model = NULL;
	}

	this->hull_model = simplified_hull.to_m(-0.001);
	hull_cd_id = g_rapid.buildModel(*this->hull_model);

	cout << "hull_cd_id=" << hull_cd_id << endl;
	cout << "simplied hull=" << simplified_hull.getFaceSize() << " hull=" << hull.getFaceSize() << endl;
	this->hull = simplified_hull;

	glDeleteLists(hull_flat,1);
	glDeleteLists(hull_wire, 1);
	hull_flat = -1;
	hull_wire = -1;
}

void glmodel::create_remeshed_hull()
{
    if (this->hull.getFaceSize() == 0)
    {
        if (m->f().size() < 4) return; //too small...
        this->hull.buildhull(m);
        if (this->hull.getFaceSize() == 0) return; //nothing get created...
        this->hull_model = this->hull.to_m(-0.01); //convert hull to model obj
        hull_cd_id = g_rapid.buildModel(*this->hull_model); //for checking overlapping between convert....
    }

    //JML: # of faces should mostly remain the same...
    // if (this->hull.getFaceSize() <= state.hull_remeshing_min_face_count) //yh-copied
    // {
    //     cout<<"(n-1) stop the process: hull_face < min_face "<<endl;
    //     return;
    // }

    int fsize = state.hull_simplification_min_face_count; //max(state.hull_simplification_min_face_count, this->hull.getFaceSize()*state.hull_simplification_face_reduction_ratio); //yh-cmt-max_fsize = max(min_fsize, reducded_fsize)

    auto vol_inc = hull.computeHullVol()*state.hull_simplification_max_vol_inc_ratio; //yh-cmt-max_vol_inc(increased_vol)

    cd_hull remeshed_hull; //yh-cmt-the result hull

    this->hull.remesh(vol_inc, this->hull_trim_planes, remeshed_hull); //yh-cmt-(min_fsize, max_vol_inc, enclosing constraints, result), yh-hyll_trim_planes?

    if (this->hull_model != NULL)
    {
        this->hull_model->destroy();
        delete this->hull_model;
        this->hull_model = NULL;
    }

    this->hull_model = remeshed_hull.to_m(-0.001);

    hull_cd_id = g_rapid.buildModel(*this->hull_model);

    this->hull = remeshed_hull;

    glDeleteLists(hull_flat,1);
    glDeleteLists(hull_wire, 1);
    hull_flat = -1;
    hull_wire = -1;
}

void glmodel::build_solid_flat()
{
	if (solid_flat<0)
    {
        solid_flat=glGenLists(1);
        glNewList(solid_flat,GL_COMPILE );
            //glPushAttrib(GL_CURRENT_BIT);
            const FL& fl=m->f();
            double v[3];
            glEnable( GL_POLYGON_OFFSET_FILL );
            glPolygonOffset( 2.0, 2.0 );
            glBegin(GL_TRIANGLES);
            FOREACH_CF(fl)
                f->n.normalize().get(v);
                glNormal3dv(v);
                for( int D=0;D<3;D++ ){
                    f->v[D]->pos.get(v);
                    glVertex3dv(v);
                }//end D
            FOREACH_END
            glEnd();
            //glPopAttrib();
            glDisable( GL_POLYGON_OFFSET_FILL );
        glEndList();
    }
}


void glmodel::build_solid_smooth()
{
    if(solid_smooth<0)
    {
        VL& vl=m->v();
        double * normals=new double[vl.size()*3];
        int nid=0;
        FOREACH_V(vl)
            Vector3d n;
            FOREACH_E(v->edges)
                for(int i=0;i<2;i++){
                    if( e->f[i]==NULL ) continue;
                    n[0]+=e->f[i]->n[0];
                    n[1]+=e->f[i]->n[1];
                    n[2]+=e->f[i]->n[2];
                }
            FOREACH_END
            n=n.normalize();
            normals[nid*3]=n[0];
            normals[nid*3+1]=n[1];
            normals[nid*3+2]=n[2];
            nid++;
        FOREACH_END

        solid_smooth=glGenLists(1);
        glNewList(solid_smooth,GL_COMPILE );
        {
            const FL& fl=m->f();
            double v[3];
            glEnable( GL_POLYGON_OFFSET_FILL );
            glPolygonOffset( 2.0, 2.0 );
            glBegin(GL_TRIANGLES);
            FOREACH_CF(fl)
                for( int D=0;D<3;D++ ){
                    glNormal3dv(&normals[f->v[D]->id*3]);
                    f->v[D]->pos.get(v);
                    glVertex3dv(v);
                }//end D
            FOREACH_END
            glEnd();
            glDisable( GL_POLYGON_OFFSET_FILL );
        }
        glEndList();

    delete [] normals;
    }
}

void glmodel::build_wire()
{
    if(wire<0)
    {
        wire=glGenLists(1);
        glNewList(wire,GL_COMPILE );
        {
            glLineWidth(1);
            glBegin(GL_LINES);
            FOREACH_E(m->e())
                const Point3d p1=e->v[0]->pos; const Point3d p2=e->v[1]->pos;
                glColor3f(color[0]/4,color[1]/4,color[2]/4);
                glVertex3d(p1[0],p1[1],p1[2]);
                glVertex3d(p2[0],p2[1],p2[2]);
            FOREACH_END
            glEnd();
        }
        glEndList();
    }
}

void glmodel::build_weight()
{
    if(weight<0)
    {
        //find the max conavity
        double max_concavity=0;
        const VL& vl=m->v();
        double * w=new double[vl.size()];
        FOREACH_CV(vl)
            if( v->concavity>max_concavity )
               max_concavity=v->concavity;
        FOREACH_END

        //draw
        weight=glGenLists(1);
        glNewList(weight,GL_COMPILE );
        {
            glDisable(GL_LIGHTING);
            double v[3];
            glEnable( GL_POLYGON_OFFSET_FILL );
            glPolygonOffset( 2.0, 2.0 );
            glBegin(GL_TRIANGLES);
            FOREACH_CF(m->f())
                f->n.get(v);
                for( int D=0;D<3;D++ ){
                    f->v[D]->pos.get(v);
                    double w=(f->v[D]->concavity)/max_concavity;
                    //if( w>0.3 ) w=1;
                    //else w=0;
                    w=(1-w)*1;
                    glColor3d(w,w,w);
                    glVertex3dv(v);
                }//end D
            FOREACH_END
            glEnd();
            glDisable( GL_POLYGON_OFFSET_FILL );
        }
        glEndList();

        delete [] w;
    }
}


void glmodel::build_openings()
{
    if(openings<0)
    {
        //opening boundary
        openings=glGenLists(1);
        glNewList(openings,GL_COMPILE );
        glDisable(GL_LIGHTING);
        glLineWidth(2);
        glBegin(GL_LINES);
        FOREACH_E(m->e())
            if( e->f[1]!=NULL ) continue;
            const Point3d& p1=e->v[0]->pos;
            const Point3d& p2=e->v[1]->pos;
            glVertex3d(p1[0],p1[1],p1[2]);
            glVertex3d(p2[0],p2[1],p2[2]);
        FOREACH_END
        glEnd();
        glEndList();
    }
}

void glmodel::build()
{
    color[0]=drand48();
    color[1]=drand48();
    color[2]=drand48();

    ///////////////////////////////////////////////////////////////////////////
    //draw solid
    build_solid_flat();
    if( solid_smooth!=-1 ){
        glDeleteLists(solid_smooth,1);
        solid_smooth=-1;
    }
    if( wire!=-1 ){
        glDeleteLists(wire,1);
        wire=-1;
    }
    if( weight!=-1 ){
        glDeleteLists(weight,1);
        weight=-1;
    }
}
void glmodel::showVertexNumber()
{
    //if( state.show_hull || state.show_pbridge ) return;

    if(vertex_number_gid==-1)
    {
        vertex_number_gid=glGenLists(1);
        glNewList(vertex_number_gid,GL_COMPILE );

        //typedef list<cd_f*>::const_iterator FIT;
        char value[10]; //string for storing value
        glColor3d(1,0,0);
        double gap=state.model_radus/50;
        glLineWidth(0.5);

        FOREACH_V(m->v())
            sprintf(value,"%d",v->id);
            Point3d p=v->pos;
            //normal direction
            Vector3d n;
            FOREACH_E(v->edges)
                n=n+e->f[0]->n;
                if(e->f[1]!=NULL) n=n+e->f[1]->n;
                //cd_v * ov=e->otherv(v);
                //n=n+(v->pos-ov->pos).normalize();
            FOREACH_END
            n=n.normalize();
            //draw
            glBegin(GL_LINES);
            glVertex3d(p[0],p[1],p[2]);
            p=p+(n*gap);
            glVertex3d(p[0],p[1],p[2]);
            glEnd();
            p=p+n*(gap/3);
            drawstr(p[0],p[1],p[2],value);
        FOREACH_END

        glEndList();
    }

    glCallList(vertex_number_gid);
}

void glmodel::showFaceNumber()
{
    if( state.show_hull ) return;

    if(face_number_gid==-1)
    {
        face_number_gid=glGenLists(1);
        glNewList(face_number_gid,GL_COMPILE );

        //typedef list<cd_f*>::const_iterator FIT;
        char value[10]; //string for storing value
        glColor3d(0,0,0);
        double gap=state.model_radus/50;
        glLineWidth(0.5);

        FOREACH_F(m->f())
            sprintf(value,"%d",f->id);
            Point3d p=compCenter(f);
            glBegin(GL_LINES);
            glVertex3d(p[0],p[1],p[2]);
            p=p+(f->n*gap);
            glVertex3d(p[0],p[1],p[2]);
            glEnd();
            p=p+(f->n*(gap/3));
            drawstr(p[0],p[1],p[2],value);
        FOREACH_END

        glEndList();
    }

    glCallList(face_number_gid);
}

void glmodel::separate()
{
    Vector3d v=(com-state.model_com);
    double dist=v.norm();
    if( dist<=1e-5 ) return;
    Vector3d pv=v/dist; //normalize
    translate=pv*0.1;//*radius/2;
}

//save hull model to file
void glmodel::save_hull(const string& filename)
{
	ofstream fout(filename.c_str());
	if (fout.good())
	{
		if (hull.getFaceSize() != 0)
		{
			cd_m * m = hull.to_m();
			saveObj(m, fout);
			m->destroy();
			delete m;
		}
	}
	fout.close();
}

void glmodel::getMissclassifiedSamples(const hull_plane& plane, vector<Point3d>& missed) const
{
  int count=0;
  Point3d missed_sample;
  for (auto& pt : hull.getSampledPoints())
  {
    if( (pt-plane.o)*plane.n<0){ missed.push_back(pt); count++; missed_sample=pt; }
  }

  if(count>0)
  {
    cd_hull residual_hull;
    //compute the convex object
    hull_plane tmp_plane=plane;
    tmp_plane.n=-tmp_plane.n;
    //compute intersection
    list<hull_plane> planes;
    planes.push_back(tmp_plane);
    trim_hull_by_planes(this->hull, planes, residual_hull,missed_sample);
    vector<Point3d> pts;
    residual_hull.to_pvector(pts);
    for(auto& pt:pts)
    {
      missed.push_back(pt);
    }
  }

}

///////////////////////////////////////////////////////////////////////////////
//
//
// functions that operate on the glmodel
//
//
//
///////////////////////////////////////////////////////////////////////////////
void rebuild_glmodel_hulls()
{
	g_rapid.destroy(); //clean up

	for (auto ig = g_glmodels.begin(); ig != g_glmodels.end(); ig++)
	{
		if (ig->hull_completely_covered) continue;

		ig->hull_trim_planes.clear();
    ig->trim_hull_failed=true;

		if (ig->hull_model!=NULL)
		{
			ig->hull_model->destroy();
			delete  ig->hull_model;
			ig->hull_model = NULL;
		}

		cd_hull newhull;
		newhull.buildhull(ig->m);

		//cout << "ig->m vsize=" << ig->m->v().size() << " newhull fsize=" << newhull.getFaceSize() << endl;
		ig->hull = newhull;
		ig->hull_model = ig->hull.to_m(-0.001);

		//cout << "ig->hull_model vsize=" << ig->hull_model->v().size() << endl;

		ig->hull_cd_id = g_rapid.buildModel(*ig->hull_model);

		glDeleteLists(ig->hull_flat, 1);
		glDeleteLists(ig->hull_wire, 1);
		ig->hull_flat = -1;
		ig->hull_wire = -1;
	}
}

void print_fatness()
{
    double tmp[]={0,0,DBL_MAX,0,0,0}; // fsize, total_area, min, max, sum, weigthed sum
    FATNESS_STAT curr_fatness(tmp);

    for (auto ig = g_glmodels.begin(); ig != g_glmodels.end(); ig++)
    {
        curr_fatness = getFatness(ig->hull_model, curr_fatness);
    }

    cout<<setw(3)<<"#i"<<setw(15)<<"Min"<<setw(15)<<"Max"<<setw(15)<<"Avg"<<setw(15)<<" Weigted_Avg"<<endl;
    cout<<setw(3)<<0<<setw(15)<<curr_fatness[2]<<setw(15)<<curr_fatness[3]<<setw(15)<<curr_fatness[4]/curr_fatness[0]<<setw(15)<<curr_fatness[5]/curr_fatness[1]<<endl;
}

void simplify_glmodel_hulls()
{
  	cout << "------------simplify_glmodel_hulls-----------------------" << endl;
	for (auto ig = g_glmodels.begin(); ig != g_glmodels.end(); ig++)
	{
        ig->create_simplified_hull();
	}
  cout << "-------------------- done -------------------------------" << endl;
}

void remesh_glmodel_hulls(int i)
{
    double tmp[]={0,0,DBL_MAX,0,0,0}; // fsize, total_area, min, max, sum, weigthed sum
    FATNESS_STAT init_fatness(tmp);
    FATNESS_STAT curr_fatness(tmp);

    for (auto ig = g_glmodels.begin(); ig != g_glmodels.end(); ig++)
    {
        if (ig->hull_completely_covered) continue;

        init_fatness = getFatness(ig->hull_model, init_fatness);    //before
        ig->create_remeshed_hull(); //remesh
        curr_fatness = getFatness(ig->hull_model, curr_fatness);    //after
    }


    if(true)
    {
        if(i==0)    cout<< "saving fatness data"<<endl;
        save_fatness(i, init_fatness, curr_fatness);
    }

    //if(state.b_show_hull_fatness){
    if(false)
    {
        show_fatness(i, init_fatness, curr_fatness);
    }
    else    cout<<"."; //loading
}

//r/R = cosA + cosB + cosC -1
inline double two_balls_fatness(cd_f * f)
{
    const Point3d& p1 = f->v[0]->pos;
    const Point3d& p2 = f->v[1]->pos;
    const Point3d& p3 = f->v[2]->pos;

    Vector3d v1 = p2 - p1;
    Vector3d v2 = p3 - p1;
    Vector3d v3 = p3 - p2;

    double cosA = v1*v2 / (v1.norm()*v2.norm());
    double cosB = (-v1)*v3 / (v1.norm()*v3.norm());
    double cosC = (-v2)*(-v3) / (v2.norm()*v3.norm());

    return (cosA+cosB+cosC-1);  //fatness, max fatness=0.5
}

//return fatness = fsize, area, min, max, sum, weighted sum (NOT avg!!)
  inline FATNESS_STAT computeFatness(cd_m * m)
{
    double tmp[]={0,0,DBL_MAX,0,0,0};
    FATNESS_STAT result(tmp);

    for(cd_f * f: m->f())
    {
        double fatness = two_balls_fatness(f);
        result[1] += f->area(); //area
        result[4] += fatness;   //sum
        result[5] += f->area()*fatness; //weighted sum
        if(fatness < result[2])  result[2] = fatness;   //min
        if(fatness > result[3])  result[3] = fatness;   //max
    }

    result[0] = m->f().size();  //# of face

    return result;
}

//return accumulated fatness = fsize, area, min, max, sum, weighted sum (NOT avg!!)
FATNESS_STAT getFatness(cd_m * m, FATNESS_STAT tot_fatness)
{
    FATNESS_STAT ig_fatness = computeFatness(m);

    tot_fatness[0] += ig_fatness[0];    //# of face
    tot_fatness[1] += ig_fatness[1];    //area
    if(ig_fatness[2] < tot_fatness[2]) tot_fatness[2] = ig_fatness[2];    //min
    if(ig_fatness[3] > tot_fatness[3]) tot_fatness[3] = ig_fatness[3];    //max
    tot_fatness[4] += ig_fatness[4];    //sum
    tot_fatness[5] += ig_fatness[5];    //weighted sum

    return tot_fatness;
}

void show_fatness(int i, FATNESS_STAT init_fatness, FATNESS_STAT curr_fatness){
    if(i==0)
    {
        cout<<"#iterations: "<<state.hull_simplification_iteration<<" volume_inc_ratio :" <<state.hull_simplification_max_vol_inc_ratio<<endl;
        cout<<setw(3)<<"#i"<<setw(15)<<"Min"<<setw(15)<<"Max"<<setw(15)<<"Avg"<<setw(15)<<" Weigted_Avg"<<endl;
        cout<<setw(3)<<i<<setw(15)<<init_fatness[2]<<setw(15)<<init_fatness[3]<<setw(15)<<init_fatness[4]/init_fatness[0]<<setw(15)<<init_fatness[5]/init_fatness[1]<<endl;
    }
    cout<<setw(3)<<i+1<<setw(15)<< curr_fatness[2]<<setw(15)<<curr_fatness[3]<<setw(15)<<curr_fatness[4]/curr_fatness[0]<<setw(15)<<curr_fatness[5]/curr_fatness[1]<<endl;
}

void save_fatness(int i, FATNESS_STAT init_fatness, FATNESS_STAT curr_fatness)
{
    stringstream ss;
    auto name = state.str_input.front().substr(0, state.str_input.front().find_last_of("."));
    ss << name << "_fatness_all" << ".txt";
    auto file = ss.str().c_str();

    ofstream fout;
    fout.open(ss.str().c_str(), ios::app);

    if(fout.good())
    {
        if(i==0)
        {
            cout<< ss.str() <<endl;
            fout<<"#iterations: "<<state.hull_simplification_iteration<<" volume_inc_ratio :" <<state.hull_simplification_max_vol_inc_ratio<<"\n";
            fout<<setw(3)<<"#i"<<setw(15)<<"Min"<<setw(15)<<"Max"<<setw(15)<<"Avg"<<setw(15)<<" Weigted_Avg\n";
            fout<<setw(3)<<i<<setw(15)<<init_fatness[2]<<setw(15)<<init_fatness[3]<<setw(15)<<init_fatness[4]/init_fatness[0]<<setw(15)<<init_fatness[5]/init_fatness[1]<<"\n";
        }
        fout<<setw(3)<<i+1<<setw(15)<< curr_fatness[2]<<setw(15)<<curr_fatness[3]<<setw(15)<<curr_fatness[4]/curr_fatness[0]<<setw(15)<<curr_fatness[5]/curr_fatness[1]<<"\n";
        if(i==state.hull_simplification_iteration-1)    fout<<"\n\n";
    }
    fout.close();
}

void save_gmodel_hulls()
{
	auto name = state.str_input.front().substr(0, state.str_input.front().find_last_of("."));
	int id = 0;
	for (auto ig = g_glmodels.begin(); ig != g_glmodels.end(); ig++)
	{
		if (ig->hull_completely_covered) continue; //no need to anything here
    //if (ig->trim_hull_failed) continue; //failed, continue
		if (ig->hull.getFaceSize() == 0) continue; //no hull to save...

		stringstream ss;
		ss << name << "_hull_" << setfill('0') << setw(3)<< id++ << ".obj";
		cout << "- Save file: " << ss.str() << endl;
		ig->save_hull(ss.str());
	}
}

#include "hull_vol.h"
list<hull_plane>  all_trim_planes;
extern vector<Point3d> all_sur_pts;
extern vector<Point3d> all_sv_A;
extern vector<Point3d> all_sv_B;

//use csvm to trim convex hulls
void trim_glmodel_hulls(TRIM_METHOD method)
{
  all_trim_planes.clear();
  all_sur_pts.clear();
  all_sv_A.clear();
  all_sv_B.clear();

  switch(method)
  {
    case SVM_Trim:       cout<<"- Trim Convex Hulls using SVM"<<endl; break;
    case Heuristic_Trim: cout<<"- Trim Convex Hulls using Heuristics"<<endl; break;
    case Dualspace_Trim: cout<<"- Trim Convex Hulls using Dual space arrangment"<<endl; break;
    case Exact_Trim:     cout<<"- Trim Convex Hulls using Exact Translation"<<endl; break;
    default: cerr<<"! Error: Unknown Trim method"<<endl;
  }

	//map< glmodel *, list<hull_plane> > hull2trim_planes;
	double hull_vol_sum = 0;
  double hull_vol_min = FLT_MAX;
	int total_sample_size = state.hull_trim_vol_sample_size;

	//
	// build the convex hull for each decomposed model
	//
	for (auto ig = g_glmodels.begin(); ig != g_glmodels.end(); ig++)
	{
		if (ig->hull_completely_covered) continue;
		if (ig->hull.getFaceSize() == 0)
		{
			ig->hull.buildhull(ig->m);
			ig->hull_model = ig->hull.to_m(-0.001);
			ig->hull_cd_id = g_rapid.buildModel(*ig->hull_model);
		}
    ig->trim_hull_failed=false;
    double vol=ig->hull.computeHullVol();
		hull_vol_sum += vol;
    if(vol<hull_vol_min) hull_vol_min=vol;
	}


	//
	// sample points inside the convex hulls
	//
	for (auto ig = g_glmodels.begin(); ig != g_glmodels.end(); ig++)
	{
		if (ig->hull_completely_covered) continue;
    ig->miss_classified_points.clear();
    ig->hull.clearSamplePoints();

		//the number of points is proportion to the volume of the component
		int k = (int)ceil(ig->hull.computeHullVol() *total_sample_size *1.0f/ hull_vol_sum);
		if (k < 3) k = 3;
		ig->hull.samplePoints(k, state.hull_trim_surface_sample_coverage);
	}

	//check if this hull is completely inside the union of the other hulls
	for (auto ig = g_glmodels.begin(); ig != g_glmodels.end(); ig++)
	{
		if (ig->hull_completely_covered) continue;

		ig->hull_completely_covered = true;

		for (const auto& pt : ig->hull.getSampledPoints())
		{
			bool pt_is_outside_all_hulls=true;
			for (auto ig2 = g_glmodels.begin(); ig2 != g_glmodels.end(); ig2++)
			{
				if (ig == ig2) continue;
				if (ig2->hull_completely_covered) continue;

				if (ig2->hull.isInsideHull(pt) )
				{
					pt_is_outside_all_hulls = false;
					break;
				}
			}//end for ig2

			if (pt_is_outside_all_hulls)
			{
				ig->hull_completely_covered = false;
				break;
			}
		}//end for pt

	}//end for ig

	//
	// find intersecting convex hull pairs and, for each pair, find a separating plane
	//

	for (auto ig = g_glmodels.begin(); ig != g_glmodels.end(); ig++)
	{
		if (ig->hull_completely_covered) //no need to trim this hull
		{
			cd_hull empty_hull;
			ig->hull = empty_hull;
			//reset the display list ids for the hull
			glDeleteLists(ig->hull_flat, 1);
			glDeleteLists(ig->hull_wire, 1);
			ig->hull_flat = -1;
			ig->hull_wire = -1;
			continue;
		}

		//get pairs of intersecting hulls
		auto ig2 = ig;
		ig2++;
		for (; ig2 != g_glmodels.end(); ig2++)
		{
			if (ig2->hull_completely_covered) continue; //no need to trim this hull

			//check if model and this model are in collision.
			//cout << "ig->hull_cd_id=" << ig->hull_cd_id << ", ig2->hull_cd_id=" << ig2->hull_cd_id << endl;

			if (g_rapid.isInCollision(ig->hull_cd_id, ig2->hull_cd_id) == false) continue; //no overlapping, good

			//overlapping convex hulls. need to be trimmed.
			hull_plane trim_plane;
      bool result=false;
      switch(method)
      {
        case SVM_Trim:       result=trim_overlapping_hulls_svm(ig->hull, ig2->hull, trim_plane, state.hull_trim_svm_C); break;
        case Heuristic_Trim: result=trim_overlapping_hulls_heuristic(ig->hull, ig2->hull, trim_plane); break;
        case Dualspace_Trim: cout<<"! Error: Not implemented: Trim Convex Hulls using Dual space arrangment"<<endl; break;
        case Exact_Trim:     result=trim_overlapping_hulls_svm_opt(ig->hull, ig2->hull, trim_plane, state.hull_trim_svm_C); break;
        default: cerr<<"! Error: Unknown Trim method"<<endl;
      }

      if(!result)
      {
        ig->trim_hull_failed = true;
        ig2->trim_hull_failed = true;
        continue;
      }

			all_trim_planes.push_back(trim_plane);
			ig->hull_trim_planes.push_back(trim_plane);
      vector<Point3d> ig_missed;

      // cout<<"!!!!!!!!!!!!!!!!"<<endl;
      // cout<<"compute_volume="<<compute_volume(ig->hull,trim_plane)<<endl;
      // cout<<"compute_volume_explicit="<<compute_volume_explicit(ig->hull,trim_plane)<<endl;
      // cout<<"!!!!!!!!!!!!!!!!"<<endl;

      ig->getMissclassifiedSamples(trim_plane, ig_missed);
      //cout<<"ig_missed size="<<ig_missed.size()<<endl;

			//flip the normal and add to the other hull
			trim_plane.n = -trim_plane.n;
			ig2->hull_trim_planes.push_back(trim_plane);
      vector<Point3d> ig2_missed;
      ig2->getMissclassifiedSamples(trim_plane, ig2_missed);
      //cout<<"ig2_missed size="<<ig2_missed.size()<<endl;

      //remember the missclassified samples
      ig->miss_classified_points.push_back(ig2_missed);
      ig2->miss_classified_points.push_back(ig_missed);


			//return;
		}//end for ig2

	}//end for ig

	///
	/// reconstruct the convex hull including the separating planes
	///
	//for (auto& trims : hull2trim_planes)
	for (auto ig = g_glmodels.begin(); ig != g_glmodels.end(); ig++)
	{
    if (ig->hull_completely_covered) continue;

		//reset the display list ids for the hull
		glDeleteLists(ig->hull_flat, 1);
		glDeleteLists(ig->hull_wire, 1);
		ig->hull_flat = -1;
		ig->hull_wire = -1;

		if (ig->hull_trim_planes.empty() && ig->trim_hull_failed){
      cd_hull empty_hull;
      ig->hull=empty_hull;
      continue;
    }

    //trim the convex hull of this gmodel by the trim_planes
		cd_hull trimed_hull;
		bool r=trim_hull_by_planes(ig->hull, ig->hull_trim_planes, trimed_hull, state.hull_trim_vertex_collapse_dist);

    //consider missclassified points
    //BUG: this will expand to other hulls that were not in collision before
    //cd_hull expanded_hull;
    //auto vol_inc = trimed_hull.computeHullVol()*0.3;
    //r=expand_hull_by_points(trimed_hull, ig->hull_trim_planes, ig->miss_classified_points, expanded_hull, vol_inc);

    //save the result
    //ig->hull = expanded_hull;
    ig->hull = trimed_hull;
    ig->trim_hull_failed = false; //done
	}
}

///////////////////////////////////////////////////////////////////////////////
//
//
//
//
//
//
///////////////////////////////////////////////////////////////////////////////

void draw(cd_hull& hull)
{
    glBegin(GL_TRIANGLES);

    int tsize=hull.getFaceSize();

    for(int i=0;i<tsize;i++)
    {
        glNormal3dv(hull.getTriNormal(i).get());
        for(int d=2;d>=0;d--) glVertex3dv(hull.getVertex(3*i+d).get());
    }
    glEnd();
}

void drawEnvBox()
{
    const float r_scale=10;
    const double r=state.model_radus*r_scale;
    const double left_x=state.model_com[0]-r;
    const double right_x=state.model_com[0]+r;
    const double min_y=state.bbox[2]; //-r/10; //-state.model_radus;
    const double y_plus=state.model_com[1]+r;
    const double z=state.model_com[2]-r;
    const double z_plus=state.model_com[2]+r;

    if(state.b_lighting)
	{
		glEnable(GL_LIGHTING);
		glDisable(GL_LIGHT1);
		glDisable(GL_LIGHT2);
	}

    //draw box...
    glColor3f(1,1,1);

    glBegin(GL_QUADS);

    //floor
    glNormal3f(0,1,0);
    glVertex3d(left_x,min_y,z);
    glVertex3d(left_x,min_y,z_plus);
    glVertex3d(right_x,min_y,z_plus);
    glVertex3d(right_x,min_y,z);

    //back
    glNormal3f(0,0,1);
    glVertex3d(left_x,min_y,z);
    glVertex3d(right_x,min_y,z);
    glVertex3d(right_x,y_plus,z);
    glVertex3d(left_x,y_plus,z);

    //left
    glNormal3f(1,0,0);
    glVertex3d(left_x,min_y,z);
    glVertex3d(left_x,y_plus,z);
    glVertex3d(left_x,y_plus,z_plus);
    glVertex3d(left_x,min_y,z_plus);

    //right
    glNormal3f(-1,0,0);
    glVertex3d(right_x,min_y,z);
    glVertex3d(right_x,min_y,z_plus);
    glVertex3d(right_x,y_plus,z_plus);
    glVertex3d(right_x,y_plus,z);

    //front
    glNormal3f(0,0,-1);
    glVertex3d(left_x,min_y,z_plus);
    glVertex3d(left_x,y_plus,z_plus);
    glVertex3d(right_x,y_plus,z_plus);
    glVertex3d(right_x,min_y,z_plus);

    //top
    glNormal3f(0,-1,0);
    glVertex3d(left_x,y_plus,z);
    glVertex3d(right_x,y_plus,z);
    glVertex3d(right_x,y_plus,z_plus);
    glVertex3d(left_x,y_plus,z_plus);

    //
    glEnd();
    //glFlush();

	glEnable(GL_LIGHT1);
	glEnable(GL_LIGHT2);
}


void buildglmodels(const ML& ms)
{
    //double min_d=1e10; //for compute g_model_center
    FOREACH_CM(ms)
        glmodel model(m);
        g_glmodels.push_back(model);
        /*
        double d=(model.com-state.model_com).normsqr();
        if( d<min_d ){
            g_model_center=model.com;
            d=min_d;
        }
        */
    FOREACH_END

    rebuild_glmodel_hulls();
}

//<--- tmp
//void DrawLevelSets();


void draw()
{
    glPushMatrix();

    glTranslatef(state.usrModelTranslate[0],state.usrModelTranslate[1],state.usrModelTranslate[2]);

	typedef list<glmodel>::iterator GIT;

	for( GIT ig=g_glmodels.begin();ig!=g_glmodels.end();ig++ )
	{
		ig->translate.set(0,0,0);
		if( state.show_apart ) ig->separate();
	}

	if( state.show_smooth )
	{
		for( GIT ig=g_glmodels.begin();ig!=g_glmodels.end();ig++ ) ig->draw_solid_smooth();
    }
    else if( state.show_solid ){
		for( GIT ig=g_glmodels.begin();ig!=g_glmodels.end();ig++ ) ig->draw_solid_flat();
    }

	if (state.show_wire && !state.show_hull && !state.show_b_hull)
	{
		for( GIT ig=g_glmodels.begin();ig!=g_glmodels.end();ig++ ) ig->draw_wire();
    }

    if( state.show_b_hull ){
		for( GIT ig=g_glmodels.begin();ig!=g_glmodels.end();ig++ ) ig->draw_hull();
    }

    if( state.show_b_hull && state.show_wire ){
		for( GIT ig=g_glmodels.begin();ig!=g_glmodels.end();ig++ ) ig->draw_hull_wire();
    }

	//show samples
	if (state.b_show_hull_samples)
	{
		for (GIT ig = g_glmodels.begin(); ig != g_glmodels.end(); ig++) ig->draw_samples();
	}

    //
    if( state.show_face_number ){
        glDisable(GL_LIGHTING);
        for( GIT ig=g_glmodels.begin();ig!=g_glmodels.end();ig++ ) ig->showFaceNumber();
    }

    //
    if( state.show_vertex_number ){
        glDisable(GL_LIGHTING);
        for( GIT ig=g_glmodels.begin();ig!=g_glmodels.end();ig++ ) ig->showVertexNumber();
    }

    //drawDebug();

    glPopMatrix();


    //return;
    // ================== temp =========================

	/*
    typedef list<list<pair<Point3d,Point3d> > >::iterator LLPIT;
    typedef list<pair<Point3d,Point3d> >::iterator LPIT;

    glColor3f(1, 0, 1);
    for(LLPIT llpit = g_minCutEdges.begin(); llpit!=g_minCutEdges.end(); ++llpit)
    {
    	glBegin(GL_LINES);
    	for(LPIT lpit = llpit->begin(); lpit != llpit->end(); ++lpit)
    	{
    		glVertex3dv(lpit->first.get());
    		glVertex3dv(lpit->second.get());
    	}
    	glEnd();
    }
    */




    //if(state.show_solid)
    //{
    //	glPopMatrix();
    //	return;
    //}

    //cout<<"after return"<<endl;


    // =================================================


    //glPopMatrix();
}

inline void rebuild(const ML& models)
{
    g_glmodels.clear();
	buildglmodels(models);
}

///////////////////////////////////////////////////////////////////////

inline void draw_selected()
{
    glDisable(GL_LIGHTING);

    ////////////////////////////////////////////////////////////////////
    if( state.selectedVID>0 )
    {
        bool update=false;
        if(state.selectedV==NULL) update=true;
        else if(state.selectedV->id!=state.selectedVID) update=true;

        if(update) //find the vertex...
        {
            for(list<glmodel>::iterator i=g_glmodels.begin();i!=g_glmodels.end();i++)
            {
                cd_m * m=i->m;
                FOREACH_V(m->v())
                    if(v->id==state.selectedVID){
                        state.selectedV=v;
                        update=false;
                        break;
                    }
                FOREACH_END

                if(update==false) break;
            }
        }

        if(state.selectedV!=NULL)
        {
            glDisable( GL_DEPTH_TEST);
            glColor3d(1,1,0);
            Point3d& pos=state.selectedV->pos;
            drawBall_GL(pos,state.model_radus/60);
            glEnable( GL_DEPTH_TEST);
        }
    }

    ////////////////////////////////////////////////////////////////////
//    glDisable(GL_LIGHTING);
//    if( !state.selectedF.empty() ){
//        double gap=0.1;
//        char value[256];
//        glColor3d(0,0,1);
//        FOREACH_F(state.selectedF)
//            sprintf(value,"%d",f->id);
//            Point3d p=compCenter(f);
//            glLineWidth(1);
//            //draw triangle outline
//            glBegin(GL_LINE_LOOP);
//            for(int i=0;i<3;i++){
//                Point3d& pos=f->v[i]->pos;
//                glVertex3d(pos[0],pos[1],pos[2]);
//            }
//            glEnd();
//
//            //draw index
//            glBegin(GL_LINES);
//            glVertex3d(p[0],p[1],p[2]);
//            p=p+(f->n*gap);
//            glVertex3d(p[0],p[1],p[2]);
//            glEnd();
//            p=p+(f->n*(gap/3));
//            drawstr(p[0],p[1],p[2],value);
//        FOREACH_END
//    }
    ////////////////////////////////////////////////////////////////////
}

extern vector<Point3d> all_sur_pts;
extern vector<Point3d> all_sv_A;
extern vector<Point3d> all_sv_B;

///////////////////////////////////////////////////////////////////////
void draw(const list<cd_m*>& models)
{

    if( state.rebuild )
	{
        if(state.rebuild) state.rebuild=false;
		rebuild(models);
    }

    {
		gli::ApplyCameraTrasnform();
		glTranslatef(-state.model_com[0], -state.model_com[1], -state.model_com[2]);

        if(state.b_draw_bbox) drawEnvBox();
        draw();
    }

    draw_selected();

	if (state.b_show_hull_support_vectors)
	{
		glPointSize(5);

		glBegin(GL_POINTS);
		glColor3d(0, 1, 0);
		for (auto &pts : all_sv_A)
		{
			glVertex3dv(pts.get());
		}
		glEnd();

		glBegin(GL_POINTS);
		glColor3d(0, 0, 1);
		for (auto &pts : all_sv_B)
		{
			glVertex3dv(pts.get());
		}
		glEnd();
	}


	//draw trim planes
	if (state.b_show_hull_separators)
	{
		glPointSize(5);
		glBegin(GL_POINTS);
		glColor3d(1, 0, 1);
		for (auto &pts : all_sur_pts)
		{
			glVertex3dv(pts.get());
		}
		glEnd();

    glDisable(GL_CULL_FACE);
		for (auto& plane : all_trim_planes)
		{
			Vector3d u(0.5, 0.5, 0.5);//(drand48() - 0.5, drand48() - 0.5, drand48() - 0.5);
			u = u.normalize();
			auto v = u%plane.n;
			u = plane.n%v;
			glBegin(GL_POLYGON);
			glColor3d(1, 1, 0);
			for (int i = 0; i < 20; i++)
			{
				Point3d pt = plane.o + (u*cos(i*PI2 / 20)*plane.error + v*sin(i*PI2 / 20)*plane.error);
				glVertex3dv(pt.get());
			}
			glEnd();

			glBegin(GL_LINES);
			glColor3d(1, 0, 0);
			glVertex3dv(plane.o.get());
			glVertex3dv((plane.o + plane.n*plane.error).get());
			glEnd();
		}
    glEnable(GL_CULL_FACE);

    for (auto& plane : all_trim_planes)
    {
      Vector3d u(0.5, 0.5, 0.5);//(drand48() - 0.5, drand48() - 0.5, drand48() - 0.5);
      u = u.normalize();
      auto v = u%plane.n;
      u = plane.n%v;
      glBegin(GL_LINE_LOOP);
      glColor3d(0, 0, 0);
      for (int i = 0; i < 20; i++)
      {
        Point3d pt = plane.o + (u*cos(i*PI2 / 20)*plane.error + v*sin(i*PI2 / 20)*(plane.error*1.01));
        glVertex3dv(pt.get());
      }
      glEnd();
    }

	}

    if(state.show_txt) drawTextInfo();
}

void drawReset()
{
    g_glmodels.clear();
}

void randColor()
{
    typedef list<glmodel>::iterator GIT;
    for( GIT ig=g_glmodels.begin();ig!=g_glmodels.end();ig++ ){
        ig->setColor(drand48(),drand48(),drand48());
    }

	//glutPostRedisplay();
}

///////////////////////////////////////////////////////////////////////////////

//This will draw an ellipsoid centered at the origin, and having the "radius"
//around X, Y, and Z given by fA, fB, and fC respectively
//written by "terminate"
void DrawEllipsoid(unsigned int uiStacks, unsigned int uiSlices, float fA, float fB, float fC)
{
    glEnable(GL_NORMALIZE);
    glPushMatrix();
    glScalef(fA, fB, fC);
    glutSolidSphere(1, uiStacks, uiSlices);
    glPopMatrix();
    glDisable(GL_NORMALIZE);
}

void DrawBox(float x1, float x2, float y1, float y2, float z1, float z2)
{
    glBegin(GL_POLYGON);
        glVertex3f(x1, y1, z2);
        glVertex3f(x2, y1, z2);
        glVertex3f(x2, y1, z1);
        glVertex3f(x1, y1, z1);
    glEnd();

    glBegin(GL_POLYGON);
        glVertex3f(x1, y2, z2);
        glVertex3f(x2, y2, z2);
        glVertex3f(x2, y2, z1);
        glVertex3f(x1, y2, z1);
    glEnd();

    glBegin(GL_POLYGON);
        glVertex3f(x1, y1, z2);
        glVertex3f(x2, y1, z2);
        glVertex3f(x2, y2, z2);
        glVertex3f(x1, y2, z2);
    glEnd();

    glBegin(GL_POLYGON);
        glVertex3f(x1, y1, z1);
        glVertex3f(x2, y1, z1);
        glVertex3f(x2, y2, z1);
        glVertex3f(x1, y2, z1);
    glEnd();

    glBegin(GL_POLYGON);
        glVertex3f(x1, y1, z2);
        glVertex3f(x1, y2, z2);
        glVertex3f(x1, y2, z1);
        glVertex3f(x1, y1, z1);
    glEnd();

    glBegin(GL_POLYGON);
        glVertex3f(x2, y1, z2);
        glVertex3f(x2, y2, z2);
        glVertex3f(x2, y2, z1);
        glVertex3f(x2, y1, z1);
    glEnd();
}

///////////////////////////////////////////////////////////////////////////////

void drawE_GL(cd_e * e)
{
    glBegin(GL_LINES);
    Point3d p1=e->v[0]->pos;
    Point3d p2=e->v[1]->pos;
    glVertex3d(p1[0],p1[1],p1[2]);
    glVertex3d(p2[0],p2[1],p2[2]);
    glEnd();
}


void drawVL_GL( const VL& vl, GLenum type)
{
    //draw path
    glBegin(type);
    FOREACH_CV(vl)
        const Point3d& p=v->pos;
        glVertex3d(p[0],p[1],p[2]);
    FOREACH_END
    glEnd();
}

void drawBall_GL( const Point3d& p, double size)
{
    glPushMatrix();
    glTranslated(p[0],p[1],p[2]);
    glutSolidSphere(size,10,10);
    glPopMatrix();
}

void drawPlane_GL
(const Point3d& p, const Vector3d& n, const Vector3d& v1, const Vector3d& v2)
{
    //Draw The Plane
    float scale=state.model_radus/2;

    Point3d p1=p+v1*scale+v2*scale;
    Point3d p2=p+v1*(-scale)+v2*scale;
    Point3d p3=p+v1*(-scale)+v2*(-scale);
    Point3d p4=p+v1*scale+v2*(-scale);

    glBegin(GL_LINE_LOOP);
    glNormal3d(n[0],n[1],n[2]);
    glVertex3d(p1[0],p1[1],p1[2]);
    glVertex3d(p2[0],p2[1],p2[2]);
    glVertex3d(p3[0],p3[1],p3[2]);
    glVertex3d(p4[0],p4[1],p4[2]);
    glEnd();

    glColor4f(0.5f,0.8f,0.3f,0.4f);
    glBegin(GL_QUADS);
    glNormal3d(n[0],n[1],n[2]);
    glVertex3d(p1[0],p1[1],p1[2]);
    glVertex3d(p2[0],p2[1],p2[2]);
    glVertex3d(p3[0],p3[1],p3[2]);
    glVertex3d(p4[0],p4[1],p4[2]);
    glEnd();



    //Draw the Axis
    /*
    Point3d x=p+v1*scale;
    Point3d y=p+v2*scale;
    Point3d z=p+n*scale;
    glLineWidth(1);
    glBegin(GL_LINES);
    glColor3d(1,0,0);
    glVertex3d(p[0],p[1],p[2]);
    glVertex3d(x[0],x[1],x[2]);
    glColor3d(0,1,0);
    glVertex3d(p[0],p[1],p[2]);
    glVertex3d(y[0],y[1],y[2]);
    glColor3d(0,0,1);
    glVertex3d(p[0],p[1],p[2]);
    glVertex3d(z[0],z[1],z[2]);
    glEnd();

    glColor3d(1,0,0);
    drawstr(x[0],x[1],x[2],"u");
    glColor3d(0,1,0);
    drawstr(y[0],y[1],y[2],"v");
    glColor3d(0,0,1);
    drawstr(z[0],z[1],z[2],"n");
*/
}

void drawVinBall_GL( const cd_v * v, double size)
{
    const Point3d& p=v->pos;
    drawBall_GL(p,size);
}

void drawVLinBalls_GL( const VL& vl, double size)
{
    glEnable(GL_LIGHTING);
    FOREACH_CV(vl) drawVinBall_GL(v,size); FOREACH_END
}


void drawEL_GL( const EL& el )
{
    //draw path
    glBegin(GL_LINES);
    FOREACH_CE(el);
        for(int i=0;i<2;i++){
            const Point3d& p=e->v[i]->pos;
            glVertex3d(p[0],p[1],p[2]);
        }
    FOREACH_END
    glEnd();
}

void drawEdges_GL(EL& el)
{
    glDisable(GL_LIGHTING);
    glLineWidth(1);
    glBegin(GL_LINES);
    FOREACH_E(el)
        const Point3d& p1=e->v[0]->pos;
        const Point3d& p2=e->v[1]->pos;
        glVertex3d(p1[0],p1[1],p1[2]);
        glVertex3d(p2[0],p2[1],p2[2]);
    FOREACH_END
    glEnd();
}

void drawFacets_GL(FL& fl)
{
    glEnable(GL_LIGHTING);
    glEnable( GL_POLYGON_OFFSET_FILL );
    glPolygonOffset( 1.0, 1.0 );
    glBegin(GL_TRIANGLES);
    FOREACH_F(fl)
        const Vector3d& n=f->n;
        glNormal3d(n[0],n[1],n[2]);
        for(int j=0;j<3;j++){
            const Point3d& p1=f->v[j]->pos;
            glVertex3d(p1[0],p1[1],p1[2]);
        }
    FOREACH_END
    glEnd();
    glDisable( GL_POLYGON_OFFSET_FILL );
}


void drawCubes(list<Point3d>& cubes, double xsize, double ysize, double zsize)
{
	int csize = cubes.size();
	//cout<<"the cube size is "<<csize<<endl;
	glPushMatrix();
	glColor3f(1.0f, 0.0f, 0.0f);
	for(list<Point3d>::iterator pit = cubes.begin(); pit != cubes.end(); ++pit)
	{
		glPushMatrix();
		Point3d& pd = *pit;
		glTranslated(pd[0], pd[1], pd[2]);
		glScalef(xsize, ysize, zsize);
		glutSolidCube(1);
		glPopMatrix();
	}
	glPopMatrix();

//    glPushMatrix();
//    glTranslated(x, y, z);
//    glutSolidCube(size);
//    //glutSolidSphere(size,10,10);
//    glPopMatrix();

}
