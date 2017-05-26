#ifndef _CD3D_DRAW_GL_H_
#define _CD3D_DRAW_GL_H_

//#include "util/model/IPolygonal.h"
#include "GL/gli.h"
#include "CD3d_model.h"
#include "measure/CD3d_hull.h"

class cd3d;
class cd_m;

void draw(const list<cd_m*>& models);
void draw(cd_hull& hull);
void drawReset();
void randColor();
void randTipColor();
void setSelectedColor(float R,float G, float B);
void getInfo(const list<cd_m*>& todo,const list<cd_m*>& done);
void drawPath( const VL& path, GLenum type=GL_LINE_STRIP );
void print_fatness();
void buildglmodels(const ML& ms);
void simplify_glmodel_hulls();
void remesh_glmodel_hulls(int i);

///////////////////////////////////////////////////////////////////////////////
//fuctions to draw primitives
void DrawEllipsoid(unsigned int uiStacks, unsigned int uiSlices, float fA, float fB, float fC);
void DrawBox(float x1, float x2, float y1, float y2, float z1, float z2);
void drawE_GL(cd_e * e);
void drawVL_GL( const VL& vl, GLenum type=GL_LINE_STRIP );
void drawBall_GL( const Point3d& p, double size);
void drawPlane_GL
(const Point3d& p, const Vector3d& n, const Vector3d& v1, const Vector3d& v2);
void drawVinBall_GL( const cd_v * v, double size);
void drawVLinBalls_GL( const VL& vl, double size);
void drawEL_GL( const EL& el );
void drawEdges_GL(EL& el);
void drawFacets_GL(FL& fl);

///////////////////////////////////////////////////////////////////////////////
//functions operate on the glmodels
enum TRIM_METHOD { SVM_Trim, Heuristic_Trim, Dualspace_Trim, Exact_Trim };
void simplify_glmodel_hulls();
void remesh_glmodel_hulls(int i);
void trim_glmodel_hulls(TRIM_METHOD method);
void rebuild_glmodel_hulls();

//////////////////////////////////////////////////////////////////////////
struct plane
{
	Point3d o;
	Vector3d n;
};

///////////////////////////////////////////////////////////////////////////////
class glmodel{
public:

    glmodel(cd_m * m);

    ~glmodel()
	{
        glDeleteLists(solid_flat,1);
        glDeleteLists(wire,1);
        glDeleteLists(weight,1);
        glDeleteLists(solid_smooth,1);
        glDeleteLists(face_number_gid,1);
		glDeleteLists(hull_flat,1);
		glDeleteLists(hull_wire,1);

		if (hull_model != NULL)
		{
			hull_model->destroy();
			delete hull_model;
		}
    }

    void draw();
	void draw_solid_flat();
    void draw_solid_smooth();
    void draw_wire();
    void draw_weight();
    void draw_openings();
	void draw_hull();
	void draw_hull_wire();
	void draw_samples();

	void rm_draw(ostream& out);  //renderman draw
	void rm_weight(ostream& out);

	void rm_solid_flat(ostream& out);
	void rm_hull(ostream& out);
	void rm_wire(ostream& out);


    void setColor(float r,float g,float b)
	{
        color[0]=r; color[1]=g; color[2]=b;
    }

    cd_m * m;

    void build();
    void build_solid_flat();
    void build_solid_smooth();
    void build_wire();
    void build_weight();
    void build_openings();
	void build_hull();
	void build_hull_wire();
	void build_simplified_hull();
	void build_simplified_hull_wire();

	void create_simplified_hull();
	void create_remeshed_hull();

    void showFaceNumber();
    void showVertexNumber();
    void separate();

    Point3d  com;
    Vector3d translate;
    double radius;


	//GL display ids
    int solid_flat;
    int solid_smooth;
    int face_number_gid;
    int vertex_number_gid;
    int wire;
    int weight;       // paint the vertices with weight
    int openings;
	int hull_flat;
	int hull_wire;
	int hull_sample_gid;

    float color[3];

	//convex hull related
	cd_hull hull;
	cd_m* hull_model;
	int hull_cd_id;     //id for checking overlapping between convex hulls
	list<hull_plane> hull_trim_planes; //list of planes that trims the convex hull
	bool hull_completely_covered;
	bool trim_hull_failed;

	//list of points from other hulls missclassified to this hull
	//one list<Point3d> per hull
	vector< vector<Point3d> > miss_classified_points;

	void save_hull(const string& filename);
	void getMissclassifiedSamples(const hull_plane& plane, vector<Point3d>& missed) const;
};


///////////////////////////////////////////////////////////////////////////////
inline int BuildGLModel
(const PtVector& ptV, const TriVector& triV)
{
    //build GL model
    int GLID=glGenLists(1);
    glNewList(GLID,GL_COMPILE);

    glBegin(GL_TRIANGLES);
    int size=triV.size();
    double p[3][3]; double nv[3];
    for( int i=0;i<size;i++ ){
        const Tri & tri=triV[i];
        int id1=tri[0]; int id2=tri[1]; int id3=tri[2];
        Vector3d n=(ptV[id2]-ptV[id1])%(ptV[id3]-ptV[id1]);
        n.normalize().get(nv);
        glNormal3dv(nv);
        for( int D=0;D<3;D++ ){
            (ptV[tri[D]]).get(p[D]);
            glVertex3dv(p[D]);
        }//end D
    }//end i

    glEnd();
    glEndList();

    return GLID;
}

#endif //_CD3D_DRAW_GL_H_
