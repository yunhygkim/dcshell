#ifndef _CD_MODEL_H_
#define _CD_MODEL_H_

///////////////////////////////////////////////////////////////////////////////
#include "util/mathtool/Point.h"
#include "util/mathtool/Vector.h"
#include "util/model/IPolygonal.h"
using namespace mathtool;

///////////////////////////////////////////////////////////////////////////////
#include <list>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <set>
#include <assert.h>

using namespace std;

///////////////////////////////////////////////////////////////////////////////
//#include "util/ManageTip/TipPair.h"
//#include "util/CVI/CVIIndex.h"
///////////////////////////////////////////////////////////////////////////////

#include "io/CD3d_msg.h"

///////////////////////////////////////////////////////////////////////////////
#ifdef WIN32
#pragma warning(disable : 4786)
#endif

///////////////////////////////////////////////////////////////////////////////
// Defines

const double SMLL_D=1e-10; // a small double value, used for check refexity

#ifdef flagT
#undef flagT
#endif
#define flagT unsigned int

const flagT INIT_FLAG=0;

///////////////////////////////////////////////////////////////////////////////
//
//  FLAG ID for marking
//
///////////////////////////////////////////////////////////////////////////////

flagT generateID();

///////////////////////////////////////////////////////////////////////////////
//Relative Position

enum CP_POS {CP_ON, CP_ABOVE, CP_BELOW};

///////////////////////////////////////////////////////////////////////////////

struct cd_v;
struct cd_e;
struct cd_f;
class  cd_m;

///////////////////////////////////////////////////////////////////////////////

typedef list<cd_v*>          VL;
typedef vector<cd_v*>		 VV;
typedef list<cd_e*>          EL;
typedef list<cd_f*>          FL;
typedef list<cd_m*>          ML;
typedef vector<Point3d>      PV;
typedef VL::iterator         VIT;
typedef VV::iterator		 VVIT;
typedef EL::iterator         EIT;
typedef FL::iterator         FIT;
typedef ML::iterator         MIT;
typedef VL::const_iterator   CVIT;
typedef VV::const_iterator	 CVVIT;
typedef EL::const_iterator   CEIT;
typedef FL::const_iterator   CFIT;
typedef ML::const_iterator   CMIT;
typedef VL::reverse_iterator RVIT;
typedef EL::reverse_iterator REIT;
typedef FL::reverse_iterator RFIT;
typedef ML::reverse_iterator RMIT;
typedef list<VL>             VLL;
typedef VLL::iterator        VLIT;
typedef VLL::const_iterator  CVLIT;
typedef PV::iterator         PVIT;
typedef PV::const_iterator   CPVIT;

///////////////////////////////////////////////////////////////////////////////

#define FOREACH_V(vl)  {for(VIT iv=vl.begin();iv!=vl.end();iv++)  { cd_v* v=*iv;
#define FOREACH_VV(vv) {for(VVIT iv=vv.begin();iv!=vv.end();++iv)  { cd_v* v=*iv;
#define FOREACH_E(el)  {for(EIT ie=el.begin();ie!=el.end();ie++)  { cd_e* e=*ie;
#define FOREACH_F(fl)  {for(FIT fi=fl.begin();fi!=fl.end();fi++)  { cd_f* f=*fi;
#define FOREACH_M(ml)  {for(MIT im=ml.begin();im!=ml.end();im++)  { cd_m* m=*im;
#define FOREACH_Pt(pv) {for(PVIT ip=pv.begin();ip!=pv.end();ip++){Point3d& p=*ip;
#define FOREACH_Vi(_s,_e) {for(VIT iv=_s;iv!=_e;iv++) { cd_v* v=*iv;
#define FOREACH_CV(vl) {for(CVIT iv=vl.begin();iv!=vl.end();iv++){ cd_v* v=*iv;
#define FOREACH_CVV(vv) {for(CVVIT iv=vv.begin();iv!=vv.end();++iv)  { cd_v* v=*iv;
#define FOREACH_CVi(_s,_e) {for(CVIT iv=_s;iv!=_e;iv++) { const cd_v* v=*iv;
#define FOREACH_CE(el) {for(CEIT ie=el.begin();ie!=el.end();ie++){ cd_e* e=*ie;
#define FOREACH_CF(fl) {for(CFIT fi=fl.begin();fi!=fl.end();fi++){ cd_f* f=*fi;
#define FOREACH_CM(ml) {for(CMIT im=ml.begin();im!=ml.end();im++){ cd_m* m=*im;
#define FOREACH_CPt(pv){for(CPVIT ip=pv.begin();ip!=pv.end();ip++){ const Point3d& p=*ip;
#define FOREACH_RV(vl) {for(RVIT iv=vl.rbegin();iv!=vl.rend();iv++){ cd_v* v=*iv;
#define FOREACH_RE(el) {for(REIT ie=el.rbegin();ie!=el.rend();ie++){ cd_e* e=*ie;
#define FOREACH_RF(fl) {for(RFIT fi=fl.rbegin();fi!=fl.rend();fi++){ cd_f* f=*fi;
#define FOREACH_RM(ml) {for(RMIT im=ml.rbegin();im!=ml.rend();im++){ cd_m* m=*im;
#define FOREACH_VL(vll){for(VLIT ivl=vll.begin();ivl!=vll.end();ivl++){VL& vl=*ivl;
#define FOREACH_CVL(vll){for(CVLIT ivl=vll.begin();ivl!=vll.end();ivl++){const VL& vl=*ivl;
#define FOREACH_END   }}

///////////////////////////////////////////////////////////////////////////////
//vertex for model
struct cd_v {

    cd_v();
    ~cd_v()
    {
    	//assert(false);
    }
    
    //delete the edge from the list. O(n).
    void del(cd_e * e)
    { 
        EL::iterator d=std::find(edges.begin(), edges.end(),(const cd_e*)e); 
        if( d!=edges.end() ) edges.erase(d);
        else cd_err(cd_msg("cd_v::del"));
    }
    
    //find the edge with given v in the list. O(n).
    cd_e * find(cd_v * v) const;

    // get message from this vertex
    cd_msg to_msg();

    ///////////////////////////////////////////////////////////////////////////
    int id;
    int hull_cid; //convex hull vertex cluster ID
    Point3d pos;
    EL edges;

    CP_POS cp_pos;        // on, below or above the cutting plane
    float concavity;      // concavity measurement
    flagT flag;
    cd_v * parent;        // used for path search
    bool b_feature_point; // true if it's feature point(knot)
    int component_id;	  // id of decomposed component

	//////////////////////////////////////////////
	int round;
};

///////////////////////////////////////////////////////////////////////////////
//edge for model
struct cd_f; //face
struct cd_e {
    cd_e();

    //compute if this notch is notch or not. stored in reflex.
    void compRef();
    void compLength(){ length=(float)(v[0]->pos-v[1]->pos).norm(); }
    Vector3d n();    //return normal
    void update();

    ///////////////////////////////////////////////////////////////////////////
    //find the other vertex
    cd_v * otherv(const cd_v * v1) const{
        if( v[0]==v1 ) return v[1];
        if( v[1]==v1 ) return v[0];
        cd_err(cd_msg("cd_v::del"));
        return NULL;
    }

    //find the other face
    cd_f * otherf(const cd_f * f1) const{
        if( f[0]==f1 ) return f[1];
        if( f[1]==f1 ) return f[0];
        cd_err(cd_msg("cd_e::otherf"));
        return NULL;
    }

    //delete the face df
    void del(cd_f * df){ 
        if( f[0]==df ){ f[0]=f[1]; f[1]=NULL; return; }
        if( f[1]==df ){ f[1]=NULL; return;}
        cd_err(cd_msg("cd_e::del"));
    }
    
    //repace the face of by nf
    void replace(cd_f * of, cd_f * nf){ 
        if( of==nf ) return; //same, no need to replace
        {for(int i=0;i<2;i++){ if(f[i]==of){ f[i]=nf; return; } }}
        {for(int i=0;i<2;i++){ if(f[i]==NULL){ f[i]=nf; return; } }}
        cd_err(cd_msg("cd_e::replace f"));
    }

    //replace ov by nv
    void replace(cd_v * ov, cd_v * nv){
        if( v[0]==ov ){ v[0]=nv; return;}
        if( v[1]==ov ){ v[1]=nv; return;}
        cd_err(cd_msg("cd_e::replace v"));
    }

    //check if this edge is colliding with the cutting plane
    bool isColliding(){
        if( v[0]->cp_pos==CP_ON && v[1]->cp_pos==CP_ON )
            return false; //this is on the plane
        if( v[0]->cp_pos==v[1]->cp_pos ) 
            return false;
        return true;
    }

    // get message from this edge
    cd_msg to_msg();

    ///////////////////////////////////////////////////////////////////////////
    cd_v * v[2];    //vertices
    cd_f * f[2];    //faces
    bool reflex;    //is a reflex? (notch, concave)
    float length;   //length
    flagT flag;
};

///////////////////////////////////////////////////////////////////////////////
//face for model
struct cd_f{
    cd_f();
    
    //get i-th neighbor
    cd_f * nei(int i){ return e[i]->otherf(this); }

    void compNormal();

    float area(){
        const Point3d& p1=v[0]->pos;
        const Point3d& p2=v[1]->pos;
        const Point3d& p3=v[2]->pos;
        return (float)((p2-p1)%(p3-p1)).norm();      
    }

    float perimeter(){ return e[0]->length+e[1]->length+e[2]->length; }

    ///////////////////////////////////////////////////////////////////////////
    int find3(cd_e * e1, cd_e * e2){
        for(int i=0;i<3;i++){ if(e[i]!=e1&&e[i]!=e2) return i;}
        cd_err(cd_msg("cd_f::find3"));
        return -1;
    }

    void add(cd_e * ne){
        for(int i=0;i<3;i++){ if(e[i]==NULL){ e[i]=ne; return; } } 
        cd_err(cd_msg("cd_f::add e"));
    }

    int find(cd_e * e1,bool outerr=true){ 
        for(int i=0;i<3;i++ ){ if(e[i]==e1) return i;} 
        if(outerr) cd_err(cd_msg("cd_f::find e"));
        return -1;
    }

    void replace( cd_v * v1, cd_v * v2, bool outerr=true ){
        for(int i=0;i<3;i++){ if(v[i]==v1){ v[i]=v2; return; } } 
        if(outerr) cd_err(cd_msg("cd_f::replace v"));
    }

    void replace( cd_e * e1, cd_e * e2, bool outerr=true ){
        for(int i=0;i<3;i++){ if(e[i]==e1){ e[i]=e2; return; } } 
        if(outerr) cd_err(cd_msg("cd_f::replace e"));
    }

    //find the vertex that is not in e
    cd_v * otherv(const cd_e * e) const{
        for(int i=0;i<3;i++){ 
            if(v[i]!=e->v[0]&&v[i]!=e->v[1]){ return v[i]; } 
        } 
        cd_err(cd_msg("cd_f::otherv"));
        return NULL;
    }

    // get message from this facet
    cd_msg to_msg();

    ///////////////////////////////////////////////////////////////////////////
    int org_id;		//face_id in mesh file
    int id;
    Vector3d n;  //normal   
    cd_v * v[3]; //vertices
    cd_e * e[3]; //edges
    flagT flag;
};

///////////////////////////////////////////////////////////////////////////////
//hole
struct cd_h{
    cd_h(){ length=0; }
    bool operator>(const cd_h& other)const{return length>other.length;}
    VL vl;
    float length;
    Point3d com;
};

typedef list<cd_h>           HL;
typedef HL::iterator         HIT;
typedef HL::const_iterator   CHIT;
typedef HL::reverse_iterator RHIT;
#define FOREACH_H(hl) {for(HIT ih=hl.begin();ih!=hl.end();ih++)  { cd_h& h=*ih;
#define FOREACH_Hi(_s,_e) {for(HIT ih=_s;ih!=_e;ih++)  { cd_h& h=*ih;

///////////////////////////////////////////////////////////////////////////////
class cd_m{
public:
    cd_m();
    void destroy();
    cd_m * clone();
    
    VL& v() { return m_v; }
    EL& e() { return m_e; }
    FL& f() { return m_f; }
    HL& holes() { return m_h; }

    const VL& v() const { return m_v; }
    const EL& e() const { return m_e; }
    const FL& f() const { return m_f; }
    const HL& holes() const { return m_h; }

    //in some cases we know that the model is convex
    //enough so we don't want to decompose it!
    bool isDecomposable() const { return m_decompose_me; }
    void setDecomposable(bool f){ m_decompose_me=f; }

    float getMaxConcavity() const { return m_max_c; }
    void setMaxConcavity(float v) { m_max_c = v; }

private:
    VL  m_v;       //vertices
    EL  m_e;       //edges
    FL  m_f;       //faces
    HL  m_h;       //hole boundaries
    bool m_decompose_me; //default is true
    // added by Zhonghua 10/1/2013
    float m_max_c;   // max_concavity
};

///////////////////////////////////////////////////////////////////////////////
// build models

// create models from points and triangles
void buildModels(const PartVector&, const PtVector&, const TriVector&, ML&);
// create models from Face List
void buildModels(const list<FL>&, ML&);

///////////////////////////////////////////////////////////////////////////////
// identify hole boundaries
void identifyHoles( cd_m * m );

#endif //_CD_MODEL_H_


