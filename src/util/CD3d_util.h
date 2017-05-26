#ifndef _MODEL_UTILITY_H_
#define _MODEL_UTILITY_H_
///////////////////////////////////////////////////////////////////////////////
//Global Functions
///////////////////////////////////////////////////////////////////////////////

#include "CD3d_model.h"
#include <functional>
#include <cassert>
#include <algorithm>
#include <float.h>
using namespace std;

///////////////////////////////////////////////////////////////////////////////

//count the number of notches from a given edge list (EL)
inline unsigned int countNotch(EL& el)
{
    unsigned int count=0;
    FOREACH_E(el)
        if(e->reflex) count++;
    FOREACH_END
    return count;
}

///////////////////////////////////////////////////////////////////////////////
//
// fake VE
//
///////////////////////////////////////////////////////////////////////////////

inline void addFakeVE(VL& vl, cd_v * _v)
{
    FOREACH_V(vl)
        cd_e * e=new cd_e();
        e->v[0]=_v; e->v[1]=v;
        e->length=(float)(_v->pos-v->pos).norm();
        _v->edges.push_back(e);
        v->edges.push_back(e);
        _v->concavity+=v->concavity;
    FOREACH_END
    _v->concavity/=vl.size();
}

inline void removeFakeVE(cd_v * fakeV)
{
    FOREACH_E(fakeV->edges)
        cd_v * ov=e->otherv(fakeV);
        EL& edges=ov->edges;
        EIT pos=find(edges.begin(),edges.end(),e);
        ov->edges.erase(pos);
        delete e;
    FOREACH_END
    delete fakeV;
}

inline void addFakeVE(HL& holes)
{
    FOREACH_H(holes)
        cd_v * fakev=new cd_v();
        fakev->pos=h.com;
        ::addFakeVE(h.vl,fakev);
    FOREACH_END
}

inline void removeFakeVE(HL& holes)
{
    FOREACH_H(holes)
        //find fake v
        cd_v * fakeV=NULL;
        cd_v * v=h.vl.front();
        FOREACH_E(v->edges)
            cd_v * ov=e->otherv(v);
            if( ov->id<0 ){
                fakeV=ov;
                break;
            }
        FOREACH_END
        ::removeFakeVE(fakeV);
    FOREACH_END
}

///////////////////////////////////////////////////////////////////////////////
inline void markFlagAs(VL& vl, flagT flag)
{
    FOREACH_V(vl)
        if( v->id==-1 ) continue;
        v->flag=flag;
    FOREACH_END
}

inline void markFlagAs_E(VL& vl, flagT flag)
{
    FOREACH_V(vl)
        VIT nv=iv; nv++;
        if(nv==vl.end()) break;
        cd_e * e=v->find(*nv);
        if( e==NULL ) continue;
        e->flag=flag;
    FOREACH_END
}


inline void markFlagAs(EL& el, flagT flag)
{
    FOREACH_E(el)
        e->flag=flag;
    FOREACH_END
}

inline void markFlagAs(FL& fl, flagT flag)
{
    FOREACH_F(fl)
        f->flag=flag;
    FOREACH_END
}

inline bool isMarkedAs(const VL& vl, flagT flag)
{
    FOREACH_CV(vl)
        if( v->id==-1 ) continue;
        if( v->flag==flag ) return true;
    FOREACH_END
    return false;
}

inline bool isMarkedAs_E(const VL& vl, flagT flag)
{
    FOREACH_CV(vl)
        CVIT nv=iv; nv++;
        if(nv==vl.end()) break;
        cd_e * e=v->find(*nv);
        if( e==NULL ) continue;
        if( e->flag==flag ) return true;
    FOREACH_END
    return false;
}


inline bool isMarkedAs(const EL& el, flagT flag)
{
    FOREACH_CE(el)
        if( e->flag==flag ) return true;
    FOREACH_END
    return false;
}

inline bool isMarkedAs(const FL& fl, flagT flag)
{
    FOREACH_CF(fl)
        if( f->flag==flag ) return true;
    FOREACH_END
    return false;
}

///////////////////////////////////////////////////////////////////////////////
inline cd_v * getRootV(cd_v * v)
{
    while(v->parent!=v) v=v->parent;
    return v;
}

///////////////////////////////////////////////////////////////////////////////
//return right,left facet
inline pair<cd_f*,cd_f*> orderFacets(cd_v * v1,cd_v *v2)
{
    cd_e * e=v1->find(v2);
    cd_f * f1=e->f[0];
    cd_v * v3=f1->otherv(e);
    Vector3d vec1=v2->pos-v1->pos;
    Vector3d vec2=v3->pos-v1->pos;
    double dir=(vec1%vec2)*f1->n;
    if( dir>0 )
        return pair<cd_f*,cd_f*>(f1,e->f[1]);
    else
        return pair<cd_f*,cd_f*>(e->f[1],f1);
}

//return right,left facet
inline pair<cd_f*,cd_f*> orderFacets(cd_v * v1,cd_v * v2, cd_e * e)
{
    cd_f * f1=e->f[0];
    cd_v * v3=f1->otherv(e);
    Vector3d vec1=v2->pos-v1->pos;
    Vector3d vec2=v3->pos-v1->pos;
    double dir=(vec1%vec2)*f1->n;
    if( dir>0 )
        return pair<cd_f*,cd_f*>(f1,e->f[1]);
    else
        return pair<cd_f*,cd_f*>(e->f[1],f1);
}

///////////////////////////////////////////////////////////////////////////////
//
//
//  Compute center of a list of vertices
//
//
///////////////////////////////////////////////////////////////////////////////

inline Point3d computeCenter(const PV& pv)
{
	Point3d com;
	FOREACH_CPt(pv)
		for (int i = 0; i<3; i++)  com[i] += p[i];
	FOREACH_END
	int size = pv.size();
	for (int i = 0; i<3; i++) com[i] /= size;
	return com;
}


inline Point3d computeCenter( const VL& vl )
{
    Point3d com;
    FOREACH_CV(vl)
        for( int i=0;i<3;i++ )  com[i]+=v->pos[i];
    FOREACH_END
    int size=vl.size();
    for( int i=0;i<3;i++ ) com[i]/=size;
    return com;
}

inline float computeRadius( const VL& vl, const Point3d& com )
{
    float max_R=-1e10f;
    FOREACH_CV(vl)
        float d=(float)(v->pos-com).normsqr();
        if(d>max_R) max_R=d;
    FOREACH_END
    return sqrt(max_R);
}

inline float computeRadius( const VL& vl )
{
    if(vl.empty()) return 0;
    Point3d com=computeCenter(vl);
    return computeRadius(vl,com);
}

inline void computeBoundingBox( const VL& vl, float bbox[6] )
{
    bbox[0]=bbox[2]=bbox[4]= FLT_MAX;
    bbox[1]=bbox[3]=bbox[5]=-FLT_MAX;

    FOREACH_CV(vl)
        if(v->pos[0]<bbox[0]) bbox[0]=v->pos[0];
        if(v->pos[0]>bbox[1]) bbox[1]=v->pos[0];
        if(v->pos[1]<bbox[2]) bbox[2]=v->pos[1];
        if(v->pos[1]>bbox[3]) bbox[3]=v->pos[1];
        if(v->pos[2]<bbox[4]) bbox[4]=v->pos[2];
        if(v->pos[2]>bbox[5]) bbox[5]=v->pos[2];
    FOREACH_END
}

///////////////////////////////////////////////////////////////////////////////
//
//
//  Set Parent to
//
//
///////////////////////////////////////////////////////////////////////////////

inline void parent_V( const VL& vl, cd_v * p )
{
    for(CVIT iv=vl.begin();iv!=vl.end();iv++ )
        (*iv)->parent=p;
}

///////////////////////////////////////////////////////////////////////////////
//
//
//  VL functions (VL, vertex list)
//
//
///////////////////////////////////////////////////////////////////////////////

inline float pathLength( CVIT s, CVIT g )
{
    if( s==g ) return 0;
    float L=0;
    CVIT end=g; end--;
    for(CVIT iv=s;iv!=end;iv++){
        CVIT niv=iv; niv++;
        cd_e * e=(*iv)->find(*niv);
        if(e!=NULL) L+=e->length;
        else L+=(float)((*iv)->pos-(*niv)->pos).norm();
    }//end for
    return L;
}

inline float pathLength( const VL& path )
{
    return pathLength(path.begin(),path.end());
}

//compute the distance between points
inline float vlLength( VIT s, VIT g )
{
    if( s==g ) return 0;

    float L=0;
    VIT end=g; end--;
    for(VIT iv=s;iv!=end;iv++){
        VIT niv=iv; niv++;
        float d=(float)((*iv)->pos-(*niv)->pos).norm();
        L+=d;
    }//end for
    return L;
}

inline float pathConcavity( VIT s, VIT g )
{
    float C=0;
    VIT end=g; end--;
    for(VIT iv=s;iv!=end;iv++){
        VIT niv=iv; niv++;
        cd_e * e=(*iv)->find(*niv);
        if( !e->reflex ) continue;
        C+=(e->length);
    }//end for
    return C;
}

inline void printPath(VIT s, VIT g)
{
    VIT e=g; e--;
    for(VIT iv=s;iv!=e;iv++)
        cout<<(*iv)->id<<"->";
    cout<<(*e)->id;
}

//return true if success
inline bool vl2el(VL& vl, EL& el)
{
    bool find_NUL=false;
    FOREACH_Vi(++vl.begin(),vl.end())
        VIT pre=iv; pre--;
        cd_e * e=v->find(*pre);
        if(e==NULL){
            cd_wrn(cd_msg("! can't find edge betweem v1=")+v->id+" v2="+(*pre)->id);
            //return false;
            find_NUL=true;
            continue;
        }
        el.push_back(e);
    FOREACH_END
    return !find_NUL;
    //return true; //done
}


///////////////////////////////////////////////////////////////////////////////
// Change the vertex flag of those vertices in the same connected components
// of v to the type of v

inline void propagateFlag( cd_v * v )
{
    flagT flag=v->flag;
    VL cur;
    cur.push_back(v);

    while(cur.empty()==false ){
        cd_v * v=cur.front();
        cur.pop_front();
        FOREACH_E(v->edges)
            cd_v * o=(e->v[0]!=v)?e->v[0]:e->v[1];
            if( o->flag!=flag ){
                o->flag=flag;
                cur.push_back(o);
            }
        FOREACH_END
    }//end while
}


inline float signed_dist(const Point3d& s, const Point3d& t, const Point3d& q)
{
    Vector3d vec1=t-s;
    Vector3d vec2=q-s;
    Vector3d vec3=q-t;

    float norm1=vec1.normsqr();
    float norm2=vec2.normsqr();
    float norm3=vec3.normsqr();

    if(norm2==0 || norm3==0) return 0;
    if(norm1==0) return sqrt(norm2);

    float fSign = vec1[0]*vec2[1] - vec1[1]*vec2[0];

    vec1=vec1/sqrt(norm1);
    float dot=vec1*vec2;
    float diff=norm2-dot*dot;
    if(diff<0) diff=0;
    float d=sqrt(diff);

    if(fSign >= 0)
    	return d;
    else
    	return d*(-1);
}

inline float dist2Line(const Point3d& s, const Point3d& t, const Point3d& q)
{
    Vector3d vec1=t-s;
    Vector3d vec2=q-s;
    Vector3d vec3=q-t;

    float norm1=vec1.normsqr();
    float norm2=vec2.normsqr();
    float norm3=vec3.normsqr();

    if(norm2==0 || norm3==0) return 0;
//    if(norm1==0) return sqrt(norm2);


    vec1=vec1/sqrt(norm1);
    float dot=vec1*vec2;
    float diff=norm2-dot*dot;
    if(diff<0) diff=0;
    float d=sqrt(diff);

    return d;
}

inline float dist(const Point3d& s, const Point3d& t, const Point3d& q)
{
    Vector3d vec1=t-s;
    Vector3d vec2=q-s;
    Vector3d vec3=q-t;

    float norm1=vec1.normsqr();
    float norm2=vec2.normsqr();
    float norm3=vec3.normsqr();

    if(norm2==0 || norm3==0) return 0;
    if(norm1==0) return sqrt(norm2);


    vec1=vec1/sqrt(norm1);
    float dot=vec1*vec2;
    float diff=norm2-dot*dot;
    if(diff<0) diff=0;
    float d=sqrt(diff);

    return d;
}


//distance from a point q to a segment (st)
inline float dist(cd_v * s, cd_v * t, cd_v *q)
{
    return dist(s->pos,t->pos,q->pos);
}

//volume of a tetra
inline double vol(const Point3d& p1,const Point3d& p2,const Point3d& p3,const Point3d& p4)
{
    return
    -(p1[2]-p4[2])*(p2[1]-p4[1])*(p3[0]-p4[0])
    +(p1[1]-p4[1])*(p2[2]-p4[2])*(p3[0]-p4[0])
    +(p1[2]-p4[2])*(p2[0]-p4[0])*(p3[1]-p4[1])
    -(p1[0]-p4[0])*(p2[2]-p4[2])*(p3[1]-p4[1])
    -(p1[1]-p4[1])*(p2[0]-p4[0])*(p3[2]-p4[2])
    +(p1[0]-p4[0])*(p2[1]-p4[1])*(p3[2]-p4[2]);
}

//compute the volume of this model
//assume it is closed....
inline double vol(cd_m * m)
{
  double v=0;
  const Point3d O;
  for(cd_f* f : m->f())
  {
    v+=vol(f->v[0]->pos,f->v[2]->pos,f->v[2]->pos,O);
  }

  return v;
}

inline double area(cd_m * m)
{
  double a=0;
  for(cd_f* f : m->f())
  {
    a+=f->area();
  }

  cout<<"m->f().size()="<<m->f().size()<<endl;

  return a;
}

inline bool isflat(cd_m * m)
{
  if(m->f().empty()) return true; //nothing??
  Point3d& p=m->v().front()->pos;
  Vector3d& n=m->f().front()->n;
  for(cd_v * v:m->v())
  {
    if( fabs((v->pos-p)*n)>1e-10) return false;
  }
  return true;
}
///////////////////////////////////////////////////////////////////////////////
//
//
//  Compute Left/Right edges/facets of a path
//
//
///////////////////////////////////////////////////////////////////////////////

//note: in order to get all edges/facet on the side, this model
//need to be manifold
void getSideEF(cd_v *v, cd_e * e[2], cd_f *rf[2], cd_f *lf[2], EL& rE, EL& lE, FL& rF, FL& lF);
void getSideEF(cd_v * v1,cd_v * v2,cd_v * v3,EL& rE, EL& lE, FL& rF, FL& lF);
void getSideEF(VIT begin, VIT end, EL& rE, EL& lE, FL& rF, FL& lF);
void getSideEdges(VIT begin, VIT end, EL& rE, EL& lE);



#endif //_MODEL_UTILITY_H_
