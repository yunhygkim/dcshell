#include "hull_vol.h"
#include "CD3d_model.h"

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include "intersection.h"
inline void getOrderedNormals(cd_v * v, vector<Vector3d>& normals)
{
  flagT flag=generateID();
  cd_f * f=NULL;

  //find a face with non-zero area
  for(  cd_e * e  : v->edges)
  {
    f=e->f[0];
    if(f->area()>0) break;
    f=e->f[1];
    if(f->area()>0) break;
  }

  if(f==NULL) return;
  normals.push_back(f->n);

  while(f->flag!=flag)
  {
    f->flag=flag;
    //cout<<"f->area()="<<f->area()<<endl;
    for(int i=0;i<3;i++)
    {
      if( f->e[i]->v[0]!=v && f->e[i]->v[1]!=v ) continue; //this e dose not contain v
      cd_f * of=f->e[i]->otherf(f); //of also contains v

      if(of->flag!=flag)
      {
        if( AlmostEqual3( of->n.get(),normals.back().get())==false )
        {
          //if(!(of->n==normals.back()))
          if(of->area()>0)
             normals.push_back(of->n);
        }
        f=of;
        break;
      }
    }//end for i
  }//end while

  //if(normals.front()==normals.back()) normals.pop_back();
  //cout<<"normals size="<<normals.size()<<endl;
  if( AlmostEqual3(normals.front().get(), normals.back().get()) ) normals.pop_back();
}

//
// for computing volume
//

inline double compute_cone_volume(const Point3d& p, const hull_plane& hp, const Vector3d& a0, const Vector3d& a1, const Vector3d& a2, bool print=false)
{
  Vector3d n=hp.n;
  Eigen::Matrix3d A;
	Eigen::Vector3d c;
  //cout<<"a0="<<a0<<" a1="<<a1<<" a2="<<a2<<endl;

	A << a0[0], a1[0], a2[0], a0[1], a1[1], a2[1], a0[2], a1[2], a2[2];
	c << n[0], n[1], n[2];
  Eigen::Vector3d r = A.colPivHouseholderQr().solve(c);; // A.fullPivHouseholderQr().solve(c); //A.colPivHouseholderQr().solve(c);

  if(r[0]==0 || r[1]==0 || r[2]==0)
  {
    cout<<"! Error: Plane normal is parallet to face normal..."<<endl;
    return 0;
  }

  double delta_v = fabs(A.determinant());
  double d=-(hp.o[0]*n[0]+hp.o[1]*n[1]+hp.o[2]*n[2]);
  //double f_v = max(0.0,pow(d+(p[0]*n[0]+p[1]*n[1]+p[2]*n[2]),3));
  double f_v = max(0.0,d+(p[0]*n[0]+p[1]*n[1]+p[2]*n[2]));

  if(print)
  {
    cout<<"------"<<endl;
    cout<<"A="<<A<<" \nc="<<c<<" \nr="<<r<<"\ndelta="<<delta_v<<" d="<<d<<" f_v="<<f_v<<" vol="<<
    ((f_v/r[0])*(f_v/r[1])*(f_v/r[2]))/(6*delta_v)<<endl;;
    //f_v/(6*delta_v*r[0]*r[1]*r[2])<<endl;
  }

  return ((f_v/r[0])*(f_v/r[1])*(f_v/r[2]))/(6*delta_v);
  //return f_v/(6*delta_v*r[0]*r[1]*r[2]);
}


double compute_cone_volume(cd_v * v, hull_plane& hp, bool print=false)
{
  //get incident faces
  vector<Vector3d> normals;
  getOrderedNormals(v, normals);
  //if(print) cout<<"normals size="<<normals.size()<<endl;

  if(normals.size()==3)
  {
    double vol=compute_cone_volume(v->pos, hp, normals[0], normals[1], normals[2], false);
    if(print) cout<<"! regular case vol="<<vol<<endl;
    return vol;
  }

  //compute tangent vector
  Vector3d an;
  for(auto& n:normals){
    an=an+n;
    //cout<<"n="<<n<<endl;
  }
  an=an.normalize();
  //cout<<"an="<<an<<endl;

  double vol=0;
  for(auto it=normals.begin();it!=normals.end();it++)
  {
    auto nit=(it+1==normals.end())?normals.begin():it+1;
    auto tmp=compute_cone_volume(v->pos, hp, an, *it, *nit, false);
    if(print) cout<<"case["<<normals.size()<<"] V="<<tmp<<endl;
    vol+=tmp;
  }

  if(print) cout<<"case["<<normals.size()<<"] vol="<<vol<<endl;
  //cout<<"------"<<endl;
  return vol;
}


//compute the volume of intersection between hull and p
//without computing the intersection explicitly
double compute_volume(const cd_hull& hull, hull_plane& hp, bool print)
{
  double vol=0;
  cd_m * m = hull.to_m();
  for(cd_v * v: m->v())
  {
    double vv=compute_cone_volume(v, hp, print);
    //cout<<"v:"<<v->id<<"/"<<m->v().size()<<endl;
    vol+=vv;
  }

  m->destroy();
  delete m;
  return vol;
}

//
// for computing volume derivatives
//

inline void compute_cone_volume_derivative_coeeficient
(const Point3d& p, hull_plane& hp,
 const Vector3d& a0, const Vector3d& a1, const Vector3d& a2,
 double& c2, double& c1, double& c0)
{
  Vector3d n=hp.n;
  Eigen::Matrix3d A;
  Eigen::Vector3d c;
  //cout<<"a0="<<a0<<" a1="<<a1<<" a2="<<a2<<endl;

  A << a0[0], a1[0], a2[0], a0[1], a1[1], a2[1], a0[2], a1[2], a2[2];
  c << n[0], n[1], n[2];
  Eigen::Vector3d r = A.colPivHouseholderQr().solve(c);; // A.fullPivHouseholderQr().solve(c); //A.colPivHouseholderQr().solve(c);
  double delta_v = fabs(A.determinant());
  //double d=-(hp.o[0]*n[0]+hp.o[1]*n[1]+hp.o[2]*n[2]);
  double denom = (6*delta_v*r[0]*r[1]*r[2]);
  double ct = (p[0]*n[0]+p[1]*n[1]+p[2]*n[2]);

  double myc2=3.0/denom;
  double myc1=(6.0*ct)/denom;
  double myc0=(3.0*ct*ct)/denom;

  if(!isfinite(myc2) || !isfinite(myc1) || !isfinite(myc0) )
  {
    cout<<"A="<<A<<endl;
    cout<<"c="<<c<<endl;
    cout<<"r="<<r<<endl;
    cout<<"myc2="<<myc2<<" myc1="<<myc1<<" myc0="<<myc0<<" delta_v"<<delta_v<<" denom="<<denom<<" ct="<<ct<<endl;
    exit(1);
  }

  c2+=myc2;
  c1+=myc1;
  c0+=myc0;
}

void compute_cone_volume_derivative_coeeficient
(cd_v * v, hull_plane& hp, double& c2, double& c1, double& c0)
{
  //get incident faces
  vector<Vector3d> normals;
  getOrderedNormals(v, normals);

  if(normals.size()==3)
  {
    compute_cone_volume_derivative_coeeficient(v->pos, hp, normals[0], normals[1], normals[2], c2, c1, c0);
    return;
  }

  Vector3d an;
  for(auto& n:normals) an=an+n;

  an=an.normalize();
  for(auto it=normals.begin();it!=normals.end();it++)
  {
    auto nit=(it+1==normals.end())?normals.begin():it+1;
    compute_cone_volume_derivative_coeeficient(v->pos, hp, an, *it, *nit, c2, c1, c0);
  }
}

//
//
// test only
//
//

#include "CD3d_trim_hull.h"
double compute_volume_explicit(const cd_hull& hull, hull_plane& hp)
{
  //compute the trimed hull and then
  cd_hull residual_hull;
  list<hull_plane> planes;
  planes.push_back(hp);
  trim_hull_by_planes(hull, planes, residual_hull,0);
  return residual_hull.computeHullVol();
}
