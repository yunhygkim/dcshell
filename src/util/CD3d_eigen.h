#ifndef _CD3d_EIGEN_H_
#define _CD3d_EIGEN_H_

/*******************************************************************************

  Copied from RAPID

 *******************************************************************************/
//compute eigenvector and eigenvalue.
//vout is eigenvectors and dout is eigenvalues and a is matrix

#define rfabs(x) ((x < 0) ? -x : x)
#define ROT(a,i,j,k,l) g=a[i][j]; h=a[k][l]; a[i][j]=g-s*(h+g*tau); a[k][l]=h+s*(g-h*tau);

inline void
McM(double Mr[3][3], double M[3][3])
{
  Mr[0][0] = M[0][0];  Mr[0][1] = M[0][1];  Mr[0][2] = M[0][2];
  Mr[1][0] = M[1][0];  Mr[1][1] = M[1][1];  Mr[1][2] = M[1][2];
  Mr[2][0] = M[2][0];  Mr[2][1] = M[2][1];  Mr[2][2] = M[2][2];
}

inline void
VcV(double Vr[3], double V[3])
{
  Vr[0] = V[0];  Vr[1] = V[1];  Vr[2] = V[2];
}

int inline
Meigen(double vout[3][3], double dout[3], double a[3][3])
{
    int i;
    double tresh,theta,tau,t,sm,s,h,g,c;
    int nrot;
    double b[3];
    double z[3];
    double v[3][3];
    double d[3];
    
    v[0][0] = v[1][1] = v[2][2] = 1.0;
    v[0][1] = v[1][2] = v[2][0] = 0.0;
    v[0][2] = v[1][0] = v[2][1] = 0.0;
    
    b[0] = a[0][0]; d[0] = a[0][0]; z[0] = 0.0;
    b[1] = a[1][1]; d[1] = a[1][1]; z[1] = 0.0;
    b[2] = a[2][2]; d[2] = a[2][2]; z[2] = 0.0;
    
    nrot = 0;
    
    for(i=0; i<100; i++)
    {
        
        sm=0.0; sm+=fabs(a[0][1]); sm+=fabs(a[0][2]); sm+=fabs(a[1][2]);
        if (sm == 0.0) { McM(vout,v); VcV(dout,d); return i; }
        
        if (i < 3) tresh=0.2*sm/(3*3); else tresh=0.0;
        
        {
            g = 100.0*rfabs(a[0][1]);  
            if (i>3 && rfabs(d[0])+g==rfabs(d[0]) && rfabs(d[1])+g==rfabs(d[1]))
                a[0][1]=0.0;
            else if (rfabs(a[0][1])>tresh)
            {
                h = d[1]-d[0];
                if (rfabs(h)+g == rfabs(h)) t=(a[0][1])/h;
                else
                {
                    theta=0.5*h/(a[0][1]);
                    t=1.0/(rfabs(theta)+sqrt(1.0+theta*theta));
                    if (theta < 0.0) t = -t;
                }
                c=1.0/sqrt(1+t*t); s=t*c; tau=s/(1.0+c); h=t*a[0][1];
                z[0] -= h; z[1] += h; d[0] -= h; d[1] += h;
                a[0][1]=0.0;
                ROT(a,0,2,1,2); ROT(v,0,0,0,1); ROT(v,1,0,1,1); ROT(v,2,0,2,1); 
                nrot++;
            }
        }
        
        {
            g = 100.0*rfabs(a[0][2]);
            if (i>3 && rfabs(d[0])+g==rfabs(d[0]) && rfabs(d[2])+g==rfabs(d[2]))
                a[0][2]=0.0;
            else if (rfabs(a[0][2])>tresh)
            {
                h = d[2]-d[0];
                if (rfabs(h)+g == rfabs(h)) t=(a[0][2])/h;
                else
                {
                    theta=0.5*h/(a[0][2]);
                    t=1.0/(rfabs(theta)+sqrt(1.0+theta*theta));
                    if (theta < 0.0) t = -t;
                }
                c=1.0/sqrt(1+t*t); s=t*c; tau=s/(1.0+c); h=t*a[0][2];
                z[0] -= h; z[2] += h; d[0] -= h; d[2] += h;
                a[0][2]=0.0;
                ROT(a,0,1,1,2); ROT(v,0,0,0,2); ROT(v,1,0,1,2); ROT(v,2,0,2,2); 
                nrot++;
            }
        }
        
        
        {
            g = 100.0*rfabs(a[1][2]);
            if (i>3 && rfabs(d[1])+g==rfabs(d[1]) && rfabs(d[2])+g==rfabs(d[2]))
                a[1][2]=0.0;
            else if (rfabs(a[1][2])>tresh)
            {
                h = d[2]-d[1];
                if (rfabs(h)+g == rfabs(h)) t=(a[1][2])/h;
                else
                {
                    theta=0.5*h/(a[1][2]);
                    t=1.0/(rfabs(theta)+sqrt(1.0+theta*theta));
                    if (theta < 0.0) t = -t;
                }
                c=1.0/sqrt(1+t*t); s=t*c; tau=s/(1.0+c); h=t*a[1][2];
                z[1] -= h; z[2] += h; d[1] -= h; d[2] += h;
                a[1][2]=0.0;
                ROT(a,0,1,0,2); ROT(v,0,1,0,2); ROT(v,1,1,1,2); ROT(v,2,1,2,2); 
                nrot++;
            }
        }
        
        b[0] += z[0]; d[0] = b[0]; z[0] = 0.0;
        b[1] += z[1]; d[1] = b[1]; z[1] = 0.0;
        b[2] += z[2]; d[2] = b[2]; z[2] = 0.0;
        
    }
    
    fprintf(stderr, "eigen: too many iterations in Jacobi transform (%d).\n", i);
    
    return i;
}

/*******************************************************************************
 *              My Stuff                                                       *
 *******************************************************************************/

///////////////////////////////////////////////////////////////////////////////
#include "util/mathtool/Point.h"
#include "util/mathtool/Vector.h"
using namespace mathtool;

///////////////////////////////////////////////////////////////////////////////
inline
Point3d ComputeCOM(const PtVector& vl)
{
    typedef PtVector::const_iterator IT;

    Point3d com;
    for(IT i=vl.begin();i!=vl.end();i++ ){
        const Point3d& pt=*i;
        for( int D=0;D<3;D++ ) com[D]+=pt[D];
    }
    int size=vl.size();
    for( int D=0;D<3;D++ ) com[D]/=size;
    return com;
}

inline
double Covariance( const PtVector& pts, const Point3d& com, int j, int k )
{
    typedef PtVector::const_iterator IT;

    double r=0;
    for(IT i=pts.begin();i!=pts.end();i++ ){
        const Point3d& pt=*i;
        r+=((pt[j]-com[j])*(pt[k]-com[k]));
    }

    return r/pts.size();
}

inline
void Covariance( double c[3][3], const PtVector& pts, const Point3d& com )
{
    for( int j=0;j<3;j++ ){
        for( int k=0;k<3;k++ ){
            c[j][k]=Covariance(pts,com,j,k);
        }
    }
}

inline
void EigenVectors(double c[3][3], Vector3d& v1, Vector3d& v2, Vector3d& v3)
{
    double evec[3][3]; double eval[3];
    Meigen(evec,eval,c);
    if( eval[0]>eval[1] && eval[0]>eval[2] ) //eval[0] is the largest
    { 
        v1[0]=evec[0][0]; v1[1]=evec[1][0]; v1[2]=evec[2][0]; 
        if( eval[1]>eval[2] ){ 
            v2[0]=evec[0][1]; v2[1]=evec[1][1]; v2[2]=evec[2][1]; 
            v3[0]=evec[0][2]; v3[1]=evec[1][2]; v3[2]=evec[2][2]; 
        }
        else{
            v2[0]=evec[0][2]; v2[1]=evec[1][2]; v2[2]=evec[2][2]; 
            v3[0]=evec[0][1]; v3[1]=evec[1][1]; v3[2]=evec[2][1]; 
        }
    }
    else if( eval[1]>eval[0] && eval[1]>eval[2] )
    { 
        v1[0]=evec[0][1]; v1[1]=evec[1][1]; v1[2]=evec[2][1]; 
        if( eval[0]>eval[2] ){  //0>2
            v2[0]=evec[0][0]; v2[1]=evec[1][0]; v2[2]=evec[2][0]; 
            v3[0]=evec[0][2]; v3[1]=evec[1][2]; v3[2]=evec[2][2]; 
        }
        else{ //2>0
            v2[0]=evec[0][2]; v2[1]=evec[1][2]; v2[2]=evec[2][2]; 
            v3[0]=evec[0][0]; v3[1]=evec[1][0]; v3[2]=evec[2][0]; 
        }
    }
    else //eval[2] is the largest
    { 
        v1[0]=evec[0][2]; v1[1]=evec[1][2]; v1[2]=evec[2][2]; 
        if( eval[0]>eval[1] ){  //0>1
            v2[0]=evec[0][0]; v2[1]=evec[1][0]; v2[2]=evec[2][0]; 
            v3[0]=evec[0][1]; v3[1]=evec[1][1]; v3[2]=evec[2][1]; 
        }
        else{ //1>0
            v2[0]=evec[0][1]; v2[1]=evec[1][1]; v2[2]=evec[2][1]; 
            v3[0]=evec[0][0]; v3[1]=evec[1][0]; v3[2]=evec[2][0]; 
        }   
    }
}

inline Point3d EigenVectors
(const PtVector& pts, Vector3d& v1, Vector3d& v2, Vector3d& v3)
{
    double c[3][3];
    Point3d com=ComputeCOM(pts);
    Covariance(c,pts,com);
    EigenVectors(c,v1,v2,v3);
    return com;
}

#endif //_CD3d_EIGEN_H_
