/********************************************************************

    Quaternion.h    Header File

    Jyh-Ming Lien 03/30/2002
    Computer Science.
    Texas A&M University

*********************************************************************/

#if !defined( _QUATERNION_H_NPRM_ )
#define _QUATERNION_H_NPRM_

#include "util/mathtool/Vector.h"
#include "util/mathtool/Matrix.h"

namespace mathtool{

    class Quaternion {
    public:
        Quaternion(){ m_s=1; m_v.set(0,0,0); }
        Quaternion( float s, Vector3d v ){ m_v=v; m_s=s; }
        Quaternion( const Quaternion & q){ *this=q; }
        Quaternion( float r[3]){  //compute rotational quaternion
            //convert rot to quaternion
            float sin_r=(float)sin(r[0]/2);
            float cos_r=(float)cos(r[0]/2);
            Quaternion qx(cos_r,sin_r*Vector3d(1,0,0));
            sin_r=(float)sin(r[1]/2);
            cos_r=(float)cos(r[1]/2);
            Quaternion qy(cos_r,sin_r*Vector3d(0,1,0));
            sin_r=(float)sin(r[2]/2);
            cos_r=(float)cos(r[2]/2);
            Quaternion qz(cos_r,sin_r*Vector3d(0,0,1));
            *this=(qz*qy*qx).normalize();
        }

        ////////////////////////////////////////////////////////////////////////
        // Operations for Quaternion
        Quaternion operator*(const Quaternion & q) const {
            float s=q.m_s*m_s-q.m_v*m_v;
            Vector3d v=q.m_v*m_s+m_v*q.m_s+m_v%q.m_v;
            return Quaternion(s,v);
        }
        Quaternion operator*(const Vector3d & v) const { return *this*Quaternion(0,v); }
        Quaternion operator/(float s) const { return Quaternion(m_s/s,m_v/s); }
        Quaternion & operator=(const Quaternion & q){ set(q.m_s,q.m_v); return *this; }
        Quaternion operator+(const Quaternion & q) const { return Quaternion(m_s+q.m_s,m_v+q.m_v); }
        Quaternion operator-(const Quaternion & q) const { return Quaternion(m_s-q.m_s,m_v-q.m_v); }
        Quaternion operator-() const { return Quaternion(m_s,-m_v); }
        friend Quaternion operator*(const Vector3d & v, const Quaternion & q);
        friend istream& operator>>(istream & in, Quaternion & q );
        friend ostream& operator<<(ostream & out, const Quaternion & q );

        //////////////////////////////////////////////////////////////////////////
        //Normalization
        Quaternion normalize(){ 
            Quaternion q(*this);
            float l=q.norm();
            q=q/l;
            return q;
        }

        float norm(){ return sqrt(normsqr()); }
        float normsqr(){ return m_v.normsqr()+sqr(m_s); }

        //////////////////////////////////////////////////////////////////////////
        //Access

        Matrix3x3 getMatrix(){
            float x_2=2*sqr(m_v[0]); float y_2=2*sqr(m_v[1]); float z_2=2*sqr(m_v[2]);
            float xy=2*m_v[0]*m_v[1]; float yz=2*m_v[1]*m_v[2]; float zx=2*m_v[2]*m_v[0]; 
            float sx=2*m_s*m_v[0]; float sy=2*m_s*m_v[1]; float sz=2*m_s*m_v[2]; 
            return Matrix3x3(1-y_2-z_2, xy-sz, zx+sy,
                             xy+sz, 1-x_2-z_2, yz-sx,
                             zx-sy, yz+sx, 1-x_2-y_2);
        }
    

        void set(float s,const Vector3d & v){ m_v=v; m_s=s; }
        void set(float q1, float q2, float q3, float q4){ m_s=q1; m_v.set(q2,q3,q4); }
        const Vector3d& getComplex() const { return m_v; }
        float getReal() const { return m_s; }

    private:
        Vector3d m_v;
        float m_s;
    };

}//emd of nprmlib

#endif //#if !defined( _QUATERNION_H_NPRM_ ) 

