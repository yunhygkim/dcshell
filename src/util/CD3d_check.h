#ifndef _CD3D_CHECK_H_
#define _CD3D_CHECK_H_

#include "CD3d_model.h"
#include <functional>
#include <cassert>
using namespace std;

///////////////////////////////////////////////////////////////////////////////
//
// Check and output the information for
// the models
//
///////////////////////////////////////////////////////////////////////////////

inline ostream& operator<<(ostream& out, const cd_f& f)
{
    out<<"f id="<<f.id<<" v=[";
    for(int iv=0;iv<3;iv++){
        if(f.v[iv]!=NULL) out<<f.v[iv]->id;
        else out<<"NULL";
        if( iv<2 ) out<<",";
    }
    out<<"]";

    return out;
}

inline ostream& operator<<(ostream& out, const cd_v& v)
{
    out<<"v "<<v.id<<" pos="<<v.pos;
    return out;
}

inline ostream& operator<<(ostream& out, const cd_e& e)
{
    out<<"e ";
    if( e.reflex ) out<<"r v=[";
    else out<<"n v=[";

    if( e.v[0]!=NULL && e.v[1]!=NULL ){
        out<<e.v[0]->id<<","<<e.v[1]->id;
    }
    else if( e.v[0]!=NULL ){
        out<<e.v[0]->id<<",NULL";
    }
    else if( e.v[1]!=NULL ){
        out<<"NULL,"<<e.v[1]->id;
    }
    else{
        out<<"NULL,NULL";
    }

    out<<"] f=[";

    if( e.f[0]!=NULL && e.f[1]!=NULL ){
        out<<e.f[0]->id<<","<<e.f[1]->id;
    }
    else if( e.f[0]!=NULL ){
        out<<e.f[0]->id<<",NULL";
    }
    else if( e.f[1]!=NULL ){
        out<<"NULL,"<<e.f[1]->id;
    }
    else{
        out<<"NULL,NULL";
    }
    out<<"]";

    return out;
}

inline bool compare_e( cd_e * e1, cd_e * e2 )
{
    double e1_max=(e1->v[0]->id>e1->v[1]->id)?e1->v[0]->id:e1->v[1]->id;
    double e1_min=(e1->v[0]->id<e1->v[1]->id)?e1->v[0]->id:e1->v[1]->id;
    double e2_max=(e2->v[0]->id>e2->v[1]->id)?e2->v[0]->id:e2->v[1]->id;
    double e2_min=(e2->v[0]->id<e2->v[1]->id)?e2->v[0]->id:e2->v[1]->id;

    if( e1_max>e2_max ) return true;
    if( e1_max<e2_max ) return false;
    return e1_min>e2_min;
}

inline ostream& operator<<(ostream& out, const cd_m& m )
{
    //Check Edges
    {
        list<cd_e*> el=m.e();
        el.sort(compare_e);
        for(EIT ie=el.begin();ie!=el.end();ie++){
            cd_e * e=*ie;
            out<<*e<<endl;
        }//end for
    }

    return out;
}// output model

///////////////////////////////////////////////////////////////////////////////
// Check if model has holes
// return true is hole is found

inline bool checkHoles( cd_m * m ){
    typedef list<cd_e*>::iterator EIT;
    list<cd_e*>& el=m->e();
    for(EIT ie=el.begin();ie!=el.end();ie++){
        cd_e * e=*ie;
        if( e->f[0]==NULL || e->f[1]==NULL ) return true;
    }
    return false;
}

inline void checkV(VL& vl,EL& el)
{
    ///////////////////////////////////////////////////////////////////////////
    //
    // Check Vertices
    //
    ///////////////////////////////////////////////////////////////////////////
    cerr<<"- Check Vertices"<<endl;

    FOREACH_V(vl)
        FOREACH_E(v->edges)
            cd_v * o = (e->v[0]!=v)?e->v[0]:e->v[1];
            if(e->v[0]!=v && e->v[1]!=v){
                cerr<<"! CheckModel Error : v:"<<*v<<" Not in the Edge:"<<*e<<endl;
            }
            else if( o->find(v)==NULL ){ //v is in this edge e
                cerr<<"! CheckModel Error : v: "<<*o<<" is not neighbor of"<<*v<<endl;
            }
            if( e->f[0]==NULL ) //|| e->f[1]==NULL)
            {
                cerr<<"! CheckModel Error : v: "<<*v<<" has bad Edge: "<<*e<<endl;
            }

            EIT pos=find(el.begin(),el.end(),e);
            if( pos==el.end() ){
                cerr<<"! CheckModel Error : v: "<<*v<<" has edge "
                    <<*e<<" not in the global list"<<endl;
            }
        FOREACH_END
    FOREACH_END
}

inline void checkE(EL& el,VL& vl,FL& fl)
{
    ///////////////////////////////////////////////////////////////////////////
    //
    // Check Edges
    //
    ///////////////////////////////////////////////////////////////////////////
    cerr<<"- Check Edge"<<endl;
    FOREACH_E(el)
        {for(int i=0;i<2;i++){
            if( e->v[i]!=NULL ){
                VIT pos=find(vl.begin(),vl.end(),e->v[i]);
                if( pos==vl.end() )
                    cerr<<"! CheckModel Error : edge="<<*e<<" has vertex="<<*(e->v[i])
                        <<" not in the global list"<<endl;

                //make sure this edges is in v1 and v2
                EIT epos=find(e->v[i]->edges.begin(),e->v[i]->edges.end(),e);
                if( epos==e->v[i]->edges.end() ){ 
                    cerr<<"! CheckModel Error : edge="<<*e
                        <<" is not in the his vertex "<<*e->v[i]<<endl;
                }
            }
            else
                cerr<<"! CheckModel Error : edge="<<*e<<" has NULL vetex="
                    <<i<<endl;
        }}//end for

        {for(int i=0;i<2;i++){
            cd_f * f=e->f[i];
            if( f!=NULL ){
                bool found=false;
                for( int j=0;j<3;j++ ){
                    if( f->e[j]==e ){ found=true; break; }
                }
                if( !found ) 
                    cerr<<"! CheckModel Error : edge "<<*e
                        <<" is not in face "<< f->id <<endl;
                FIT pos=find(fl.begin(),fl.end(),f);
                if( pos==fl.end() )
                    cerr<<"! CheckModel Error : edge "<<*e
                        <<" has face "<<f->id 
                        <<" not in the global list"<<endl;
            }//end if
            else if(i==0)
                cerr<<"! CheckModel Error : edge="<<*e<<" has NULL face="
                    <<i<<endl;
        }}//end for

        bool r=e->reflex;
        e->compRef();
        if( r!=e->reflex ){ 
            cerr<<"! CheckModel Error: edge "<<*e<<" reflectivity is wrong"
                <<" was "<<r<<" and is "<<e->reflex<<endl;
            for(int i=0;i<2;i++)
                if(e->f[i]!=NULL) 
                    cd_err(cd_msg("f=")+e->f[i]->id+" has n="+e->f[i]->n);
        }
    FOREACH_END
}

inline void checkF(FL& fl,EL& el,VL& vl)
{
    ///////////////////////////////////////////////////////////////////////////
    //
    // Check Faces
    //
    ///////////////////////////////////////////////////////////////////////////
    cerr<<"- Check Faces"<<endl;
    FOREACH_F(fl)
        if( f->area()==0 )
            cerr<<"! CheckModel Error : face="<<*f<<" has no area"<<endl;
        for( int i=0;i<3;i++ ){
            if( f->v[i]==NULL ){
                cerr<<"! CheckModel Error : face="<<*f<<" has a NULL vertex="
                    <<i<<endl;
            }
            else{
                VIT pos=find(vl.begin(),vl.end(),f->v[i]);
                if( pos==vl.end() )
                    cerr<<"! CheckModel Error : face="<<*f<<" has vertex="
                        <<*(f->v[i])<<" not in the global list"<<endl;
            }

            int nvid=(i+1)%3;
            if( f->v[nvid]==f->v[i] ){
                cerr<<"! CheckModel Error : face="<<*f<<" has vertex="
                    <<i<<" and "<<nvid<<" the same"<<endl;
            }

            if( f->e[i]==NULL ){
                cerr<<"! CheckModel Error : face="<<*f<<" has a NULL edge="
                    <<i<<endl;
                continue;
            }

            cd_e * e=f->e[i];

            EIT pos=find(el.begin(),el.end(),e);
            if( pos==el.end() ){
                cerr<<"! CheckModel Error : f: "<<*f<<" has edge "
                    <<*e<<" not in the global list"<<endl;
            }

            //edge i must contain vertex i and i+1
            int cur=i; int next=(cur+1)%3;
            if( (e->v[0]!=f->v[cur]) && (e->v[1]!=f->v[cur]) )
                cerr<<"! CheckModel Error : face "<<*f
                <<" whose edge="<<i<<" does not have vertex="<<cur<<endl;
            if( (e->v[0]!=f->v[next]) && (e->v[1]!=f->v[next]) )
                cerr<<"! CheckModel Error : face "<<*f
                <<" whose edge="<<i<<" does not have vertex="<<next<<endl;
            
            //face f must in its edge
            if( f!=e->f[0] && f!=e->f[1] ) 
                cerr<<"! CheckModel Error : face "<<*f
                    <<" is not in its edge"<< *e <<endl;
        }
    FOREACH_END
}

///////////////////////////////////////////////////////////////////////////////
// Check everything....
inline void checkModel( cd_m * m )
{
    checkV(m->v(),m->e());
    checkE(m->e(),m->v(),m->f());
    checkF(m->f(),m->e(),m->v());

    //check topological properties here
}


#endif //_CD3D_CHECK_H_
