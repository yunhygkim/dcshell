#include <cassert>
#include <climits>
#include <map>

#include "CD3d_model.h"
#include "util/CD3d_util.h"
//#include "util/helper/CD3dIdentifier.h"

#ifdef WIN32
#pragma warning(disable : 4786)
#endif

///////////////////////////////////////////////////////////////////////////////
//
//  FLAG ID for Marking
//
///////////////////////////////////////////////////////////////////////////////

flagT g_flag=INIT_FLAG+1;

///////////////////////////////////////////////////////////////////////////////
//
//  An unique id for propagating on the mesh
//
///////////////////////////////////////////////////////////////////////////////

//flags for marking prapagating
flagT generateID() {
    if( g_flag==UINT_MAX )
        cd_err("generateID: Reach the maximum value");
    return g_flag++;
} ///< Generate an unique id

///////////////////////////////////////////////////////////////////////////////
//
// cd_v: vertex
//
///////////////////////////////////////////////////////////////////////////////

cd_v::cd_v()
{
    cp_pos=CP_ON;
    concavity=0;
	parent=NULL;
	hull_cid=id=-1;
    flag=INIT_FLAG;
}

//find the edge with given v in the list. O(n).
cd_e * cd_v::find(cd_v * v) const{
    FOREACH_CE(edges)
        if( e->v[0]==v || e->v[1]==v ) return e;
    FOREACH_END
    return NULL;
}

///////////////////////////////////////////////////////////////////////////////
//
// cd_e: edge
//
///////////////////////////////////////////////////////////////////////////////

cd_e::cd_e()
{
    f[0]=f[1]=NULL;
    v[0]=v[1]=NULL;
    reflex=false;
    length=0;
    flag=INIT_FLAG;
}

void cd_e::compRef()
{
    reflex=false;
    if( f[1]==NULL ) return;
    cd_v * f0_v=f[0]->otherv(this);
    if(f0_v==NULL){ cd_err("cd_e::compRef"); return;}
    Vector3d vec1=(f0_v->pos-f[1]->v[1]->pos);//.normalize();
    float reflexity=(float)(vec1*f[1]->n);
    reflex=reflexity>SMLL_D;
}

Vector3d cd_e::n()
{
    if( f[1]==NULL ) return f[0]->n;
    Vector3d n=f[0]->n+f[1]->n;
    double nsqr=n.normsqr();
    if( nsqr==0 ){ //how can this happend???
        return Vector3d(0,0,0);
    }
    return n/sqrt(nsqr);
}

void cd_e::update()
{
    compRef();
    compLength();
    //compNormal();
    //flag=INIT_FLAG;
}

///////////////////////////////////////////////////////////////////////////////
//
// cd_f: face
//
///////////////////////////////////////////////////////////////////////////////

cd_f::cd_f()
{
    id=-1;
    v[0]=v[1]=v[2]=NULL;
    e[0]=e[1]=e[2]=NULL;
    flag=INIT_FLAG;
}

void cd_f::compNormal()
{
    Vector3d v1=(v[1]->pos-v[0]->pos);
    Vector3d v2=(v[2]->pos-v[0]->pos);
    n=(v1%v2);
    double l=n.norm();
    if(l==0){
		// cerr << "! Warning: F=" << id << " has zero area: " << endl;;
    //     if( v[1]->pos.almost_equ(v[0]->pos ) )
		// 	cerr << "\t v=" << v[1]->id << "(" << v[1]->pos << ") and " << v[0]->id << "(" << v[0]->pos<<") are too close"<<endl;
    //     if( v[2]->pos.almost_equ(v[0]->pos ) )
		// 	cerr << "\t v=" << v[2]->id << "(" << v[2]->pos << ") and " << v[0]->id << "(" << v[0]->pos << ") are too close" << endl;
    //     if( v[2]->pos.almost_equ(v[1]->pos ) )
		// 	cerr << "\t v=" << v[2]->id << "(" << v[2]->pos << ") and " << v[1]->id << "(" << v[1]->pos << ") are too close" << endl;
        n.set(1,0,0);
    }
    else n=n/l;
}

///////////////////////////////////////////////////////////////////////////////
//
// cd_m: model itself
//
///////////////////////////////////////////////////////////////////////////////

cd_m::cd_m(){ m_decompose_me=true; }

void cd_m::destroy()
{
    FOREACH_V(m_v) delete v; FOREACH_END
    FOREACH_E(m_e) delete e; FOREACH_END
    FOREACH_F(m_f) delete f; FOREACH_END
    m_v.clear(); m_e.clear(); m_f.clear();
}

cd_m * cd_m::clone()
{
    cd_m * other=new cd_m();
    assert(other);

    //create new vertices, edges, faces and holes
    map<cd_v*,cd_v*> v2v; //map from original v to new v
    map<cd_e*,cd_e*> e2e; //map from original e to new e
    map<cd_f*,cd_f*> f2f; //map from original f to new f

    FOREACH_V(m_v) cd_v* nv=new cd_v(); assert(nv); *nv=*v; v2v[v]=nv; other->m_v.push_back(nv); FOREACH_END
    FOREACH_E(m_e) cd_e* ne=new cd_e(); assert(ne); *ne=*e; e2e[e]=ne; other->m_e.push_back(ne); FOREACH_END
    FOREACH_F(m_f) cd_f* nf=new cd_f(); assert(nf); *nf=*f; f2f[f]=nf; other->m_f.push_back(nf); FOREACH_END
    FOREACH_H(m_h) cd_h  nh; nh=h; other->m_h.push_back(nh); FOREACH_END

    //map edges in v to new edges
    FOREACH_V(other->m_v)
        EL tmp;
        FOREACH_E(v->edges)
            tmp.push_back(e2e[e]);
        FOREACH_END
        v->edges.swap(tmp);
        if(v->parent!=NULL) v->parent=v2v[v->parent];
    FOREACH_END

    //map incident vertices and faces to the new ones
    FOREACH_E(other->m_e)
        e->v[0]=v2v[e->v[0]]; e->v[1]=v2v[e->v[1]];
        e->f[0]=f2f[e->f[0]]; e->f[1]=f2f[e->f[1]];
    FOREACH_END

    //map incident vertices and edge to the new ones
    FOREACH_F(other->m_f)
        for(short i=0;i<3;i++){
            f->v[i]=v2v[f->v[i]];
            f->e[i]=e2e[f->e[i]];
        }
    FOREACH_END

    //map vertices in the hole to new vertices
    FOREACH_H(other->m_h)
        VL tmp;
        FOREACH_V(h.vl)
            tmp.push_back(v2v[v]);
        FOREACH_END
        h.vl.swap(tmp);
    FOREACH_END

    //done
    return other;
}

///////////////////////////////////////////////////////////////////////////////
//
// build models from the points and triangles
//
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//init v and push into the list
inline void comp_v(const PtVector& pt, cd_v ** vl,VL& _v, VV* origVs = NULL)
{
    int size=pt.size();
    for( int i=0;i<size;i++ ){
        cd_v * v=new cd_v();
        assert(v!=NULL);
        v->pos=pt[i];
        v->id=i;
        _v.push_back(v);
        vl[i]=v;
    }
}

///////////////////////////////////////////////////////////////////////////////
// find faces and relation between faces, i.e., edges..

typedef list< pair<int,cd_e*> > edge;
inline cd_e * findedge( const edge& e, int id )
{
    typedef edge::const_iterator EIT;
    for( EIT ie=e.begin();ie!=e.end();ie++ )
        if( ie->first==id ) return ie->second;
    return NULL;
}

inline void comp_fe
(const PtVector& pt, const TriVector& tri, cd_v ** v,
 FL& _f, EL& _e, int* org_face_id = NULL, FL* origFaces = NULL)
{
    //init
    edge * edges=new edge[pt.size()]; //edges for each vertex
    if( edges==NULL ){
        cd_err("comp_fe: Not Enough Memory");
        exit(1);
    }

    //loop through each triangle
    FL::const_iterator fit;
    if(origFaces != NULL)
    {
    	fit = origFaces->begin();
    }

    typedef TriVector::const_iterator TIT;
    int fid=0;
    for( TIT it=tri.begin();it!=tri.end();it++,fid++){
        ///////////////////////////////////////////////////////////////////////
        const Tri& t=*it;
        cd_f * f=new cd_f();
        assert(f!=NULL);
        f->id=fid;
        f->org_id=fid;
        // added by Zhonghua 9/25/2013
        if(org_face_id!=NULL) {
        	f->org_id=*org_face_id;
        	org_face_id++;
        }
        for(int iv=0;iv<3;iv++){ f->v[iv]=v[t[iv]]; }
        f->compNormal();

        //Loop through each vertex
        for( int ie=0;ie<3;ie++ ){
            int c_e=t[ie]; int n_e=t[(ie+1)%3];

            //find edge
            cd_e * e=findedge(edges[n_e],c_e);
            if( e==NULL ){ //not created yet. now create it.
                e=new cd_e();
                assert(e!=NULL);
                e->f[0]=f; e->v[0]=v[c_e]; e->v[1]=v[n_e];
                v[c_e]->edges.push_back(e);
                v[n_e]->edges.push_back(e);
                edges[c_e].push_back(pair<int,cd_e*>(n_e,e));
                _e.push_back(e);
                f->e[ie]=e;
            }
            else{
                if( e->f[1]==NULL ){
                    e->f[1]=f;
                    f->e[ie]=e;
                }
                else{
                    //non-manifold...
                    f->e[ie]=e;
                }
            }
        }//end for ie

        //add face
        _f.push_back(f);

        if(origFaces != NULL)
        {
        	++fit;
        }

    }//end for it

    delete [] edges;
}

///////////////////////////////////////////////////////////////////////////////
//Find the type of edges, notch or not
inline void comp_e(EL& _e)
{
    FOREACH_E(_e)
        e->update();
    FOREACH_END
}

///////////////////////////////////////////////////////////////////////////////
inline cd_m *
extractModel(VL& _v,EL& _e,FL& _f)
{
    cd_m * m=new cd_m();

    if(m==NULL){
        cd_err("extractModel: Not Enough Memory");
        exit(1);
    }

    //blacken_V();
	flagT flag=_e.front()->v[0]->flag=generateID();
	propagateFlag(_e.front()->v[0]);

    //collect vertex
    int vid=0;
    //FOREACH_V(_v)
	for(VIT iv = _v.begin();iv!= _v.end();)
	{
		cd_v* v=*iv;
        if( v->flag==flag ){ // same cc
            //VIT cur=iv;//iv--;
            //(*cur)->id=vid++;
            //m->v().push_back(*cur);
            //_v.erase(cur);
			 //VIT cur=iv;//iv--;
            (*iv)->id=vid++;
            m->v().push_back(*iv);
            iv = _v.erase(iv);
        }
		else{
			++iv;
		}
	}
    //FOREACH_END//end for

    //collect edges
    //FOREACH_E(_e)
	for(EIT ie = _e.begin();ie!=_e.end();)
	{
		 cd_e* e=*ie;
        if( e->v[0]->flag==flag ){// same cc
            //EIT cur=ie; //ie--;
            //m->e().push_back(*cur);
            //_e.erase(cur);
            m->e().push_back(*ie);
            ie = _e.erase(ie);
        }
		else
		{
			++ie;
		}
	}
    //FOREACH_END//end for

    //collect faces
    int fid=0;
    //FOREACH_F(_f)
	for(FIT fi = _f.begin();fi!=_f.end();)
	{
		cd_f* f=*fi;
        if( f->v[0]->flag == flag ){// same cc
            //FIT cur=fi; //--;
            //(*cur)->id=fid++;
            //m->f().push_back(*cur);
            //_f.erase(cur);
            (*fi)->id=fid++;
            m->f().push_back(*fi);
           fi = _f.erase(fi);
        }
		else{
			++fi;
		}
	}
    //FOREACH_END//end for

    return m;
}

///////////////////////////////////////////////////////////////////////////////

//use to check if there exists non-manifold vertices.
inline void checkNon2ManifoldV(VL& vl)
{
    FOREACH_V(vl)
        int e_count=0;
        FOREACH_E(v->edges)
            if(e->f[1]==NULL)
                e_count++;
        FOREACH_END
        if( e_count>2 ) cout<<"v="<<v->id<<" is not 2 manifold: has "<<e_count<<" B edges"<<endl;
    FOREACH_END
}

void buildModels(const list<FL>& _FL, ML& ms)
{
	for(list<FL>::const_iterator it = _FL.begin(); it != _FL.end(); ++it)
	{
		flagT FLAG = generateID();
		PtVector pt;
		TriVector tri;

		FL origFs;//add by Guilin 10/10/2013
		VV origVs;//add by Guilin 10/10/2013

		int* org_face_ids = new int[it->size()];

		int vid = 0;
		int face_index = 0;
		FOREACH_CF((*it))
			Tri triangle;
			for(int i=0;i<3;i++)
			{
				if(f->v[i]->flag == FLAG) {
					triangle[i] = f->v[i]->id;
				}
				else
				{
					f->v[i]->flag = FLAG;
					f->v[i]->id = vid++;
					pt.push_back(f->v[i]->pos);
					triangle[i] = f->v[i]->id;
				}
			}

			org_face_ids[face_index++] = f->org_id;
			tri.push_back(triangle);
		FOREACH_END

//		//remeshing add by Guilin 9/26
//		vector<int> vecfaceids;
//		int initSize = tri.size();
//		remesh(*it, tri);
//		int afterSize = tri.size();
//		int * org_face_ids = new int[afterSize];
//		for(int k = 0; k < initSize; k++)
//		{
//			org_face_ids[k] = org_face_idsX[k];
//		}
//		for(int k = initSize; k < afterSize; k++)
//		{
//			org_face_ids[k] = -100;
//		}

		VL _v; //vertices
		int ptsize=pt.size();
		cd_v ** v=(cd_v**)calloc(ptsize,sizeof(cd_v*));
		if( v==NULL ){
			cd_err("buildModels: Not Enough Memory");
			exit(1);
		}
		memset(v,0,ptsize);
		comp_v(pt,v,_v, &origVs);

		EL _e; //edges
		FL _f; //faces

		comp_fe(pt,tri,v,_f,_e, org_face_ids, &origFs);
		comp_e(_e);

		//classify into vertices
		while( !_e.empty() && !_f.empty() ){

			cd_m * m=extractModel(_v,_e,_f);
			if( m!=NULL ){
        if(isflat(m)){ m->destroy(); delete m; continue; }
				identifyHoles(m);
				ms.push_back(m);
				//identifyBound(m);
			}
		}//end while

	   delete [] v;
	   delete [] org_face_ids;
	   //delete [] org_face_idsX;
	}
}


void buildModels
(cd_v ** v, VL& _v, const PtVector& pt, const TriVector& tri, ML& ms)
{
    EL _e; //edges
    FL _f; //faces

    comp_fe(pt,tri,v,_f,_e);
    comp_e(_e);

    //classify into vertices
    while( !_e.empty() && !_f.empty() ){
        cd_m * m=extractModel(_v,_e,_f);
        if( m!=NULL ){
            if(isflat(m)){ m->destroy(); delete m; continue; }
            identifyHoles(m);
            ms.push_back(m);
        }
    }//end while
}

void buildModels
(const PartVector& part, const PtVector& pt, const TriVector& tri, ML& ms)
{
    ///////////////////////////////////////////////////////////////////////////
    VL _v; //vertices
    int ptsize=pt.size();
    cd_v ** v=(cd_v**)calloc(ptsize,sizeof(cd_v*));
    if( v==NULL ){
        cd_err("buildModels: Not Enough Memory");
        exit(1);
    }
    memset(v,0,ptsize);
    comp_v(pt,v,_v);

    ///////////////////////////////////////////////////////////////////////////
    if( part.size()==1 ){
        buildModels(v,_v,pt,tri,ms);
    }
    else{
        typedef PartVector::const_iterator CPIT;
        int beginID=0;
        for( CPIT ip=part.begin();ip!=part.end();ip++ ){
            //cd_out(cd_msg("Create Model (")+(ip-part.begin())+")");
            int endID=beginID+ip->second;
            TriVector part_tri(tri.begin()+beginID,tri.begin()+endID);
            buildModels(v,_v,pt,part_tri,ms);
            beginID=endID;
        }
    }
    delete [] v;
}

///////////////////////////////////////////////////////////////////////////////

inline void gethole(cd_v * v, cd_h& h)
{
    VL open;
    EL el;
    open.push_back(v);
    flagT flag=v->flag;
    while(!open.empty()){
        cd_v * v=open.front();
        open.pop_front();
        h.vl.push_back(v);
        for(int i=0;i<3;i++) h.com[i]+=v->pos[i];
        //explore
        FOREACH_E(v->edges)
            if( e->flag!=flag ) continue; //not a hole bd
            cd_v * o=e->otherv(v);
            if( o->flag==flag ) continue; //visited...
            h.length+=e->length;
            o->flag=flag;
            open.push_back(o);
            el.push_back(e);
        FOREACH_END
    }//end while
    //compute com of hole
    int vsize=h.vl.size();
    for(int i=0;i<3;i++) h.com[i]/=vsize;
}

//go through all edges in m and identify connected edges
//as the boundaries of holes
void identifyHoles( cd_m * m )
{
    flagT holeFlag=generateID();

    VL hole_v;
    FOREACH_E(m->e())
        if(e->f[1]==NULL){
            e->flag=holeFlag;
            hole_v.push_back(e->v[0]);
            hole_v.push_back(e->v[1]);
        }
    FOREACH_END
    if( hole_v.empty() ) return;

    while(!hole_v.empty()){
        cd_v * v=hole_v.front();
        hole_v.pop_front();
        if(v->flag==holeFlag) continue;
        //explore
        v->flag=holeFlag;
        m->holes().push_back(cd_h());
        gethole(v,m->holes().back());
    }
    m->holes().sort(greater<cd_h>());
}

///////////////////////////////////////////////////////////////////////////////

cd_msg cd_f::to_msg()
{
    cd_msg msg=cd_msg("f id=")+id+" v=[";
    for(int iv=0;iv<3;iv++){
        if(v[iv]!=NULL) msg=msg+v[iv]->id;
        else msg+="NULL";
        if( iv<2 ) msg+=",";
    }
    msg+="]";

    return msg;
}

cd_msg cd_e::to_msg()
{
    cd_msg msg("e ");
    if( reflex ) msg+="r v=[";
    else msg+="n v=[";

    if( v[0]!=NULL && v[1]!=NULL ){
        msg=msg+v[0]->id+","+v[1]->id;
    }
    else if( v[0]!=NULL ){
        msg=msg+v[0]->id+",NULL";
    }
    else if( v[1]!=NULL ){
        msg=msg+"NULL,"+v[1]->id;
    }
    else{
        msg+="NULL,NULL";
    }

    msg+="] f=[";

    if( f[0]!=NULL && f[1]!=NULL ){
        msg=msg+f[0]->id+","+f[1]->id;
    }
    else if( f[0]!=NULL ){
        msg=msg+f[0]->id+",NULL";
    }
    else if( f[1]!=NULL ){
        msg=msg+"NULL,"+f[1]->id;
    }
    else{
        msg+="NULL,NULL";
    }
    msg+="]";

    return msg;
}

cd_msg cd_v::to_msg()
{
    return cd_msg("v ")+id+" pos="+pos;
}
