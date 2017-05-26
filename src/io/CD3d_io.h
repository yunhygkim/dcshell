#ifndef _CD3D_IO_H_
#define _CD3D_IO_H_

#include <iostream>
#include "CD3d_model.h"

#define F_OUT_PRECISION 12

///////////////////////////////////////////////////////////////////////////////
//
//
// Save the model into BYU files
//
//
///////////////////////////////////////////////////////////////////////////////

inline
bool SaveBYU(ostream& out, cd_m * m)
{
    out<<1<<" "<<m->v().size()<<" "<<m->f().size()<<" "<<m->f().size()*3<<"\n";
    out<<1<<" "<<m->f().size()<<"\n";

    FOREACH_V(m->v())
        Point3d& pos=v->pos;
        out<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<"\n";
    FOREACH_END

    FOREACH_F(m->f())
        out<<(1+f->v[0]->id)<<" "<<(1+f->v[1]->id)<<" "<<-(1+f->v[2]->id)<<"\n";
    FOREACH_END
    return true;
}


inline 
bool SaveBYU(const string& filename, cd_m * m)
{
    //open file
    ofstream fout(filename.c_str());
    fout.precision(F_OUT_PRECISION);
    if( fout.good()==false ){
        cd_err(cd_msg("Can't Open file : ")+filename);
        return false;
    }
    
    bool r=SaveBYU(fout,m);
    fout.close();
    return r;
}

inline
bool SaveBYU(ostream& out, ML& models)
{
    long total_v=0, total_f=0;
    //count
    FOREACH_M(models)
        total_v+=m->v().size();
        total_f+=m->f().size();
    FOREACH_END
    out<<models.size()
       <<" "<<total_v<<" "<<total_f<<" "<<total_f*3<<"\n";
    
    int id=1;
    FOREACH_M(models)
        out<<id++<<" "<<m->f().size()<<"\n";
    FOREACH_END

    FOREACH_M(models)
        FOREACH_V(m->v())
            Point3d& pos=v->pos;
            out<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<"\n";
        FOREACH_END
    FOREACH_END

    int vid=0;
    FOREACH_M(models)
        FOREACH_F(m->f())
            out<< (1+f->v[0]->id+vid)<<" "
               << (1+f->v[1]->id+vid)<<" "
               <<-(1+f->v[2]->id+vid)<<"\n";
        FOREACH_END
        vid+=m->v().size();
    FOREACH_END

    return true;
}

inline 
bool SaveBYU(const string& filename, ML& models)
{
    //open file
    ofstream fout(filename.c_str());
    fout.precision(F_OUT_PRECISION);
    if( fout.good()==false ){
        cd_err(cd_msg("Can't Open file : ")+filename);
        return false;
    }
    
    bool r=SaveBYU(fout,models);
    fout.close();
    return r;
}

///////////////////////////////////////////////////////////////////////////////
//
//
// Save CD_CM_DATA
//
//
///////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////////////////
////
//// READ FROM A ISTREAM
////
/////////////////////////////////////////////////////////////////////////////////
//inline bool readCMData(cdm_data * data, istream& in, cd_cr * cr)
//{
//    int size=0;
//    int pid=0;
//    int vid=0;
//    double tmp;
//    //check signature
//    in>>tmp; if( fabs(tmp-get_cd3d()->m_tau)>1e-10 ) return false;
//    in>>size; if(size!=qh num_facets) return false;
//    ///////////////////////////////////////////////////////////////////////////
//    // read handle cut data
//    in>>size;
//    for(int ihc=0;ihc<size;ihc++){
//        int vsize;
//        in>>vsize;
//        data->handle_cuts.push_back(cd_cut());
//        VL& vl=data->handle_cuts.back().vl;
//        for(int iv=0;iv<vsize;iv++){
//            in>>vid; vl.push_back(data->V(vid).v);
//        }//end for iv
//    }//enf for ihc
//
//    if( !data->handle_cuts.empty() )
//	{
//        cr->resolveConcavity(data->m,data->handle_cuts);
//        data->resize();
//    }
//
//    ///////////////////////////////////////////////////////////////////////////
//    // read v_data
//    in>>size;
//    if( size!=data->m_vsize ) return false;
//    for(int iv=0;iv<size;iv++) in>>data->vdata[iv].flags;
//    ///////////////////////////////////////////////////////////////////////////
//    // read patches data
//    in>>size;
//    cd_patch ** phs=(cd_patch **)calloc(size,sizeof(cd_patch*));
//    for(int ip=0;ip<size;ip++){
//        phs[ip]=new cd_patch();
//        phs[ip]->id=ip;
//        in>>phs[ip]->o[0]>>phs[ip]->o[1]>>phs[ip]->o[2]>>
//            phs[ip]->n[0]>>phs[ip]->n[1]>>phs[ip]->n[2];
//        data->patches.push_back(phs[ip]);
//    }
//    ///////////////////////////////////////////////////////////////////////////
//    // read f2p data
//    in>>size;
//    if( size!=qh num_facets ) return false;
//    facetT ** fs=(facetT **)calloc(size,sizeof(facetT*));
//    facetT *facet;
//    FORALLfacets {
//        in>>pid;
//        data->f2p[facet->id].patch=phs[pid];
//        phs[pid]->facets.push_back(facet);
//        fs[facet->id]=facet;
//    }
//    ///////////////////////////////////////////////////////////////////////////
//    // read bedges data
//    in>>size;
//    for(int ib=0;ib<size;ib++){
//        b_edge * be=new b_edge(0);
//        int f1,f2,v1,v2,vsize;
//        in>>f1>>f2>>v1>>v2>>be->length>>vsize;
//        be->f1=fs[f1]; be->f2=fs[f2];
//        be->v1=data->V(v1).v; be->v2=data->V(v2).v;
//        for(int iv=0;iv<vsize;iv++){
//            in>>vid; be->path.push_back(data->V(vid).v);
//        }
//        data->getPatch(fs[f1])->bedges.push_back(be);
//        data->getPatch(fs[f2])->bedges.push_back(be);
//        data->bedges.push_back(be);
//    }
//    ///////////////////////////////////////////////////////////////////////////
//    // free mem
//    free(fs);
//    free(phs);
//
//    return true;
//}
//
//inline bool saveCMData(cdm_data * data, ostream& out)
//{
//    //signature...
//    out<<get_cd3d()->m_tau<<" "<<qh num_facets<<"\n";
//    ///////////////////////////////////////////////////////////////////////////
//    // save handle cut data
//    out<<data->handle_cuts.size()<<"\n";
//    FOREACH_CUT(data->handle_cuts)
//        out<<cut.vl.size()<<" ";
//        FOREACH_V(cut.vl)
//            out<<v->id<<" ";
//        FOREACH_END
//        out<<"\n";
//    FOREACH_END
//    ///////////////////////////////////////////////////////////////////////////
//    // save v_data
//    out<<data->m_vsize<<"\n";
//    for(int iv=0;iv<data->m_vsize;iv++) out<<data->vdata[iv].flags<<" ";
//    out<<"\n";
//    ///////////////////////////////////////////////////////////////////////////
//    // save patches data
//    out<<data->patches.size()<<"\n";
//    FOREACH_P(data->patches)
//        const Point3d& o=p->o;
//        const Vector3d& n=p->n;
//        out<<o[0]<<" "<<o[1]<<" "<<o[2]<<" "<<n[0]<<" "<<n[1]<<" "<<n[2]<<"\n";
//    FOREACH_END
//    ///////////////////////////////////////////////////////////////////////////
//    // save f2p data
//    out<<qh num_facets<<"\n";
//    for(int it=0;it<qh num_facets;it++) out<<data->f2p[it].patch->id<<" ";
//    out<<"\n";
//    ///////////////////////////////////////////////////////////////////////////
//    // save bedges data
//    out<<data->bedges.size()<<"\n";
//    FOREACH_B(data->bedges)
//        out<<be->f1->id<<" "<<be->f2->id<<" "<<be->v1->id<<" "<<be->v2->id
//           <<" "<<be->length<<" "<<be->path.size()<<" ";
//        FOREACH_V(be->path)
//            out<<v->id<<" ";
//        FOREACH_END
//        out<<"\n";
//    FOREACH_END
//    //done
//    return true;
//}
//
//inline bool readCMData(cdm_data * data, const string& name,cd_cr * cr)
//{
//    ifstream fin(name.c_str());
//    if( !fin.good() )
//        return false;
//    bool r=readCMData(data,fin,cr);
//    fin.close();
//    return r;
//}
//
//inline bool saveCMData(cdm_data * data, const string& name)
//{
//    ofstream fout(name.c_str());
//    fout.precision(F_OUT_PRECISION);
//    if( !fout.good() )
//        return false;
//    bool r=saveCMData(data,fout);
//    fout.close();
//    return r;
//}
//
//inline bool saveSegData(cd_m * model, ML components, ostream& out)
//{
//	int cid = 0;
//
//	FOREACH_M(components)
//		FOREACH_CV(m->v())
//			v->component_id = cid;
//		FOREACH_END
//		cid++;
//	FOREACH_END
//
//	FOREACH_CV(model->v())
//		out<<v->component_id<<endl;
//	FOREACH_END
//	
//	return true;
//}
//
//inline bool saveSegData(cd_m * model, ML components, const string& name)
//{
//	ofstream fout(name.c_str());
//	if( !fout.good() )
//		return false;
//	bool r=saveSegData(model, components, fout);
//	fout.close();
//	return r;
//}

#endif//_CD3D_IO_H_
