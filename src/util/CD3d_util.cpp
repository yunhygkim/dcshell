#include "util/CD3d_util.h"

///////////////////////////////////////////////////////////////////////////////
//
//
//  Compute Left/Right edges/facets of a path
//
//
///////////////////////////////////////////////////////////////////////////////

//return the last facet visited
inline cd_f * traceNextEdge
(cd_f * f, cd_e * e, cd_v * v, cd_f * endf, EL& edges, FL& facets /*return*/)
{
    if( f==NULL || f==endf ) return f;
    facets.push_back(f);
    //find the edge of f that contains v but not e
    //e but contains v
    cd_e * ne=NULL;
    for(int i=0;i<3;i++){ //check each edge of f
        if( f->e[i]==e ) continue;
        if( f->e[i]->v[0]==v || f->e[i]->v[1]==v ){ //found the next edge
            ne=f->e[i];
            break;
        }
    }
    assert(ne); //make sure be is not null
    edges.push_back(ne);
    return traceNextEdge(ne->otherf(f),ne,v,endf,edges,facets);
}

void getSideEF
(cd_v *v, cd_e * e[2], cd_f *rf[2], cd_f *lf[2], 
 EL& rE, EL& lE, FL& rF, FL& lF)
{
    cd_f * lastL=NULL, *lastR=NULL;
    if(lf[0]!=NULL){ lastL=traceNextEdge(lf[0],e[0],v,lf[1],lE,lF); }
    if(lf[1]!=NULL){ 
        if(lastL==NULL) //there is an empty edge between lf[0] and lf[1]
            traceNextEdge(lf[1],e[1],v,lf[0],lE,lF);
        else{ assert(lastL==lf[1]); lF.push_back(lastL); }
    }
    if(rf[0]!=NULL){ lastR=traceNextEdge(rf[0],e[0],v,rf[1],rE,rF); }
    if(rf[1]!=NULL){ 
        if(lastR==NULL) //there is an empty edge between lf[0] and lf[1]
            traceNextEdge(rf[1],e[1],v,rf[0],rE,rF);
        else{ assert(lastR==rf[1]); rF.push_back(lastR); }
    }
}

void getSideEF
(cd_v * v1,cd_v * v2,cd_v * v3,EL& rE, EL& lE, FL& rF, FL& lF)
{
    cd_e * e[2]={v1->find(v2),v2->find(v3)};
    pair<cd_f*,cd_f*> f1=orderFacets(v1,v2,e[0]);
    pair<cd_f*,cd_f*> f2=orderFacets(v2,v3,e[1]);
    cd_f *lf[2]={f1.first,f2.first};
    cd_f *rf[2]={f1.second,f2.second};
    getSideEF(v2,e,rf,lf,rE,lE,rF,lF);
}


//note: in order to get all edges/facet on the side, this model
//need to be manifold

void getSideEF
(VIT begin, VIT end, EL& rE, EL& lE, FL& rF, FL& lF)
{
    VIT s=begin; s++;
    VIT e=end;   e--;
    FOREACH_Vi(s,e)
        VIT pre=iv;  pre--;
        VIT next=iv; next++;
        getSideEF(*pre,v,*next,rE,lE,rF,lF);
    FOREACH_END
}

void getSideEdges(VIT begin, VIT end, EL& rE, EL& lE)
{
    FL tmp1,tmp2;
    getSideEF(begin,end,rE,lE,tmp1,tmp2);
}

