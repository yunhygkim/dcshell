#include "util/CD3d_convexity.h"
#include "measure/CD3d_hull.h"

double vol
(const Point3d& p1,const Point3d& p2,const Point3d& p3,const Point3d& p4)
{
    return
    -(p1[2]-p4[2])*(p2[1]-p4[1])*(p3[0]-p4[0])
    +(p1[1]-p4[1])*(p2[2]-p4[2])*(p3[0]-p4[0])
    +(p1[2]-p4[2])*(p2[0]-p4[0])*(p3[1]-p4[1])
    -(p1[0]-p4[0])*(p2[2]-p4[2])*(p3[1]-p4[1])
    -(p1[1]-p4[1])*(p2[0]-p4[0])*(p3[2]-p4[2])
    +(p1[0]-p4[0])*(p2[1]-p4[1])*(p3[2]-p4[2]);
}

double vol(cd_h& h)
{
    //collected edges
    EL el; //edges extracted fromt h.vl
    flagT flag=generateID();
    FOREACH_V(h.vl)
        FOREACH_E(v->edges)
            if(e->f[1]==NULL&&e->flag!=flag)
            {e->flag=flag; el.push_back(e);}
        FOREACH_END
    FOREACH_END

    //compute vol
    double volume=0;
    const Point3d O(0,0,0);

    FOREACH_E(el)
        cd_v * ov=e->f[0]->otherv(e);
        const Point3d& p1=e->v[0]->pos;
        const Point3d& p2=e->v[1]->pos;
        Vector3d v=(p2-p1)%(ov->pos-p1);

        if(v*e->f[0]->n>0)
            volume+=vol(p2,p1,h.com,O);
        else
            volume+=vol(p1,p2,h.com,O);
    FOREACH_END

    return volume;
}

double vol(HL& hl)
{
    double volume=0;
    FOREACH_H(hl)
        volume+=vol(h);
    FOREACH_END
    return volume;
}

double vol(FL& fl)
{
    double volume=0;
    const Point3d O(0,0,0);
    FOREACH_F(fl)
        volume+=vol(f->v[0]->pos,f->v[1]->pos,f->v[2]->pos,O);
    FOREACH_END
    return volume;
}

double vol(cd_hull& hull)
{
    //global varible for qull
    facetT *facet;
    vertexT *vertex, **vertexp;
    setT *vertices;
    const Point3d O(0,0,0);
    double volum=0;
    FORALLfacets {
        vertices= qh_facet3vertex (facet);
        PtVector ptV;
        FOREACHvertex_(vertices) {
            Point3d pt; 
            pt.set(vertex->point);
            ptV.push_back(pt);
        }
        qh_settempfree(&vertices);
        volum+=vol(ptV[0],ptV[1],ptV[2],O);
    }//FORALLfacets
    return volum;
}

double convexity(cd_m * m)
{
    if(m->f().size()==1) return 1;
    //
    cd_hull hull;
    hull.buildhull(m);
    //
    double m_vol=vol(m->f())+vol(m->holes());
    double h_vol=vol(hull);
    //
    return m_vol/h_vol;
}
