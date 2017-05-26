#ifndef _CD_HULL2D_H_
#define _CD_HULL2D_H_

///////////////////////////////////////////////////////////////////////////////
// QHull Headers
extern "C"{
#include <stdio.h>
#include "qhull.h"
#include "poly.h"
#include "qset.h"
}

///////////////////////////////////////////////////////////////////////////////

#include <cassert>
using namespace std;

#include "util/mathtool/Point.h"
using namespace mathtool;

#include <list>
using namespace std;

///////////////////////////////////////////////////////////////////////////////
class hull_2d{
public:
    hull_2d(){ ptsize=0; pt=NULL; }
    ~hull_2d(){ free(); }
    
    void reset(){ 
        ptsize=0; //fsize=0; 
    }
    
    bool alloc(int size){
        size=size*2; //for 3 coord
        pt=new coordT[size];
        return (pt!=NULL);
    }
    
    void free(){
        ptsize=0;
        delete [] pt; pt=NULL; 
    }
    
    void addpt(const Point2d& pt){ 
        for( int i=0;i<2;i++ ) this->pt[ptsize*2+i]=pt[i]; 
        ptsize++;
    }
    
    void buildhull(){
        //global varible for qt
        int curlong, totlong;
        
        //using qhull
        static char options[]="qhull Qs Qx QJ i Tcv C-0 Pp";
        qh_init_A(stdin, stdout, stderr, 0, NULL);
        qh_initflags(options);
        qh_init_B (pt, ptsize, 2, false);
        qh_qhull();
        qh_check_output();
        
        //find bridge
        findBridges();
        
        //free mem
        qh_freeqhull(!qh_ALL);
        qh_memfreeshort (&curlong, &totlong);
    }
    
    /**
     * Get Bridge
     */
    const list< pair<int,int> > & GetBridge() const { return bridge; }
       
    ///////////////////////////////////////////////////////////////////////////
    //Data
private:
    
    ///////////////////////////////////////////////////////////////////////////
    void findBridges()
    {
        //global variable for Qhull
        facetT *facet;
        vertexT *vertex;
        vertexT **vertexp;
        //setT *vertices;

        FORALLfacets {

            list<int> ids;
            FOREACHvertex_(facet->vertices)
            {
                ids.push_back(qh_pointid(vertex->point));
            }

            assert(ids.size()==2);

            if( abs(ids.front()-ids.back())!=1 )
            {
                if(facet->toporient)
                    bridge.push_back(make_pair(ids.front(),ids.back()));
                else
                    bridge.push_back(make_pair(ids.back(),ids.front()));
            }
        }//end FORALLfacets

    }//end findBridges

    ///////////////////////////////////////////////////////////////////////////
    //Input
    int ptsize;  //number of point
    coordT * pt; //points in the patch (*3 for each coord)
    
    //Output
    list< pair<int,int> > bridge;
};
    
#endif //_CD_HULL2D_H_
    

