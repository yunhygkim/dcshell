// IPolygonal.h: interface for the IDisplay class.
//
//////////////////////////////////////////////////////////////////////

#ifndef _IPOLYGONAL_H_
#define _IPOLYGONAL_H_

//////////////////////////////////////////////////////////////////////
//  This is an interface for all classes which constians
//  polygonal data. For example, deformable object, sculpture deformer
//  and evene model loader
//////////////////////////////////////////////////////////////////////

#include "util/mathtool/Point.h"
#include "util/mathtool/Vector.h"
using namespace mathtool;

#include <vector>
using namespace std;

typedef vector<Point3d> PtVector;
typedef Vector<int, 3> Tri;
typedef vector<Tri> TriVector;
typedef pair<int,int> Part;
typedef vector<Part>  PartVector;

class IPolygonal  
{
public:

//////////////////////////////////////////////////////////////////////
//  Interface for Retrive general polgons
//////////////////////////////////////////////////////////////////////

	virtual int GetVertices(PtVector& v) const =0;
	virtual int GetTriangles(TriVector& v) const =0;
	virtual const PtVector & GetVertices() const =0;
	virtual const TriVector & GetTriangles() const =0;
};

#endif // !defined(_IPOLYGONAL_H_)
