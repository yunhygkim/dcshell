//------------------------------------------------------------------------------
//  Copyright 2007-2009 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#ifndef _OBJ_READER_H_
#define _OBJ_READER_H_

#include <string>
#include <fstream>
#include <iostream>
#include <list>
#include <vector>
#include <ctype.h>
#include <math.h>

using namespace std;

#include "util/mathtool/Point.h"
#include "util/mathtool/Vector.h"
using namespace mathtool;

#include "util/model/ILoadable.h"

class objReader : public ILoadable
{
public:

    //////////////////////////////////////////////////////////////////////////////////////
    // Implemetation of ILoadable interface
    //////////////////////////////////////////////////////////////////////////////////////
    virtual bool ParseFile()
    {
        if( CheckCurrentStatus()==false ) return false;
        ifstream in(m_strFileName);
        bool r=Read(in);
        in.close();
        return r;
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // Implemetation of IPolygonal interface
    //////////////////////////////////////////////////////////////////////////////////////
    int GetVertices(PtVector& v) const{ v=points; return v.size(); }
    int GetTriangles(TriVector& v) const{ v=triangles; return v.size(); }
    const PtVector & GetVertices() const{ return points; }
    const TriVector & GetTriangles() const{ return triangles; }
    const PartVector& GetParts() const { return parts; }

    //////////////////////////////////////////////////////////////////////////////////////
    //
    //
    //  Protected Member Functions
    //
    //
    //////////////////////////////////////////////////////////////////////////////////////


protected:

    bool CheckCurrentStatus()
    {
        if( m_strFileName==NULL )
            return false;

        //Check if file exist
        ifstream fin(m_strFileName, ios::in);
        if( !fin.good() )
        {   //not good. file not found
            cerr<<"! Error: File "<< m_strFileName <<" not found"<<endl;
            return false;
        }

        fin.close();

        return true;
    }

    bool Read(istream& in){

        string tmp;
        
        //read pts
        while(true)
        {
            in>>tmp;
            if( tmp=="f" ) break;
            if( tmp=="v" ){
                Point3d pt;
                in>>pt[0]>>pt[1]>>pt[2];
                points.push_back(pt);
            }

            getline(in,tmp);
        }

        //read faces
        list<int> poly;
        do{

            in>>tmp;
            if( in.eof() ) break;
            if( isdigit(tmp[0]) ){ //this defines a vetex
                
                int pos1=tmp.find('/');
                int pos2=tmp.rfind('/');

                int id_v=atoi(tmp.substr(0,pos1).c_str())-1;
                int id_n=atoi(tmp.substr(pos2+1).c_str())-1;

                poly.push_back(id_v);
            }
            else if( tmp=="f" )
            {
                if(!poly.empty()){
                    //create a triangle
                    if(poly.size()!=3){
                        cerr<<"! Warning: non-triangle face found and ignored"<<endl;
                    }
                    else{
                        triangles.push_back(list2Tri(poly));
                    }
                }
                poly.clear();
            }
            else 
                getline(in,tmp);

        }
        while( !in.eof() );

        if(!poly.empty()){
            //create a triangle
            if(poly.size()!=3){
                cerr<<"! Warning: non-triangle face found and ignored"<<endl;
            }
            else{
                triangles.push_back(list2Tri(poly));
            }
        }

        //create part (only one part anyway...)
        Part part(1,triangles.size());
        parts.push_back(part);

        return true;
    }

    Tri list2Tri(const list<int>& poly)
    {
        int id=0;
        Tri tri;
        for(list<int>::const_iterator i=poly.begin();i!=poly.end();i++, id++)
        {
            tri[id]=*i;
        }//end i
        return tri;
    }

    string m_filename;

    PtVector points;
    TriVector triangles;
    PartVector parts;

};

#endif //_OBJ_READER_H_



