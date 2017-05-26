/*
 * CD3dexport.cpp
 *
 *  Created on: Oct 16, 2013
 *      Author: guilin
 */

#include "io/CD3dexport.h"
#include "measure/CD3d_hull.h"
#include <iostream>
#include <fstream>
#include <map>
#include <sstream>
#include <assert.h>

ExportOption::ExportOption()
{
	this->b_simply_convex_hull = 0;
	this->simply_convex_hull_ratio = 0.02;
	this->simply_convex_hull_offset = 0.01;
}


    Material::Material(void)
    {
        m_diffuseColor[0]  = 0.5;
        m_diffuseColor[1]  = 0.5;
        m_diffuseColor[2]  = 0.5;
        m_specularColor[0] = 0.5;
        m_specularColor[1] = 0.5;
        m_specularColor[2] = 0.5;
        m_ambientIntensity  = 0.4;
        m_emissiveColor[0] = 0.0;
        m_emissiveColor[1] = 0.0;
        m_emissiveColor[2] = 0.0;
        m_shininess         = 0.4;
        m_transparency      = 0.0;
    }
    void generateRandMaterial(Material& mat)
    {
    	while (mat.m_diffuseColor[0] == mat.m_diffuseColor[1] ||
           mat.m_diffuseColor[2] == mat.m_diffuseColor[1] ||
           mat.m_diffuseColor[2] == mat.m_diffuseColor[0]  )
    	{
    		mat.m_diffuseColor[0] = (rand()%100) / 100.0;
    		mat.m_diffuseColor[1] = (rand()%100) / 100.0;
    		mat.m_diffuseColor[2] = (rand()%100) / 100.0;
    	}
    }
    void savePrepare(int nV, int nT, ofstream& fout, const Material& material)
    {
        fout <<"#VRML V2.0 utf8" << std::endl;
        fout <<"" << std::endl;
        fout <<"# Vertices: " << nV << std::endl;
        fout <<"# Triangles: " << nT << std::endl;
        fout <<"" << std::endl;
        fout <<"Group {" << std::endl;
        fout <<"	children [" << std::endl;
        fout <<"		Shape {" << std::endl;
        fout <<"			appearance Appearance {" << std::endl;
        fout <<"				material Material {" << std::endl;
        fout <<"					diffuseColor "      << material.m_diffuseColor[0]      << " "
                                                        << material.m_diffuseColor[1]      << " "
                                                        << material.m_diffuseColor[2]      << std::endl;
        fout <<"					ambientIntensity "  << material.m_ambientIntensity      << std::endl;
        fout <<"					specularColor "     << material.m_specularColor[0]     << " "
                                                        << material.m_specularColor[1]     << " "
                                                        << material.m_specularColor[2]     << std::endl;
        fout <<"					emissiveColor "     << material.m_emissiveColor[0]     << " "
                                                        << material.m_emissiveColor[1]     << " "
                                                        << material.m_emissiveColor[2]     << std::endl;
        fout <<"					shininess "         << material.m_shininess             << std::endl;
        fout <<"					transparency "      << material.m_transparency          << std::endl;
        fout <<"				}" << std::endl;
        fout <<"			}" << std::endl;
        fout <<"			geometry IndexedFaceSet {" << std::endl;
        fout <<"				ccw TRUE" << std::endl;
        fout <<"				solid TRUE" << std::endl;
        fout <<"				convex FALSE" << std::endl;
    }

    void saveOne(cd_m* m,ofstream& fout, const Material& material)
    {
    	 //if (fout.is_open())
    	        //{
    	            size_t nV = m->v().size();//m_vertices.GetSize();
    	            size_t nT = m->f().size();//m_triangles.GetSize();

    	            savePrepare(nV, nT, fout, material);

    	            //save vertices
	                fout <<"				coord DEF co Coordinate {" << std::endl;
	                fout <<"					point [" << std::endl;
	                FOREACH_V(m->v())
	                		fout <<"						" <<v->pos[0]<< " "
	                    	     	 	 	 	 	 	 	 << v->pos[1]<< " "
	                    	     	 	 	 	 	 	 	 << v->pos[2] << "," << std::endl;
	                FOREACH_END

	                fout<<"					]"<<std::endl;
	                fout<<"					}"<<std::endl;

	                //save faces
	                fout <<"				coordIndex [ " << std::endl;
	                FOREACH_F(m->f())
	                fout <<"						" << f->v[0]->id << ", "
	                    	                          << f->v[1]->id << ", "
	                    	                          << f->v[2]->id << ", -1," << std::endl;
	                FOREACH_END


	                fout<<"				]" << std::endl;
    	            fout <<"			}" << std::endl;
    	            fout <<"		}" << std::endl;
    	            fout <<"	]" << std::endl;
    	            fout <<"}" << std::endl;

//
//    	           // if (GetNVertices() > 0) {
//    	                fout <<"				coord DEF co Coordinate {" << std::endl;
//    	                fout <<"					point [" << std::endl;
//    	                for(size_t v = 0; v < nV; v++)
//    	                {
//    	                    TMMVertex & currentVertex = m_vertices.GetData();
//    	                    fout <<"						" << currentVertex.m_pos.X() << " "
//    	                                                      << currentVertex.m_pos.Y() << " "
//    	                                                      << currentVertex.m_pos.Z() << "," << std::endl;
//    	                    currentVertex.m_id = v;
//    	                    m_vertices.Next();
//    	                }
//    	                fout <<"					]" << std::endl;
//    	                fout <<"				}" << std::endl;
//    	            //}
//    	            //if (GetNTriangles() > 0) {
//    	                fout <<"				coordIndex [ " << std::endl;
//    	                for(size_t f = 0; f < nT; f++)
//    	                {
//    	                    TMMTriangle & currentTriangle = m_triangles.GetData();
//    	                    fout <<"						" << currentTriangle.m_vertices[0]->GetData().m_id << ", "
//    	                                                      << currentTriangle.m_vertices[1]->GetData().m_id << ", "
//    	                                                      << currentTriangle.m_vertices[2]->GetData().m_id << ", -1," << std::endl;
//    	                    m_triangles.Next();
//    	                }
//    	                fout <<"				]" << std::endl;
//    	           //}
//    	            fout <<"			}" << std::endl;
//    	            fout <<"		}" << std::endl;
//    	            fout <<"	]" << std::endl;
//    	            fout <<"}" << std::endl;
    	       // }
    }
    
	void saveObj(cd_m* m, ofstream& fout)
	{
		for (auto& v : m->v())
		{
			fout << std::setprecision(10) << "v " << v->pos << "\n";
		}

		for (auto& f : m->f())
		{
			auto area = f->area();

			if (area < 1e-3) 
			{
				cout << "tiny triangle has area=" << area << endl;
				auto v1 = (f->v[1]->pos - f->v[0]->pos).normalize();
				auto v2 = (f->v[2]->pos - f->v[1]->pos).normalize();
				auto v3 = (f->v[0]->pos - f->v[2]->pos).normalize();
				 
			}

			if (area != 0) //ignore triangle with 0 area...
				fout << "f " << f->v[0]->id + 1 << " " << f->v[1]->id + 1 << " " << f->v[2]->id + 1 << "\n";
			else
				cout << "! Warning: Face has 0 area; ignored" << endl;
		}
		
	}


    void getHullInfo(map<int,int>& fv2v, map<int,int>& v2fv, cd_hull& hull, vector<cd_v*>& vs)
    {
    	int initvn = vs.size();

    	int f = hull.getFaceSize();
    	int vnum = hull.getFaceSize() * 3;
    	bool* b = new bool[initvn];
    	memset(b, false, sizeof(bool)*initvn);

    	for(int i = 0; i < vnum; i++)
    	{
    		assert(hull.getHullVID()[i] <= initvn);
    		b[ hull.getHullVID()[i] ] = true;
    	}

    	int k = 0;
    	for(int i = 0; i < initvn; i++)
    	{
    		if(b[i])
    		{
    			fv2v.insert(make_pair(i, k));
    			v2fv.insert(make_pair(k, i));
    			k++;
    		}
    	}
    	delete[] b;
    }


    void saveHullOff(map<int,int>& fv2v, map<int,int>& v2fv,
    		vector<cd_v*>& vs, cd_m* m, cd_hull& hull, const string& offile, float offset=0.0f)
    {
        int nV = fv2v.size();
    	int nT = hull.getFaceSize();
    	
    	Point3d com;
		for(int k = 0; k < nV; k++)
		{
			int oid = v2fv[k];
			for(int d=0;d<3;d++) com[d]+=vs[oid]->pos[d];
		}
		for(int d=0;d<3;d++) com[d]/=nV; 


    	cout<<"writing off file directly"<<endl;
    	ofstream fout;
    	fout.open(offile.c_str());

		fout<<"OFF\n";
        fout<<nV<<" "<<nT<<" "<<nT*3/2<<"\n";
    	for(int k = 0; k < nV; k++)
    	{
    		int oid = v2fv[k];
    		Point3d offsetpt= vs[oid]->pos+(vs[oid]->pos-com).normalize()*offset;
    		fout<<offsetpt[0]<<" "<<offsetpt[1]<<" "<<offsetpt[2]<<"\n";
    	}
    	
    	for(int i = 0; i < nT; i++)
    	{
    		fout<<"3 ";
    		for(int j = 0; j < 3; j++)
    		{
    			fout<<" "<<fv2v[hull.vid(i,j)];// + 1;
    		}
    		fout<<endl;
    	}
    	fout.close();
    }
    
    void saveHullObj(map<int,int>& fv2v, map<int,int>& v2fv,
    		vector<cd_v*>& vs, cd_m* m, cd_hull& hull, string& objfile)
    {
    	cout<<"writing obj file "<<objfile.c_str()<<endl;
    	cout<<"init v size "<<vs.size()<<endl;
    	cout<<"hull v size "<<fv2v.size()<<endl;
    	cout<<"hull f size "<<hull.getFaceSize()<<endl;
    	    	int nV = fv2v.size();
    	    	int nT = hull.getFaceSize();
    	ofstream fout;
    	fout.open(objfile.c_str());
    	//save vertices
//    	for(map<int,int>::iterator mit = v2fv.begin(); mit != v2fv.end(); ++mit)
//    	{
//    		int oid = mit->second;
//    		fout<<"v "<<vs[oid]->pos[0]<<" "<<vs[oid]->pos[1]<<" "<<vs[oid]->pos[2]<<"\n";
//    	}
    	for(int k = 0; k < nV; k++)
    	{
    		int oid = v2fv[k];
    		fout<<"v "<<vs[oid]->pos[0]<<" "<<vs[oid]->pos[1]<<" "<<vs[oid]->pos[2]<<"\n";
    	}
    	//save faces
    	for(int i = 0; i < nT; i++)
    	{
    		fout<<"f";
    		for(int j = 0; j < 3; j++)
    		{
    			fout<<" "<<fv2v[hull.vid(i,j)] + 1;
    		}
    		if(i != nT - 1)
    			fout<<"\n";
    	}
    	fout.close();
    }

    void saveHull(map<int,int>& fv2v, map<int,int>& v2fv,
    		vector<cd_v*>& vs, cd_m* m, cd_hull& hull,ofstream& fout, const Material& material )
    {
    	int f = hull.getFaceSize();
    	int nV = fv2v.size();
    	int nT = hull.getFaceSize();

    	/*******************************************************/
    	//saveHull2Off(nV, nT, )

    	/*******************************************************/
    	//save
    	savePrepare(nV, nT, fout, material);

    	//save v
//        vector<cd_v*> vs;//(m->v().begin(), m->v().end());
        int mvsize = m->v().size();
        cout<<"the model's v size is "<<mvsize<<endl<<endl;


        fout <<"				coord DEF co Coordinate {" << std::endl;
        fout <<"					point [" << std::endl;
        for(map<int,int>::iterator mit = v2fv.begin(); mit != v2fv.end(); ++mit)
        {
        	int oid = mit->second;
        	fout <<"						" <<vs[oid]->pos[0]<< " "
        	            	     	 	 	  << vs[oid]->pos[1]<< " "
        	            	     	 	 	  << vs[oid]->pos[2] << "," << std::endl;
        }

        fout<<"					]"<<std::endl;
        fout<<"					}"<<std::endl;

        //save faces
        fout <<"				coordIndex [ " << std::endl;
        for(int i = 0; i < f; i++)
        {
        	int vid1 = hull.vid(i, 0);
        	int vid2 = hull.vid(i, 1);
        	int vid3 = hull.vid(i, 2);
            fout <<"						" << fv2v[vid1] << ", "
                	                          << fv2v[vid2]<< ", "
                	                          << fv2v[vid3] << ", -1," << std::endl;
        }

        fout<<"				]" << std::endl;
        fout <<"			}" << std::endl;
        fout <<"		}" << std::endl;
        fout <<"	]" << std::endl;
        fout <<"}" << std::endl;


        //
        //fv2v.clear();
        //v2fv.clear();

        fout.flush();
    }

void saveWRL(const string& wrl_filename, ML& ml)
{
    int max_fid = -1;
    int mid = 0;

	ofstream fout(wrl_filename.c_str());
	Material mat;

	FOREACH_M(ml)
    	    	//save
    	mat.m_diffuseColor[0] = mat.m_diffuseColor[1] = mat.m_diffuseColor[2] = 0.0;
		generateRandMaterial(mat);

    	saveOne(m, fout,mat );
	FOREACH_END
	fout.close();
}

void saveCHull(const string& offname, cd_m * m, float scale)
{
	vector<cd_v*> vs(m->v().begin(), m->v().end());

    //save a convex hull
    cd_hull hull;
    hull.buildhull(m);

    //loop over
    map<int, int> fv2v;
    map<int, int> v2fv;
    getHullInfo(fv2v, v2fv, hull, vs);//get hull information
    saveHullOff(fv2v, v2fv, vs, m, hull, offname, scale);
}

void saveCompHull(const string& filename, ML& ml, bool objformat)
{
	ExportOption& eo = ExportOption::getInstance();

	//vector<>
	vector<string> offiles;
	int mid = 0;
	FOREACH_M(ml)

		cd_m* hull_m = m;

    	//save a convex hull
    	cd_hull hull;
    	hull.buildhull(hull_m);

    	if(eo.b_simply_convex_hull)
    	{
    		int k = m->f().size() * eo.simply_convex_hull_ratio;
    		hull_m = hull.simplify(k<10?10:k, eo.simply_convex_hull_offset);
    		hull.buildhull(hull_m);
    	}

    	vector<cd_v*> vs(hull_m->v().begin(), hull_m->v().end());

    	//loop over
    	map<int, int> fv2v;
    	map<int, int> v2fv;
    	getHullInfo(fv2v, v2fv, hull, vs);//get hull information
//
//    	//save the convex hull to VRML
//    	saveHull(fv2v, v2fv, vs, m, hull, fout, mat);

//    	//save the convex hull to off
//    	/*************************************/
    	stringstream ss;
    	ss<<filename<<"_hull"<<"_"<<mid;
    	ss<<(objformat ? ".obj" : ".off");

    	string output_filename = ss.str();

    	//saveHull2Off(fv2v, v2fv, vs, m, hull, objname, offname);

    	if(objformat)
    		saveHullObj(fv2v, v2fv, vs, hull_m, hull, output_filename);
    	else
    		saveHullOff(fv2v, v2fv, vs, hull_m, hull, output_filename);

    	offiles.push_back(output_filename);

    	mid++;
	FOREACH_END

	//unionOff(offiles);
}

void saveCHullWRL(const string & wrl_filename, ML& ml)
{
	//save the convex hull of each model(cd_m)
	ofstream fout(wrl_filename.c_str());
	Material mat;

	//vector<>
	vector<string> offiles;
	int mid = 0;
	FOREACH_M(ml)
    	    	//save
    	mat.m_diffuseColor[0] = mat.m_diffuseColor[1] = mat.m_diffuseColor[2] = 0.0;
		generateRandMaterial(mat);

		vector<cd_v*> vs(m->v().begin(), m->v().end());

    	//save a convex hull
    	cd_hull hull;
    	hull.buildhull(m);

    	//loop over
    	map<int, int> fv2v;
    	map<int, int> v2fv;
    	getHullInfo(fv2v, v2fv, hull, vs);//get hull information

    	//save the convex hull to VRML
    	saveHull(fv2v, v2fv, vs, m, hull, fout, mat);

//    	//save the convex hull to off
//    	/*************************************/
    	stringstream ssobj;
    	ssobj<<wrl_filename.substr(0, wrl_filename.find_last_of("."))<<"_"<<mid<<".obj";
    	string objname = ssobj.str();
    	stringstream ssoff;
    	ssoff<<wrl_filename.substr(0, wrl_filename.find_last_of("."))<<"_"<<mid<<".off";
    	string offname = ssoff.str();

    	//saveHull2Off(fv2v, v2fv, vs, m, hull, objname, offname);
    	saveHullOff(fv2v, v2fv, vs, m, hull, offname);

    	offiles.push_back(offname);

    	mid++;
	FOREACH_END
	fout.close();

	//unionOff(offiles);
}
//
//void saveCHullOff(string& file, )














