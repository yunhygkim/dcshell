/*
 * CD3dexport.h
 *
 *  Created on: Oct 16, 2013
 *      Author: guilin
 */

#ifndef CD3DEXPORT_H_
#define CD3DEXPORT_H_

#include <fstream>
#include <iostream>

#include "CD3d_model.h"
#include "util/mathtool/Vector.h"
//#include "CD3d_polyUnion.h"

using namespace std;
using namespace mathtool;

    class Material
    {
    public:
                                                                Material(void);
                                                                ~Material(void){}
//    private:
        Vector3d                                            	m_diffuseColor;
        double                                                  m_ambientIntensity;
        Vector3d                                            	m_specularColor;
        Vector3d                                            	m_emissiveColor;
        double                                                  m_shininess;
        double                                                  m_transparency;

        //friend class TMMesh;
		//friend class HACD;
    };

    class ExportOption
    {
    public:
    	static ExportOption& getInstance()
    	{
    		static ExportOption instance;

    		return instance;
    	}

    	int b_simply_convex_hull; 		// export simplified convex hulls
		float simply_convex_hull_ratio;	// % facets remained in the simplified convex hulls [default = 0.02]
		float simply_convex_hull_offset;
    private:
		ExportOption();
    };

    //save one model
    void saveOne(cd_m* m,ofstream& fout, const Material& material);

	void saveObj(cd_m* m, ofstream& fout);
	
	void saveWRL(const string& wrl_filename, ML& ml);

    void saveCompHull(const string& offname, ML& ml, bool objformat = false);
    
	void saveCHullWRL(const string & wrl_filename, ML& ml);

    //scale scales the chull of m
    void saveCHull(const string& offname, cd_m * m, float offset=0.0f);

#endif /* CD3DEXPORT_H_ */
