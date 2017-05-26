#ifndef _CD3d_STATE_H_
#define _CD3d_STATE_H_

#pragma once

#include "CD3d_model.h"

class cd_state
{
public:

    cd_state();
    void reset();

    int windowID; //glut window id

	int b_Verbose;
	int b_showGL;

    int show_hull;       //show convex hull of the entire input model
	int show_b_hull;     //show bounding hull of all components

    int show_wire;              //show wire frame
    int show_solid;             //solid the model
    int show_smooth;            //smooth solid model
    int show_face_number;       //show face id on models/hull
    int show_vertex_number;     //show vertex id on models/hull
    int show_apart;             //seperate decomposed models
    int show_axis;
    int show_txt;

    int b_grey;                 //make model grey
    int b_rotate_x;		        //rotate 90 degree around x axis

	//hull related flags
    //int b_show_hull_fatness;
	int b_show_hull_samples;
	int b_show_hull_support_vectors;
	int b_show_hull_separators;
	int rebuild;

	//data
	double  model_radus;  // radius of the model
	Point3d model_com;    // center of mass of the model

	float opacity;

	//enclosing box
	int b_draw_bbox;
	float bbox[6];

	//light/shadow
	int b_lighting;
	int b_cast_shadow;

	int selectedVID;
	int selectedFID;

	cd_v * selectedV;
	cd_f * selectedF;

    list<string> str_input;          //input filenames
    string str_output;         //output filename

    Vector3d usrModelTranslate;

	///hull operation parameters
	float hull_simplification_max_vol_inc_ratio;
	float hull_simplification_face_reduction_ratio;
	float hull_simplification_min_face_count;
  int hull_simplification_iteration;

  //JML we will use the same values above
  //  float hull_remeshing_max_vol_inc_ratio;
  //  float hull_remeshing_face_reduction_ratio;
  //  float hull_remeshing_min_face_count;
  float hull_trim_vol_sample_size;

	float hull_trim_surface_sample_coverage;
	float hull_trim_svm_C;
	float hull_trim_vertex_collapse_dist;

  //
  string trim_method;

	float scale;
};

#endif //_CD3d_STATE_H_
