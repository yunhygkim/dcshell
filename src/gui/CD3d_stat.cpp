#include "CD3d_stat.h"
#include <cfloat>

#ifdef WIN32
#pragma warning(disable : 4786)
#endif

cd_state::cd_state()
{
	b_Verbose=false;
	b_showGL=true;
	b_cast_shadow=false;
	b_lighting=true;
	b_grey=false;
	b_rotate_x=false;
	windowID=-1;

    show_hull=false;
    show_b_hull=false;
	b_show_hull_samples = false;
	b_show_hull_support_vectors = false;
	b_show_hull_separators = false;

    show_wire=false;
    show_solid=true;
    show_face_number=false;
    show_vertex_number=false;
    show_apart=false;
    show_smooth=false;
    show_axis=true;
    show_txt=true;
    b_draw_bbox=false;


    model_radus=-1;  // radius of the model

	//for making movies
	rebuild=true;

    opacity=1; // opacity of components
    scale=1;   // scale of components

	hull_simplification_max_vol_inc_ratio = 0.01;
	hull_simplification_min_face_count = 10;
	hull_simplification_face_reduction_ratio = 0.1;
	hull_simplification_iteration=1;

	hull_trim_vol_sample_size = 10000;
	hull_trim_surface_sample_coverage = 5;
	hull_trim_svm_C = 100;
	hull_trim_vertex_collapse_dist = 0.01;

    //
    bbox[0]=bbox[2]=bbox[4]= FLT_MAX;
    bbox[1]=bbox[3]=bbox[5]=-FLT_MAX;

    selectedVID=-1;
    selectedFID=-1;

    selectedV=NULL;
    selectedF=NULL;
}

void cd_state::reset()
{
//    selectedV.clear();
//    selectedF.clear();
//    selectedP.clear();
//    selectedM.clear();
}
