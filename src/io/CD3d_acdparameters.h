#ifndef _CD3D_ACDPARAMS_H_
#define _CD3D_ACDPARAMS_H_

#include <string>
#include <list>
using namespace std;


struct cd_params
{

    cd_params()
    {
        tau=cp_min=cp_cur=cut_cur=convx=0;
		use_convx=false;
        verbosity = false;
        base_tau = 0;

		hull_simplification_max_vol_inc_ratio  = 0.1;
		hull_simplification_min_face_count = 100;
		hull_simplification_face_reduction_ratio=0.1;
        hull_trim_vol_sample_size=10000;
		hull_trim_surface_sample_coverage=0.5;
		hull_trim_svm_C=1e4;
		hull_trim_vertex_collapse_dist = 0.001;
    }

    list<string> in_filename;
    string out_filename;
    float tau, cp_min, cp_cur, cut_cur, convx;
	bool use_convx;

    //display option
    bool verbosity;

    ///base tolerance
    float base_tau;

	///hull operation parameters
	float hull_simplification_max_vol_inc_ratio;
	float hull_simplification_face_reduction_ratio;
	float hull_simplification_min_face_count;
    float hull_trim_vol_sample_size;
	float hull_trim_surface_sample_coverage;
	float hull_trim_svm_C;
	float hull_trim_vertex_collapse_dist;

};

#endif //_CD3D_ACDPARAMS_H_
