#include "cshell.h"
#include "io/CD3dexport.h"
#include "util/CD3d_convexity.h"

///////////////////////////////////////////////////////////////////////////////
//States
cd_state state;            // displaying state
ML models;                 // all models
//////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
#include <vector>
#ifndef _CORISELIB_
vector<string> filenames;
#endif
#ifndef _CORISELIB_
///////////////////////////////////////////////////////////////////////////////

#include <CD3d_hull.h>
#include <hull_vol.h>
#include "CD3d_draw_GL.h"

void random_hull_vol_test();
void random_hull_vol_test(cd_m * m);

int main( int argc, char ** argv )
{
    //initialized
    if( parseCmd(argc,argv)==false )
    {
      cout<<"----------- run random hull vol test -----------"<<endl;
      random_hull_vol_test();
      return 1;
    }

 	  if (state.str_input.empty())
    {
        cerr<<"! Error: No input file"<<endl;
        return 1;
    }

    //create models
  	createModels(state.str_input, models);
	  int msize = models.size();

    if(msize==1)
    {
      cout<<"----------- only one model is loaded: compute vol of the model -----------"<<endl;
      random_hull_vol_test(models.front());
      return 0; //no need to do anything
    }

    /////////////////////////////////////////////////////////////////
    //setup glut/gli
    if(state.b_showGL)
    {
		    /////////////////////////////////////////////////////////////
		    glutInit( &argc, argv );
		    CreateGLUI(); //see dude3d_draw.h/cpp

        /////////////////////////////////////////////////////////////
        gli::gliMainLoop();
    }
    else
    {
        //
        glutInit( &argc, argv );
        CreateGLUI(); //see dude3d_draw.h/cpp
        buildglmodels(models);

        for(int i=0;i<state.hull_simplification_iteration;i++){
          simplify_glmodel_hulls();
          remesh_glmodel_hulls(i);
        }

        if(state.trim_method=="svm")
          trim_glmodel_hulls(SVM_Trim);
        else if(state.trim_method=="exact")
          trim_glmodel_hulls(Exact_Trim);
        else trim_glmodel_hulls(Heuristic_Trim);
    }


    return 0;
}

#endif



///////////////////////////////////////////////////////////////////////////////
bool parseCmd(int argc,char ** argv)
{
    if( argc==1 ){
        printUsage(argv[0]);
        return false;
    }

    for(int i=1;i<argc;i++){
        if( strlen(argv[i])<2 ){ printUsage(argv[0]); return false; }
        if( argv[i][0]=='-' ){//an option
            switch(argv[i][1]){
              case 'v': state.b_Verbose=true;
                    break;
      				case 'o':
      					state.str_output=argv[++i];
      					break;
      				case 'g': state.b_showGL=false;
      				    break;
              case 'm': //method
                  state.trim_method = argv[++i];
                  break;
              case 'C': //penalty
                  state.hull_trim_svm_C = atof(argv[++i]);
                  break;
              case 'i': //iterations
                  state.hull_simplification_iteration = atoi(argv[++i]);
                  break;
              case 'e': //max volume expansion allowed
                      state.hull_simplification_max_vol_inc_ratio = atof(argv[++i]);
                      break;
              default:
                  printUsage(argv[0]);
                  return false;
            }//end switch
        }
        else{// not an option
			      state.str_input.push_back(argv[i]);
            filenames.push_back(argv[i]);
        }
    }//end for
    return true;

}



///////////////////////////////////////////////////////////////////////////////

void createModels(const list<string>& filenames, ML& ml)
{
  for(auto i=filenames.begin();i!=filenames.end();i++)
	{
      ///////////////////////////////////////////////////////////////////////
      if(state.b_Verbose) cout<<"- Read File, "<<*i<<endl;
      objReader loader;
      loader.SetDataFileName(i->c_str());
      loader.ParseFile();
      ///////////////////////////////////////////////////////////////////////
      if(state.b_Verbose) cout<<"- Create Model "<<"Vertices = "<<loader.GetVertices().size()<<" Triangles = "<<loader.GetTriangles().size()<< endl;
      buildModels(loader.GetParts(),loader.GetVertices(),loader.GetTriangles(),ml);
  }//for end

  if(state.b_Verbose) cout<<"- Found "<<ml.size()<<" model(s)"<<endl;

	if (ml.empty())
	{
		cerr << "! Error: Can't find any models" << endl;
		exit(1);
	}

	//compute the com and radius of the models
	model_COM_R(ml);
}

///////////////////////////////////////////////////////////////////////////////
void printUsage(char * arg0)
{
    cerr<<"Usage: "<<arg0<<" [option] filename \n"
        <<"Options: \n"
        <<"\t-v \t\t Verbose mode. Output statisitic data. \n"
		<<"\t-o \t\t Output file.\n"
        <<"Report bugs to jmlien@cs.gmu.edu\n"
        <<endl;
}

void model_COM_R(ML& models)
{
	VL allv; //all points in the model...
    FOREACH_CM(models)
        VL& vl=m->v();
        FOREACH_V(vl) allv.push_back(v); FOREACH_END
    FOREACH_END
    state.model_com=computeCenter(allv);
    state.model_radus=computeRadius(allv,state.model_com);
    computeBoundingBox(allv,state.bbox);
    //cd_out(cd_msg("Compute COM=")+state.model_com+"/Radius="+state.model_radus);
}


//----------------------------------------------------------------
void random_hull_vol_test()
{
  srand48(time(NULL));
  cout.precision(20);
  double max_diff=0;

  for(int d=0;d<100;d++)
  {
    cout<<"----------"<<endl;
    cd_hull tmp;
    vector<Point3d> pts;
    for(uint i=0;i<300;i++)
    {
      Point3d pt=Point3d(0,0,0)+Vector3d(drand48(),drand48(),drand48())*100;
      if(drand48()<0.5) pt[0]=-pt[0];
      if(drand48()<0.5) pt[1]=-pt[1];
      if(drand48()<0.5) pt[2]=-pt[2];
      pts.push_back(pt);
    }

    tmp.buildhull(pts);

    hull_plane hp;
    hp.o.set(0,0,0);
    hp.n.set(0,-1,0);
    //hp.n = hp.n+Vector3d(drand48()*1e-5,drand48()*1e-5,drand48()*1e-5);
    hp.n=hp.n.normalize();

    double a=compute_volume_explicit(tmp,hp);
    double b=compute_volume(tmp,hp);
    cout<<"original vol="<<tmp.computeHullVol()<<endl;
    cout<<"compute_volume_explicit="<<a<<endl;
    cout<<"compute_volume="<<b<<endl;
    cout<<"diff="<<fabs(a-b)<<endl;
    if(fabs(a-b)>max_diff)
    {
      max_diff=fabs(a-b);
    }
  }

  cout<<"max_diff="<<max_diff<<endl;
}


void random_hull_vol_test(cd_m * m)
{
  cd_hull hull;
  hull.buildhull(m);

  hull_plane hp;
  //hp.o=-0.156393 -0.463652 2.38332  hp.n=0.0642792 0.190566 -0.979568
  hp.o.set(-0.156393, -0.463652, 3.38332);
  hp.n.set(0.0642792, 0.190566, -0.979568 );

  double a=compute_volume_explicit(hull,hp);
  double b=compute_volume(hull,hp, true);

  cout<<"original vol="<<hull.computeHullVol()<<endl;
  cout<<"compute_volume_explicit="<<a<<endl;
  cout<<"compute_volume="<<b<<endl;
  cout<<"diff="<<fabs(a-b)<<endl;
}
