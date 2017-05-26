#include <stdlib.h>
#define GLEW_STATIC


#include <GL/glew.h>
#include "GL/glui.h"
#include "GL/gliLight.h"
#include "GL/gliDump.h"

#include <sstream>
#include "GL/shadowMap_glsl.h"

#include "dude3d_draw.h"
#include "CD3d_stat.h"
#include "CD3d_draw_GL.h"
#include "io/CD3d_acdparameters.h"
#include "util/mathtool/Point.h"

using namespace mathtool;


extern cd_state state;            // displaying state
extern ML models;

//for random color
int g_color_seed=31;


//-----------------------------------------------------------------------------
void Idle( void )
{
  // According to the GLUT specification, the current window is
  //   undefined during an idle callback.  So we need to explicitly change
  //   it if necessary
  if ( glutGetWindow() != state.windowID )
    glutSetWindow(state.windowID);

  glutPostRedisplay();
}


void CreateGLUI()
{
    //
    // must call glutInit before this...
    //
    int Hsize=600;
    int Wsize=Hsize*1.5;
    int x_pos=20;
    int y_pos=50;

    glutInitDisplayMode( GLUT_RGBA|GLUT_DOUBLE|GLUT_DEPTH );
    glutInitWindowSize( Wsize, Hsize);
    glutInitWindowPosition( x_pos, y_pos );

    stringstream ss;
    string title = ss.str();
    state.windowID = glutCreateWindow((char*)title.c_str());

    InitGL();
    gli::gliInit();
    gli::gliDisplayFunc(Display);
    glutReshapeFunc(Reshape);
    gli::gliKeyboardFunc(Keyboard);
    //gli::gliMouseFunc(Mouse);
    GLUI_Master.set_glutIdleFunc(Idle);

    GLUI *glui = GLUI_Master.create_glui( title.c_str(), 0, x_pos+40+Wsize, y_pos );

#if 1
    // SETTING UP THE CONTROL PANEL:
    GLUI_Panel * top_panel=glui->add_panel("CShell Controls");

    glui->add_checkbox_to_panel(top_panel,"turn-on light", &state.b_lighting)->set_int_val(state.b_lighting);
    glui->add_checkbox_to_panel(top_panel,"show solid (S)", &state.show_solid)->set_int_val(state.show_solid);
    glui->add_checkbox_to_panel(top_panel,"show convex hull (h)", &state.show_b_hull)->set_int_val(state.show_b_hull);
    glui->add_checkbox_to_panel(top_panel,"smooth (s)", &state.show_smooth)->set_int_val(state.show_smooth);
    glui->add_checkbox_to_panel(top_panel,"grey (no color)", &state.b_grey)->set_int_val(state.b_grey);

    glui->add_checkbox_to_panel(top_panel,"show face number (N)", &state.show_face_number)->set_int_val(state.show_face_number);
    glui->add_checkbox_to_panel(top_panel,"show vertex number (V)", &state.show_vertex_number)->set_int_val(state.show_vertex_number);
    glui->add_checkbox_to_panel(top_panel,"show bounding box", &state.b_draw_bbox)->set_int_val(state.b_draw_bbox);
    glui->add_checkbox_to_panel(top_panel,"show text", &state.show_txt)->set_int_val(state.show_txt);
    glui->add_checkbox_to_panel(top_panel,"rotate around X", &state.b_rotate_x)->set_int_val(state.b_rotate_x);

    glui->add_edittext_to_panel(top_panel, "select vertex:", GLUI_EDITTEXT_INT, &state.selectedVID)->set_float_val(state.selectedVID);
    glui->add_edittext_to_panel(top_panel, "select face:", GLUI_EDITTEXT_INT, &state.selectedFID)->set_float_val(state.selectedFID);

//    glui->add_button_to_panel(top_panel, "Move Model Up", -1, (GLUI_Update_CB)moveModelUp);
//    glui->add_button_to_panel(top_panel, "Move Model Down", -1, (GLUI_Update_CB)moveModelDown);
//    glui->add_button_to_panel(top_panel, "Move Model Front", -1, (GLUI_Update_CB)moveModelFront);
//    glui->add_button_to_panel(top_panel, "Move Model Back", -1, (GLUI_Update_CB)moveModelBack);

    glui->add_button_to_panel(top_panel, "Random Color", -1, (GLUI_Update_CB)randColor);
    glui->add_button_to_panel(top_panel, "Reset Camera", -1, (GLUI_Update_CB)resetCamera);

#endif

#if 1 //GUIs for convex hull operations
	GLUI_Panel * hull_panel = glui->add_panel("Hulls");
	glui->add_button_to_panel(hull_panel, "Rebuild Hulls", -1, (GLUI_Update_CB)rebuild_glmodel_hulls);
    glui->add_button_to_panel(hull_panel, "Print Fatness", -1, (GLUI_Update_CB)print_fatness);

    GLUI_Panel * simp_hull_panel = glui->add_panel_to_panel(hull_panel,"Simplify&Remesh Hulls");
    //glui->add_checkbox_to_panel(simp_hull_panel, "show fatness", &state.b_show_hull_fatness)->set_int_val(state.b_show_hull_fatness);
    glui->add_button_to_panel(simp_hull_panel, "Simplify Hulls", -1, (GLUI_Update_CB)simplify_hull);
    glui->add_button_to_panel(simp_hull_panel, "Remesh Hulls", -1, (GLUI_Update_CB)remesh_hull);
    glui->add_edittext_to_panel(simp_hull_panel, "iteration:", GLUI_EDITTEXT_INT, &state.hull_simplification_iteration)->set_float_val(state.hull_simplification_iteration);
    glui->add_edittext_to_panel(simp_hull_panel, "max vol. incease %:", GLUI_EDITTEXT_FLOAT, &state.hull_simplification_max_vol_inc_ratio)->set_float_val(state.hull_simplification_max_vol_inc_ratio);
	//glui->add_edittext_to_panel(simp_hull_panel, "min face count:", GLUI_EDITTEXT_FLOAT, &state.hull_simplification_min_face_count)->set_float_val(state.hull_simplification_min_face_count);

    //common trim options
    GLUI_Panel * trim_hull_panel = glui->add_panel_to_panel(hull_panel,"Trim Hulls");
	glui->add_edittext_to_panel(trim_hull_panel, "collapse vertex dist:", GLUI_EDITTEXT_FLOAT, &state.hull_trim_vertex_collapse_dist)->set_float_val(state.hull_trim_vertex_collapse_dist);
	glui->add_edittext_to_panel(trim_hull_panel, "volume sample size:", GLUI_EDITTEXT_FLOAT, &state.hull_trim_vol_sample_size)->set_float_val(state.hull_trim_vol_sample_size);
	//glui->add_edittext_to_panel(trim_hull_panel, "surface sample density:", GLUI_EDITTEXT_FLOAT, &state.hull_trim_surface_sample_coverage)->set_float_val(state.hull_trim_surface_sample_coverage);
glui->add_button_to_panel(trim_hull_panel, "Use Heuristic", -1, (GLUI_Update_CB)trim_hull_heuristic);
    //svm
    glui->add_button_to_panel(trim_hull_panel, "Use Exact Volume", -1, (GLUI_Update_CB)trim_hull_exact);
    glui->add_button_to_panel(trim_hull_panel, "Use SVM", -1, (GLUI_Update_CB)trim_hull_svm);
    glui->add_edittext_to_panel(trim_hull_panel, "c-svm C:", GLUI_EDITTEXT_FLOAT, &state.hull_trim_svm_C)->set_float_val(state.hull_trim_svm_C);
    //glui->add_button_to_panel(trim_hull_panel, "Use Duality", -1, (GLUI_Update_CB)trim_hull_duality);
	glui->add_button_to_panel(trim_hull_panel, "Save Hulls", -1, (GLUI_Update_CB)save_gmodel_hulls);

    glui->add_checkbox_to_panel(hull_panel, "show samples", &state.b_show_hull_samples)->set_int_val(state.b_show_hull_samples);
	glui->add_checkbox_to_panel(hull_panel, "show support vectors", &state.b_show_hull_support_vectors)->set_int_val(state.b_show_hull_support_vectors);
	glui->add_checkbox_to_panel(hull_panel, "show separating planes", &state.b_show_hull_separators)->set_int_val(state.b_show_hull_separators);
#endif
}

bool InitGL()
{
    srand48(g_color_seed);

	//
    // disable cast shadow if openGL 2.0 is not available
	// JML: glewInit() must be called before enable some properties below...
	//
    glewInit();
    GLboolean GL2=glewIsSupported("GL_VERSION_2_0");

    if(GL2==false) cout<<"! Warning: GL2 is Not supported"<<endl;

    /////////////////////////////////////////////////////////////////////
    //update view
    state.reset();
    drawReset();
    resetCamera();

    // Antialias and others
	glEnable( GL_LINE_SMOOTH );
	glEnable( GL_BLEND );
	glShadeModel(GL_SMOOTH);
	glEnable(GL_DEPTH_TEST);
	glClearColor(1,1,1,1.0f);
	glEnable(GL_CULL_FACE);
	glHint(GL_PERSPECTIVE_CORRECTION_HINT,GL_NICEST);
	glEnable(GL_AUTO_NORMAL);
	glEnable(GL_NORMALIZE);
    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
	glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);

    //let's have light
    setupLight();

    //create shadow map
    if(GL2)
    {
        int vp[4];
        glGetIntegerv(GL_VIEWPORT, vp);
        ShadowMap_GLSL::generateShadowFBO(vp[2], vp[3], light0_position[0], light0_position[1], light0_position[2]);
        ShadowMap_GLSL::loadShadowShader();
    }
    //

    return true;
}

void Display( void )
{
    //Init Draw
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glMatrixMode(GL_MODELVIEW);

	glLoadIdentity();

    ///////////////////////////////////////////////////////////////////////////

    if(state.b_lighting)
    {
    	glEnable(GL_LIGHTING);
    	//glPushMatrix();
    	////make the lights rotate with user control (this fix the relative transform with the env)
    	//glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
    	//glLightfv(GL_LIGHT1, GL_POSITION, light1_position);
    	//glLightfv(GL_LIGHT2, GL_POSITION, light2_position);
    	//glLightfv(GL_LIGHT3, GL_POSITION, light3_position);
    	//glPopMatrix();
    }
    else
    {
    	glDisable(GL_LIGHTING);
    }

    glPushMatrix();
	//move everything to the original of the world
    //glTranslated(-state.model_com[0],-state.model_com[1],-state.model_com[2]);
    drawAll();
    glPopMatrix();
}

///////////////////////////////////////////////////////////////////////////////
void drawAll()
{
    //draw model
    draw(models);
}

void Reshape( int w, int h)
{
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();

    gluPerspective(60, w*1.0f/h, 0.1, state.model_radus*100 );

    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();

    // regenerate shadow map
    ShadowMap_GLSL::clean();
    ShadowMap_GLSL::generateShadowFBO(w, h, light0_position[0], light0_position[1], light0_position[2]);
    ShadowMap_GLSL::loadShadowShader();
}



void Keyboard( unsigned char key, int x, int y )
{
    // find closest colorPt3D if ctrl is pressed...
    switch( key ){
        case 'w': state.show_wire=!state.show_wire; break;
        case 'S': state.show_solid=!state.show_solid; break;
        case 's': state.show_smooth=!state.show_smooth; break;
        case 'N': state.show_face_number=!state.show_face_number; break;
        case 'V': state.show_vertex_number=!state.show_vertex_number; break;
        case 'E': state.show_apart=!state.show_apart; break; //explode
        case 'C': randColor(); break;

        case 'h': state.show_b_hull=!state.show_b_hull; break;
		case 'r': resetCamera(); break;

		//
        case 'L' : state.b_lighting=!state.b_lighting; break;
    }

    GLUI_Master.sync_live_all();
}

void resetCamera()
{
    //reset camera
    gli::setScale(1.25);
    gli::setCameraPosZ(2*state.model_radus);
    gli::setCameraPosX(0);
    gli::setCameraPosY(0);
    gli::setAzim(0);
    gli::setElev(0);

	glutPostRedisplay();
}


void setupLight()
{
    //Let's have light!
    light0_position[0]=3 * state.model_radus;
    light0_position[1]=state.bbox[3]+10*state.model_radus;
    light0_position[2]=3 * state.model_radus;

    light1_position[0]=-2*state.model_radus;
    light1_position[1]=state.bbox[3] + state.model_radus;//*2;
    light1_position[2]=-2*state.model_radus;//*2;

    light2_position[0] = -state.model_radus*4;;//-state.model_radus*2;
    light2_position[1] = state.bbox[3]+state.model_radus*4;
    light2_position[2] = state.model_radus *4;

    light3_position[0] = 4*state.model_radus;
    light3_position[1] = -state.model_radus *4;
    light3_position[2] = state.model_radus *4;


    glLightfv(GL_LIGHT0, GL_DIFFUSE,  light0_diffuse);
    glLightfv(GL_LIGHT0, GL_AMBIENT,  light0_ambient);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light0_specular);
    glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
    glEnable(GL_LIGHT0);

    glLightfv(GL_LIGHT1, GL_DIFFUSE,  light1_diffuse);
    glLightfv(GL_LIGHT1, GL_AMBIENT,  light1_ambient);
    glLightfv(GL_LIGHT1, GL_SPECULAR, light1_specular);
    glLightfv(GL_LIGHT1, GL_POSITION, light1_position);
    glEnable(GL_LIGHT1);

    glLightfv(GL_LIGHT2, GL_DIFFUSE, light2_diffuse);
    glLightfv(GL_LIGHT2, GL_AMBIENT,  light2_ambient);
    glLightfv(GL_LIGHT2, GL_SPECULAR, light2_specular);
    glLightfv(GL_LIGHT2, GL_POSITION, light2_position);
    glEnable(GL_LIGHT2);

	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT, GL_SHININESS, mat_low_shininess);
}

void simplify_hull()
{
    if(state.hull_simplification_iteration<=0) state.hull_simplification_iteration=1;
    for(int i=0;i<state.hull_simplification_iteration;i++) simplify_glmodel_hulls();
}

void remesh_hull()
{
    cout << "------------ remesh_glmodel_hulls ------------" << endl;

    if(state.hull_simplification_iteration<=0) state.hull_simplification_iteration=1;
    for(int i=0;i<state.hull_simplification_iteration;i++) remesh_glmodel_hulls(i);

    cout << "\n------------ done ------------" << endl;
}

void trim_hull_svm()
{
	trim_glmodel_hulls(SVM_Trim);
}

void trim_hull_heuristic()
{
  trim_glmodel_hulls(Heuristic_Trim);
}

void trim_hull_duality()
{
  trim_glmodel_hulls(Dualspace_Trim);
}

void trim_hull_exact()
{
  trim_glmodel_hulls(Exact_Trim);
}

void moveModelUp()
{
    state.usrModelTranslate[1]+=state.model_radus/50;
}

void moveModelDown()
{
    state.usrModelTranslate[1]-=state.model_radus/50;
}

void moveModelFront()
{
	state.usrModelTranslate[2]+=state.model_radus/50;
}

void moveModelBack()
{
	state.usrModelTranslate[2]-=state.model_radus/50;
}
