/**
  * Interaction class for opengl.
  * OpenGL Interactor. (gliInteractor)
  *
  * 2001/3/19 Jyh-Ming
  */

#ifdef _WIN32
#pragma warning(disable : 4786 4244 4305)
#endif

#include "gli.h"
#include <iostream>
#include <math.h>
#include <float.h>

using namespace std;

GLfloat gli::m_CameraPos[3] = {0,0,100};
GLfloat gli::m_deltaDis[3]  = {0,0,0};

bool gli::m_DisalbeMouseControl=false;

GLfloat gli::m_currentAzim = 0;
GLfloat gli::m_deltaAzim   = 0;
GLfloat gli::m_currentElev = 0; 
GLfloat gli::m_deltaElev   = 0;

int gli::m_StartX        =0; 
int gli::m_StartY        =0; 
int gli::m_PressedButton =0;

gli::GLI_DISPLAY_FUNC   gli::m_dislpayFunc  = NULL;
gli::GLI_MOUSE_FUNC     gli::m_mouseFunc    = NULL;
gli::GLI_MOTION_FUNC    gli::m_motionFunc   = NULL;
gli::GLI_SPECIAL_FUNC   gli::m_specialFunc  = NULL;
gli::GLI_KEYBOARD_FUNC  gli::m_keyboardFunc  = NULL; 

float gli::m_WindowX[3]={1,0,0};
float gli::m_WindowY[3]={0,1,0};

float gli::m_Scale=1;
float gli::m_deltaScale=1;

bool gli::m_2D=false;

void gli::gliInit()
{
    glutDisplayFunc( gliDisplay );
    glutMouseFunc( gliMouse );
    glutMotionFunc( gliMotion );
    glutSpecialFunc( gliSpecialFunc );
    glutKeyboardFunc( gliKeyboardFunc );
    glutSpecialFunc( gliSpecialFunc );
}

void gli::gliRotate()
{
    glRotated(m_currentElev+m_deltaElev, 1.0, 0.0, 0.0);
    glRotated(m_currentAzim+m_deltaAzim, 0.0, 1.0, 0.0);

//
//    cout<<"rot= "<<(m_currentElev+m_deltaElev)<<","<<
//            (m_currentAzim+m_deltaAzim)<<endl;
}

void gli::gliTranslate()
{

    //cout<<"pos= "<<(m_CameraPos[0]+m_deltaDis[0])<<","<<
    //        (m_CameraPos[1]+m_deltaDis[1])<<","<<
    //        (m_CameraPos[2]+m_deltaDis[2])<<endl;

    glTranslatef( (m_CameraPos[0]+m_deltaDis[0]),
                 -(m_CameraPos[1]+m_deltaDis[1]),
                 -(m_CameraPos[2]+m_deltaDis[2]) );
}

void  gli::ApplyCameraTrasnform()
{
    glMatrixMode( GL_MODELVIEW );
    if( !m_DisalbeMouseControl ){ //mouse control is enabled
        gliTranslate();
        gliRotate();		
		float s=getScale();
		glScalef(s,s,s);
    }
}

void gli::drawAxes( void )
{
    //draw reference axis
    glMatrixMode(GL_PROJECTION); //change to Ortho view
    glPushMatrix(); 
    glLoadIdentity();
    gluOrtho2D(0,20,0,20);

    glMatrixMode( GL_MODELVIEW );
    glPushMatrix();
    glLoadIdentity();

	glDisable(GL_LIGHTING);
	glDisable(GL_DEPTH_TEST);
    glTranslatef(1.1f,1,0);
    glRotated(m_currentElev+m_deltaElev, 1.0, 0.0, 0.0);
    glRotated(m_currentAzim+m_deltaAzim, 0.0, 1.0, 0.0);
    glBegin(GL_LINES);
    glColor3f(1,0,0);
    glVertex3f(0,0,0);
    glVertex3f(1,0,0);
    glColor3f(0,1,0);
    glVertex3f(0,0,0);
    glVertex3f(0,1,0);
    glColor3f(0,0,1);
    glVertex3f(0,0,0);
    glVertex3f(0,0,1);
    glEnd();
    glPopMatrix();
	glEnable(GL_DEPTH_TEST);
    
    //pop GL_PROJECTION
    glMatrixMode(GL_PROJECTION); //change to Pers view
    glPopMatrix();
    glMatrixMode( GL_MODELVIEW );
}

void gli::gliDisplay( void )
{
    //Draw scene
    glLoadIdentity();
    ApplyCameraTrasnform();

    if( m_dislpayFunc!=NULL ){
        m_dislpayFunc();
    }
    
    drawAxes();
    glutSwapBuffers();
}

void gli::
gliMouse( int button, int state, int x, int y )
{

    if( !m_DisalbeMouseControl ){
        //if( glutGetModifiers()!=GLUT_ACTIVE_CTRL )
        {
            if( state == GLUT_UP ) //reset every thing
            {
                for( int iD=0;iD<3;iD++ ){
                    m_CameraPos[iD] += m_deltaDis[iD];
                    m_deltaDis[iD]=0;
                }


                if(m_2D){
                    m_Scale = m_Scale*m_deltaScale;
                    m_deltaScale= 1.0f;
                }
                else{
                    m_currentElev += m_deltaElev;
                    m_currentAzim += m_deltaAzim;

                    m_deltaElev = 0.0;
                    m_deltaAzim = 0.0;
                }
            }
            else { //press down

                m_StartY = y;
                m_StartX = x;
                m_StartY = y;

                m_PressedButton = button;
            }
        }
        //else m_PressedButton=-1;
    }

    //call it when alt is not pressed
    if( m_mouseFunc!=NULL )
        m_mouseFunc(button, state, x, y);

    //redraw any way
    glutPostRedisplay();
}

void gli::
gliMotion( int x, int y )
{

    static float initZ=-FLT_MAX;
    if(initZ==-FLT_MAX) initZ=m_CameraPos[2];

    if( !m_DisalbeMouseControl ){

        if( m_PressedButton == GLUT_RIGHT_BUTTON )
        {
            if(m_2D){
                double d = (y-m_StartY)*1.0/100;
                if(d>0) m_deltaScale = 1.0/(1+fabs(d));
                if(d<0) m_deltaScale = (1+fabs(d));
            }
            else{
                m_deltaDis[2] = ((m_CameraPos[2]>10)?m_CameraPos[2]:10) * ((GLfloat)(y - m_StartY))/100.0;
            }
        }
        else if( m_PressedButton == GLUT_MIDDLE_BUTTON )
        {
            if(!m_2D)
            {
                float zoom_diff;
                if(m_CameraPos[2]>initZ){ //zoom out
                    zoom_diff=(m_CameraPos[2]-initZ)/initZ+1;
                }
                else{ //zoom in
                    zoom_diff=fabs(initZ)/(initZ-m_CameraPos[2]+fabs(initZ));
                }

                m_deltaDis[0] = ((GLfloat)(x - m_StartX))*40*zoom_diff; //((m_CameraPos[0]>5)?m_CameraPos[0]:5) * ((GLfloat)(x - m_StartX))/20.0;
                m_deltaDis[1] = ((GLfloat)(y - m_StartY))*40*zoom_diff; //((m_CameraPos[1]>5)?m_CameraPos[1]:5) * ((GLfloat)(y - m_StartY))/20.0;
            }
        }
        else if(m_PressedButton == GLUT_LEFT_BUTTON)
        {
            if(m_2D){
            	float z=(m_CameraPos[2]/600)/m_Scale;
                m_deltaDis[0] = ((GLfloat)(x - m_StartX))*z;
                m_deltaDis[1] = ((GLfloat)(y - m_StartY))*z;
            }
            else{
                m_deltaAzim = ((GLfloat)(x - m_StartX))/5.0;
                m_deltaElev = ((GLfloat)(y - m_StartY))/5.0;
                //compute window x, y dir
                m_WindowX[0]=1; m_WindowX[1]=0; m_WindowX[2]=0;
                m_WindowY[0]=0; m_WindowY[1]=1; m_WindowY[2]=0;
                gliRotateX(m_WindowX, m_currentElev+m_deltaElev);
                gliRotateY(m_WindowX, m_currentAzim+m_deltaAzim);
                m_WindowX[2]=-m_WindowX[2];
                gliRotateX(m_WindowY, m_currentElev+m_deltaElev);
                gliRotateY(m_WindowY, m_currentAzim+m_deltaAzim);
                m_WindowY[2]=-m_WindowY[2];
            }
        }
    }

    if( m_motionFunc!=NULL )
        m_motionFunc(x, y);

    //redraw any way
    glutPostRedisplay();
}

void gli::gliRotateX(float v[3], float degree){
    float c=cos(3.1415926*degree/180);
    float s=sin(3.1415926*degree/180);
    float v1=v[1]*c-v[2]*s;
    float v2=v[1]*s+v[2]*c;
    v[1]=v1; v[2]=v2;
}

void gli::gliRotateY(float v[3], float degree){
    float c=cos(3.1415926*degree/180);
    float s=sin(3.1415926*degree/180);
    float v0=v[0]*c+v[2]*s;
    float v2=-v[0]*s+v[2]*c;
    v[0]=v0; v[2]=v2;
}

void gli::gliKeyboardFunc(unsigned char key, int x, int y)
{
    switch(key)
    {
		case 27: exit(0);
	    case '5' : gli::setScale(gli::getScale()*0.95);
				break;
	    case '6' : gli::setScale(gli::getScale()*1.05);
				break;
	}

    if( m_keyboardFunc!=NULL)
        m_keyboardFunc(key,x,y);

    //redraw any way
    glutPostRedisplay();
}

void gli::gliSpecialFunc(int key, int x, int y)
{
	float z=m_CameraPos[2]/100;
	
    switch(key){
        case GLUT_KEY_LEFT: m_CameraPos[0]-=z/m_Scale; break;
        case GLUT_KEY_RIGHT: m_CameraPos[0]+=z/m_Scale; break;
        case GLUT_KEY_UP: m_CameraPos[1]-=z/m_Scale; break;
        case GLUT_KEY_DOWN: m_CameraPos[1]+=z/m_Scale; break;
    }

    if( m_specialFunc!=NULL)
        m_specialFunc(key,x,y);

    //redraw any way
    glutPostRedisplay();
}
