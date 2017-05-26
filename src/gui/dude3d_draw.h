#pragma once
#ifndef _DUDE3D_DRAW_H_
#define _DUDE3D_DRAW_H_

#define GLEW_STATIC

#include "GL/gli.h"

//OpenGL Related Functions
void CreateGLUI();
bool InitGL();
void Display( void );
void Reshape( int w, int h);
void Keyboard( unsigned char key, int x, int y );
void TimerCallback(int value);


void drawAll();
void setupLight();
void resetCamera();
void randColor(); //defined/implemented in CD3d_draw_GL
void animate();
void moveModelUp();
void moveModelDown();
void moveModelFront();
void moveModelBack();

//operations on the hull
void trim_hull_svm();
void trim_hull_heuristic();
void trim_hull_duality();
void trim_hull_exact();
void simplify_hull();
void remesh_hull();
void save_gmodel_hulls();

//?
#endif //_DUDE3D_DRAW_H_
