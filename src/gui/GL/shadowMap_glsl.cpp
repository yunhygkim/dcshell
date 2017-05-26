#include <cstdio>
#include <cstdlib>
#include <string.h>

#ifdef _WIN32
#pragma warning(disable : 4786 4244 4305)
#endif

#include "shadowMap_glsl.h"
#include "shadowMap.frag"
#include "shadowMap.vert"

// static init
// Hold id of the framebuffer for light POV rendering
GLuint ShadowMap_GLSL::fboId;

// Z values will be rendered to this texture when using fboId framebuffer
GLuint ShadowMap_GLSL::depthTextureId;

// Use to activate/disable shadowShader
GLuint ShadowMap_GLSL::shadowShaderId;
GLuint ShadowMap_GLSL::shadowMapUniform;
GLuint ShadowMap_GLSL::shadowMapStepXUniform;
GLuint ShadowMap_GLSL::shadowMapStepYUniform;
GLuint ShadowMap_GLSL::lightShaderId;
float  ShadowMap_GLSL::lightPos[3];
float ShadowMap_GLSL::screenWidth;
float ShadowMap_GLSL::screenHeight;


ShadowMap_GLSL::ShadowMap_GLSL() {}
ShadowMap_GLSL::~ShadowMap_GLSL() {}

// Loading shader function
GLuint ShadowMap_GLSL::loadShader(unsigned char* buffer, unsigned int type)
{
	//FILE *pfile;
	GLuint handle;
	const GLcharARB* files[1];
	
	// shader Compilation variable
	GLint result;				// Compilation code result
	GLint errorLoglength ;
	char* errorLogText;
	GLsizei actualErrorLogLength;
	
	handle = glCreateShaderObjectARB(type);
	if (!handle)
	{
		//We have failed creating the vertex shader object.
		printf("! ERROR: ShadowMap_GLSL: Failed creating vertex shader object\n");
		exit(0);
	}
	
	files[0] = (const GLcharARB*)buffer;
	glShaderSourceARB(
					  handle, //The handle to our shader
					  1, //The number of files.
					  files, //An array of const char * data, which represents the source code of theshaders
					  NULL);
	
	glCompileShaderARB(handle);
	
	//Compilation checking.
	glGetObjectParameterivARB(handle, GL_OBJECT_COMPILE_STATUS_ARB, &result);
	
	// If an error was detected.
	if (!result)
	{
		//We failed to compile.
		printf("! ERROR: ShadowMap_GLSL: Shader failed compilation.\n");
		
		//Attempt to get the length of our error log.
		glGetObjectParameterivARB(handle, GL_OBJECT_INFO_LOG_LENGTH_ARB, &errorLoglength);
		
		//Create a buffer to read compilation error message
		errorLogText = (char*)malloc(sizeof(char) * errorLoglength);
		
		//Used to get the final length of the log.
		glGetInfoLogARB(handle, errorLoglength, &actualErrorLogLength, errorLogText);
		
		// Display errors.
		printf("! ERROR: ShadowMap_GLSL: %s\n",errorLogText);
		
		// Free the buffer malloced earlier
		free(errorLogText);
	}
	
	return handle;
}

void ShadowMap_GLSL::loadShadowShader()
{
	GLuint vertexShaderHandle;
	GLuint fragmentShaderHandle;
	//GLuint vertexShaderHandle2;
	//GLuint fragmentShaderHandle2;
	
	vertexShaderHandle   = loadShader(shadowMap_vert,GL_VERTEX_SHADER);
	fragmentShaderHandle = loadShader(shadowMap_frag,GL_FRAGMENT_SHADER);
	
	shadowShaderId = glCreateProgramObjectARB();
	
	glAttachObjectARB(shadowShaderId,vertexShaderHandle);
	glAttachObjectARB(shadowShaderId,fragmentShaderHandle);
	glLinkProgramARB(shadowShaderId);
	
	shadowMapUniform = glGetUniformLocationARB(shadowShaderId,"ShadowMap");
	shadowMapStepXUniform = glGetUniformLocationARB(shadowShaderId,"xPixelOffset");
	shadowMapStepYUniform = glGetUniformLocationARB(shadowShaderId,"yPixelOffset");


	// light shader
	//vertexShaderHandle2   = loadShader("lightShader.vert",GL_VERTEX_SHADER);
	//fragmentShaderHandle2 = loadShader("lightShader.frag",GL_FRAGMENT_SHADER);

	//lightShaderId = glCreateProgramObjectARB();
	//glAttachObjectARB(lightShaderId,vertexShaderHandle2);
	//glAttachObjectARB(lightShaderId,fragmentShaderHandle2);
	//glLinkProgramARB(lightShaderId);
}

void ShadowMap_GLSL::generateShadowFBO(int w, int h, float lightPx, float lightPy, float lightPz)
{
	int shadowMapWidth = w * SHADOW_MAP_RATIO;
	int shadowMapHeight = h * SHADOW_MAP_RATIO;
	screenWidth = w;
	screenHeight = h;
	lightPos[0] = lightPx;
	lightPos[1] = lightPy;
	lightPos[2] = lightPz;

	
	GLenum FBOstatus;
	
	// Try to use a texture depth component
	glGenTextures(1, &depthTextureId);
	glBindTexture(GL_TEXTURE_2D, depthTextureId);
	
	// GL_LINEAR does not make sense for depth texture. However, next tutorial shows usage of GL_LINEAR and PCF. Using GL_NEAREST
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	
	// Remove artefact on the edges of the shadowmap
	glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP );
	glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP );
	
	// This is to allow usage of shadow2DProj function in the shader
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_COMPARE_R_TO_TEXTURE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_FUNC, GL_LEQUAL);
	glTexParameteri(GL_TEXTURE_2D, GL_DEPTH_TEXTURE_MODE, GL_INTENSITY); 
	
	
	
	// No need to force GL_DEPTH_COMPONENT24, drivers usually give you the max precision if available 
	glTexImage2D( GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, shadowMapWidth, shadowMapHeight, 0, GL_DEPTH_COMPONENT, GL_UNSIGNED_BYTE, 0);
	glBindTexture(GL_TEXTURE_2D, 0);
	
	// create a framebuffer object
	glGenFramebuffersEXT(1, &fboId);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fboId);
	
	// Instruct openGL that we won't bind a color texture with the currently binded FBO
	glDrawBuffer(GL_NONE);
	glReadBuffer(GL_NONE);
	
	// attach the texture to FBO depth attachment point
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT,GL_TEXTURE_2D, depthTextureId, 0);
	
	// check FBO status
	FBOstatus = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
	if(FBOstatus != GL_FRAMEBUFFER_COMPLETE_EXT)
		printf("! ERROR: ShadowMap_GLSL: GL_FRAMEBUFFER_COMPLETE_EXT failed, CANNOT use FBO\n");
	
	// switch back to window-system-provided framebuffer (id = 0)
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
}

void ShadowMap_GLSL::setupMatrices(float position_x,float position_y,float position_z,float lookAt_x,float lookAt_y,float lookAt_z)
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60, (float)screenWidth/screenHeight,0.1,10000);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(position_x,position_y,position_z,lookAt_x,lookAt_y,lookAt_z,0,1,0);
}

void ShadowMap_GLSL::setTextureMatrix(void)
{
	static double modelView[16];
	static double projection[16];
	
	// This is matrix transform every coordinate x,y,z
	// x = x* 0.5 + 0.5 
	// y = y* 0.5 + 0.5 
	// z = z* 0.5 + 0.5 
	// Moving from unit cube [-1,1] to [0,1]  
	const GLdouble bias[16] = {	
		0.5, 0.0, 0.0, 0.0, 
		0.0, 0.5, 0.0, 0.0,
		0.0, 0.0, 0.5, 0.0,
	0.5, 0.5, 0.5, 1.0};
	
	// Grab modelview and transformation matrices
	glGetDoublev(GL_MODELVIEW_MATRIX, modelView);
	glGetDoublev(GL_PROJECTION_MATRIX, projection);
	
	
	glMatrixMode(GL_TEXTURE);
	glActiveTextureARB(GL_TEXTURE0);
	
	glLoadIdentity();	
	glLoadMatrixd(bias);
	
	// concatating all matrice into one.
	glMultMatrixd (projection);
	glMultMatrixd (modelView);
	

	// Go back to normal matrix mode
	glMatrixMode(GL_MODELVIEW);
}

void ShadowMap_GLSL::clean()
{
	glDeleteFramebuffersEXT(1, &fboId);
	glDeleteRenderbuffersEXT(1, &depthTextureId);
	glUseProgramObjectARB(0);
}

void ShadowMap_GLSL::printLog(GLuint obj)
{
	int infologLength = 0;
	int maxLength;
	
	if(glIsShader(obj))
		glGetShaderiv(obj,GL_INFO_LOG_LENGTH,&maxLength);
	else
		glGetProgramiv(obj,GL_INFO_LOG_LENGTH,&maxLength);
			
	char infoLog[4096];
 
	if (glIsShader(obj))
		glGetShaderInfoLog(obj, maxLength, &infologLength, infoLog);
	else
		glGetProgramInfoLog(obj, maxLength, &infologLength, infoLog);
 
	if (infologLength > 0)
		printf("- ShadowMap_GLSL: %s\n",infoLog);
}

void ShadowMap_GLSL::beginShadowMap()
{
	//First step: Render from the light POV to a FBO, story depth values only
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,ShadowMap_GLSL::fboId);	//Rendering offscreen
	
	//Using the fixed pipeline to render to the depthbuffer
	glUseProgramObjectARB(0);
	
	// In the case we render the shadowmap to a higher resolution, the viewport must be modified accordingly.
	glViewport(0,0,screenWidth * SHADOW_MAP_RATIO, screenHeight* SHADOW_MAP_RATIO);
	
	// Clear previous frame values
	glClear( GL_DEPTH_BUFFER_BIT);
	
	//Disable color rendering, we only want to write to the Z-Buffer
	glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE); 

	//make sure that setup matrix function does not overwrite the existing matrices
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();

	ShadowMap_GLSL::setupMatrices(lightPos[0], lightPos[1], lightPos[2], 0, 0, 0);

	// Culling switching, rendering only backface, this is done to avoid self-shadowing
	glCullFace(GL_FRONT);
	//glDisable(GL_LIGHTING);
}

void ShadowMap_GLSL::endShadowMap()
{
    //glEnable(GL_LIGHTING);
    //lightInit();

    //glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
    //glLightfv(GL_LIGHT1, GL_POSITION, light1_position);

    //Save modelview/projection matrice into texture0, also add a biais
    ShadowMap_GLSL::setTextureMatrix();


    // Now rendering from the camera POV, using the FBO to generate shadows
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,0);

    glViewport(0,0,screenWidth, screenHeight);

    //Enabling color write (previously disabled for light POV z-buffer rendering)
    glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);

    // Clear previous frame values
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    //Using the shadow shader
    glUseProgramObjectARB(ShadowMap_GLSL::shadowShaderId);
    glUniform1fARB(ShadowMap_GLSL::shadowMapStepXUniform, 1.0/ (screenWidth * SHADOW_MAP_RATIO));
    glUniform1fARB(ShadowMap_GLSL::shadowMapStepYUniform, 1.0/ (screenHeight *SHADOW_MAP_RATIO));
    glUniform1iARB(ShadowMap_GLSL::shadowMapUniform,0);
    glActiveTextureARB(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D,ShadowMap_GLSL::depthTextureId);
    //printLog(ShadowMap_GLSL::shadowShaderId);

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
}

void ShadowMap_GLSL::beginSceneRender()
{
   // glEnable(GL_DEPTH_TEST);	       // use depth buffering
    glCullFace(GL_BACK);	           // select backside of polygons for culling
    //glEnable(GL_CULL_FACE);	        // cull backside of polygons

}

void ShadowMap_GLSL::endSceneRender()
{
	glUseProgramObjectARB(0);
	//glDisable(GL_LIGHTING);
}
