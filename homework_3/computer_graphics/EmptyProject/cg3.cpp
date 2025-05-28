/*
  This program does 3D transforms.
*/

#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <iostream>
#include <string>
#include <ctime>
#include <chrono>

//#include "glut.h" //MSVC local library install
#include "glut.h" //system-wide install (or compiler default path)

const double circle = atan(1) * 8;
const double halfCircle = atan(1) * 4;
const double tau = circle; // 2 * PI = TAU
const double pi = halfCircle; // TAU / 2 = PI

const double deg2rad = circle / 360.0;
const double rad2deg = 360.0 / circle;

int g_w = 1000, g_h = 1000;

unsigned char g_prevKey = '1';

//For setting the scene properly:
float g_sceneDepth = 105;
float g_sceneZRadius = 5;

//For animations and frame time computing
bool g_animating = true;
float g_animationProgress = 0.0;
float g_animationSpeed = 1.0;
float g_secondsToCompleteAnimation = 10.0;


typedef std::chrono::steady_clock Clock;
typedef std::chrono::duration<float> FloatSeconds;
auto g_lastTime = Clock::now();


//----------------Utility functions----------------------

//FloatSeconds timeDelta(auto time) { //Works after c++20
template <typename T>
float durationInSeconds(T time0, T time1) {
  //Get a duration in floating-point seconds.
  FloatSeconds duration = time1 - time0;
  return duration.count();
}

float timeDelta() {
  //Get the time in floating-point seconds since the last call of this function.
  auto now = Clock::now();
  float duration = durationInSeconds(now, g_lastTime);
  g_lastTime = now;
  return duration;
}

void setRasterPos(float x, float y) {
  //We need to ensure we're back in the raster coordinate system.
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix(); //M0
  glLoadIdentity();
  glMatrixMode(GL_PROJECTION);
  glPushMatrix(); //P0
  glLoadIdentity();
  glRasterPos2f(x, y);
  glPopMatrix(); //P0
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix(); //M0
}

void bitmapString(void* font, const char* str) {
  //Draw a string, character-by-character.
  char cp;
  for(const char* c = str; *c != 0; ++c) {
    cp = *c; //to respect const
    glutBitmapCharacter(font, cp);
  }
}

void drawBitmapString(const char* str, float x = -2, float y = -2) {
  //Draw a string, optionally setting raster position.
  /*
    We define the convetion that both values -2 mean 'do not change
    raster position'.
  */
  if((-2 != x) || (-2 != y)) {
    //glRasterPos2f(x, y);
    setRasterPos(x, y);
  }
  //freeglut, not old glut: glutBitmapString(GLUT_BITMAP_8_BY_13, str);
  bitmapString(GLUT_BITMAP_8_BY_13, str);
}

template <typename Numtype>
void drawBitmapNumber(Numtype number, float x = -2, float y = -2) {
  //Convert a number to a string, then draw it.
  //We need the template so we don't display '2' as '2.000000'.
  if((-2 != x) || (-2 != y)) {
    //glRasterPos2f(x, y);
    setRasterPos(x, y);
  }
  bitmapString(GLUT_BITMAP_8_BY_13, std::to_string(number).c_str());
}
//^^^^^^^^^^^^^^^^^Utility functions^^^^^^^^^^^^^^^^^^

void drawCube(float sideLength = 1, float opacity = 1) {
  float l = sideLength;
  glBegin(GL_QUADS); {
    //front
    glColor4f(0.9, 0, 0, opacity);
    glVertex3f(0, 0, 0);
    glVertex3f(l, 0, 0);
    glVertex3f(l, l, 0);
    glVertex3f(0, l, 0);
    //back
    glColor4f(0, 0.9, 0.9, opacity);
    glVertex3f(0, 0, l);
    glVertex3f(0, l, l);
    glVertex3f(l, l, l);
    glVertex3f(l, 0, l);

    //left
    glColor4f(0, 0.9, 0, opacity);
    glVertex3f(0, 0, 0);
    glVertex3f(0, 0, l);
    glVertex3f(0, l, l);
    glVertex3f(0, l, 0);
    //right
    glColor4f(0.9, 0, 0.9, opacity);
    glVertex3f(l, 0, 0);
    glVertex3f(l, 0, l);
    glVertex3f(l, l, l);
    glVertex3f(l, l, 0);

    //bottom
    glColor4f(0, 0, 0.9, opacity);
    glVertex3f(0, 0, 0);
    glVertex3f(l, 0, 0);
    glVertex3f(l, 0, l);
    glVertex3f(0, 0, l);
    //top
    glColor4f(0.9, 0.9, 0, opacity);
    glVertex3f(0, l, 0);
    glVertex3f(l, l, 0);
    glVertex3f(l, l, l);
    glVertex3f(0, l, l);    
  } glEnd();
}

void drawAxes(float length = 1) {
  float l = length;
  glBegin(GL_LINES); {
    glColor4f(1, 0, 0, 1);
    glVertex3f(0, 0, 0);
    glVertex3f(l, 0, 0);
    
    glColor4f(0, 1, 0, 1);
    glVertex3f(0, 0, 0);
    glVertex3f(0, l, 0);

    glColor4f(0, 0, 1, 1);
    glVertex3f(0, 0, 0);
    glVertex3f(0, 0, l);
  } glEnd();
};

void advanceAnimationTime() {
  //How much has passed since the last frame?
  float dt = timeDelta();
  //g_animationSpeed works by accelerating simulation time
  dt *= g_animationSpeed;
  
  //the more seconds it takes for an animation to finish, the slower it progresses in [0, 1].
  // true = 1, false = 0;
  g_animationProgress += dt / g_secondsToCompleteAnimation * g_animating;
  if(g_animationProgress > 1)
    g_animationProgress -= 1; //decrease, rather than reset, for a smoother reset.
}

void Display1() {
  glColor3f(0, 0, 0);
  drawBitmapString("Keys: Space, r, R, -, +, x, X, y, Y, z, Z.", -0.95, -0.95);
  advanceAnimationTime();
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();//P0
  glLoadIdentity();
  //l, r, b, t, -n, -f
  glOrtho(-2, 2, -2, 2, g_sceneDepth - g_sceneZRadius, g_sceneDepth + g_sceneZRadius);
  
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();//M0
  //Translate so we're in the view volume.
  glTranslated(0, 0, -g_sceneDepth);
  //Rotate so we can see along the first diagonal
  glRotated(20, 1, 0, 0);
  glRotated(-45, 0, 1, 0);

  drawAxes();

  //glRotated receives angles in degrees
  //the rotation axis is the first diagonal of the space
  glRotated(g_animationProgress * 360, 1, 1, 1);
  /*
    Task 1: decompose the above rotation so it happens only alongside
    axes. I.e.: glRotated(angle, 0, 0, 1);
  */
  drawCube(1, 0.95);
  glPopMatrix();//M0
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();//P0
}



// Task 1 - Scene 1

void Display2 ()
{
	// Set the color of the text to black.
    glColor3f(0, 0, 0);

	// Display the text "Task 1 - Scene 1" at the bottom left corner of the window.
    drawBitmapString("Task 1 - Scene 1", -0.95f, -0.95f);

	// Advance the animation time.
    advanceAnimationTime();

	// Set up the orthographic projection matrix.
    glMatrixMode(GL_PROJECTION);

	// Push the current projection matrix onto the stack.
    glPushMatrix();

	// Load the identity matrix into the projection matrix.
    glLoadIdentity();

	// Set up the orthographic projection with the specified parameters (left, right, bottom, top, near, far).
    glOrtho(-2, 2, -2, 2, g_sceneDepth - g_sceneZRadius, g_sceneDepth + g_sceneZRadius);

	// Switch to the model view matrix.
    glMatrixMode(GL_MODELVIEW);

	// Push the current model view matrix onto the stack.
    glPushMatrix();

	// Translate the model view matrix to the camera position.
    glTranslated(0, 0, -g_sceneDepth);

	// Rotate the model view matrix by 20 degrees around the X-axis.
    glRotated(20, 1, 0, 0);

	// Rotate the model view matrix by -45 degrees around the Y-axis.
    glRotated(-45, 0, 1, 0);

	// Draw the axes in the scene.
    drawAxes();

	// Declare the angle of rotation around the Y-axis.
    const double angY = 50.0;

	// Rotate the model view matrix by the specified angle around the Y-axis.
    glRotated(angY, 0, 1, 0);

	// Declare the angle of rotation around the X-axis.
    const double angX = atan(1.0 / sqrt(2.0)) * rad2deg;

	// Rotate the model view matrix by the specified angle around the X-axis.
    glRotated(-angX, 1, 0, 0);

	// Rotate the model view matrix by the animation progress around the Z-axis.
    glRotated(g_animationProgress * 360.0, 0, 0, 1);

	// Rotate the model view matrix by the negative angle around the Y-axis.
    glRotated(angX, 1, 0, 0);

	// Rotate the model view matrix by the negative angle around the Y-axis.
    glRotated(-angY, 0, 1, 0);

	// Draw the cube at the origin.
    drawCube(1, 0.95f);

	// Pop the current model view matrix from the stack.
    glPopMatrix();

	// Set up the orthographic projection matrix.
    glMatrixMode(GL_PROJECTION);

	// Pop the projection matrix from the stack.
    glPopMatrix();
}



// Task 1 - Scene 2

void Display3 ()
{
	// Set the color of the text to black.
    glColor3f(0, 0, 0);

	// Display the text "Task 1 - Scene 2" at the bottom left corner of the window.
    drawBitmapString("Task 1 - Scene 2", -0.95f, -0.95f);

	// Advance the animation time.
    advanceAnimationTime();

	// Set up the orthographic projection matrix.
    glMatrixMode(GL_PROJECTION);

	// Push the current projection matrix onto the stack.
    glPushMatrix();

	// Load the identity matrix into the projection matrix.
    glLoadIdentity();

	// Set up the orthographic projection with the specified parameters (left, right, bottom, top, near, far).
    glOrtho(-2, 2, -2, 2, g_sceneDepth - g_sceneZRadius, g_sceneDepth + g_sceneZRadius);

	// Switch to the model view matrix.
    glMatrixMode(GL_MODELVIEW);

	// Push the current model view matrix onto the stack.
    glPushMatrix();

	// Translate the model view matrix to the camera position.
    glTranslated(0, 0, -g_sceneDepth);

	// Rotate the model view matrix by 20 degrees around the X-axis.
    glRotated(20, 1, 0, 0);

	// Rotate the model view matrix by -45 degrees around the Y-axis.
    glRotated(-45, 0, 1, 0);

	// Declare the angles of rotation around the Y and X axes.
    const double angY = 45.0;
    const double angX = atan(1.0 / sqrt(2.0)) * rad2deg;

	// Declare the full rotation angle.
    const double full = 360.0;

	// Declare the animation progress variable.
    double p = fmod(g_animationProgress, 1.0);

	// Define a lambda function to calculate the phase based on the animation progress.
    auto phase = [&](double start)
    {
        double t = (p - start) / 0.2;
        
        if (t < 0)
            t = 0;
        
        if (t > 1)
            t = 1;
        
        return t;
    };

	// Declare the rotation angles based on the animation progress.
    double y1 = angY * phase(0.0);
    double x1 = -angX * phase(0.2);
    double z = full * phase(0.4);
    double x2 = angX * phase(0.6);
    double y2 = -angY * phase(0.8);

	// Apply the rotations to the model view matrix.
    glRotated(y1, 0, 1, 0);
    glRotated(x1, 1, 0, 0);
    glRotated(z, 0, 0, 1);
    glRotated(x2, 1, 0, 0);
    glRotated(y2, 0, 1, 0);

	// Draw the cube at the origin.
    drawCube(1, 0.95f);

	// Disable depth testing to draw the axes on top of the cube.
    glDisable(GL_DEPTH_TEST);

	// Set the line width for the axes.
    glLineWidth(4);

	// Draw the axes in the scene (length = 1.5).
    drawAxes(1.5f);

	// Set the width of the lines for the axes.
    glLineWidth(3);

	// Enable depth testing again.
    glEnable(GL_DEPTH_TEST);

	// Set the color of the axes to black.
    glPopMatrix();

	// Set the color of the text to black.
    glMatrixMode(GL_PROJECTION);

	// Pop the projection matrix from the stack.
    glPopMatrix();
}



// Task 2 - Scene 1

void Display4 ()
{
	// Set the color of the text to black.
    glColor3f(0, 0, 0);

	// Display the text "Task 2 - Scene 1" at the bottom left corner of the window.
    drawBitmapString("Task 2 - Scene 1", -0.95f, -0.95f);

	// Advance the animation time.
    advanceAnimationTime();

	// Set up the orthographic projection matrix.
    glMatrixMode(GL_PROJECTION);

	// Push the current projection matrix onto the stack.
    glPushMatrix();

	// Load the identity matrix into the projection matrix.
    glLoadIdentity();

	// Set up the orthographic projection with the specified parameters (left, right, bottom, top, near, far).
    glOrtho(-2, 2, -2, 2, g_sceneDepth - g_sceneZRadius, g_sceneDepth + g_sceneZRadius);

	// Switch to the model view matrix.
    glMatrixMode(GL_MODELVIEW);

	// Push the current model view matrix onto the stack.
    glPushMatrix();

	// Translate the model view matrix to the camera position.
    glTranslated(0, 0, -g_sceneDepth);

	// Draw the axes in the scene (length = 2).
    drawAxes(2.0f);

	// Draw the first cube at the origin.
    drawCube(1.0f, 1.0f);

	// Push the current model view matrix onto the stack.
    glPushMatrix();

	// Translate the model view matrix to the right by 1.1 units.
    glTranslatef(1.1f, 0.0f, 0.0f);

	// Rotate the model view matrix by 45 degrees around the Y-axis.
    drawCube(1.0f, 1.0f);

	// Pop the current model view matrix from the stack.
    glPopMatrix();

	// Pop the current model view matrix from the stack.
    glPopMatrix();

	// Set up the orthographic projection matrix.
    glMatrixMode(GL_PROJECTION);

	// Pop the model view matrix from the stack.
    glPopMatrix();
}



// Task  2 - Scene 2

void Display5 ()
{
	// Set the color of the text to black.
    glColor3f(0, 0, 0);

	// Display the text "Task 2 - Scene 2" at the bottom left corner of the window.
    drawBitmapString("Task 2 - Scene 2", -0.95f, -0.95f);

	// Advance the animation time.
    advanceAnimationTime();

    // Set up the orthographic projection matrix.
    glMatrixMode(GL_PROJECTION);

	// Push the current projection matrix onto the stack.
    glPushMatrix();

	// Load the identity matrix into the projection matrix.
    glLoadIdentity();

	// Set up the orthographic projection with the specified parameters (left, right, bottom, top, near, far).
    glOrtho(-2, 2, -2, 2, g_sceneDepth - g_sceneZRadius, g_sceneDepth + g_sceneZRadius);

	// Switch to the model view matrix.
    glMatrixMode(GL_MODELVIEW);

	// Push the current model view matrix onto the stack.
    glPushMatrix();

	// Translate the model view matrix to the camera position.
    glTranslated(0, 0, -g_sceneDepth);

	// Rotate the model view matrix by -10 degrees around the Y-axis.
    glRotated(-10.0, 0, 1, 0);

	// Draw the axes in the scene (length = 2).
    drawAxes(2.0f);

	// Draw the first cube at the origin.
    drawCube(1.0f, 1.0f);

	// Push the current model view matrix onto the stack.
    glPushMatrix();

	// Translate the model view matrix to the right by 1.1 units.
    glTranslatef(1.1f, 0.0f, 0.0f);

	// Rotate the model view matrix by 45 degrees around the Y-axis.
    drawCube(1.0f, 1.0f);

	// Push the current model view matrix onto the stack.
    glPopMatrix();

	// Pop the current model view matrix from the stack.
    glPopMatrix();

	// Set up the orthographic projection matrix.
    glMatrixMode(GL_PROJECTION);

	// Pop the model view matrix from the stack.
    glPopMatrix();
}



// Task 3 - Scene 1

void Display6 ()
{
	// Set the color of the text to black.
    glColor3f(0, 0, 0);

	// Display the text "Task  3 - Scene 1" at the bottom left corner of the window.
    drawBitmapString("Task 3 - Scene 1", -0.95f, -0.95f);

	// Advance the animation time.
    advanceAnimationTime();

	// Set up the orthographic projection matrix.
    glMatrixMode(GL_PROJECTION);

	// Push the current projection matrix onto the stack.
    glPushMatrix();

	// Load the identity matrix into the projection matrix.
    glLoadIdentity();

	// Set up the orthographic projection with the specified parameters (left, right, bottom, top, near, far).
    glOrtho(-2, 2, -2, 2, g_sceneDepth - g_sceneZRadius, g_sceneDepth + g_sceneZRadius);

	// Switch to the model view matrix.
    glMatrixMode(GL_MODELVIEW);

	// Push the current model view matrix onto the stack.
    glPushMatrix();

	// Translate the model view matrix to the camera position.
    glTranslated(0, 0, -g_sceneDepth);

	// Draw the axes in the scene (length = 2).
    drawAxes(2.0f);

	// Draw the first cube at the origin.
    drawCube(1.0f, 1.0f);

	// Push the current model view matrix onto the stack.
    glPushMatrix();

	// Translate the model view matrix to the right by 1.1 units.
    glTranslatef(1.1f, 0.0f, 0.0f);

	// Rotate the model view matrix by 45 degrees around the Y-axis.
    glScalef(1.0f, 0.8f, 1.0f);

	// Draw the second cube at the translated position.
    drawCube(1.0f, 1.0f);

	// Push the current model view matrix onto the stack.
    glPopMatrix();

	// Pop the current model view matrix from the stack.
    glPopMatrix();

	// Set up the orthographic projection matrix.
    glMatrixMode(GL_PROJECTION);

	// Pop the model view matrix from the stack.
    glPopMatrix();
}



// Task 3 - Scene 2

void Display7 ()
{
	// Set the color of the text to black.
    glColor3f(0, 0, 0);

	// Display the text "Task  3 - Scene 2" at the bottom left corner of the window.
    drawBitmapString("Task  3 - Scene 2", -0.95f, -0.95f);

	// Advance the animation time.
    advanceAnimationTime();

	// Set up the orthographic projection matrix.
    glMatrixMode(GL_PROJECTION);

	// Push the current projection matrix onto the stack.
    glPushMatrix();

	// Load the identity matrix into the projection matrix.
    glLoadIdentity();

	// Set up the orthographic projection with the specified parameters (left, right, bottom, top, near, far).
    glOrtho(-2, 2, -2, 2, g_sceneDepth - g_sceneZRadius, g_sceneDepth + g_sceneZRadius);

	// Switch to the model view matrix.
    glMatrixMode(GL_MODELVIEW);

	// Push the current model view matrix onto the stack.
    glPushMatrix();

	// Translate the model view matrix to the camera position.
    glTranslated(0, 0, -g_sceneDepth);

	// Rotate the model view matrix by -10 degrees around the Y-axis.
    glRotated(-10.0, 0, 1, 0);

	// Draw the axes in the scene (length = 2).
    drawAxes(2.0f);

	// Draw the first cube at the origin.
    drawCube(1.0f, 1.0f);

	// Push the current model view matrix onto the stack.
    glPushMatrix();

	// Translate the model view matrix to the right by 1.1 units.
    glTranslatef(1.1f, 0.0f, 0.0f);

	// Rotate the model view matrix by 45 degrees around the Y-axis.
    glScalef(1.0f, 0.8f, 1.0f);

	// Draw the second cube at the translated position.
    drawCube(1.0f, 1.0f);

	// Push the current model view matrix onto the stack.
    glPopMatrix();

	// Pop the current model view matrix from the stack.
    glPopMatrix();

	// Set up the orthographic projection matrix.
    glMatrixMode(GL_PROJECTION);

	// Pop the model view matrix from the stack.
    glPopMatrix();
}



// Task 4 - Scene 1

void Display8 ()
{
	// Set the color of the text to black.
    glColor3f(0, 0, 0);

	// Display the text "Task 4 - Scene 1t" at the bottom left corner of the window.
    drawBitmapString("Task 4 - Scene 1", -0.95f, -0.95f);

	// Advance the animation time.
    advanceAnimationTime();

	// Set up the orthographic projection matrix.
    glMatrixMode(GL_PROJECTION);

	// Push the current projection matrix onto the stack.
    glPushMatrix();

	// Load the identity matrix into the projection matrix.
    glLoadIdentity();

	// Set up the orthographic projection with the specified parameters (left, right, bottom, top, near, far).
    glOrtho(-2, 2, -2, 2, g_sceneDepth - g_sceneZRadius, g_sceneDepth + g_sceneZRadius);

	// Switch to the model view matrix.
    glMatrixMode(GL_MODELVIEW);

	// Push the current model view matrix onto the stack.
    glPushMatrix();

	// Translate the model view matrix to the camera position.
    glTranslated(0, 0, -g_sceneDepth);

	// Rotate the model view matrix by -30 degrees around the X-axis.
    glRotated(-30.0, 1, 0, 0);

	// Rotate the model view matrix by 30 degrees around the Y-axis.
    glRotated(30.0, 0, 1, 0);

	// Draw the axes in the scene (length = 1.5).
    drawAxes(1.5f);

	// Draw the first cube at the origin.
    drawCube(1.0f, 1.0f);

	// Pop the current model view matrix from the stack.
    glPopMatrix();

	// Set up the orthographic projection matrix.
    glMatrixMode(GL_PROJECTION);

	// Pop the projection matrix from the stack.
    glPopMatrix();
}



// Task 4 - Scene 2

void Display9 ()
{
	// Set the color of the text to black.
    glColor3f(0, 0, 0);

	// Display the text "Task 4 - Scene 2" at the bottom left corner of the window.
    drawBitmapString("Task 4 - Scene 2", -0.95f, -0.95f);

	// Advance the animation time.
    advanceAnimationTime();

	// Set up the orthographic projection matrix.
    glMatrixMode(GL_PROJECTION);

	// Push the current projection matrix onto the stack.
    glPushMatrix();

	// Load the identity matrix into the projection matrix.
    glLoadIdentity();

	// Set up the orthographic projection with the specified parameters (left, right, bottom, top, near, far).
    glOrtho(-2, 2, -2, 2, g_sceneDepth - g_sceneZRadius, g_sceneDepth + g_sceneZRadius);

	// Switch to the model view matrix.
    glMatrixMode(GL_MODELVIEW);

	// Push the current model view matrix onto the stack.
    glPushMatrix();

	// Translate the model view matrix to the camera position.
    glTranslated(0, 0, -g_sceneDepth);

	// Rotate the model view matrix by -30 degrees around the X-axis.
    glRotated(30.0, 1, 0, 0);

	// Rotate the model view matrix by -30 degrees around the Y-axis.
    glRotated(-30.0, 0, 1, 0);

	// Draw the axes in the scene (length = 1.5).
    drawAxes(1.5f);

	// Draw the first cube at the origin.
    drawCube(1.0f, 1.0f);

	// Pop the current model view matrix from the stack.
    glPopMatrix();

	// Set up the orthographic projection matrix.
    glMatrixMode(GL_PROJECTION);

	// Pop the projection matrix from the stack.
    glPopMatrix();
}



void Display10() {
}

void init(void) {
  glColor3f(1, 0, 0); //Just a starting default drawing colour.
  glClearColor(1.0,1.0,1.0,1.0);
  glLineWidth(3);
  glPointSize(4);
  //glPolygonMode(GL_FRONT, GL_LINE);
  //Smoothing - don't enable it without multisampling, which Freeglut has, but Glut doesn't.
  /*glEnable(GL_SMOOTH);
  glEnable(GL_POINT_SMOOTH);
  glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
  glEnable(GL_LINE_SMOOTH);
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
  glEnable(GL_POLYGON_SMOOTH);
  glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);*/

  //Alpha-blending
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  //glBlendFunc(GL_ONE_MINUS_DST_ALPHA,GL_DST_ALPHA); 
  //Actually enable depth-testing (the init function just requests the display capability).
  glEnable(GL_DEPTH_TEST);
  /*
    Display a pixel (well, fragment, but we haven't talked about them yet) only if its
    depth is smaller than that of other fragments at the same (x, y) position.
    Smaller meaning closer to the eye position / camera.
    You may notice drawing errors when mixed with transparency.
  */
  glDepthFunc(GL_LESS);
  
 
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

void Display(void) {
  // Clear the colour buffer. See init();
  // Also, clear the depth buffer, which stores depth information for pixels.
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  switch(g_prevKey) {
  case '1':
    Display1();
    break;
  case '2':
    Display2();
    break;
  case '3':
    Display3();
    break;
  case '4':
    Display4();
    break;
  case '5':
    Display5();
    break;
  case '6':
    Display6();
    break;
  case '7':
    Display7();
    break;
  case '8':
    Display8();
    break;
  case '9':
    Display9();
    break;
  case '0':
    Display10();
    break;
  default:
    break;
  }
  //glFlush();
  //As we're working double-buffered, we need to swap the buffers once we've finished drawing an image.
  glutSwapBuffers();
}

void Reshape(int w, int h) {
  g_w = w;
  g_h = h;
  glViewport(0, 0, (GLsizei) w, (GLsizei) h);
}
void KeyboardFunc(unsigned char key, int x, int y) {
  switch(key) {
  case 27: // escape
    exit(0);
    break;
  case 32: //space: pause / unpause animation
    g_animating = !g_animating;
    break;
  case 'r': //reset animation to first frame
    g_animationProgress = 0;
    g_animationSpeed = 1.0;
    break;
  case 'R': //reset transforms
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
  case '-': //decrease animation speed
    g_animationSpeed /= 1.5;
    break;
  case '+': //increase animation speed
    g_animationSpeed *= 1.5;
    break;
  case 'x': //rotate along x, trig-
    glMatrixMode(GL_MODELVIEW);
    glTranslated(0, 0, -g_sceneDepth);
    glRotated(-10 * g_animationSpeed, 1, 0, 0);
    glTranslated(0, 0, g_sceneDepth);
    break;
  case 'X': //rotate along x, trig+
    glMatrixMode(GL_MODELVIEW);
    glTranslated(0, 0, -g_sceneDepth);
    glRotated(10 * g_animationSpeed, 1, 0, 0);
    glTranslated(0, 0, g_sceneDepth);
    break;
  case 'y':
    glMatrixMode(GL_MODELVIEW);
    glTranslated(0, 0, -g_sceneDepth);
    glRotated(-10 * g_animationSpeed, 0, 1, 0);
    glTranslated(0, 0, g_sceneDepth);
    break;
  case 'Y':
    glMatrixMode(GL_MODELVIEW);
    glTranslated(0, 0, -g_sceneDepth);
    glRotated(10 * g_animationSpeed, 0, 1, 0);
    glTranslated(0, 0, g_sceneDepth);
    break;
  case 'z':
    glMatrixMode(GL_MODELVIEW);
    glTranslated(0, 0, -g_sceneDepth);
    glRotated(-10 * g_animationSpeed, 0, 0, 1);
    glTranslated(0, 0, g_sceneDepth);
    break;
  case 'Z': 
    glMatrixMode(GL_MODELVIEW);
    glTranslated(0, 0, -g_sceneDepth);
    glRotated(10 * g_animationSpeed, 0, 0, 1);
    glTranslated(0, 0, g_sceneDepth);
    break;
  default:
    //Only change the image if a 'special' key wasn't pressed.
    g_prevKey = key;
  }
  //The proper way to ask glut to redraw the window.
  glutPostRedisplay();
}

/*
  Callback upon mouse press or release.
  The button can be:
  GLUT_LEFT_BUTTON, GLUT_MIDDLE_BUTTON, GLUT_RIGHT_BUTTON
  (and further for mousewheel and other mouse buttons)
  The state can be either GLUT_DOWN or  GLUT_UP, for
  a pressed or released button.
  (x, y) are the coordinates of the mouse.
*/
void MouseFunc(int button, int state, int x, int y) {
  std::cout<< "Mouse button ";
  std::cout<<( (button == GLUT_LEFT_BUTTON) ? "left" : ((button == GLUT_RIGHT_BUTTON) ? "right": "middle") ) << " ";
  std::cout<< ( (state == GLUT_DOWN) ? "pressed" : "released" );
  std::cout<< " at coordinates: " << x <<" x " << y << std::endl;
}

int main(int argc, char** argv) {
  glutInit(&argc, argv);
  glutInitWindowSize(g_w, g_h);
  glutInitWindowPosition(-1, -1);
  glutInitDisplayMode (GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGBA); //Double buffering and depth testing / '3D'
  glutCreateWindow (argv[0]);
  init();
  glutReshapeFunc(Reshape);
  glutKeyboardFunc(KeyboardFunc);
  glutMouseFunc(MouseFunc);
  glutDisplayFunc(Display);
  //For animations
  glutIdleFunc(Display);
  glutMainLoop();

  return 0;
}
