// 
//  main.cc
//  convexHullTest/rectangularVolume
//  
//  Created by Robert Snapp on 2008-07-21.
//  Copyright 2008 Robert R. Snapp. All rights reserved.
// 

#include "env.h" // OpenGL and Glut headers
#include "polyhedron.h"
#include "convexhull.h"
#include "ntuple.h"
#include <iostream>
#include <list>

using namespace std;

// Model parameters
int const npts = 1000;  // number of points generated 
int const r3[3] = {100, 50, 75}; // Bounds for the rectangular region in test 3.


// GLUT parameters
const int DefaultWindowWidth  = 1024*1.618;
const int DefaultWindowHeight = 1024;
const int DefaultWindowX      =  10;
const int DefaultWindowY      =  50;

// OpenGL constants:
const GLdouble rho_def   = 6.0;
const GLdouble theta_def = 45.;
const GLdouble phi_def   = 30.;
const double radiansPerDegree = M_PI/180.0;

const GLdouble centerX_def = 0.0;
const GLdouble centerY_def = 0.0;
const GLdouble centerZ_def = 0.0;

const GLdouble upX_def = 0.0;
const GLdouble upY_def = 0.0;
const GLdouble upZ_def = 1.0;

const GLdouble fieldOfViewY_def = 45.;
const GLdouble aspect_def       = 1.0;
const GLdouble near_def         = 0.0001;
const GLdouble far_def          = 10000;

const int defaultWidth = 512;
const int defaultHeight = 512;

// OpenGL viewing variables:
static int windowWidth  = defaultWidth; 
static int windowHeight = defaultHeight;

static GLdouble rho  = rho_def;
static GLdouble theta = theta_def;
static GLdouble phi   = phi_def;

static GLdouble centerX = centerX_def;
static GLdouble centerY = centerY_def;
static GLdouble centerZ = centerZ_def;

static GLdouble upX     = upX_def;
static GLdouble upY     = upY_def;
static GLdouble upZ     = upZ_def;

static GLdouble fieldOfViewY = fieldOfViewY_def;
static GLdouble aspect = aspect_def;
static GLdouble near   = near_def;
static GLdouble far    = far_def;

// OpenGL lighting variables
const int n = 4;
static GLfloat light_position[n][4] = {{ 1.0,  1.0,  1.0,  0.0},
									   {-1.0, -1.0,  1.0,  0.0},
                                       {-1.0,  1.0,  1.0,  0.0},
                                       { 0.0,  0.0,  1.0,  0.0}};
static GLfloat light_diffuse[n][4]  = {{1.0, 1.0, 1.0, 1.0},   // white
									   {1.0, 0.0, 0.0, 1.0},
                                       {0.0, 1.0, 0.0, 1.0},
                                       {0.0, 0.0, 1.0, 1.0}};  // red
static GLfloat dark[] = {0.0, 0.0, 0.0, 1.0};
static bool light_switch[n] = {true, false, false, false};
static GLint light[n] = {GL_LIGHT0, GL_LIGHT1, GL_LIGHT2, GL_LIGHT3};

// OpenGL materials
enum Material {BRASS = 0, BRONZE, CHROME, COPPER, 
							 GOLD, PEWTER, SILVER, POLISHED_SILVER};
void setMaterial(Material m);
const int nMaterials = 8;
static GLfloat rhoAmbient[nMaterials][4] = {
	{0.329412, 0.223529, 0.027451, 1.0},
	{0.2125,   0.1275,   0.054,    0.9},
	{0.25,     0.25,     0.25,     1.0},
	{0.19125,  0.0735,   0.0225,   1.0},
	{0.24725,  0.1995,   0.0745,   1.0},
	{0.10588,  0.058824, 0.113725, 1.0},
    {0.19225,  0.19225,  0.19225,  0.6},
	{0.23125,  0.23125,  0.23125,  1.0}};

static GLfloat rhoDiffuse[nMaterials][4] = {
	{0.780392, 0.568627, 0.113725, 1.0},
	{0.714,    0.4284,   0.18144,  0.9},
	{0.4,      0.4,      0.4,      1.0},
	{0.7038,   0.27048,  0.0828,   1.0},
	{0.75164,  0.60648,  0.22648,  1.0},
	{0.427451, 0.470588, 0.541176, 1.0},
	{0.50754,	 0.50754,  0.50754,  0.6},
	{0.2775,   0.2755,   0.2755,   1.0}};

static GLfloat rhoSpecular[nMaterials][4] = {
	{0.992157, 0.941176, 0.807843, 1.0},
	{0.393548, 0.271906, 0.166721, 0.9},
	{0.774597, 0.774597, 0.774597, 1.0},  // Chrome
	{0.256777, 0.137622, 0.086014, 1.0},
	{0.628281, 0.555802, 0.366065, 1.0},
	{0.3333,   0.3333,   0.521569, 1.0},  // Pewter
	{0.508273, 0.508273, 0.508273, 1.0},
	{0.773911, 0.773911, 0.773911, 1.0}};

static GLfloat specularExponent[nMaterials] = {
	27.8972, // Brass
	25.6,    // Bronze
	76.8,    // Chrome
	12.8,	   // Copper
	51.2,		 // Gold
	 9.84615, // Pewter
    51.2,    // Silver
    89.6};   // Polished Silver



// Global polyhedral models
ConvexHull<int> hull;
vector<ntuple<int,3> > ipoints;

// OpenGL display list indices:
GLuint theFaces;
GLuint thePoints;
GLuint theAxes;

void initLights(void);
void toggleLight(int i);
void updateLights(void);


// Initialization

// display an arrow with a conical head, of unit length, based at the origin,
// oriented along the positive z-axis.
void displayArrow() {
  GLfloat origin[3] = {0.0, 0.0, 0.0};
  GLfloat zHat[3]   = {0.0, 0.0, 1.0};

  glBegin(GL_LINES);
    glVertex3fv(origin);
    glVertex3fv(zHat);
  glEnd();

  glPushMatrix();
    glTranslatef(0.0, 0.0, 1.0);
    glutSolidCone(0.02, 0.08, 10, 10);
  glPopMatrix();

  return;
}



// displayAxes draws a cartesian triple
void displayAxes(float length, float alpha) {
  GLfloat red[4]   = {1.0, 0.0, 0.0, 1.0};
  GLfloat green[4] = {0.0, 0.8, 0.0, 1.0};
  GLfloat blue[4]  = {0.2, 0.2, 1.0, 1.0};
  // GLfloat black[4]  = {0.0, 0.0, 0.0, 1.0};


  // update alphas
  red[3] = green[3] = blue[3] = alpha;

#ifdef COMMENT
  GLuint zArrow = glGenLists(1);
  glNewList(zArrow, GL_COMPILE);
  glPushMatrix();
    glScalef(length, length, length);
    displayArrow();
  glPopMatrix();
  glEndList();
#endif

  glScalef(length, length, length);
  glColor4fv(blue);
  displayArrow();             // z axis

  glPushMatrix();
	glRotated(90, -1, 0, 0);
	glColor4fv(green);
	displayArrow();			  // y axis
  glPopMatrix();

  glPushMatrix();
	glRotated(90, 0, 1, 0);
	glColor4fv(red);
	displayArrow();			  // x axis
  glPopMatrix();

  return;
}

void displayPoints() {
  glBegin(GL_POINTS);
     for(size_t i = 0; i < ipoints.size(); i++) {
	   ntuple<int,3> v = ipoints[i];
	   glVertex3f(static_cast<GLfloat>(v[0]),
				  static_cast<GLfloat>(v[1]), 
				  static_cast<GLfloat>(v[2]));
	 }
  glEnd();
  return;
}

void createConvexHull(bool verbose) {

  for(int i = 0; i < npts; i++) {
	ntuple<int, 3> x;
	for(int j = 0; j < 3; j++) {
	  x[j] = rand() % r3[j];
	}
	ipoints.push_back(x);
  }

  hull.compute(ipoints);
  if (verbose) {
	hull.describe(cout);
  }
}

void init() {

  
  createConvexHull(true);


  glClearColor( 0.0, 0.0, 0.0, 0.0);  // Black
  glShadeModel(GL_SMOOTH);            // Gouraud Shading
  setMaterial(SILVER);                // see material.cc
  glEnable(GL_DEPTH_TEST);            // Z-buffer mode
  initLights();                       // see light.cc

  ntuple<int,3> hmax = hull.getMax();
  ntuple<int,3> hmin = hull.getMin();
  ntuple<int,3> center = (hmax + hmin)/2;
  ntuple<int,3> diff   = hmax - hmin;
  double diameter = max(max(abs(static_cast<double>(diff[0])), 
							abs(static_cast<double>(diff[1]))), 
						abs(static_cast<double>(diff[2])));

  cout << "hmin = " << hmin << endl;
  cout << "hmax = " << hmax << endl;
  cout << "center = " << center << endl;
  cout << "diameter = " << diameter << endl;

  GLdouble s = 2.0/static_cast<double>(diameter);

  glScaled(s, s, s);
  glTranslatef(static_cast<GLfloat>(-center[0]), 
			   static_cast<GLfloat>(-center[1]), 
			   static_cast<GLfloat>(-center[2]));

  // Create top level display lists that render the points and the hull.
  theFaces = glGenLists(1);
  glNewList(theFaces, GL_COMPILE);
    glPushMatrix();
      glScaled(s, s, s);
      glTranslatef(static_cast<GLfloat>(-center[0]), 
				   static_cast<GLfloat>(-center[1]), 
				   static_cast<GLfloat>(-center[2]));
      hull.display();
    glPopMatrix();
  glEndList();

  thePoints = glGenLists(1);
  glNewList(thePoints, GL_COMPILE);
    glPushMatrix();
      glScaled(s, s, s);
      glTranslatef(static_cast<GLfloat>(-center[0]), 
				   static_cast<GLfloat>(-center[1]), 
				   static_cast<GLfloat>(-center[2]));
      displayPoints();
    glPopMatrix();
  glEndList();

  theAxes = glGenLists(1);
  glNewList(theAxes, GL_COMPILE);
     displayAxes(2.0, 0.8);
  glEndList();

  return;
}


//////////// OpenGL callbacks:

// Called in response to a redisplay event.
void display() {

  GLdouble theta_rad = theta*radiansPerDegree;
  GLdouble phi_rad   = phi*radiansPerDegree;

  GLdouble eyeX = rho*sin(theta_rad)*cos(phi_rad);
  GLdouble eyeY = rho*sin(theta_rad)*sin(phi_rad);
  GLdouble eyeZ = rho*cos(theta_rad);

  glClearColor(1.0, 1.0, 1.0, 1.0); // a white background saves printer toner.
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();
  gluLookAt(eyeX, eyeY, eyeZ,
			centerX, centerY, centerZ,
			upX, upY, upZ);

  //tetrahedron.display();
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_BLEND);

  glDisable(GL_LIGHTING);
  glLineWidth(2.0);
  glCallList(theAxes);


  glPointSize(6);
  glEnable(GL_POINT_SMOOTH);  // antialias the points
  glColor4f(0.8, 0.1, 0.2, 0.6);
  glCallList(thePoints);

  glEnable(GL_LIGHTING);
  glEnable(GL_NORMALIZE);
  glPolygonMode(GL_FRONT, GL_FILL);
  glCullFace(GL_BACK);
  glEnable(GL_CULL_FACE);
 

  // Polygon offset is an OpenGL kludge that helps render lines and faces simultaneously.
  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(2.0, 0.02);

  glCallList(theFaces);
  glDisable(GL_POLYGON_OFFSET_FILL);

  glDisable(GL_LIGHTING);
  glPolygonMode(GL_FRONT, GL_LINE);
  glLineWidth(1.0);
  glColor4f(0.0, 0.0, 0.0, 1.0);
  glCallList(theFaces);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

 

  glFlush();  /* Single buffered, so needs a flush. */
}

// Called after a keyboard event is received.
void keyboard(unsigned char c, int x, int y) {
	switch(c) {

	case 'h': // home
		theta = theta_def;
		phi   = phi_def;
		rho   = rho_def;
		break;
	case 'n':											// North
		theta = max(theta - 5.0, 5.0);
		break;
	case 's':											// South
		theta = min(theta + 5.0, 175.0);
		break;
	case 'e':											// East
		phi += 5.0;
		break;
	case 'w':											// West
		phi -= 5.0;
		break;
	case 'u':											// Up
		rho += 0.1;
		break;
	case 'd':											// Down
		rho = max(rho - 0.1, 0.1);
		break;
	case 'q':											// Quit
		exit(0);
		break;

	case '!':
		setMaterial(BRASS);
		break;
	case '@':
		setMaterial(BRONZE);
		break;
	case '#':
		setMaterial(CHROME);
		break;
	case '$':
		setMaterial(COPPER);
		break;
	case '%':
		setMaterial(GOLD);
		break;
	case '^':
		setMaterial(PEWTER);
		break;
	case '&':
		setMaterial(SILVER);
		break;
	case '*':
		setMaterial(POLISHED_SILVER);
		break;

	case '0':
		toggleLight(0);
		break;

	case '1':
		toggleLight(1);
		break;

	case '2':
		toggleLight(2);
		break;

	case '3':
		toggleLight(3);
		break;

	}
	glutPostRedisplay();
	return;
}


// Called when the user changes the dimensions of the GLUT window
void reshape(int w, int h) {  
  int x = 0;
  int y = 0;
  int s;		// size of the viewport side
  windowWidth  = w;
  windowHeight = h;

  if (w > h) {
	s = h;
	x = (w - h)/2;
  } else {
	s = w;
	y = (h - w)/2;
  }
  aspect = double(windowWidth)/double(windowHeight);
  //  glViewport(x, y, s, s); 
  
  glViewport(0, 0, windowWidth, windowHeight);

  glMatrixMode(GL_PROJECTION);     /* Start modifying the projection matrix. */
  glLoadIdentity();                /* Reset project matrix. */
  gluPerspective(fieldOfViewY, aspect, near, far);
  glMatrixMode(GL_MODELVIEW);
}


// OpenGL lighting functions:

void turnOnLight(int i) {
	glLightfv(light[i], GL_AMBIENT,  light_diffuse[i]);
	glLightfv(light[i], GL_DIFFUSE,  light_diffuse[i]);
	//	glLightfv(light[i], GL_SPECULAR, light_diffuse[i]);
	return;
}

void turnOffLight(int i) {
	glLightfv(light[i], GL_AMBIENT,  dark);
	glLightfv(light[i], GL_DIFFUSE,  dark);
	glLightfv(light[i], GL_SPECULAR, dark);
	return;
}

void setLight(int i) {
	if (light_switch[i]) {
		turnOnLight(i);
		cout << "Enable light " << i << "." << endl;
	} else {
		//		glDisable(light[i]);
		turnOffLight(i);

		cout << "Disable light " << i << "." << endl;
	}
	return;
}
	
void toggleLight(int i) {
	light_switch[i] = ! light_switch[i];
	setLight(i);
	glutPostRedisplay();
	return;
}

void setLights(void) {
	for (int i = 0; i < n; i++) {
		glEnable(light[i]);
		setLight(i);
	}	
}

void initLights() {
	for(int i = 0; i < n; i++) {
		glLightfv(light[i], GL_POSITION, light_position[i]);
	}
	glEnable(GL_LIGHTING);
	setLights();
}		


// OpenGL material

void setMaterial(Material m) {
	glMaterialfv(GL_FRONT, GL_AMBIENT, rhoAmbient[m]);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, rhoDiffuse[m]);
	glMaterialfv(GL_FRONT, GL_SPECULAR, rhoSpecular[m]);
	glMaterialf(GL_FRONT, GL_SHININESS, specularExponent[m]);
  return;
}



int main (int argc, char *argv[]) {
  glutInitWindowSize(DefaultWindowWidth, DefaultWindowHeight);
  glutInitWindowPosition(DefaultWindowX, DefaultWindowY);
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH);
	
  glutCreateWindow("polyhedronTest");

  init();

  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(keyboard);
  glutMainLoop();


  return 0;
}

