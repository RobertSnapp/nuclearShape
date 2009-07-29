// 
//  main.cc
//  polyhedronTest
//  
//  Created by Robert Snapp on 2007-02-17.
//  Copyright 2007 Robert R. Snapp. All rights reserved.
// 

#include "env.h" // OpenGL and Glut headers
#include "polyhedron.h"
#include "ntuple.h"
#include <iostream>
#include <list>

using namespace std;

// GLUT parameters
const int DefaultWindowWidth  = 512;
const int DefaultWindowHeight = 512;
const int DefaultWindowX      =  50;
const int DefaultWindowY      =  50;

// OpenGL constants:
const GLdouble rho_def   = 2.0;
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
	{0.2125,   0.1275,   0.054,    1.0},
	{0.25,     0.25,     0.25,     1.0},
	{0.19125,  0.0735,   0.0225,   1.0},
	{0.24725,  0.1995,   0.0745,   1.0},
	{0.10588,  0.058824, 0.113725, 1.0},
  {0.19225,  0.19225,  0.19225,  1.0},
	{0.23125,  0.23125,  0.23125,  1.0}};

static GLfloat rhoDiffuse[nMaterials][4] = {
	{0.780392, 0.568627, 0.113725, 1.0},
	{0.714,    0.4284,   0.18144,  1.0},
	{0.4,      0.4,      0.4,      1.0},
	{0.7038,   0.27048,  0.0828,   1.0},
	{0.75164,  0.60648,  0.22648,  1.0},
	{0.427451, 0.470588, 0.541176, 1.0},
	{0.50754,	 0.50754,  0.50754,  1.0},
	{0.2775,   0.2755,   0.2755,   1.0}};

static GLfloat rhoSpecular[nMaterials][4] = {
	{0.992157, 0.941176, 0.807843, 1.0},
	{0.393548, 0.271906, 0.166721, 1.0},
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
Polyhedron<double> tetrahedron;
Polyhedron<int> cube;


void initLights(void);
void toggleLight(int i);
void updateLights(void);


// Initialization

void createTetrahedron(bool verbose) {
  ntuple<double,3> t0(0, 0, 0);        // Each vertex is an ntuple
  ntuple<double,3> t1(1, 1, 0);
  ntuple<double,3> t2(1, 0, 1);
  ntuple<double,3> t3(0, 1, 1);
	
  list< ntuple<double,3> > face012;    // Each face is specified as a list of vertices
  list< ntuple<double,3> > face023;    // in counterclockwise order.
  list< ntuple<double,3> > face031;
  list< ntuple<double,3> > face132;
	
  face012.push_back(t0);
  face012.push_back(t1);
  face012.push_back(t2);
	
  face023.push_back(t0);
  face023.push_back(t2);
  face023.push_back(t3);
		
  face031.push_back(t0);
  face031.push_back(t3);
  face031.push_back(t1);
	
  face132.push_back(t1);
  face132.push_back(t3);
  face132.push_back(t2);
	
   index_t index;
   bool status;

   status = tetrahedron.addFace(face012, &index);
   if (! status) {
	 cerr << "Error when adding face012 to tetrahedron." << endl;
	 tetrahedron.describe(cerr);
	 exit(1);
   }

   if (verbose) {
	 cout << "Adding face012 yields status = " << status << ", index = " << index << endl;
	 tetrahedron.describe(cout);
   }

   status = tetrahedron.addFace(face023, &index);
   if (! status) {
	 cerr << "Error when adding face023 to tetrahedron." << endl;
	 tetrahedron.describe(cerr);
	 exit(1);
   }

   if (verbose) {
	 cout << "Adding face023 yields status = " << status << ", index = " << index << endl;
	 tetrahedron.describe(cout);
   }

   status = tetrahedron.addFace(face031, &index);
   if (! status) {
	 cerr << "Error when adding face031 to tetrahedron." << endl;
	 tetrahedron.describe(cerr);
	 exit(1);
   }

   if (verbose) {
	 cout << "Adding face031 yields status = " << status << ", index = " << index << endl;
	 tetrahedron.describe(cout);
   }

   status = tetrahedron.addFace(face132, &index);
   if (! status) {
	 cerr << "Error when adding face132 to tetrahedron." << endl;
	 tetrahedron.describe(cerr);
	 exit(1);
   }

   if (verbose) {
	 cout << "Adding face132 yields status = " 
			   << status << ", index = " << index << endl << endl;	
	 cout << "Output of tetrahedron.describe(cout):\n\n";
	 tetrahedron.describe(cout);
   }
 
   status = tetrahedron.checkTables(true);
   if (! status) {
	 cerr << "Consistency error detected while checking the tables for the tetrahedron." << endl;
   }

   if (verbose) {
	 cout << "Output of tetrahedron.checkTables(true) yields status = " << status << endl;
   }

   status = tetrahedron.eulerTest(true);
   if (! status) {
	 cerr << "Euler test error for tetrahedron." << endl;
   }

   if (verbose) {
	 cout << "Output of tetrahedron.eulerTest(true) yields status = " << status << endl;
   }
   return;
}


void createCube(bool verbose) {
  bool status;
  index_t index;

  ntuple<int,3> c0(0,0,0);   // Each vertex is an ntuple
  ntuple<int,3> c1(1,0,0);
  ntuple<int,3> c2(1,1,0);
  ntuple<int,3> c3(0,1,0);
  ntuple<int,3> c4(0,0,1);
  ntuple<int,3> c5(1,0,1);
  ntuple<int,3> c6(1,1,1);
  ntuple<int,3> c7(0,1,1);
	
  list< ntuple<int,3> > face0321;  // Each face is specified as a list of vertices
  list< ntuple<int,3> > face0154;
  list< ntuple<int,3> > face0473;
  list< ntuple<int,3> > face2376;
  list< ntuple<int,3> > face1265;
  list< ntuple<int,3> > face4567;
	
  face0321.push_back(c0);
  face0321.push_back(c3);
  face0321.push_back(c2);
  face0321.push_back(c1);
	
  face0154.push_back(c0);
  face0154.push_back(c1);
  face0154.push_back(c5);
  face0154.push_back(c4);

  face0473.push_back(c0);
  face0473.push_back(c4);
  face0473.push_back(c7);
  face0473.push_back(c3);
	
  face2376.push_back(c2);
  face2376.push_back(c3);
  face2376.push_back(c7);
  face2376.push_back(c6);
	
  face1265.push_back(c1);
  face1265.push_back(c2);
  face1265.push_back(c6);
  face1265.push_back(c5);
	
  face4567.push_back(c4);
  face4567.push_back(c5);
  face4567.push_back(c6);
  face4567.push_back(c7);
	
  status = cube.addFace(face0321, &index);
  if ((! status) || verbose) {
	cout << "Adding face0321 yields status = " << status << ", index = " << index << endl;
	cube.describe(cout);
  }
 
  status = cube.addFace(face0154, &index);
  if ((! status) || verbose) {
	cout << "Adding face0154 yields status = " << status << ", index = " << index << endl;
	cube.describe(cout);
  }

  status = cube.addFace(face0473, &index);
  if ((! status) || verbose) {
	cout << "Adding face0473 yields status = " << status << ", index = " << index << endl;
	cube.describe(cout);
  }

  status = cube.addFace(face2376, &index);
  if ((! status) || verbose) {
   cout << "Adding face2376 yields status = " << status << ", index = " << index << endl;
   cube.describe(cout);
  }

  status = cube.addFace(face1265, &index);
  if ((! status) || verbose) {
	cout << "Adding face1265 yields status = " << status << ", index = " << index << endl;
	cube.describe(cout);
  }

  status = cube.addFace(face4567, &index);
  if ((! status) || verbose) {
	cout << "Adding face4567 yields status = " << status << ", index = " << index << endl << endl;
	cube.describe(cout);
  }

  return;
}


void destroyCube(bool verbose) {
  bool status;
  index_t index;

  if (verbose) {
	cout << "Deleting faces 1, 2, 3, and 4 from cube ...\n";
  }

  for(index_t i = 1; i <= 4; i++) {
	status = cube.deleteFace(i);
	if (! status) {
	  cerr << "Error encountered while deleting face" << i << "." << endl;
	}
  }

  status = cube.purgeTables();
  if ((! status) || verbose) {
	cout << "cube.purgeTables() returned status = " << status << "." << endl <<  endl;  
	cout << "Output of cube.describe(cout) after the purge:\n\n";
	cube.describe(cout);
  }

  status = cube.checkTables(true);
  cout << "Output of cube.checkTables(true) yields status = " << status << endl;
  status = cube.eulerTest(true);
  cout << "Output of cube.eulerTest(true) yields status = " << status << endl;

  // Now we add some new faces.

  // Insert the point (2, 0, 1) into the vertex table
  index = cube.addVertex(ntuple<int,3>(2,0,1));
  cube.extendFace(index, 5);
  index = cube.addVertex(ntuple<int,3>(2,0,0));
  cube.extendFace(index, 2);

  // Construct a face that connects vertices 8, 9, 2, and 7.
  list<index_t> vertexList;
  vertexList.push_back(8);
  vertexList.push_back(9);
  vertexList.push_back(2);
  vertexList.push_back(7);
  index_t faceIndex;
  cube.addFace(vertexList, &faceIndex);

  cout << "Output of cube.describe(cout) after the additions:\n\n";
  cube.describe(cout);
  return;
}

void initializeModels(bool verbose) {
   // Create a tetrahedron
   cout << "Test One: Create a Tetrahedron with double precision vertices :\n\n";
   createTetrahedron(verbose);

   // Create a cube
   cout << "\n\nTest Two: Building a cube with integer valued vertex coordinates.\n\n";
   createCube(verbose);
   return;
}

//////////// OpenGL callbacks:

void display() {
  glMatrixMode(GL_PROJECTION);     /* Start modifying the projection matrix. */
  glLoadIdentity();                /* Reset project matrix. */
  gluPerspective(fieldOfViewY, aspect, near, far);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  GLdouble theta_rad = theta*radiansPerDegree;
  GLdouble phi_rad   = phi*radiansPerDegree;

  GLdouble eyeX = rho*sin(theta_rad)*cos(phi_rad);
  GLdouble eyeY = rho*sin(theta_rad)*sin(phi_rad);
  GLdouble eyeZ = rho*cos(theta_rad);

  gluLookAt(eyeX, eyeY, eyeZ,
			centerX, centerY, centerZ,
			upX, upY, upZ);

  glClearColor(0.0, 0.0, 0.0, 0.0); /* Black background */
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  //tetrahedron.display();
  cube.display();
  glFlush();  /* Single buffered, so needs a flush. */
}

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


void reshape(int w, int h)   /* Called when the user changes the dimensions of the X window. */
{
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
}


// OpenGL lighting functions:

void turnOnLight(int i) {
	glLightfv(light[i], GL_AMBIENT,  light_diffuse[i]);
	glLightfv(light[i], GL_DIFFUSE,  light_diffuse[i]);
	glLightfv(light[i], GL_SPECULAR, light_diffuse[i]);
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

void init(void) {
	glClearColor( 0.0, 0.0, 0.0, 0.0);  // Black
	glShadeModel(GL_SMOOTH);            // Gouraud Shading
	setMaterial(SILVER);                // see material.cc
	glEnable(GL_DEPTH_TEST);            // Z-buffer mode
	initLights();                       // see light.cc
}

int main (int argc, char *argv[]) {
  bool verbose = true;
  initializeModels(verbose);
 
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

