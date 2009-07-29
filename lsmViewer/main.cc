/* lsmViewer.cc
 * 
 * Created on 09 Jan 2007
 *
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <stdexcept>

#include "env.h"
#include "lsm.h"
#include "lsmVoxelImage.h"

const bool debug = true;

/* Initial dimensions and attributes of the GUI Window */
const int WINDOW_WIDTH = 512;
const int WINDOW_HEIGHT = 512;
const int X_OFFSET = 10;
const int Y_OFFSET = 10;

typedef struct rect_tag {
	float xmin;
	float xmax;
	float ymin;
	float ymax;} rect;
	
/* World coordinate extrema */
const float XMIN = 0.0;
const float XMAX = 1.0;
const float YMIN = 0.0;
const float YMAX = 1.0;
const rect defaultWorld= {XMIN, XMAX, YMIN, YMAX};

struct glutWindowInfo {
	int x;  // x-offset
	int y;  // y-offset
	int height;
	int width;
	std::string title;
	int id;
	rect world;
};
		
const int SEARCH_FAILED = -1;


//string& eraseDirectories(std::string &);
void resizeGlutImageWindow(glutWindowInfo &w, lsmVoxelImage &image, double sx, double sy);
void initializeGlutWindow(glutWindowInfo &w);

// Graphics State Information.
//ImageInfo *imageInfo;
//grayImage<GLubyte> inputImage;
uint32 imageIndex;
uint32 channelIndex;
//LSM_File *lsm;

lsmVoxelImage *lsmImage;

glutWindowInfo inputWindow = 
  {X_OFFSET, Y_OFFSET, WINDOW_HEIGHT, WINDOW_WIDTH, "Image", 0, defaultWorld};
	
using namespace std;

/* eraseDirectories is a string utility that removes the leading subirectories 
 * from the string referenced by path. Thus, the following code:
 *
 * string path("/home/users/smith/data.dat");
 * eraseDirectories(path);
 * cout << path << endl;
 *
 * should send the string "data.dat" to standard output.
 */

string& eraseDirectories(string &path) {
 	size_t last_slash = path.find_last_of("/");  // Remove leading directories
   	if (last_slash != string::npos)
   		(void) path.erase(0, last_slash + 1);
   	return path;
}

void reshape(int w, int h) {
	inputWindow.width = w;
	inputWindow.height = h;
	glViewport(0, 0, (GLsizei) w, (GLsizei) h);
   								/* using an orthographic projection */
   	glMatrixMode(GL_PROJECTION);
   	glLoadIdentity();
   	gluOrtho2D(inputWindow.world.xmin, inputWindow.world.xmax, 
				inputWindow.world.ymin, inputWindow.world.ymax);
   	glMatrixMode(GL_MODELVIEW);
   
   	glClearColor(1.0, 1.0, 1.0, 1.0); /* white background */
   	glColor3f(1.0, 0.0, 0.0);         /* draw in red */
}

void redraw() {
	resizeGlutImageWindow(inputWindow, *lsmImage, 0.5, 0.5);
	glutReshapeWindow(inputWindow.width, inputWindow.height);
	glutPostRedisplay();	
}

void keyboard(unsigned char c, int x, int y) {
	int candidate;
	
	switch(c) {
		case '0':
			if (channelIndex != 0) {
				channelIndex = 0;
				redraw();
			}
			break;
		case '1':
			if (channelIndex != 1 && lsmImage->bands() > 1) {
				channelIndex = 1;
				redraw();
			}
			break;	
		case 'q':
			exit(0);
			break;
		case 'n':
		    candidate = imageIndex + 1;
			if (candidate < lsmImage->layers()) {
				imageIndex = candidate;
				redraw();
		    }
		
			break;	
		case 'p':
			candidate = imageIndex - 1;
			if (candidate >= 0) {
				imageIndex = candidate;
				redraw();
			}
			break;
		default:
		    break;
	}			
}

void display() {
	glClear(GL_COLOR_BUFFER_BIT);  // clear the window (paint it with the background color)
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

	//glRasterPos2i(256, 256);
	glPixelZoom(0.5, -0.5);
	glRasterPos2i(inputWindow.world.xmin, inputWindow.world.ymax);
	lsmImage->drawPixels(imageIndex, channelIndex);
	glFlush();
}

void resizeGlutImageWindow(	glutWindowInfo &w, 
							lsmVoxelImage &image, 
							double sx, double sy) {
	w.height = image.rows()*sy;
	w.width  = image.cols()*sx;
	w.world.xmin = 0;
	w.world.xmax = image.cols();
	w.world.ymin = 0;
	w.world.ymax = image.rows();
}

void initializeGlutWindow(glutWindowInfo &w) {
	if (debug) {
		cout << "initializeGlutWindow()" << endl;
	}
	glutInitWindowSize(w.width, w.height); 
    glutInitWindowPosition(w.x, w.y); 
    w.id = glutCreateWindow(w.title.c_str()); 
    glutReshapeFunc(reshape); // register reshape as a callback function 
    glutDisplayFunc(display); // register display as a callback function
    glutKeyboardFunc(keyboard); // register keyboard as a callback function
}


int main (int argc, char **argv) {
	if (debug) {
   		cerr << "argc = " << argc << endl;
	}

	// Initialize GLUT:
	glutInit(&argc,argv);
	glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB); /* default, not needed */
 	
   	if (argc < 2) {
   		string command = argv[0];                       // Identify the command name
   		(void) eraseDirectories(command);  // Remove leading directories (if present).
   		 
   		cerr << "Usage: " << command << " <filename>" << endl;
   		return EXIT_FAILURE;
   	}
   
	lsmImage = new lsmVoxelImage(argv[1]);


	imageIndex   = 0;
	channelIndex = 0;

	resizeGlutImageWindow(inputWindow, *lsmImage, 0.5, 0.5);
	initializeGlutWindow(inputWindow);
	glutMainLoop();
    return 0;
}