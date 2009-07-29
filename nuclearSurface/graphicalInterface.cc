/*
 * graphicalInterface.cc
 * Project: nuclearSurface
 *
 * Created by Robert R. Snapp on 2008-08-24
 *
 * This file contains routines to view a specified layer of an lsm
 * data file using OpenGL, and the outcomes of various image processing
 * and analysis procedures.
 *
 * 
 */

#include "env.h" // Contains the paths to the glut.h, etc.
#include "lsm.h"
//#include "pixelArray.h"
#include "lsmVoxelImage.h"
#include "nuclearSurface.h"
#include "imageComponentArray.h"
#include "grayImage3d.h"
#include "binaryImage3d.h"


#include <iostream>
#include <iomanip>
#include <vector>
#include <stdexcept>


/* Initial dimensions and attributes of the GUI Window */
const int MARGIN = 5;
const int WINDOW_WIDTH  = 840;
const int WINDOW_HEIGHT = 840;
const int X0 = 10;
const int Y0 = 10;

template<typename T>
struct rect {
  T xmin;
  T xmax;
  T ymin;
  T ymax;};

/* World coordinate extrema */
const float XMIN = 0.0;
const float XMAX = 1.0;
const float YMIN = 0.0;
const float YMAX = 1.0;
const rect<float> defaultWorld= {XMIN, XMAX, YMIN, YMAX};


// Constants used for screen fonts
void* const glutFont = GLUT_BITMAP_8_BY_13;
int const glutFontWidth = 8;
int const glutFontLineSkip = 18;
int const BUFLEN = 128;

// Graphics state information.
enum View {LSM_SLICE, THRESHOLD, DIALATED, CLOSED, LABELS};
View currentView;
enum RegionSelectMode {DISABLED, ENABLED, ONE_CORNER_SELECTED, RECTANGLE_SELECTED};
RegionSelectMode regionSelectMode;
Rect selectionRegion;

ulong layer;

bool infoVisibility;
bool showConvexHull;
bool invVideo = false;
float hullVolume;
float clusterVolume;

struct glutWindowInfo {
  int x;  // x-offset
  int y;  // y-offset
  int height;
  int width;
  std::string title;
  int id;
  rect<float> world;
};
		
const int SEARCH_FAILED = -1;

// private functions: only should be used within this file. (Declarations of other
// functions appear in project.h)
void initializeGlutWindow(glutWindowInfo &w);
void resizeMainWindow();
void resizeWindow(glutWindowInfo &w, Rect &region, double sx, double sy);
void reshape(int w, int h);
void keyboard(unsigned char c, int x, int y);
void imageKeyboard(unsigned char c, int x, int y);
void display();
void specialKeyboard(int key, int x, int y);
void cropImage();

/* The mainWindow displays a specified layer of a specified channel of an lsm file.
 */

glutWindowInfo mainWindow =
   {X0, Y0, WINDOW_HEIGHT, WINDOW_WIDTH, "LSM Image Layer", -1, defaultWorld};

using namespace std;

//*******************************
//*  PRIVATE FUNCTIIONS
//*******************************

void cropImage() {
  cout << "Changing to subregion" 
       << subvol.col << ", "
       << subvol.row << ", "
       << subvol.width << ", "
       << subvol.height
       << endl;

  loadAndPrepareLsmFile(fileName.c_str(), subvol);
  comp.clear();
  currentView = LSM_SLICE;
  resizeMainWindow();
  return;
 }

int drawValue(string &v) {
  size_t i;

  for(i = 0; i < v.size(); i++) {
    glutBitmapCharacter(glutFont, v[i]);
  }
  return glutFontWidth*i;
}

int drawValue(char* buf, size_t n) {
  size_t i;
  for(i = 0; buf[i] !='\0' && i < n; i++) {
    glutBitmapCharacter(glutFont, buf[i]);
  }
  return glutFontWidth*i;
}

int drawValue(char* const buf) {
  int i;
  for(i = 0; buf[i] !='\0'; i++) {
    glutBitmapCharacter(glutFont, buf[i]);
  }
  return glutFontWidth*i;
}

int drawValue(int v) {
  char buffer[BUFLEN];

  sprintf(buffer, "%d", v);
  return drawValue(buffer, BUFLEN);
}

int drawValue(ulong v) {
  char buffer[BUFLEN];

  sprintf(buffer, "%ld", v);
  return drawValue(buffer, BUFLEN);
}

int drawValue(float v) {
  char buffer[BUFLEN];

  sprintf(buffer, "%f", v);
  return drawValue(buffer, BUFLEN);
}  

int drawValue(double v) {
  char buffer[BUFLEN];

  sprintf(buffer, "%lf", v);
  return drawValue(buffer, BUFLEN);
}  

void updateProcessedImage() {
  switch (currentView) {
	case LSM_SLICE:  // pass through
	case THRESHOLD:  // pass through
	  break;
	case DIALATED:
	case CLOSED:
	  closureOp();
	  break;
	case LABELS:
	  closureOp();
	  comp.clear();
	  comp.computeComponents(processedImage[1]);
	  int nLabels = comp.size();
	  comp.makeColorLabels(nLabels);
	  break;
	}
  return;
}
	  
void mouseButton(int button, int state, int x, int y) {
  switch(regionSelectMode) {
  case ENABLED:
  case RECTANGLE_SELECTED:
    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
      regionSelectMode = ONE_CORNER_SELECTED;
      selectionRegion.col = x;
      selectionRegion.row = y;
    }
    break;
  case ONE_CORNER_SELECTED:
    if (button == GLUT_LEFT_BUTTON && state == GLUT_UP) {
      selectionRegion.width = max(x, selectionRegion.col) - min(x, selectionRegion.col);
      selectionRegion.col   = min(x, selectionRegion.col);

      int yp = y;
      selectionRegion.height = max(yp, selectionRegion.row) - min(yp, selectionRegion.row);
      selectionRegion.row    = min(yp, selectionRegion.row);

      if (selectionRegion.width > 0 && selectionRegion.height > 0) {
        // scale selection window to world coordinates
        double dx = (mainWindow.world.xmax - mainWindow.world.xmin)/static_cast<double>(mainWindow.width);
        double dy = (mainWindow.world.ymax - mainWindow.world.ymin)/static_cast<double>(mainWindow.height);
        selectionRegion.col    = static_cast<int>(mainWindow.world.xmin + dx*selectionRegion.col);
        selectionRegion.width  = static_cast<int>(dx*selectionRegion.width);
        selectionRegion.row    = static_cast<int>(mainWindow.world.ymin + dy*selectionRegion.row);
        selectionRegion.height = static_cast<int>(dy*selectionRegion.height);
        
        regionSelectMode = RECTANGLE_SELECTED;
      } else {
        regionSelectMode = DISABLED;
      }
    }
    break;
  default:
    break;
  }
  glutPostRedisplay();
}

// Procedure reshape is an OpenGL callback.
void reshape(int w, int h) {
  int currentWindowID = glutGetWindow();
  glutWindowInfo* wind;

  if (currentWindowID == mainWindow.id) {
    wind = &mainWindow;
    double sx = static_cast<double>(w)/static_cast<double>(subvol.width);
    double sy = static_cast<double>(h)/static_cast<double>(subvol.height);
    // Ensure that the aspects of the window and image agree.
    if (sx < sy) {
		h = static_cast<int>(ceil(sx*subvol.height)); // Reduce the window height to preserve the aspect ratio.
    } else {
		w = static_cast<int>(ceil(sy*subvol.width)); // Reduce the window width to preserve the aspect ratio.
    }
    glutReshapeWindow(w, h);  // Adjust the window geometry.

  } else {
    cerr << "reshape(...) glut callback (file " << __FILE__ << ", line " << __LINE__ << ")"
         << "FATAL ERROR: unknown currentWindowID = " << currentWindowID
         << endl;
    abort();
  }
  
  // 
  wind->width  = w;
  wind->height = h;

  // Draw to the entire window.
  glViewport(0, 0, (GLsizei) w, (GLsizei) h);

  // Use an orthographic projection that maps the extrema of the wind->world data structure to the
  // viewport boundaries.
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D((wind->world).xmin, (wind->world).xmax, (wind->world).ymin, (wind->world).ymax);
  glMatrixMode(GL_MODELVIEW);
   
  glClearColor(1.0, 1.0, 1.0, 1.0); /* white background */
  glColor3f(1.0, 0.0, 0.0);         /* draw in red */
}

void keyboard(unsigned char c, int x, int y) {
  int currentWindowID = glutGetWindow();
  if (currentWindowID == mainWindow.id) {
    imageKeyboard(c, x, y);
  } else {
    cerr << "keyboard(...) glut callback (file " << __FILE__ << ", line " << __LINE__ << ")"
         << "FATAL ERROR: unknown currentWindowID = " << currentWindowID
         << endl;
    abort();
  }
}

void specialKeyboard(int key, int x, int y) {
#pragma unused(x)
#pragma unused(y)
  switch(key) {
  case GLUT_KEY_F1:
    if (currentView != LSM_SLICE) {
	  currentView = LSM_SLICE;
      glutPostRedisplay();
    }
    break;
  case GLUT_KEY_F2:
    if (currentView != THRESHOLD) {
      currentView = THRESHOLD;
      glutPostRedisplay();
    }
    break;
  case GLUT_KEY_F3:
    if (currentView != DIALATED) {
      currentView = DIALATED;
      glutPostRedisplay();
    }
    break;
  case GLUT_KEY_F4:
    if (currentView != CLOSED) {
      currentView = CLOSED;
      glutPostRedisplay();
    }
    break;
  case GLUT_KEY_F5:
    if (currentView != LABELS) {
      currentView = LABELS;
      glutPostRedisplay();
    }
    break;
  case GLUT_KEY_LEFT: // previous band
	regionSelectMode = DISABLED;
    if (band > 0) {
      band--;
	  applyThreshold(lsmData, band, subvol, threshold, afterThreshold);
	  updateProcessedImage();
    } else {
	  cerr << '\007' << "Command ignored." << endl;
	}
	glutPostRedisplay();
    break;
  case GLUT_KEY_RIGHT:    // next band.
    regionSelectMode = DISABLED;
    if (band < static_cast<int>(lsmData.bands()) - 1) {
      band++;
	  applyThreshold(lsmData, band, subvol, threshold, afterThreshold);
	  updateProcessedImage();
    } else {
	  cerr << '\007' << "Command ignored." << endl;
	}
	glutPostRedisplay();
    break;
  case GLUT_KEY_UP: // increase threshold: does not require recomputing the projection.
    if (threshold < 255) {
      threshold++;
	  applyThreshold(lsmData, band, subvol, threshold, afterThreshold);
	  updateProcessedImage();
      glutPostRedisplay();
    }
    break;
  case GLUT_KEY_DOWN: // decrease threshold: does not require recomputing the projection.
    if (threshold > 0) {
      threshold--;
	  applyThreshold(lsmData, band, subvol, threshold, afterThreshold);
	  updateProcessedImage();
      glutPostRedisplay();
    }
    break;
  case GLUT_KEY_PAGE_UP: // increase threshold by 5: does not require recomputing the projection.
    if (threshold < 255) {
      threshold = min(255, threshold + 5);
	  applyThreshold(lsmData, band, subvol, threshold, afterThreshold);
	  updateProcessedImage();
      glutPostRedisplay();
    }
    break;
  case GLUT_KEY_PAGE_DOWN: // decrease threshold by 5
    if (threshold > 0) {
      threshold = max(0, threshold - 5);
	  applyThreshold(lsmData, band, subvol, threshold, afterThreshold);
	  updateProcessedImage();
      glutPostRedisplay();
    }
    break;
  default:
    break;
  }
}

// Keyboard commands
void imageKeyboard(unsigned char c, int x, int y) {
#pragma unused(x)
#pragma unused(y)
  switch(c) {
  case 27: // ESC
    regionSelectMode = DISABLED;
    break;
  case '0':  // pass through
  case '1':  // pass through
  case '2':  // pass through
  case '3':  // pass through
  case '4':  // pass through
  case '5':  // pass through
  case '6':  // pass through
  case '7':  // pass through
  case '8':  // pass through
  case '9':
    regionSelectMode = DISABLED;
	ulong requestedBand = c - '0';
	if (band != requestedBand && requestedBand < lsmData.bands()) {
	  band = requestedBand;
	  applyThreshold(lsmData, band, subvol, threshold, afterThreshold);
	  updateProcessedImage();
    } else {
	  if (requestedBand >= lsmData.bands()) {
		cerr << '\007' << "Command ignored." << endl;
	  }
	}
    break;

  case 'i':
    infoVisibility = ! infoVisibility;
    break;

  case 'I':
	invVideo = ! invVideo;
	cerr << "Setting inverseVideo mode to " << invVideo << endl;
	break;

  case 'c':
    regionSelectMode = DISABLED;
	comp.clear();
	comp.computeComponents(processedImage[1]);

    int nLabels = comp.size();
    comp.makeColorLabels(nLabels);
    currentView = LABELS;
    showConvexHull = false;
    break;
  case 'h':
    regionSelectMode = DISABLED;
    if (comp.rows() > 0 && comp.cols() > 0 && comp.layers() > 0) {
	  if (verbosity > 1) {
		cerr << "Largest component has " << comp.clusterSize(0) << " voxels." << endl;
	  }
 
	  comp.computeConvexHull(0, ch, VOXEL_CENTERS);
      hullVolume = ch.volume();
	  clusterVolume = comp.clusterSize(0);
      showConvexHull = true;
    } else {
	  cerr << '\007' << "Command ignored." << endl;
	}
    break;
  case 'H':
    regionSelectMode = DISABLED;
    if (comp.rows() > 0 && comp.cols() > 0 && comp.layers() > 0) {
	  if (verbosity > 1) {
		cerr << "Largest component has " << comp.clusterSize(0) << " voxels." << endl;
	  }
 
	  comp.computeConvexHull(0, ch, VOXEL_CORNERS);
      hullVolume = ch.volume();
	  clusterVolume = comp.clusterSize(0);
      showConvexHull = true;
    } else {
	  cerr << '\007' << "Command ignored." << endl;
	}
    break;
  case 'm':
    regionSelectMode = DISABLED;
    currentView = CLOSED;
	updateProcessedImage();
	break;
  case 'q':
    exit(0);
    break;

  case 'r':
    if (regionSelectMode == DISABLED) {
      regionSelectMode = ENABLED;
    } else {
      regionSelectMode = DISABLED;
    }
    break;

  case 'x':
	showConvexHull = false;
    if (regionSelectMode == RECTANGLE_SELECTED) {
	  // The user requests to reduce the focus to the indicated subvol relative to the window origin.
	  subvol.col += selectionRegion.col;
	  subvol.row += selectionRegion.row;
	  subvol.width = selectionRegion.width;
	  subvol.height = selectionRegion.height;
      cropImage();
      regionSelectMode = DISABLED;
    } else if (regionSelectMode == ENABLED) {
	  cerr << "Restoring subvol to entire image." << endl;
      subvol.col = 0;
      subvol.row = 0;
	  subvol.layer = 0;
      subvol.width = lsmData.cols();
      subvol.height = lsmData.rows();
	  subvol.depth = lsmData.layers();
      cropImage();
      regionSelectMode = DISABLED;
    }
    break;

  case 'n':
    regionSelectMode = DISABLED;
    if (layer < min(lsmData.layers(), static_cast<uint32>(subvol.layer + subvol.depth)) - 1 ) {
      layer++;
    } else {
      cerr << '\007' << "Command ignored." << endl;
    }
		
    break;	
  case 'p':
    regionSelectMode = DISABLED;
	  if (layer > max(0, subvol.layer)) {
      layer--;
    } else {
      cerr << '\007' << "Command ignored." << endl;
    }
    break;

  case 'R': // Reset the subvolume to the full image
    regionSelectMode = DISABLED;
    showConvexHull = false;
    subvol.row = 0;
    subvol.col = 0;
	subvol.layer = 0;
    subvol.width = lsmData.cols();
    subvol.height = lsmData.rows();
	subvol.depth = lsmData.layers();
	applyThreshold(lsmData, band, subvol, threshold, afterThreshold);
	updateProcessedImage();
    resizeMainWindow();
    break;

  case 'a':  // advance fileRosterPos to a specified offset value.
    regionSelectMode = DISABLED;
    if (fileRoster.size() > 0) {
      // Prompt user for roster offset
      int offset;

      do {
        cout << "File roster offset? ";
        cin >> offset;
        if (offset < 0 || offset >= static_cast<int>(fileRoster.size())) {
          cerr << "\007";
          cout << "File roster offset should be within [0, " << fileRoster.size() << "]."
               << endl;
        }
      } while (offset < 0 || offset >= static_cast<int>(fileRoster.size()));

      // advance fileRosterPos to the requested offset value.
      fileRosterPos = fileRoster.begin(); 
      for(int i = 0; i < offset && fileRosterPos != fileRoster.end(); ++i) {
        ++fileRosterPos;
      }

      // load and display
	  subvol.row = 0;
      subvol.col = 0;
	  subvol.layer = 0;
      subvol.width = 0;
      subvol.height = 0;
	  subvol.depth = 0;

      if (parseFileLine(*fileRosterPos)) {
		loadAndPrepareLsmFile(fileName.c_str(), subvol);
		layer = subvol.layer;
		comp.clear();
		currentView = LSM_SLICE;
		showConvexHull = false;
		resizeMainWindow();
	  }
    }
    break;
  case '>':
    // advance the file index
	if (fileRoster.size() > 0) { // Ignore if there is no file roster.
	  regionSelectMode = DISABLED;
	  showConvexHull = false;

	  // seek next nonempty record in the fileroster
	  bool status = false;
	  subvol.row = 0;
	  subvol.col = 0;
	  subvol.layer = 0;
	  subvol.width = 0;
	  subvol.height = 0;
	  subvol.depth = 0;

	  while(! status && fileRosterPos < fileRoster.end() - 1) {
		++fileRosterPos;
		
		if (fileRosterPos < fileRoster.end()) {
		  status = parseFileLine(*fileRosterPos);  // returns true if record is nonempty
		}
	  }

	  if (status) {
		loadAndPrepareLsmFile(fileName.c_str(), subvol);
		layer = subvol.layer;
		comp.clear();
		currentView = LSM_SLICE;
		resizeMainWindow();
	  } else {
		cerr << "End of fileRoster: Command ignored." << endl;
	  }
    } else {
      cerr << "No fileRoster: Command ignored." << endl;
    }
    break;

  case '<':
    // reverse the file index
	if (fileRoster.size() > 0) { // Ignore if there is no file roster.
	  regionSelectMode = DISABLED;
	  showConvexHull = false;

	  // seek previous nonempty record in the fileroster
	  bool status = false;
	  subvol.row = 0;
	  subvol.col = 0;
	  subvol.layer = 0;
	  subvol.width = 0;
	  subvol.height = 0;
	  subvol.depth = 0;

	  while (!status && fileRosterPos > fileRoster.begin()) {
		--fileRosterPos;
		if (fileRosterPos >= fileRoster.begin()) {
		  status = parseFileLine(*fileRosterPos);
		}
	  }
		
	  if (status) {
		loadAndPrepareLsmFile(fileName.c_str(), subvol);
		layer = subvol.layer;
		comp.clear();
		currentView = LSM_SLICE;
		showConvexHull = false;
		resizeMainWindow();
	  } else {
		cerr << "Beginning of fileRoster: Command ignored." << endl;
	  }
	} else {
	  cerr << "No fileRoseter: command ignored." << endl;
	}
    break;

  case 'F':
    // explicit filename entry.
   
    bool status;
    regionSelectMode = DISABLED;
    subvol.row = 0;
    subvol.col = 0;
	subvol.layer = 0;
    subvol.width = 0;
    subvol.height = 0;
	subvol.depth = 0;
    do {
       string f;
       cout << "Filename? ";
       cin >> f;
       fileName = f;
       status = loadAndPrepareLsmFile(f.c_str(), subvol);
       if (! status) {
         cerr << "Could not parse file " << f << endl;
       }
    } while (status == false);
	layer = subvol.layer;
    comp.clear();
    currentView = LSM_SLICE;
    showConvexHull = false;
    resizeMainWindow();
    break;

  default:
    break;
  }
  glutPostRedisplay();
}

void display() {
  int currentWindowID = glutGetWindow();
  glClear(GL_COLOR_BUFFER_BIT); 

  if (currentWindowID == mainWindow.id) {
    glRasterPos2f(mainWindow.world.xmin, mainWindow.world.ymax);

    double sx = static_cast<double>(mainWindow.width)/static_cast<double>(subvol.width);
    double sy = static_cast<double>(mainWindow.height)/static_cast<double>(subvol.height);
    double sf = min(sx, sy);

    glPixelZoom(sf, -sf);

    switch(currentView) {
    case LSM_SLICE:
      lsmData.drawPixels(layer, band, subvol.width, subvol.height, subvol.col, subvol.row, invVideo);
      break;
    case THRESHOLD:
	  if (afterThreshold.rows() > 0 && afterThreshold.cols() > 0)
		afterThreshold.drawPixels(layer - subvol.layer, invVideo);
      break;
    case DIALATED:
	  if (processedImage[0].rows() > 0 && processedImage[0].cols() > 0)
		processedImage[0].drawPixels(layer - subvol.layer, invVideo);
      break;
    case CLOSED:
	  if (processedImage[1].rows() > 0 && processedImage[1].cols() > 0)
		processedImage[1].drawPixels(layer - subvol.layer, invVideo);
      break;
    case LABELS:
	  if (comp.rows() > 0 && comp.cols() > 0 && comp.layers() > 0)
		comp.displayLayer(layer - subvol.layer);
      break;
    }

     if (infoVisibility) {
      glColor3f(0.1, 0.6, 1.0);
      int x = static_cast<int>(10.0f/sf);
      int y = static_cast<int>(mainWindow.world.ymax - 24.0f/sf);
      glRasterPos2i(x, y);
      drawValue(fileName);
      y -= static_cast<int>(glutFontLineSkip/sf);
      glRasterPos2i(x, y);
      drawValue("band = ");
      drawValue(band);
	  drawValue(", layer = ");
	  drawValue(layer);
      
      y -= static_cast<int>(glutFontLineSkip/sf);
      glRasterPos2i(x, y);   
      drawValue("threshold = ");
      drawValue(threshold);

      if (showConvexHull) {
        glRasterPos2i(static_cast<GLint>(10/sf), static_cast<GLint>((10 + 2*glutFontLineSkip)/sf));
        drawValue("Hull volume ratio = ");
        drawValue(static_cast<int>(hullVolume));
        drawValue("/");
        drawValue(static_cast<int>(clusterVolume));
        drawValue(" = ");
        drawValue(static_cast<float>(hullVolume)/static_cast<float>(clusterVolume));
	  }
    }

    switch(regionSelectMode) {
    case ENABLED:
      glRasterPos2i(static_cast<GLint>(10.0f/sf), static_cast<GLint>(10.0f/sf));
      drawValue("Select first corner of crop rectangle with a mouse click.");
      break;
    case ONE_CORNER_SELECTED:
      glRasterPos2i(static_cast<GLint>(10.0f/sf), static_cast<GLint>((10 + glutFontLineSkip)/sf));
      drawValue("Selected (");
      drawValue(selectionRegion.col);
      drawValue(", ");
      drawValue(selectionRegion.row);
      drawValue(")");
      glRasterPos2i(static_cast<GLint>(10.0f/sf), static_cast<GLint>(10.0f/sf));
      drawValue("Select first corner with mouse, then drag and release.");
      break;
    case RECTANGLE_SELECTED:
      glPushMatrix();
      glTranslatef(0.0f,  mainWindow.world.ymax, 0.0f);
      glScalef(1.0f, -1.0f, 1.0f);
      glBegin(GL_LINE_LOOP);
        glVertex2i(selectionRegion.col, selectionRegion.row);
        glVertex2i(selectionRegion.col + selectionRegion.width, selectionRegion.row);
        glVertex2i(selectionRegion.col + selectionRegion.width, selectionRegion.row + selectionRegion.height) ;
        glVertex2i(selectionRegion.col, selectionRegion.row + selectionRegion.height);
      glEnd();
      glPopMatrix();
      glRasterPos2i(static_cast<GLint>(10.0f/sf), static_cast<GLint>(10.0f/sf));
      drawValue("Press x to crop, or ESC to cancel.");
      break;
	case DISABLED:
	  break;
    }
  } else { 
    cerr << "reshape(...) glut callback (file " << __FILE__ << ", line " << __LINE__ << ")"
         << "FATAL ERROR: unknown currentWindowID = " << currentWindowID
         << endl;
    abort();
  }

  //  glFlush();
  glutSwapBuffers();
}


void resizeWindow(glutWindowInfo &w, 
                  Subvolume &sv, 
                  double sx, double sy) {
  w.height = static_cast<int>(ceil(sv.height*sy));
  w.width  = static_cast<int>(ceil(sv.width*sx));
  w.world.xmin = 0;
  w.world.xmax = sv.width;
  w.world.ymin = 0;
  w.world.ymax = sv.height;
}

void initializeGlutWindow(glutWindowInfo &w) {
	glutInitWindowSize(w.width, w.height); 
    glutInitWindowPosition(w.x, w.y); 
    w.id = glutCreateWindow(w.title.c_str()); 

    glutReshapeFunc(reshape);   // register reshape as a callback function 
    glutDisplayFunc(display);   // register display as a callback function
    glutKeyboardFunc(keyboard); // register keyboard as a callback function
    glutMouseFunc(mouseButton);
    glutSpecialFunc(specialKeyboard);
}

void initColorTable() {
  GLubyte colorTable[256][3];
  for(int i = 0; i < 256; i++) {
    colorTable[i][0] = 255 - i;
    colorTable[i][1] = 255 - i;
    colorTable[i][2] = 0;
  }

  GLint width;
  glColorTable(GL_PROXY_COLOR_TABLE, GL_RGB, 256, GL_RGB, GL_UNSIGNED_BYTE, NULL);
  glGetColorTableParameteriv(GL_PROXY_COLOR_TABLE, GL_COLOR_TABLE_WIDTH, &width);
  if (width == 0) {
    cerr << "Can't create color table." << endl;
  } else {
    glColorTable(GL_COLOR_TABLE, GL_RGB, 256, GL_RGB, GL_UNSIGNED_BYTE, colorTable);
    glEnable(GL_COLOR_TABLE);
  }

}
    
void initializeGlut(int* argc, char** argv) {
	// Initialize GLUT:
	glutInit(argc,argv);
	glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB); /* default, not needed */
         // initColorTable();
}

void resizeMainWindow() {
    // The dimensions of the main window should conform to the subvolume.
  float const screenWidth  = 1280;
  float const screenHeight =  854 - 50;  // dont include the menu bar.
  float sx = min(1.0f, screenWidth/subvol.width);
  float sy = min(1.0f, screenHeight/subvol.height);
  float sf = min(sx, sy);

  resizeWindow(mainWindow, subvol, sf, sf);

  if (mainWindow.id > 0) {
    int win = glutGetWindow();
    glutSetWindow(mainWindow.id);
    glutReshapeWindow(mainWindow.width, mainWindow.height);  // Adjust the window geometry.
    glutSetWindow(win);
  }
}


void runGlutLoop() {
  resizeMainWindow();
  initializeGlutWindow(mainWindow);
  currentView = LSM_SLICE;
  infoVisibility = true;
  showConvexHull = false;
  layer = subvol.layer;
  glutMainLoop();
}
  
