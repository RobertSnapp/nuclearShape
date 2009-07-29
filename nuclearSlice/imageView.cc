/*
 * imageView.cc
 * Project: nuclearSlice
 *
 * Created by Robert R. Snapp on 2007-07-16
 *
 * This file contains routines to view a specified layer of an lsm
 * data file using OpenGL, and the outcomes of various image processing
 * and analysis procedures.
 *
 */

#include "env.h" // Contains the paths to the glut.h, etc.
#include "lsm.h"
#include "pixelArray.h"
#include "lsmVoxelImage.h"
#include "project.h"

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
void** const glutFont = GLUT_BITMAP_8_BY_13;
int const glutFontWidth = 8;
int const glutFontLineSkip = 18;
int const BUFLEN = 128;

// Graphics state information.
enum View {LSM_SLICE, PROJECTION, THRESHOLD, DIALATED, CLOSED, LABELS};
View currentView;
enum RegionSelectMode {DISABLED, ENABLED, ONE_CORNER_SELECTED, RECTANGLE_SELECTED};
RegionSelectMode regionSelectMode;
Rect selectionRegion;

bool infoVisibility;
bool showConvexHull;

float hullArea;

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
void resizeMainWindow(glutWindowInfo &w, Rect &region, double sx, double sy);
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

void cropImage() { }

int drawValue(string &v) {
  int i;

  for(i = 0; i < v.size(); i++) {
    glutBitmapCharacter(glutFont, v[i]);
  }
  return glutFontWidth*i;
}

int drawValue(char* buf, int n) {
  int i;
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
  int i;
  char buffer[BUFLEN];

  sprintf(buffer, "%d", v);
  return drawValue(buffer, BUFLEN);
}

int drawValue(float v) {
  int i;
  char buffer[BUFLEN];

  sprintf(buffer, "%f", v);
  return drawValue(buffer, BUFLEN);
}  

int drawValue(double v) {
  int i;
  char buffer[BUFLEN];

  sprintf(buffer, "%lf", v);
  return drawValue(buffer, BUFLEN);
}  

void mouseButton(int button, int state, int x, int y) {
  switch(regionSelectMode) {
  case ENABLED:
  case RECTANGLE_SELECTED:
    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
      regionSelectMode = ONE_CORNER_SELECTED;
      selectionRegion.col = x;
      selectionRegion.row = mainWindow.height - y;
    }
    break;
  case ONE_CORNER_SELECTED:
    if (button == GLUT_LEFT_BUTTON && state == GLUT_UP) {
      selectionRegion.width = max(x, selectionRegion.col) - min(x, selectionRegion.col);
      selectionRegion.col   = min(x, selectionRegion.col);

      int yp = mainWindow.height - y;
      selectionRegion.height = max(yp, selectionRegion.row) - min(yp, selectionRegion.row);
      selectionRegion.row    = min(yp, selectionRegion.row);

      if (selectionRegion.width > 0 && selectionRegion.height > 0) {
        // scale selection window to world coordinates
        float dx = (mainWindow.world.xmax - mainWindow.world.xmin)/mainWindow.width;
        float dy = (mainWindow.world.ymax - mainWindow.world.ymin)/mainWindow.height;
        selectionRegion.col    = mainWindow.world.xmin + dx*selectionRegion.col;
        selectionRegion.width  = dx*selectionRegion.width;
        selectionRegion.row    = mainWindow.world.ymin + dy*selectionRegion.row;
        selectionRegion.height = dy*selectionRegion.height;
        
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
    float sx = static_cast<float>(w)/static_cast<float>(subregion.width);
    float sy = static_cast<float>(h)/static_cast<float>(subregion.height);
    // Ensure that the aspects of the window and image agree.
    if (sx < sy) {
      h = sx*subregion.height;
    } else {
      w = sy*subregion.width;
    }
    glutReshapeWindow(w, h);

  } else {
    cerr << "reshape(...) glut callback (file " << __FILE__ << ", line " << __LINE__ << ")"
         << "FATAL ERROR: unknown currentWindowID = " << currentWindowID
         << endl;
    abort();
  }
    
  wind->width  = w;
  wind->height = h;

  glViewport(0, 0, (GLsizei) w, (GLsizei) h);
  /* using an orthographic projection */
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
    cerr << "reshape(...) glut callback (file " << __FILE__ << ", line " << __LINE__ << ")"
         << "FATAL ERROR: unknown currentWindowID = " << currentWindowID
         << endl;
    abort();
  }
}

void specialKeyboard(int key, int x, int y) {
  switch(key) {
  case GLUT_KEY_F1:
    if (currentView != LSM_SLICE) {
      if (projectionMode) {
        currentView = PROJECTION;
      } else {
        currentView = LSM_SLICE;
      }
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
    if (band > 0) {
      band--;
      if (projectionMode) {
        applyThresholdToScaledImage(threshold);
      } else {
        applyThresholdToLsmData(threshold, subregion);
      }
      glutPostRedisplay();
    }
    break;
  case GLUT_KEY_RIGHT:    // next band.
    if (band < lsmData.bands() - 1) {
      band++;
      if (projectionMode) {
        applyThresholdToScaledImage(threshold);
      } else {
        applyThresholdToLsmData(threshold, subregion);
      }
      glutPostRedisplay();
    }
    break;
  case GLUT_KEY_UP: // increase threshold
    if (threshold < 255) {
      threshold++;
      if (projectionMode) {
        applyThresholdToScaledImage(threshold);
      } else {
        applyThresholdToLsmData(threshold, subregion);
      }
      glutPostRedisplay();
    }
    break;
  case GLUT_KEY_DOWN: // decrease threshold
    if (threshold > 0) {
      threshold--;
      if (projectionMode) {
        applyThresholdToScaledImage(threshold);
      } else {
        applyThresholdToLsmData(threshold, subregion);
      }
      glutPostRedisplay();
    }
    break;
  case GLUT_KEY_PAGE_UP: // increase threshold by 5
    if (threshold < 255) {
      threshold = min(255, threshold + 5);
      if (projectionMode) {
        applyThresholdToScaledImage(threshold);
      } else {
        applyThresholdToLsmData(threshold, subregion);
      }
      glutPostRedisplay();
    }
    break;
  case GLUT_KEY_PAGE_DOWN: // decrease threshold by 5
    if (threshold > 0) {
      threshold = max(0, threshold - 5);
     if (projectionMode) {
        applyThresholdToScaledImage(threshold);
      } else {
        applyThresholdToLsmData(threshold, subregion);
      }
      glutPostRedisplay();
    }
    break;
  default:
    break;
  }
}

void imageKeyboard(unsigned char c, int x, int y) {
  switch(c) {
  case 27: // ESC
    regionSelectMode = DISABLED;
    break;
  case '0':
    if (band != 0) {
      band = 0;
      if (projectionMode) {
        applyThresholdToScaledImage(threshold);
      } else {
        applyThresholdToLsmData(threshold, subregion);
      }
    }
    break;
  case '1':
    if (band != 1 && lsmData.bands() > 1) {
      band = 1;
      if (projectionMode) {
        applyThresholdToScaledImage(threshold);
      } else {
        applyThresholdToLsmData(threshold, subregion);
      }
    }
    break;	
  case 'i':
    infoVisibility = ! infoVisibility;
    break;

  case 'c':
    int nLabels = collectComponents();
    makeColorLabels(nLabels);
    currentView = LABELS;
    showConvexHull = false;
    break;
  case 'h':
    if (colorLabels.rows() > 0) {
      computeComponentHull(1, ch);
      hullArea = ch.area();
      showConvexHull = true;
    }
    break;
  case 'm':
    closureOp();
    currentView = CLOSED;
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
    if (regionSelectMode == RECTANGLE_SELECTED) {
      cropImage();
      regionSelectMode = DISABLED;
    }
    break;

  case 'n':
    if (layer < lsmData.layers() - 1) {
      layer++;
      if (projectionMode) {
        applyThresholdToScaledImage(threshold);
      } else {
        applyThresholdToLsmData(threshold, subregion);
      }
    } else {
      cerr << '\007';
    }
		
    break;	
  case 'p':
    if (layer > 0) {
      layer--;
      if (projectionMode) {
        applyThresholdToScaledImage(threshold);
      } else {
        applyThresholdToLsmData(threshold, subregion);
      }
    } else {
      cerr << '\007';
    }
    break;

  case 'P':
    if (projectionMode) {
      projectionMode = false;
      if (currentView == PROJECTION) {
        currentView = LSM_SLICE;
      }
      applyThresholdToLsmData(threshold, subregion);
      if (verbosity > 1) {
        cerr << "Projection mode disabled." << endl;
      }
    } else {
      projectionMode = true;
      if (currentView == LSM_SLICE) {
        currentView = PROJECTION;
      }
      voxelProjection(lsmData, subregion, projectedImage); // compute the projection
      rescaleImage(projectedImage, scaledImage);
      applyThresholdToScaledImage(threshold);
      if (verbosity > 1) {
        cerr << "Projection mode enabled." << endl;
      }
    }
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
    float sx = static_cast<float>(mainWindow.width)/static_cast<float>(subregion.width);
    float sy = static_cast<float>(mainWindow.height)/static_cast<float>(subregion.height);
    float sf = min(sx, sy);
    glPixelZoom(sf, -sf);
     glRasterPos2i(mainWindow.world.xmin, mainWindow.world.ymax);
    switch(currentView) {
    case LSM_SLICE:
      lsmData.drawPixels(layer, band, subregion.width, subregion.height, subregion.col, subregion.row);
      break;
    case PROJECTION:
      scaledImage.drawPixels();
      break;
    case THRESHOLD:
      afterThreshold.drawPixels();
      if (showConvexHull) {
        glPushMatrix();
          glTranslatef(0, subregion.height, 0);
          glScalef(1.0, -1.0, 1.0);
          glColor3f(0.8, 0.8, 0.0);
          ch.drawPerimeter();
        glPopMatrix();
      }
      break;
    case DIALATED:
      processedImage[0].drawPixels();
      if (showConvexHull) {
        glPushMatrix();
          glTranslatef(0, subregion.height, 0);
          glScalef(1.0, -1.0, 1.0);
          glColor3f(0.8, 0.8, 0.0);
          ch.drawPerimeter();
        glPopMatrix();
      }
      break;
    case CLOSED:
      processedImage[1].drawPixels();
      if (showConvexHull) {
        glPushMatrix();
          glTranslatef(0, subregion.height, 0);
          glScalef(1.0, -1.0, 1.0);
          glColor3f(0.8, 0.8, 0.0);
          ch.drawPerimeter();
        glPopMatrix();
      }
    case LABELS:
      colorLabels.rgbPixels();
      if (showConvexHull) {
        glPushMatrix();
          glTranslatef(0, subregion.height, 0);
          glScalef(1.0, -1.0, 1.0);
          glColor3f(0.8, 0.8, 0.0);
          ch.drawPerimeter();
        glPopMatrix();
      }
      break;
    }

    if (infoVisibility) {
      glColor3f(0.1, 0.6, 1.0);
      int x = 10;
      int y = mainWindow.world.ymax - 24;
      glRasterPos2i(x, y);
      drawValue(fileName);
      y -= glutFontLineSkip;
      glRasterPos2i(x, y);
      drawValue("band = ", 7);
      drawValue(band);
      if (projectionMode) {
        drawValue(", projection");
      } else {
        drawValue(", layer = ", 10);
        drawValue(layer);
      }
      y -= glutFontLineSkip;
      glRasterPos2i(x, y);   
      drawValue("threshold = ", 12);
      drawValue(threshold);

      if (showConvexHull) {
        glRasterPos2i(10, 10 + 2*glutFontLineSkip);
        drawValue("Hull area = ");
        drawValue(static_cast<int>(hullArea));
        drawValue("/");
        drawValue(labelCount[0].second);
        drawValue(" = ");
        drawValue(static_cast<float>(hullArea/labelCount[0].second));
      }
    }

    switch(regionSelectMode) {
    case ENABLED:
      glRasterPos2i(10, 10);
      drawValue("Select first corner of crop rectangle with a mouse click.");
      break;
    case ONE_CORNER_SELECTED:
      glRasterPos2i(10, 10 + glutFontLineSkip);
      drawValue("Selected (");
      drawValue(selectionRegion.row);
      drawValue(", ");
      drawValue(selectionRegion.col);
      drawValue(")");
      glRasterPos2i(10, 10);
      drawValue("Select first corner with mouse, then drag and release.");
      break;
    case RECTANGLE_SELECTED:
      glBegin(GL_LINE_LOOP);
        glVertex2i(selectionRegion.col, selectionRegion.row);
        glVertex2i(selectionRegion.col + selectionRegion.width, selectionRegion.row);
        glVertex2i(selectionRegion.col + selectionRegion.width, selectionRegion.row + selectionRegion.height) ;
        glVertex2i(selectionRegion.col, selectionRegion.row + selectionRegion.height);
      glEnd();
      glRasterPos2i(10, 10);
      drawValue("Press x to crop, or ESC to cancel.");
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

void resizeMainWindow(	glutWindowInfo &w, 
						Rect &region, 
						double sx, double sy) {
	w.height = region.height*sy;
	w.width  = region.width*sx;
	w.world.xmin = 0;
	w.world.xmax = region.width;
	w.world.ymin = 0;
	w.world.ymax = region.height;
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

void runGlutLoop() {
	resizeMainWindow(mainWindow, subregion, 1.0, 1.0);
	initializeGlutWindow(mainWindow);
    currentView = LSM_SLICE;
    infoVisibility = true;
    showConvexHull = false;
	glutMainLoop();
}
  
