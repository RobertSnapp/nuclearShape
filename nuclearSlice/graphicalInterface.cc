/*
 * graphicalInterface.cc
 * Project: nuclearSlice
 *
 * Created by Robert R. Snapp on 2007-07-16
 *
 * This file contains routines to view a specified layer of an lsm
 * data file using OpenGL, and the outcomes of various image processing
 * and analysis procedures.
 *
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
void* const glutFont = GLUT_BITMAP_8_BY_13;
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
bool invVideo = false;
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
       << subregion.col << ", "
       << subregion.row << ", "
       << subregion.width << ", "
       << subregion.height
       << endl;

  loadAndPrepareLsmFile(fileName.c_str(), subregion);
  componentLabels.clear();
  colorLabels.clear();      
  currentView = (projectionMode ? PROJECTION : LSM_SLICE);
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
	case PROJECTION: // pass through
	case THRESHOLD:  // pass through
	  //	  applyThreshold(subregion);  // Recomputes projection if needed.
	  // do nothing
	  break;
	case DIALATED:
	case CLOSED:
	  applyThreshold(subregion);
	  closureOp();
	  break;
	case LABELS:
	  applyThreshold(subregion);
	  closureOp();
	  int nLabels = collectComponents();
	  makeColorLabels(nLabels);
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
      // selectionRegion.row = mainWindow.height - y;
      selectionRegion.row = y;
    }
    break;
  case ONE_CORNER_SELECTED:
    if (button == GLUT_LEFT_BUTTON && state == GLUT_UP) {
      selectionRegion.width = max(x, selectionRegion.col) - min(x, selectionRegion.col);
      selectionRegion.col   = min(x, selectionRegion.col);

      // int yp = mainWindow.height - y;
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
    double sx = static_cast<double>(w)/static_cast<double>(subregion.width);
    double sy = static_cast<double>(h)/static_cast<double>(subregion.height);
    // Ensure that the aspects of the window and image agree.
    if (sx < sy) {
		h = static_cast<int>(ceil(sx*subregion.height)); // Reduce the window height to preserve the aspect ratio.
    } else {
		w = static_cast<int>(ceil(sy*subregion.width)); // Reduce the window width to preserve the aspect ratio.
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
    cerr << "reshape(...) glut callback (file " << __FILE__ << ", line " << __LINE__ << ")"
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
	regionSelectMode = DISABLED;
    if (band > 0) {
      band--;
	  applyThreshold(subregion);
	  updateProcessedImage();
    } else {
	  cerr << '\007' << "Command ignored." << endl;
	}
	glutPostRedisplay();
//       if (projectionMode) {
//         applyThresholdToScaledImage(threshold);
//       } else {
//         applyThresholdToLsmData(threshold, subregion);
//       }
// 	  updateProcessedImage();
//       glutPostRedisplay();
//     } else {
// 	  cerr << '\007' << "Command ignored." << endl;
// 	}
    break;
  case GLUT_KEY_RIGHT:    // next band.
    regionSelectMode = DISABLED;
    if (band < static_cast<int>(lsmData.bands()) - 1) {
      band++;

	  applyThreshold(subregion);
	  updateProcessedImage();
    } else {
	  cerr << '\007' << "Command ignored." << endl;
	}
	glutPostRedisplay();
//       if (projectionMode) {
//         applyThresholdToScaledImage(threshold);
//       } else {
//         applyThresholdToLsmData(threshold, subregion);
//       }
// 	  updateProcessedImage();
//       glutPostRedisplay();
//     } else {
// 	  cerr << '\007' << "Command ignored." << endl;
// 	}
    break;
  case GLUT_KEY_UP: // increase threshold: does not require recomputing the projection.
    if (threshold < 255) {
      threshold++;
      if (projectionMode) {
        applyThresholdToScaledImage(threshold);
      } else {
        applyThresholdToLsmData(threshold, subregion);
      }
	  updateProcessedImage();
      glutPostRedisplay();
    }
    break;
  case GLUT_KEY_DOWN: // decrease threshold: does not require recomputing the projection.
    if (threshold > 0) {
      threshold--;
      if (projectionMode) {
        applyThresholdToScaledImage(threshold);
      } else {
        applyThresholdToLsmData(threshold, subregion);
      }
	  updateProcessedImage();
      glutPostRedisplay();
    }
    break;
  case GLUT_KEY_PAGE_UP: // increase threshold by 5: does not require recomputing the projection.
    if (threshold < 255) {
      threshold = min(255, threshold + 5);
      if (projectionMode) {
        applyThresholdToScaledImage(threshold);
      } else {
        applyThresholdToLsmData(threshold, subregion);
      }
	  updateProcessedImage();
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
	int requestedBand = c - '0';
	if (band != requestedBand && requestedBand < lsmData.bands()) {
	  band = requestedBand;
	  applyThreshold(subregion);  // may require new image projection before generating threshold map.
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
    int nLabels = collectComponents();
    makeColorLabels(nLabels);
    currentView = LABELS;
    showConvexHull = false;
    break;
  case 'h':
    regionSelectMode = DISABLED;
    if (colorLabels.rows() > 0 && colorLabels.cols() > 0) {
	  // Compute eignevalues of largest blob
	  computeEigenvalues(1, blob);
	  
	  computeComponentHull(1, ch);
      hullArea = ch.area();
      showConvexHull = true;
    }
    break;
  case 'm':
    regionSelectMode = DISABLED;
    updateProcessedImage();
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
	showConvexHull = false;
    if (regionSelectMode == RECTANGLE_SELECTED) {
	  // The user requests to reduce the focus to the indicated subregion relative to the window origin.
	  subregion.col += selectionRegion.col;
	  subregion.row += selectionRegion.row;
	  subregion.width = selectionRegion.width;
	  subregion.height = selectionRegion.height;
      cropImage();
      regionSelectMode = DISABLED;
    } else if (regionSelectMode == ENABLED) {
	  cerr << "Restoring subregion to entire image." << endl;
      subregion.col = 0;
      subregion.row = 0;
      subregion.width = lsmData.cols();
      subregion.height = lsmData.rows();
      cropImage();
      regionSelectMode = DISABLED;
    }
    break;

  case 'n':
    regionSelectMode = DISABLED;
    if (! projectionMode && layer < static_cast<int>(lsmData.layers()) - 1) {
      layer++;
      applyThresholdToLsmData(threshold, subregion);
	  updateProcessedImage();
    } else {
      cerr << '\007' << "Command ignored." << endl;
    }
		
    break;	
  case 'p':
    regionSelectMode = DISABLED;
    if (! projectionMode && layer > 0) {
      layer--;
      applyThresholdToLsmData(threshold, subregion);
	  updateProcessedImage();
    } else {
      cerr << '\007' << "Command ignored." << endl;
    }
    break;

  case 'P':
    regionSelectMode = DISABLED;
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
	updateProcessedImage();
    break;

  case 'R': // Reset the subregion to the full image
    regionSelectMode = DISABLED;
    showConvexHull = false;
    subregion.row = 0;
    subregion.col = 0;
    subregion.width = lsmData.cols();
    subregion.height = lsmData.rows();
    if (projectionMode) {
      voxelProjection(lsmData, subregion, projectedImage); // compute the projection
      rescaleImage(projectedImage, scaledImage);
      applyThresholdToScaledImage(threshold);
    } else {
    applyThresholdToLsmData(threshold, subregion);
    }
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
      subregion.row = 0;
      subregion.col = 0;
      subregion.width = 0;
      subregion.height = 0;

      if (parseFileLine(*fileRosterPos)) {
		loadAndPrepareLsmFile(fileName.c_str(), subregion);
		componentLabels.clear();
		colorLabels.clear();      
		currentView = (projectionMode ? PROJECTION : LSM_SLICE);
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
	  while(! status && fileRosterPos < fileRoster.end() - 1) {
		++fileRosterPos;
		
		if (fileRosterPos < fileRoster.end()) {
		  status = parseFileLine(*fileRosterPos);  // returns true if record is nonempty
		}
	  }

	  if (status) {
		subregion.row = 0;
		subregion.col = 0;
		subregion.width = 0;
		subregion.height = 0;
		loadAndPrepareLsmFile(fileName.c_str(), subregion);
		componentLabels.clear();
		colorLabels.clear();      
		currentView = (projectionMode ? PROJECTION : LSM_SLICE);
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
	  while (!status && fileRosterPos > fileRoster.begin()) {
		--fileRosterPos;
		if (fileRosterPos >= fileRoster.begin()) {
		  status = parseFileLine(*fileRosterPos);
		}
	  }
		
	  if (status) {
		subregion.row = 0;
		subregion.col = 0;
		subregion.width = 0;
		subregion.height = 0;
		loadAndPrepareLsmFile(fileName.c_str(), subregion);
		componentLabels.clear();
		colorLabels.clear();
		currentView = (projectionMode ? PROJECTION : LSM_SLICE);
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
    subregion.row = 0;
    subregion.col = 0;
    subregion.width = 0;
    subregion.height = 0;
    do {
       string f;
       cout << "Filename? ";
       cin >> f;
       fileName = f;
       status = loadAndPrepareLsmFile(f.c_str(), subregion);
       if (! status) {
         cerr << "Could not parse file " << f << endl;
       }
    } while (status == false);
    componentLabels.clear();
    currentView = (projectionMode ? PROJECTION : LSM_SLICE);
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

    double sx = static_cast<double>(mainWindow.width)/static_cast<double>(subregion.width);
    double sy = static_cast<double>(mainWindow.height)/static_cast<double>(subregion.height);
    double sf = min(sx, sy);

    glPixelZoom(sf, -sf);

    switch(currentView) {
    case LSM_SLICE:
      lsmData.drawPixels(layer, band, subregion.width, subregion.height, subregion.col, subregion.row, invVideo);
      break;
    case PROJECTION:
	  if (scaledImage.rows() > 0 && scaledImage.cols() > 0)
		scaledImage.drawPixels(invVideo);
      break;
    case THRESHOLD:
	  if (afterThreshold.rows() > 0 && afterThreshold.cols() > 0)
		afterThreshold.drawPixels(invVideo);
      break;
    case DIALATED:
	  if (processedImage[0].rows() > 0 && processedImage[0].cols() > 0)
		processedImage[0].drawPixels(invVideo);
      break;
    case CLOSED:
	  if (processedImage[1].rows() > 0 && processedImage[1].cols() > 0)
		processedImage[1].drawPixels(invVideo);
      break;
    case LABELS:
	  if (colorLabels.rows() > 0 && colorLabels.cols() > 0)
		colorLabels.rgbPixels();
      break;
    }

    if (showConvexHull) {
      glPushMatrix();
      glTranslatef(0.5, mainWindow.world.ymax-0.5, 0);
      glScalef(1.0f, -1.0f, 1.0f);
      glColor3f(0.8, 0.8, 0.0);
      ch.drawPerimeter();
      glPopMatrix();
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
      if (projectionMode) {
        drawValue(", projection");
      } else {
        drawValue(", layer = ");
        drawValue(layer);
      }
      y -= static_cast<int>(glutFontLineSkip/sf);
      glRasterPos2i(x, y);   
      drawValue("threshold = ");
      drawValue(threshold);

      if (showConvexHull) {
        glRasterPos2i(static_cast<GLint>(10/sf), static_cast<GLint>((10 + 3*glutFontLineSkip)/sf));
        drawValue("Hull area = ");
        drawValue(static_cast<int>(hullArea));
        drawValue("/");
        drawValue(labelCount[0].second);
        drawValue(" = ");
        drawValue(static_cast<float>(hullArea/labelCount[0].second));
		drawValue(", EVs = (");
		drawValue(blob.lam1);
		drawValue(", ");
		drawValue(blob.lam2);
		drawValue("),");
		glRasterPos2i(static_cast<GLint>(10/sf), static_cast<GLint>((10 + 2*glutFontLineSkip)/sf));
		drawValue("ecc = ");
		drawValue(sqrt(blob.lam1/blob.lam2));
		drawValue(", angle of elongation = ");
		double angle = (180./M_PI)*atan2(blob.ev1[1], blob.ev1[0]);
		drawValue(angle);
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
                  Rect &region, 
                  double sx, double sy) {
  w.height = static_cast<int>(ceil(region.height*sy));
  w.width  = static_cast<int>(ceil(region.width*sx));
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

void resizeMainWindow() {
    // The dimensions of the main window should conform to the subregion.
  float const screenWidth  = 1280;
  float const screenHeight =  854 - 50;  // dont include the menu bar.
  float sx = min(1.0f, screenWidth/subregion.width);
  float sy = min(1.0f, screenHeight/subregion.height);
  float sf = min(sx, sy);

  resizeWindow(mainWindow, subregion, sf, sf);

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
  currentView = (projectionMode ? PROJECTION : LSM_SLICE);
  infoVisibility = true;
  showConvexHull = false;
  glutMainLoop();
}
  
