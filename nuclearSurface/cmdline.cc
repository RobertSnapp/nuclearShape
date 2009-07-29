/***********
 * cmdline.cc
 * Project nuclearSurface
 *
 * Created by Robert R. Snapp on 2007-09-04
 *
 * Functions for parsing the command line. Note that all variable options
 * are communicated to the rest of the program via global variables contained
 * in the header file project.h. Current options are
 *
 * -b<int> channel index (or band).
 * -l<int> image layer index
 * -t<int> initial threshold
 * -v<int> verbosity level
 * -G      enable graphics interface
 * -H      help.
 *
 *********/


#include "nuclearSurface.h"
#include <cstdlib>
#include <unistd.h>  // getopt
//#include <string>
#include <iostream>

extern int   optind; // Used by getopt
extern int   opterr; // Used by getopt
extern char *optarg; // Used by getopt

using namespace std;
static const char* OPTION_FLAGS = "b:c:d:h:l:o:p:r:t:v:w:GH";


void describeUsage(char const *commandName) {
  cerr << "Usage: " << commandName
       << "[options] <filename>" << endl
       << "where [options] include: " << endl
       << "  -b <int>      sets the current band index to <int>; default is 0." << endl
       << "  -c <int>      sets the column of the subwindow to <int>; default is 0." << endl
       << "  -d <int>      sets the depth of the subvolume to <int>; default is the image depth." << endl
       << "  -h <int>      sets the height of the subwindow to <int>; default is image height." << endl
       << "  -l <int>      sets the current layer index to <int>; default is 0." << endl
       << "  -o <str>      print the results in the filename defined by <str>; default is stdout." << endl
       << "  -p <str>      sets the optional lsm directory to <str>; default is current directory." << endl
       << "  -r <int>      sets the row of the subwindow to <int>; default is 0" << endl
       << "  -t <int>      sets the initial threshold value to <int>; default is 30." << endl
       << "  -v <int>      sets the verbosity value to <int>. A positive verbosity" << endl
       << "                sends debugging messages to stdout. Default is 0." << endl
       << "  -w <int>      sets the width of the subwindow to <int>; default is image width." << endl
       << "  -G            enables graphics mode for interactive processing; default is disabled." << endl
       << "  -H            enables help mode, printing just this message." << endl<<endl;
  cerr << "<filename> should be the path to a valid lsm file." << endl;
  return;
}


// function parseCommandLine processes the command line options, and
// returns optind, the index of the next command line to be parsed,
// i.e., argv[optind].
int parseCommandLine(int argc, char** argv) {
  int c;

  //  Seed = init_seed();			/* Initialize the Random Number Seed */

  opterr = 0;					/* Suppress external error check */
  while( (c = getopt(argc, argv, OPTION_FLAGS)) != EOF)
	switch(c) {
    case 'G':
      graphics = true;
      break;
	case 'H':
	  describeUsage(argv[0]);
	  exit(0);
      break;
	case 'b':   // lsm channel (or band) to analyze.
	  band = atoi(optarg);
	  break;
    case 'c':
      subvol.col = atoi(optarg);
      break;
    case 'd':
      subvol.depth = atoi(optarg);
      break;
    case 'h':
      subvol.height = atoi(optarg);
      break;
	case 'l':   // lsm layer to analyze.
	  subvol.layer = atoi(optarg);		
      break;
    case 'o':
      outputFilename = optarg;
      break;
    case 'p':
      directory = optarg;
      cerr << "path = " << directory << endl;
      break;
    case 'r':
      subvol.row = atoi(optarg);
      break;
    case 't':   // theshold to apply   
      threshold = atoi(optarg);		/* size of pattern data base */
      break;
    case 'v':
      verbosity = atoi(optarg);
      break;
    case 'w':
      subvol.width = atoi(optarg);
      break;
    default:
      cerr<< "Option " << c << " is undefined.\n";
      exit(1);
    }

  return optind;
}
