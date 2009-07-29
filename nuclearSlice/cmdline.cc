/***********
 * cmdline.cc
 * Project nuclearSlice
 *
 * Created by Robert R. Snapp on 2007-07-14
 *
 * Functions for parsing the command line. Note that all variable options
 * are communicated to the rest of the program via global variables contained
 * in the header file project.h. Current options are
 *
 * -b<int> channel index (or band).
 * -l<int> image layer index
 * -t<int> initial threshold
 * -v<int> verbosity level
 * -P      projection mode (alanyze
 * -G      enable graphics interface
 * -H      help.
 *
 *********/


#include "project.h"
#include <cstdlib>
#include <unistd.h>  // getopt
//#include <string>
#include <iostream>

extern int   optind; // Used by getopt
extern int   opterr; // Used by getopt
extern char *optarg; // Used by getopt

using namespace std;
static const char* OPTION_FLAGS = "b:c:d:h:l:o:r:t:v:w:GHP";
const char directoryDelimiter = '/';

void describeUsage(char const *commandName) {
  cerr << "Usage: " << commandName
       << "[options] <filename>" << endl
       << "where [options] include: " << endl
       << "  -b <int>      sets the current band index to <int>; default is 0." << endl
       << "  -c <int>      sets the column of the subwindow to <int>; default is 0." << endl
       << "  -d <str>      sets the optional lsm directory to <str>; default is current directory." << endl
       << "  -h <int>      sets the height of the subwindow to <int>; default is image height." << endl
       << "  -l <int>      sets the current layer index to <int>; default is 0." << endl
       << "  -o <str>      print the results in the filename defined by <str>; default is stdout." << endl
       << "  -r <int>      sets the row of the subwindow to <int>; default is 0" << endl
       << "  -t <int>      sets the initial threshold value to <int>; default is 30." << endl
       << "  -v <int>      sets the verbosity value to <int>. A positive verbosity" << endl
       << "                sends debugging messages to stdout. Default is 0." << endl
       << "  -w <int>      sets the width of the subwindow to <int>; default is image width." << endl
       << "  -G            enables graphics mode for interactive processing; default is disabled." << endl
       << "  -H            enables help mode, printing just this message." << endl
       << "  -P            enables projection mode, in which the intesities of all" << endl
       << "                of the image layers are summed together. Default is disabled." << endl << endl;
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
    case 'P':
      projectionMode = true;
      break;
	case 'b':   // lsm channel (or band) to analyze.
	  band = atoi(optarg);
	  break;
    case 'c':
      subregion.col = atoi(optarg);
      break;
    case 'd':
      directory = optarg;
	  
	  int dlength = directory.length();
	  if (dlength > 0) {
		if (directory[dlength - 1] != directoryDelimiter) {
		  directory.push_back(directoryDelimiter);
		}
	  }
      cerr << "directory = " << directory << endl;
      break;
    case 'h':
      subregion.height = atoi(optarg);
      break;
	case 'l':   // lsm layer to analyze.
	  layer = atoi(optarg);		
      break;
    case 'o':
      outputFilename = optarg;
      break;
    case 'r':
      subregion.row = atoi(optarg);
      break;
    case 't':   // theshold to apply   
      threshold = atoi(optarg);		/* size of pattern data base */
      break;
    case 'v':
      verbosity = atoi(optarg);
      break;
    case 'w':
      subregion.width = atoi(optarg);
      break;
    default:
      cerr<< "Option " << c << " is undefined.\n";
      exit(1);
    }


  // srand48(Seed);				/* Initialize drand48 random number gen. */

  return optind;
}
