
#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "lulesh.h"
#include "lulesh-util.h"

using std::endl;

/* Helper function for converting strings to ints, with error checking */
int StrToInt(const char *token, int *retVal)
{
  const char *c ;
  char *endptr ;
  const int decimal_base = 10 ;

  if (token == NULL)
    return 0 ;

  c = token ;
  *retVal = (int)strtol(c, &endptr, decimal_base) ;
  if((endptr != c) && ((*endptr == ' ') || (*endptr == '\0')))
    return 1 ;
  else
    return 0 ;
}

static void ParseError(const char *message, int myRank)
{
   if (myRank == 0) {
      printf("%s\n", message);
#if USE_MPI
      //MPI_Abort(MPI_COMM_WORLD, -1);
#else
      //exit(-1);
#endif
   }
}

void CmdLineOpts::printBanner(std::ostream& os)
{
  CmdLineOpts& opts = (*this);

  os << endl;
  os << "==============================" << endl;
  os << "== DASH port of LULESH 2.03 ==" << endl;
  os << "==============================" << endl;
  os << endl;
  os << "Problem size   : " << opts.nx() << "^3 per domain" << endl;
  os << "Num processors : " << opts.numRanks() << endl;
  os << "Processor grid : " <<
    opts.px() << " x " <<
    opts.py() << " x " <<
    opts.pz() << endl;

#if _OPENMP
  os << "Num threads    : " << omp_get_max_threads() << endl;
#else
  os << "Num threads    : " << 1 << endl;
#endif
  os << "Total elements : " <<
    (long long int)(opts.numRanks()*opts.nx()*opts.nx()*opts.nx()) << endl;
  os << endl;
}

void CmdLineOpts::printHelp(std::ostream& os)
{
  os << "Options: " << endl;
  os << " [-q ]        Quiet mode, suppress all output\n";
  os << " [-i  <int> ] Number of cycles to run\n";
  os << " [-s  <int> ] Specifiy number of local elements (s^3)\n";
  os << " [-px <int> ] Number of procs in x dimension\n";
  os << " [-py <int> ] Number of procs in y dimension\n";
  os << " [-pz <int> ] Number of procs in z dimension\n";
  os << " [-r  <int> ] Number of distinct regions\n";
  os << " [-b  <int> ] Load balance between regions of a domain\n";
  os << " [-c  <int> ] Extra cost of more expensive regions\n";
  os << " [-p ]        Print out progress\n";
  os << " [-v ]        Write an output file for VisIt\n";
  os << " [-h ]        Print help message\n\n";
  os << endl << endl;
}

void CmdLineOpts::parseCommandLineOptions(int argc, char *argv[])
{
  if(argc > 1) {
    int i = 1;

    while(i < argc) {
      int ok;

      /* -i <iterations> */
      if(strcmp(argv[i], "-i") == 0) {
	if (i+1 >= argc) {
	  ParseError("Missing integer argument to -i", m_myRank);
	  m_valid=false;
	}
	ok = StrToInt(argv[i+1], &m_its);
	if(!ok) {
	  ParseError("Parse Error on option -i integer value required after argument\n", m_myRank);
	  m_valid=false;
	}
	i+=2;
      }

      /* -s <size, sidelength> */
      else if(strcmp(argv[i], "-s") == 0) {
	if (i+1 >= argc) {
	  ParseError("Missing integer argument to -s\n", m_myRank);
	}
	ok = StrToInt(argv[i+1], &m_nx);
	if(!ok) {
	  ParseError("Parse Error on option -s integer value required after argument\n", m_myRank);
	  m_valid=false;
	}
	i+=2;
      }

      /* -px <number of processes> */
      else if(strcmp(argv[i], "-px") == 0) {
	if (i+1 >= argc) {
	  ParseError("Missing integer argument to -px\n", m_myRank);
	  m_valid=false;
	}
	ok = StrToInt(argv[i+1], &m_px);
	if(!ok) {
	  ParseError("Parse Error on option -s integer value required after argument\n", m_myRank);
	  m_valid=false;
	}
	i+=2;
      }

      /* -py <number of processes> */
      else if(strcmp(argv[i], "-py") == 0) {
	if (i+1 >= argc) {
	  ParseError("Missing integer argument to -py\n", m_myRank);
	  m_valid=false;
	}
	ok = StrToInt(argv[i+1], &m_py);
	if(!ok) {
	  ParseError("Parse Error on option -s integer value required after argument\n", m_myRank);
	  m_valid=false;
	}
	i+=2;
      }

      /* -pz <number of processes> */
      else if(strcmp(argv[i], "-pz") == 0) {
	if (i+1 >= argc) {
	  ParseError("Missing integer argument to -pz\n", m_myRank);
	  m_valid=false;
	}
	ok = StrToInt(argv[i+1], &m_pz);
	if(!ok) {
	  ParseError("Parse Error on option -s integer value required after argument\n", m_myRank);
	  m_valid=false;
	}
	i+=2;
      }

      /* -r <numregions> */
      else if (strcmp(argv[i], "-r") == 0) {
	if (i+1 >= argc) {
	  ParseError("Missing integer argument to -r\n", m_myRank);
	  m_valid=false;
	}
	ok = StrToInt(argv[i+1], &m_numReg);
	if (!ok) {
	  ParseError("Parse Error on option -r integer value required after argument\n", m_myRank);
	  m_valid=false;
	}
	i+=2;
      }

      /* -f <numfilepieces> */
      else if (strcmp(argv[i], "-f") == 0) {
	if (i+1 >= argc) {
	  ParseError("Missing integer argument to -f\n", m_myRank);
	  m_valid=false;
	}
	ok = StrToInt(argv[i+1], &m_numFiles);
	if (!ok) {
	  ParseError("Parse Error on option -f integer value required after argument\n", m_myRank);
	  m_valid=false;
	}
	i+=2;
      }

      /* -p */
      else if (strcmp(argv[i], "-p") == 0) {
	m_showProg = 1;
	i++;
      }

      /* -q */
      else if (strcmp(argv[i], "-q") == 0) {
	m_quiet = 1;
	i++;
      }

      /* -b */
      else if (strcmp(argv[i], "-b") == 0) {
	if (i+1 >= argc) {
	  ParseError("Missing integer argument to -b\n", m_myRank);
	  m_valid=false;
	}
	ok = StrToInt(argv[i+1], &m_balance);
	if (!ok) {
	  ParseError("Parse Error on option -b integer value required after argument\n", m_myRank);
	  m_valid=false;
	}
	i+=2;
      }

      /* -c */
      else if (strcmp(argv[i], "-c") == 0) {
	if (i+1 >= argc) {
	  ParseError("Missing integer argument to -c\n", m_myRank);
	  m_valid=false;
	}
	ok = StrToInt(argv[i+1], &m_cost);
	if (!ok) {
	  ParseError("Parse Error on option -c integer value required after argument\n", m_myRank);
	  m_valid=false;
	}
	i+=2;
      }

      /* -v */
      else if (strcmp(argv[i], "-v") == 0) {
#if VIZ_MESH
	m_viz = 1;
#else
	ParseError("Use of -v requires compiling with -DVIZ_MESH\n", m_myRank);
	m_valid=false;
#endif
	i++;
      }

      /* -h */
      else if (strcmp(argv[i], "-h") == 0) {
	printHelp(std::cerr);
	m_valid=false;
#if USE_MPI
	//MPI_Abort(MPI_COMM_WORLD, 0);
#else
	//exit(0);
#endif
      }

      else {
	char msg[80];
	sprintf(msg, "ERROR: Unknown command line argument: %s\n", argv[i]);
	ParseError(msg, m_myRank);
	m_valid=false;
	i++;
      }
    }
  }

  // DASH: if the process configuration is not specified explicity,
  // determine it automatically for powers of 3, like the MPI version
  // of LULESH does
  if( m_px==0 && m_py==0 && m_pz==0 ) {
    Index_t side = std::cbrt(numRanks()+0.5);
    if( side*side*side == numRanks() )
      m_px = m_py = m_pz = side;
  }

  if( m_px * m_py * m_pz != numRanks()) {
    char msg[80];
    sprintf(msg, "ERROR: Invalid proc. configuration: %d x %d x %d != %d\n",
	    m_px, m_py, m_pz, numRanks());
    ParseError(msg, m_myRank);
    m_valid=false;
  }
}


