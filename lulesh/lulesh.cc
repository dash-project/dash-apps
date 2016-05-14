#include <iostream>
#include <libdash.h>

#include "lulesh.h"
#include "lulesh-util.h"
#include "lulesh-dash.h"
#ifdef USE_MPI
#include "lulesh-mpi.h"
#endif

int main(int argc, char *argv[])
{
  dash::init(&argc, &argv);

  auto myRank   = dash::myid();
  auto numRanks = dash::size();

  CmdLineOpts opts(numRanks, myRank);
  opts.parseCommandLineOptions(argc, argv);

  if( (!opts.valid()) ) {
    if( myRank==0 ) {
      opts.printHelp(std::cerr);
    }
    dash::finalize();
    return 0;
  }

  if( (myRank == 0) && (opts.quiet()==0) ) {
    opts.printBanner(std::cout);
  }

  Domain dom(opts);

  if( dash::myid()==0 ) {
    dom.print_config(std::cout);
  }
  /*
  for( int i=0; i<dom.numNode(); ++i  ) {
    dom.nodalMass(i)=100;
  }

  if( dash::myid()==0 ) {
    peek<double>(&(dom.nodalMass(0)), dom.numElem() );
  }
  */

  dash::barrier();
#ifdef USE_MPI
  Comm comm(dom);

  Domain_member fieldData;
  fieldData = &Domain::nodalMass;

  Domain *locDom = &dom;
  // Initial domain boundary communication
  CommRecv(*locDom, comm, MSG_COMM_SBN, 1,
	   locDom->sizeX() + 1, locDom->sizeY() + 1, locDom->sizeZ() + 1,
	   true, false) ;
  CommSend(*locDom, comm, MSG_COMM_SBN, 1, &fieldData,
	   locDom->sizeX() + 1, locDom->sizeY() + 1, locDom->sizeZ() +  1,
	   true, false) ;
  CommSBN(dom, comm, 1, &fieldData) ;

  // End initialization
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if( dash::myid()==0 ) {
    peek<double>(&(dom.nodalMass(0)), dom.numElem() );
  }

  dash::finalize();
}

template<typename T>
void peek(T *ptr, int len, std::ostream& os)
{
  for( int i=0; i<len; ++i ) {
    os << (*ptr++) << " ";
  }
  os << std::endl;

}

