#include <iostream>
#include <libdash.h>

#include "lulesh.h"
#include "lulesh-util.h"
#include "lulesh-dash.h"

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
  dom.InitialBoundaryExchange();

  /*
  if( dash::myid()==0 ) {
    peek<double>(&(dom.nodalMass(0)), dom.numElem() );
  }
  */

  //
  // --- main simulation loop ---
  //
  while( (dom.time() < dom.stoptime()) && (dom.cycle() < opts.its()) )
    {
      dom.TimeIncrement();
      dom.LagrangeLeapFrog();

      if( (opts.showProg()!=0) && (opts.quiet()==0) && (myRank==0) ) {
	std::cout << "cycle = " << dom.cycle() << ", ";
	std::cout << "time = " << dom.time() << ", ";
	std::cout << "dt = " << dom.deltatime() << std::endl;
      }
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

