#include <iostream>
#include <libdash.h>

#include "lulesh.h"
#include "lulesh-opts.h"
#include "lulesh-dash.h"
#include "lulesh-util.h"

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
  if( (myRank == 0) && (!opts.quiet()) ) {
    opts.printBanner(std::cout);
  }

  Domain dom(opts);

  // if(dash::myid()==0) print_config(dom, std::cout);
  // if( dash::myid()==0 ) peek( &(dom.nodalMass(0)), (size_t)20 );

  dom.InitialBoundaryExchange();

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

