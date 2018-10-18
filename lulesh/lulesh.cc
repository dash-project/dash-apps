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

  Timer::Calibrate();

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

  auto nnodes = (1+opts.nx()) * (1+opts.nx()) * (1+opts.nx());

  // if(dash::myid()==0) print_config(dom, std::cout);
  // if( dash::myid()==0 ) peek( &(dom.nodalMass(0)), (size_t)20 );
  // fill( &(dom.nodalMass(0)), nnodes, 100.0+dash::myid());

  dom.InitialBoundaryExchange();

  //  std::cout << dash::myid() << " "
  //	    << chksum(&(dom.nodalMass(0)), nnodes)  << std::endl;
  //  if( dash::myid()==0 ) peek( &(dom.nodalMass(0)), nnodes );


  //
  // --- main simulation loop ---
  //
  double start = dom.wtime();
  while( (dom.time() < dom.stoptime()) && (dom.cycle() < opts.its()) )
    {
      // if( dash::myid()==0 ) {
      // std::cout << chksum(&(dom.x(0)), dom.numNode()) << std::endl;
      // }

      dom.TimeIncrement();
      dom.LagrangeLeapFrog();

      if( (opts.showProg()!=0) && (opts.quiet()==0) && (myRank==0) ) {
        std::cout << "cycle = " << dom.cycle() << ", ";
        std::cout << "time = " << dom.time() << ", ";
        std::cout << "dt = " << dom.deltatime() << std::endl;
      }
    }
  double end = dom.wtime()-start;
  double elapsed = dom.reduce_max(end);

  if( (myRank == 0) && (!opts.quiet()) ) {
    dom.VerifyAndWriteFinalOutput(elapsed, opts.nx(), numRanks);
  }
  dash::finalize();
}

