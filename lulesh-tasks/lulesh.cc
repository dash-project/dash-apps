#include <iostream>
#include <libdash.h>

#include <mpi.h>

#include "lulesh.h"
#include "lulesh-opts.h"
#include "lulesh-dash.h"
#include "lulesh-util.h"
#include "extrae.h"

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

  MPI_Barrier(MPI_COMM_WORLD);

  // exclude thread creation from timing loop
  dash::tasks::ASYNC([&](){
    if( (myRank == 0) && (!opts.quiet()) ) {
      opts.printBanner(std::cout);
    }
  });
  dash::tasks::complete();

  Domain dom(opts);

  auto nnodes = (1+opts.nx()) * (1+opts.nx()) * (1+opts.nx());

  EXTRAE_INIT();
  EXTRAE_TASK_INIT();

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
  Int_t cycle = 0;
  while( /*(dom.time() < dom.stoptime()) &&*/ (cycle < opts.its()) )
  {
    // if( dash::myid()==0 ) {
    // std::cout << chksum(&(dom.x(0)), dom.numNode()) << std::endl;
    // }

    //dom.PrintDomain();

    dom.TimeIncrement();
    dom.LagrangeLeapFrog();

    dash::tasks::ASYNC(
      [&](){
        if( (opts.showProg()!=0) && (opts.quiet()==0) && (myRank==0) ) {
          std::cout << "cycle = " << dom.cycle() << ", ";
          std::cout << "time = " << dom.time() << ", ";
          std::cout << "dt = " << dom.deltatime() << std::endl;
        }
      },
      dash::tasks::in(dom.deltatime()),
      // make sure this tasks runs when all previous tasks of this iteration
      // are done
      dash::tasks::in(dom.dthydro())
    );

    cycle++;

#if 0
    // TODO: this check prevents us from overlapping force calculations with
    //       time constraint calculations
    dash::tasks::async(
      [&](){
        if ((dom.time() >= dom.stoptime())) {
          // cancel all remaining  tasks
          dash::tasks::cancel_barrier();
        }
      },
      dash::tasks::in(&dom.dthydro())
    );
#endif
  }
  dash::tasks::complete();
  double end = dom.wtime()-start;
  double elapsed = dom.reduce_max(end);

  if( (myRank == 0) && (!opts.quiet()) ) {
    dom.VerifyAndWriteFinalOutput(elapsed, opts.nx(), numRanks);
  }
  dash::finalize();

  return 0;
}

