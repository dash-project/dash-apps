#ifndef ALL_PAIRS_H
#define ALL_PAIRS_H

#include <libdash.h>
#include <boost/log/trivial.hpp>

#include <sstream>
#include <fstream>
#include <iomanip>
#include <ctime>

#include "logger.h"

#define VALIDATE_KERNEL 0

typedef dash::Pattern<3>                         pattern_t;
typedef dash::NArray<double, 3, long, pattern_t> array_t;
typedef dash::NArray<double, 2, long>            marray_t;
typedef std::pair<int,int>                       tupel;
// Hostname Array Type
typedef std::vector<std::string>                 harray_t;
typedef dash::util::Timer<
  dash::util::TimeMeasure::Clock>                Timer;

namespace dio = dash::io::hdf5;

class AllPairs {

private:

    int                 repeats;
    bool                make_symmetric;
    array_t             results;
    /* summarized results */
    marray_t            medians;
    marray_t            mins;
    marray_t            maxs;

    harray_t            hostnames;

    int                 sub_diags;
    int                 current_diag;
    tupel               next_pair;
    std::string         filename;

    int                 myid;

public:

    AllPairs(
        int  rep        = 50,
        int  partests   = 0,
        bool make_sym   = false
    ):
        repeats(rep),
        make_symmetric(make_sym),
        filename(generateFilename()),
        myid(dash::myid())
    {
        sub_diags = (dash::size() -1) / partests + 1;
        auto pattern = createPattern();
        results.allocate(pattern);

        // medians pattern
        dash::Pattern<2> mpat(dash::SizeSpec<2>(dash::size(), dash::size()),
                              dash::DistributionSpec<2>(dash::BLOCKED, dash::CYCLIC),
                              dash::TeamSpec<2>(dash::Team::All()));
        medians.allocate(mpat);
        mins.allocate(mpat);
        maxs.allocate(mpat);

        if(myid == 0){
          hostnames.resize(dash::size());
          gatherHostnames();
          storeHostnames();
        }
    }

    template<
        typename APKernel>
    void runKernel(APKernel &kernel)
    {
        Timer             timer;
        timer.Calibrate(0);
 
        current_diag = 0;
        int               ndiags   = dash::size();
        double            measurestart = timer.Now();
        double            kernelstart  = timer.Now();
        double            elapsed;
        bool              is_inverted = false;
        int               rounds = 2;

        LOG_UNIT(info) << "Run Kernel";

        if(myid == 0) {
            std::cout << "== Running Kernel '" << kernel.getName()
                      << "' ==" << std::endl;
        }

        // init kernel
        kernel.init(repeats);

        // clear results
        double no_measure_value = -1;
        dash::fill(results.begin(), results.end(), no_measure_value);

        if(make_symmetric) {
          rounds = 1;
        }
 
        // calculate first pair
        updatePartUnits();
       
        for(int round = 0; round<rounds; ++round) {
            // Execute only first round if make_symmetric
            LOG_UNIT(debug) << "ROUND: " << round;
            while(current_diag <= ndiags) {
                int x,y;

                LOG_UNIT(debug) << "Measure diag " << current_diag
                                         << " out of " << ndiags;

                if(!is_inverted) {
                    x = next_pair.first;
                    y = next_pair.second;
                } else {
                    y = next_pair.first;
                    x = next_pair.second;
                }

                // Test only subsection of diagonal
                for(int s=0; s<sub_diags; s++){
                  dash::barrier();
                  LOG_UNIT(trace) << "Measure subdiag " << s
                                           << " out of " << sub_diags;

                  // Skip pair
                  if((x % sub_diags) != s &&
                     ((myid == x) || (myid == y)))
                  {
                    continue;
                  }

                  if(!(is_inverted && x == y)) {
                    // Measure r times
                    for(int r=0; r<repeats; ++r) {
                        if(myid == x) {
                            measurestart = timer.Now();
                        }
                        kernel.run(x,y);
                        if(myid == x) {
                            int int_repeats = kernel.getInternalRepeats();
                            elapsed = timer.ElapsedSince(measurestart);
                            results[x][y][r] = elapsed / int_repeats;
                            #if VALIDATE_KERNEL
                            if(!results.at(x,y,r).is_local()){
                              std::cerr << "Unit " << myid << " index "
                                        << x << "," << y << "," << r
                                        << " is not local" << std::endl;
                            }
                            #endif
                        }
                    }
                  }
                }
                kernel.reset();
                updatePartUnits();
            }
            if(!make_symmetric) {
                // reset test for second half of nton plane
                is_inverted = true;
                current_diag = 0;
            }
        }

        // Calculate Statistics 
        calculateStatistics();

        // Store results
        LOG_UNIT(info) << "Store results";
        dio::HDF5OutputStream os(this->filename + ".hdf5",
                                 dio::HDF5FileOptions::Append);
        os << dio::dataset(kernel.getName())
           << results
           << dio::dataset((kernel.getName() + "_median"))
           << medians
           << dio::dataset((kernel.getName() + "_min"))
           << mins
           << dio::dataset((kernel.getName() + "_max"))
           << maxs; 

        if(myid == 0) {
            double kernElapsed = timer.ElapsedSince(kernelstart) / 1000000; // Sec
            std::cout << "== done in " << kernElapsed
                      << " seconds ==" << std::endl;
        }
    }

private:
    static std::string generateFilename()
    {
        auto        now        = std::chrono::system_clock::now();
        std::time_t now_c      = std::chrono::system_clock::to_time_t(now);
        // Avoid using put_time as not supported by icc
        //auto        timestring = std::put_time(std::localtime(&now_c),
        //                                       "%F-%H-%M");
        char tstr[32];
        strftime(tstr, sizeof(tstr), "%F-%H-%M", std::localtime(&now_c));
        auto        timestring = std::string(tstr);

        std::stringstream fname;
        fname << "all-pairs-result-";
        fname << timestring;

        return fname.str();
    }

    const pattern_t createPattern()
    {
        // number of units
        int u_size = dash::size();

        dash::SizeSpec<3> sspec(u_size, u_size, this->repeats);
        dash::DistributionSpec<3> dspec(
            dash::BLOCKED,
            dash::CYCLIC,
            dash::BLOCKED);
        dash::TeamSpec<3> tspec(dash::Team::All());

        return pattern_t(sspec, dspec, tspec);
    }

    /**
     * Calculate participating units at given k-diagonal
     */
    void updatePartUnits()
    {
        int   n    = dash::size();
        int & k    = current_diag;

        LOG_UNIT(debug) << "update participating units";

        for(int y=0; y<n; ++y) {
            for(int x=0; x<=y; ++x) {
                if(((x+y) % n) == k) {
                    if(x == myid || y == myid) {
                        next_pair.first  = x;
                        next_pair.second = y;
                        ++current_diag;
                        return;
                    }
                }
            }
        }
        ++current_diag;
    }

  /**
  * Calculates Min, Max, Median of the measurements
  */
  void calculateStatistics(){
    auto onemeasure = std::vector<double>(repeats);

    LOG_UNIT(debug) << "calculate statistics";

    for(int x=0; x<dash::size(); ++x){
      for(int r=0; r<repeats; ++r){
        onemeasure[r] = results.local[0][x][r];
      }
      std::sort(onemeasure.begin(), onemeasure.end(), std::less<double>());
      double median = onemeasure[repeats / 2];
      double min    = onemeasure[0];
      double max    = onemeasure[repeats - 1];

      #if VALIDATE_KERNEL
      if(!medians.at(myid,x).is_local()){
         std::cerr << "Unit " << myid << " medians index "
         << myid << "," << x
         << " is not local" << std::endl;
      }
      #endif
      medians[myid][x] = median;
      mins[myid][x]    = min;
      maxs[myid][x]    = max;
    }
  }

  void gatherHostnames(){
    for(int i=0; i<dash::size(); ++i){
      hostnames[i] = dash::util::Locality::Hostname(i);
    }
  }

  void storeHostnames(){
   std::ofstream os(this->filename + "-hosts.csv");
   // Write Header
   os << "id;hostname" << std::endl;
   for(int i=0; i<dash::size(); ++i){
      os << i << ";" << hostnames[i];
      if(i != dash::size()-1){
        os << std::endl;
      }
   }
  }
};

#endif // ALL_PAIRS_H
