#ifndef ALL_PAIRS_H
#define ALL_PAIRS_H

#include <libdash.h>

typedef dash::Pattern<3>                         pattern_t;
typedef dash::NArray<double, 3, long, pattern_t> array_t;
typedef std::pair<int,int>                       tupel;
typedef dash::util::Timer<
dash::util::TimeMeasure::Clock
> Timer;

namespace dio = dash::io::hdf5;

class AllPairs {

private:

    int                 repeats;
    bool                make_symmetric;
    int                 return_node;
    array_t             results;

    int                 current_diag;
    tupel               next_pair;
    std::string         filename;

    int                 myid;

public:

    AllPairs(
        int  rep        = 50,
        bool make_sym   = false,
        int  ret_node   = 0
    ):
        repeats(rep),
        make_symmetric(make_sym),
        return_node(ret_node),
        myid(dash::myid())
    {
        auto pattern = createPattern();
        results = array_t();
        results.allocate(pattern);

        filename = "all-pairs-result.hdf5";
    }

    template<
        typename APKernel>
    void runKernel(APKernel &kernel)
    {
        current_diag = 0;
        int               ndiags   = dash::size();
        Timer             timer;
        double            measurestart = timer.Now();
        double            kernelstart  = timer.Now();
        double            elapsed;
        bool              is_inverted = false;

        if(myid == 0){
          std::cout << "== Running Kernel '" << kernel.getName()
                    << "' ==" << std::endl;
        }

        // init kernel
        kernel.init(repeats);
        timer.Calibrate(0);

        // clear results
        double no_measure_value = -1;
        dash::fill(results.begin(), results.end(), no_measure_value);

        for(int round = 0; round<2; ++round){
          // Execute only first round if make_symmetric
          if(make_symmetric){
            round = 2;
          }

          // calculate first pair
          updatePartUnits();

        while(current_diag <= ndiags) {
            int x,y;

            if(!is_inverted){
              x = next_pair.first;
              y = next_pair.second;
            } else {
              y = next_pair.first;
              x = next_pair.second;

              if(x == y){
                updatePartUnits();
                 continue; // Skip reflexive in second round

              }
            }

            dash::barrier();
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
                }
            }

            updatePartUnits();
        }
        if(make_symmetric){
          continue;
        } else {
          // reset test for second half of nton plane
          kernel.reset();
          is_inverted = true;
          current_diag = 0;
        }
        }

        // Store results
        dio::HDF5OutputStream os(this->filename,
                                  dio::HDF5FileOptions::Append);
        os << dio::dataset(kernel.getName())
           << results;
        if(myid == 0){
          double kernElapsed = timer.ElapsedSince(kernelstart) / 1000000; // Sec
          std::cout << "== done in " << kernElapsed
                    << " seconds ==" << std::endl;
        }
    }

private:

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

};

#endif // ALL_PAIRS_H
