#ifndef ALL_PAIRS_H
#define ALL_PAIRS_H

#include <libdash.h>
#include <queue>
#include "all-pairs-kernel.h"

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
        bool make_sym   = true,
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
        double            elapsed;

        // calculate first pair
        updatePartUnits();

        while(current_diag <= ndiags) {
            int x = next_pair.first;
            int y = next_pair.second;

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
                    // TODO Assert if it is really local
                    results[x][y][r] = elapsed / int_repeats;
                }
            }

            updatePartUnits();
        }

        // kernel.name();
        dio::HDF5OutputStream os(this->filename,
                                  dio::HDF5FileOptions::Append);
        os << dio::dataset(kernel.name())
           << results;
        // Store Kernel Results
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
