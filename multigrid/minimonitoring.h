#ifndef MINIMONITORING_H
#define MINIMONITORING_H

#include <iostream>
#include <fstream>
#include <vector>
#include <deque>
#include <stack>
#include <chrono>
#include <map>
#include <tuple>

#include <libdash.h>

using time_point_t = std::chrono::time_point<std::chrono::high_resolution_clock>;
using time_diff_t  = std::chrono::duration<double>;

struct MiniMonValue {

    time_diff_t runtime_sum;
    time_diff_t runtime_min;
    time_diff_t runtime_max;
    uint32_t num;
 //std::numeric_limits<int>::max()
    MiniMonValue( ) : runtime_sum(0.0), runtime_min(1.0e300), runtime_max(0.0), num(0) {}

    void apply( time_diff_t value ) {

        runtime_sum += value;
        runtime_min= std::min( runtime_min, value );
        runtime_max= std::max( runtime_max, value );
        num += 1;
    }
};


class MiniMon {

public:

   MiniMon() = default;

   void start() {

      _entries.push(std::chrono::high_resolution_clock::now());
   }

   void stop( const char* n, uint32_t p, uint64_t e = 1,
         uint64_t f = 0, uint64_t r = 0, uint64_t w = 0 ) {

      auto& top = _entries.top();
      _store[ {n,p,e,f} ].apply( std::chrono::high_resolution_clock::now() - top );
      _entries.pop();
   }

   void print(uint32_t id) {
      /* print out log to individual files */

      {
      std::ofstream file;
      std::ostringstream file_name;
      file_name << "overview_" << std::setw(5) << std::setfill('0') << id << ".csv";
      file.open(file_name.str());

      file << "# function_name;par;elements;flops;num_calls;avg_runtime;min_runtime;max_runtime"
           << std::endl;

      for( auto& e : _store) {
         file << std::get<0>(e.first) << ";" << 
            std::get<1>(e.first) << ";" << 
            std::get<2>(e.first) << ";" <<
            std::get<3>(e.first) << ";" <<
            e.second.num << ";" <<
            e.second.runtime_sum.count() / e.second.num << ";" <<
            e.second.runtime_min.count() << ";" <<
            e.second.runtime_max.count() << std::endl;
      }
      file.close();
      }
   }



private:
   std::map<std::tuple<const char*,uint32_t,uint64_t,uint64_t>,MiniMonValue> _store;
   std::stack<time_point_t, std::vector<time_point_t>> _entries;
};

#endif /* MINIMONITORING_H */
