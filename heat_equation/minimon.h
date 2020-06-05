#ifndef MINIMONITORING_H
#define MINIMONITORING_H

#include <iostream>
#include <fstream>
#include <vector>
#include <stack>
#include <chrono>
#include <map>

using time_point_t = std::chrono::time_point<std::chrono::high_resolution_clock>;
using time_diff_t  = std::chrono::duration<double>;

struct MiniMonValue {

    time_diff_t runtime_sum;
    time_diff_t runtime_min;
    time_diff_t runtime_max;
    uint32_t num;
 //std::numeric_limits<int>::max()
    MiniMonValue( ) : runtime_sum(0.0), runtime_min(std::numeric_limits<double>::max()), runtime_max(0.0), num(0) {}

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

   void enter() {

      _entries.push(std::chrono::high_resolution_clock::now());
   }

   void leave( const char* n) {

      auto& top = _entries.top();
      _store[ n ].apply( std::chrono::high_resolution_clock::now() - top );
      _entries.pop();
   }
#if 0
   void print(std::string filename, uint32_t id) {
      /* print out log to individual files */

      {
      std::ofstream file;
      std::ostringstream file_name;
      file_name << "filename_" << std::setw(2) << std::setfill('0') << id << ".csv";
      file.open(file_name.str());

      file << "# function_name;;num_calls;avg_runtime;min_runtime;max_runtime"
           << std::endl;

      for( auto& e : _store) {
         file << e.first << ";" <<
            e.second.num << ";" <<
            e.second.runtime_sum.count() / e.second.num << ";" <<
            e.second.runtime_min.count() << ";" <<
            e.second.runtime_max.count() << std::endl;
      }
      file.close();
      }
   }
#endif
   void print(uint32_t id) {
      /* print out log to individual files */

      for( auto& e : _store) {
        std::cout << id << ";" << e.first << ";" <<
            e.second.num << ";" <<
            e.second.runtime_sum.count() / e.second.num << ";" <<
            e.second.runtime_min.count() << ";" <<
            e.second.runtime_max.count() << std::endl;
      }
   }


private:
   std::map<const char*,MiniMonValue> _store;
   std::stack<time_point_t, std::vector<time_point_t>> _entries;
};

#endif /* MINIMONITORING_H */
