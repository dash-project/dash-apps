#ifndef MINIMONITORING_H
#define MINIMONITORING_H

#include <iostream>
#include <fstream>
#include <vector>
#include <deque>
#include <stack>
#include <chrono>

#include <libdash.h>

using time_point_t = std::chrono::time_point<std::chrono::high_resolution_clock>;
using time_diff_t  = std::chrono::duration<double>;

const size_t NAMELEN= 30;

struct Entry {
   time_point_t time;
   const char* name;
   uint32_t par; /* number of parallel units working on this, even though the following numbers are per process */
   uint64_t elements; /* number of grid elements */
   uint64_t flops; /* number of floating point operations */
   uint64_t loads; /* number of elements (double values) read */
   uint64_t stores; /* number of elements (double values) written */


};

struct ResultEntry {
   time_diff_t diff;
   const char* name;
   uint32_t    par; /* number of parallel units working on this, even though the following numbers are per process */
   uint64_t    elements; /* number of grid elements */
   uint64_t    flops; /* number of floating point operations */
   uint64_t    loads; /* number of elements (double values) read */
   uint64_t    stores; /* number of elements (double values) written */

};

class MiniMon {

public:

   MiniMon() = default;

   void start() {

      _entries.push(std::chrono::high_resolution_clock::now());
   }

   void stop( const char n[NAMELEN], uint32_t p, uint64_t e = 1,
         uint64_t f = 0, uint64_t r = 0, uint64_t w = 0 ) {

      auto& top = _entries.top();
      _results.push_back({std::chrono::high_resolution_clock::now() - top,
          n, p, e, f, r, w });
      _entries.pop();
   }

   void print(uint32_t id) {
      /* print out log to individual files */

      std::ofstream file;
      std::ostringstream file_name;
      file_name << "trace_" << std::setw(5) << std::setfill('0') << id << ".csv";
      file.open(file_name.str());

      file << "# Unit;Function_name;Par;Duration;Elements;Flops;Loads;Stores"
           << std::endl;

      for(const auto& result : _results) {
          file << id << ";" << result.name << ";" << result.par << ";"
               << std::setprecision(12) << result.diff.count() << ";"
               << result.elements << ";" << result.flops << ";" << result.loads
               << ";" << result.stores << std::endl;
      }

      file.close();
   }



private:
   std::deque<ResultEntry>                             _results;
   std::stack<time_point_t, std::vector<time_point_t>> _entries;
};

#endif /* MINIMONITORING_H */
