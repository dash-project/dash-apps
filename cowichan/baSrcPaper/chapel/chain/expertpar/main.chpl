/* chain: chain all problems
 * 
 * input:
 *   nelts: the number of elements
 *   randmat_seed: random number generator seed
 *   thresh_percent: percentage of cells to retain
 *   winnow_nelts: the number of points to select
 *
 * output:
 *   result: a real vector, whose values are the result of the final product
 */

use Config, Randmat, Thresh, Winnow, Outer, Product;
config const is_bench = false;
require "time.h", "stdio.h";
proc main() {
  extern "struct _IO_FILE" record FILE{}
  extern "struct timespec" record chpl_timespec{
    var tv_sec : c_long;
    var tv_nsec: c_long;
  }
  extern const CLOCK_MONOTONIC_RAW: c_long;
  extern proc clock_gettime(clk_id : c_long, ref structPtr : chpl_timespec) :int;
  extern proc fprintf( fPtr: c_ptr(FILE), fmt: c_string, vals...?numvals): int;
  extern proc printf( fmt: c_string, vals...?numvals): int;
  extern proc fopen(  strFile: c_string, op: c_string): c_ptr(FILE);
  extern proc fclose( ptr: c_ptr(FILE));
  
  var FILE_PTR : c_ptr(FILE);
  var accum : c_double;
  var start,stop : chpl_timespec;
  
  if( clock_gettime( CLOCK_MONOTONIC_RAW, start) == -1 ) {
    writeln("an error in clock_gettime start has occured!");
  }
  
  randmat(nelts, nelts, randmat_seed);
  thresh(nelts, nelts, thresh_percent);
  winnow(nelts, nelts, winnow_nelts);
  outer(winnow_nelts);
  product(winnow_nelts);

  if( clock_gettime( CLOCK_MONOTONIC_RAW, stop) == -1 ) {
    writeln("an error in clock_gettime stop has occured!");
  }
  
  accum = ( stop.tv_sec - start.tv_sec ) + ( stop.tv_nsec - start.tv_nsec ) / 1e9;
  
  if(is_bench){
    FILE_PTR = fopen("./measurements.txt", "a");
    
    if( is_c_nil(FILE_PTR) ){writeln("File opening for benchmark results failed");}
    
    // Lang, Problem, rows, cols, thresh, winnow_nelts, jobs, time
    fprintf( FILE_PTR, "Chapel,Chain, %u, %u, %u, %u, %u, %.9lf\n", nelts, nelts, thresh_percent, winnow_nelts, dataParTasksPerLocale, accum ); //, locale.totalThreads()
    fclose ( FILE_PTR );
  }
  
  if (!is_bench) {
    writeln(winnow_nelts);
    for i in 1..winnow_nelts do {
      writef("%.4dr ", result[i]);
    }
    writeln();
  }
}
