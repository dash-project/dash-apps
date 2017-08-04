/* randmat: random number generation
 * 
 * input:
 *   nrows, ncols: the number of rows and columns
 *   s: the seed
 *
 * output:
 *   matrix: an nrows x ncols integer matrix
 */

use Random;

config const is_bench = false;
config const nrows = read (int(32)),
             ncols = read (int(32)),
             s     = read (int(32));
             
require "time.h", "stdio.h";

var matrix: [1..nrows, 1..ncols] int(32);

proc randmat() {
  const LCG_A: uint(32) = 1664525,
        LCG_C: uint(32) = 1013904223;
  const cols = 1 .. ncols;
  const rows = 1 .. nrows;
  forall i in rows {
    var seed : uint(32) = (s + i - 1) : uint (32);
    for j in cols {
      seed = LCG_A * seed + LCG_C;
      matrix[i, j] = (seed % 100) : int(32);
    }
  }
}

proc main() {
  extern "struct _IO_FILE" record FILE{}
  extern "struct timespec" record chpl_timespec{
    var tv_sec : c_long;
    var tv_nsec: c_long;
  }
  extern const CLOCK_MONOTONIC_RAW: c_long;
  extern proc clock_gettime(clk_id : c_long, ref structPtr : chpl_timespec) :int;
  extern proc fprintf( fPtr: c_ptr(FILE), fmt: c_string, vals...?numvals): int;
  extern proc fopen(  strFile: c_string, op: c_string): c_ptr(FILE);
  extern proc fclose( ptr: c_ptr(FILE));
  
  var FILE_PTR : c_ptr(FILE);
  var accum : c_double;
  var start,stop : chpl_timespec;
  
  if( clock_gettime( CLOCK_MONOTONIC_RAW, start) == -1 ) {
    writeln("an error in clock_gettime start has occured!");
  }
  
  randmat();

  if( clock_gettime( CLOCK_MONOTONIC_RAW, stop) == -1 ) {
    writeln("an error in clock_gettime stop has occured!");
  }
  
  accum = ( stop.tv_sec - start.tv_sec ) + ( stop.tv_nsec - start.tv_nsec ) / 1e9;
  
  if(is_bench){
    FILE_PTR = fopen("./measurements.txt", "a");
    
    if( is_c_nil(FILE_PTR) ){writeln("File opening for benchmark results failed");}
    
    // Lang, Problem, rows, cols, thresh, winnow_nelts, jobs, time
    fprintf( FILE_PTR, "Chapel,Randmat,%u, %u, , , %u, %.9lf\n", nrows, ncols, dataParTasksPerLocale, accum ); //, locale.totalThreads()
    fclose ( FILE_PTR );
  }
  
  if (!is_bench) {
    // writeln(nrows, " ", ncols);

    for i in 1..nrows do {
      for j in 1..ncols do {
        write(matrix[i, j], " ");
      }
      writeln();
    }
    writeln(); 
    
  }
    // writeln("numlocales:", numLocales, " runningTasks:", here.runningTasks(), " totalThreads():", here.totalThreads() );
    // writeln("queuedTasks():", here.queuedTasks(), " blockedTasks():", here.blockedTasks(), " idleThreads():", here.idleThreads());
    // writeln("numLocales:", numLocales, "  numPUs:", here.numPUs(true,true), " maxTaskPar:", here.maxTaskPar);
    // writeln("dataParTasksPerLocale:",dataParTasksPerLocale );
}
