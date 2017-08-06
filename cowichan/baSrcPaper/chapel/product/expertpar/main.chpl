/* product: matrix-vector product
 * 
 * input:
 *   nelts: the number of elements
 *   matrix: a real matrix
 *   vector: a real vector
 *
 * output:
 *    result: a real vector, whose values are the result of the product
 */

config const is_bench = false;
config const nelts = read(int);

require "time.h", "stdio.h";

const Space = {1..nelts, 1..nelts};
var matrix: [Space]real;
var vector: [1..nelts]real;
var result: [1..nelts]real;

proc product(nelts: int) {
  const NeltSpace = {1..nelts};
  forall i in NeltSpace {
    var sum: real = 0;
    for j in NeltSpace {
      sum += matrix[i, j] * vector[j];
    }
    result[i] = sum;
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
  extern proc printf( fmt: c_string, vals...?numvals): int;
  extern proc fopen(  strFile: c_string, op: c_string): c_ptr(FILE);
  extern proc fclose( ptr: c_ptr(FILE));
  
  var FILE_PTR : c_ptr(FILE);
  var accum : c_double;
  var start,stop : chpl_timespec;
  if (!is_bench) {
    for i in 1..nelts do {
      for j in 1..nelts do {
        read(matrix[i, j]);
      }
    }
 
    for i in 1..nelts do {
      read(vector[i]);
    }
  }

  if( clock_gettime( CLOCK_MONOTONIC_RAW, start) == -1 ) {
    writeln("an error in clock_gettime start has occured!");
  }
  
  product(nelts);
  
  if( clock_gettime( CLOCK_MONOTONIC_RAW, stop) == -1 ) {
    writeln("an error in clock_gettime stop has occured!");
  }
  
  accum = ( stop.tv_sec - start.tv_sec ) + ( stop.tv_nsec - start.tv_nsec ) / 1e9;
  
  FILE_PTR = fopen("./measurements.txt", "a");
  
  if( is_c_nil(FILE_PTR) ){writeln("File opening for benchmark results failed");}
  
  // Lang, Problem, rows, cols, thresh, winnow_nelts, jobs, time
  fprintf( FILE_PTR, "Chapel,Product, , , , %u, %u, %.9lf,isBench:%d\n", nelts, dataParTasksPerLocale, is_bench, accum ); //, locale.totalThreads()
  fclose ( FILE_PTR );


  if (!is_bench) {
    writeln(nelts);
    for i in 1..nelts do {
      writef("%.4dr ", result[i]);
    }
    writeln();
  }
}
