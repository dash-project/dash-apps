/* outer: outer product
 * 
 * input:
 *   points: a vector of (x, y) points
 *   nelts: the number of points
 *
 * output:
 *   matrix: a real matrix, whose values are filled with inter-point
 *     distances
 *   vector: a real vector, whose values are filled with origin-to-point
 *     distances
 */

config const is_bench = false;
config const nelts = read(int);

require "time.h", "stdio.h";

var matrix: [1..nelts, 1..nelts]real;
var vector: [1..nelts]real;
var points: [1..nelts](int, int);

inline
proc sqr(x: real): real {
  return x ** 2;
}

inline
proc distance(l, r: (int, int)): real {
  var lx, ly, rx, ry: real;
  (lx, ly) = l;
  (rx, ry) = r;
  return sqrt(sqr(lx - rx) + sqr(ly - ry));
}

proc outer(nelts: int) {
  const NeltSpace = {1..nelts};
  forall i in NeltSpace {
    var nmax: real = 0;
    for j in NeltSpace {
      if (i != j) {
        matrix[i, j] = distance(points[i], points[j]);
        nmax = max(nmax, matrix[i, j]);
      }
    }
    matrix[i, i] = nmax * nelts;
    vector[i] = distance((0, 0), points[i]);
  }
}

proc read_vector_of_points(nelts: int) {
  var a, b: int;
  for i in 1..nelts do {
    read(a, b);
    points[i] = (a, b);
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
    read_vector_of_points(nelts);
  } 

  if( clock_gettime( CLOCK_MONOTONIC_RAW, start) == -1 ) {
    writeln("an error in clock_gettime start has occured!");
  }
  
  outer(nelts);
  
  if( clock_gettime( CLOCK_MONOTONIC_RAW, stop) == -1 ) {
    writeln("an error in clock_gettime stop has occured!");
  }
  
  accum = ( stop.tv_sec - start.tv_sec ) + ( stop.tv_nsec - start.tv_nsec ) / 1e9;
  
  FILE_PTR = fopen("./measurements.txt", "a");
  
  if( is_c_nil(FILE_PTR) ){writeln("File opening for benchmark results failed");}
  
  // Lang, Problem, rows, cols, thresh, winnow_nelts, jobs, time
  fprintf( FILE_PTR, "Chapel,Outer  ,     ,     ,   ,%5u,%2u,%.9lf,isBench:%d\n", nelts, dataParTasksPerLocale, is_bench, accum ); //, locale.totalThreads()
  fclose ( FILE_PTR );

  if (!is_bench) {
    writeln(nelts);
    for i in 1..nelts do {
      for j in 1..nelts do {
        writef("%.4dr ", matrix[i, j]);
      }
      writeln();
    }
    writeln();

    for i in 1..nelts do {
      writef("%.4dr ",vector[i]);
    }
    writeln();
  }
}
