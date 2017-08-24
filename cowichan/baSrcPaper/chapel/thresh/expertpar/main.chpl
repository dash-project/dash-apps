/* thresh: histogram thresholding
 * 
 * input:
 *   matrix: the integer matrix to be thresholded
 *   nrows, ncols: the number of rows and columns
 *   percent: the percentage of cells to retain
 *
 * output:
 *   mask: a boolean matrix filled with true for the cells to be kept
 */

config const is_bench = false;
config const nrows = read(int),
             ncols = read(int);
             
require "time.h", "stdio.h";

const ProbSpace = {1..nrows, 1..ncols},
      HistSpace = {1..nrows, 0..100};
const RowSpace = {1..nrows};
const ColSpace = {1..ncols};

var matrix: [ProbSpace] int; 
var mask: [ProbSpace] int;
var histogram: [HistSpace] int;

proc thresh(nrows: int, ncols: int, percent: int) {
  var nmax = max reduce matrix;

  forall i in RowSpace {
    for j in 1..ncols {
      histogram[i, matrix[i, j]] += 1;
    }
  }

  const RowSpace2 = {2..nrows};

  forall j in 0..(nmax) {
    for i in RowSpace2 {
      histogram[1, j] += histogram[i, j];
    }
  }

  var count: int = (nrows * ncols * percent) / 100;
  // writeln("thresh:",count);

  var prefixsum: int = 0;
  var threshold: int = nmax;

  for i in 0..100 {
    if (prefixsum > count) then break;
    prefixsum += histogram[1,100 - i];
    threshold = 100 - i ;
  }
 
  forall (i, j) in ProbSpace {
    mask[i, j] = matrix[i, j] >= threshold;
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
  
  var percent: int;

  if (!is_bench) {
    for i in 1..nrows {
      for j in 1..ncols {
        read(matrix[i,j]);
      }
    }
  }

  read(percent);
   
  
  if( clock_gettime( CLOCK_MONOTONIC_RAW, start) == -1 ) {
    writeln("an error in clock_gettime start has occured!");
  }

  thresh(nrows, ncols, percent);
  
  if( clock_gettime( CLOCK_MONOTONIC_RAW, stop) == -1 ) {
    writeln("an error in clock_gettime stop has occured!");
  }
  
  accum = ( stop.tv_sec - start.tv_sec ) + ( stop.tv_nsec - start.tv_nsec ) / 1e9;
  
  FILE_PTR = fopen("./measurements.txt", "a");
  
  if( is_c_nil(FILE_PTR) ){writeln("File opening for benchmark results failed");}
  
  // Lang, Problem, rows, cols, thresh, winnow_nelts, jobs, time
  fprintf( FILE_PTR, "Chapel,Thresh,%u, %u, %u, , %u, %.9lf,isBench:%d\n", nrows, ncols, percent, dataParTasksPerLocale, is_bench, accum ); //, locale.totalThreads()
  fclose ( FILE_PTR );

  if (!is_bench) {
    //writeln(nrows, " ", ncols);

    for i in 1..nrows {
      for j in 1..ncols {
        if (mask[i, j]) {
          write("1", " ");
        } else {
          write("0", " ");
        }
      }
      writeln();
    }
    writeln();
  }
}
