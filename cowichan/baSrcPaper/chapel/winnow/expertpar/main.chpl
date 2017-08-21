/* winnow: weighted point selection
 * 
 * input:
 *   matrix: an integer matrix, whose values are used as masses
 *   mask: a boolean matrix showing which points are eligible for
 *     consideration
 *   nrows, ncols: the number of rows and columns
 *   nelts: the number of points to select
 *
 * output:
 *   points: a vector of (x, y) points
 */

use Sort;

config const is_bench = false;
config const nrows = read(int),
             ncols = read(int);

require "time.h", "stdio.h";

const MatrixSpace = {0..nrows-1, 0..ncols-1};

var matrix: [MatrixSpace] int;
var mask: [MatrixSpace] bool;

var values: [0..nrows*ncols] (int, (int, int)); // (value, i, j))

var count_per_line: [1..nrows+1] int;

proc winnow(nelts: int) {
  var n: int = 0;
  var points: [1..nelts] (int, int);
  const ColSpace = {0..ncols-1};

  forall i in 0..nrows-1 {
    count_per_line[i + 2] = 0;
    for j in ColSpace {
      // if (is_bench) {
        // mask[i, j] = (((i - 1) * (j - 1)) % (ncols + 1)) == 1;
      // }
      count_per_line[i + 2] += mask[i, j];
    }
  }
   
  var total = + scan count_per_line;
  n = total[nrows + 1];

  
  forall i in 0..nrows-1 {
    var count = total[i + 1];
    for j in ColSpace {
      if (mask[i, j]) {
        values[count] = (matrix[i, j], (i, j));
        count += 1;
      }
    }
  }
  
  quickSort(values[0..n]);

  var chunk: int = n / nelts;

  forall i in 1..nelts do {
    var ind: int;
    ind = (i - 1) * chunk + 1;
    (_, points[i]) = values[ind];
  }
  
  return points;
}

proc read_matrix() {
  for i in 0..nrows-1 do {
    for j in 0..ncols-1 do {
      read(matrix[i, j]);
    }
  }
}

proc read_mask() {
  for i in 0..nrows-1 do {
    for j in 0..ncols-1 do {
      var v: int;
      read(v);
      mask[i, j] = v == 1;
    }
  }
}

proc FillOnTheFly() {
  
  forall i in 0..nrows-1 {
  var b:bool = true;
  var c:int = 0;
    for j in 0..ncols-1 {
        c+=1;
        if( c == 1 ){
          mask[i, j]   = true;
          c = 0;
        }else{
          mask[i, j]   = false;
        }
        matrix[i, j] = (i*ncols + j) % 100;
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
  extern proc printf( fmt: c_string, vals...?numvals): int;
  extern proc fopen(  strFile: c_string, op: c_string): c_ptr(FILE);
  extern proc fclose( ptr: c_ptr(FILE));
  
  var FILE_PTR : c_ptr(FILE);
  var accum : c_double;
  var start,stop : chpl_timespec;
  
  var nelts: int;

  if (!is_bench) {
    read_matrix();
    read_mask();
  }else{
    FillOnTheFly();
  }

  read(nelts);
  
  var points: [1..nelts] (int, int);
  
  if( clock_gettime( CLOCK_MONOTONIC_RAW, start) == -1 ) {
    writeln("an error in clock_gettime start has occured!");
  }

  points = winnow(nelts);
  
  if( clock_gettime( CLOCK_MONOTONIC_RAW, stop) == -1 ) {
    writeln("an error in clock_gettime stop has occured!");
  }
  
  accum = ( stop.tv_sec - start.tv_sec ) + ( stop.tv_nsec - start.tv_nsec ) / 1e9;
  
  FILE_PTR = fopen("./measurements.txt", "a");
  
  if( is_c_nil(FILE_PTR) ){writeln("File opening for benchmark results failed");}
  
  // Lang, Problem, rows, cols, thresh, winnow_nelts, jobs, time
  fprintf( FILE_PTR, "Chapel,Winnow,%u, %u, , %u, %u, %.9lf,isBench:%d\n", nrows, ncols, nelts, dataParTasksPerLocale, is_bench, accum ); //, locale.totalThreads()
  fclose ( FILE_PTR );

  if (!is_bench) {
    writeln(nelts);

    for i in 1..nelts do {
      writeln(points[i](1)," ", points[i](2));
    }

    writeln();
  }
}
