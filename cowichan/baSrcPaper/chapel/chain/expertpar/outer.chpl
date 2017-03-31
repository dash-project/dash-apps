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

module Outer {

use Config;

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

proc outer(nelts: int)
{
  const NeltSpace = {1..nelts};
  forall i in NeltSpace {
    var nmax: real = -1;

    for j in NeltSpace {
      if (i != j) {
        dists[i, j] = distance(points[i], points[j]);
        nmax = max(nmax, dists[i, j]);
      }
    }

    dists[i, i] = nmax * nelts;
    vector[i] = distance((1, 1), points[i]);
  }
}
}
