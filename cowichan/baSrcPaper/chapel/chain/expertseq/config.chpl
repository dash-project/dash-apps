module Config {
  config const nelts = read (int),
               randmat_seed = read (uint(32)),
               thresh_percent = read (int),
               winnow_nelts = read (int);

  const neltSpace = {1..nelts};
  const randSpace = {1..nelts, 1..nelts};
  const histSpace = {0..99};
  const pointSpace = {1 ..winnow_nelts};
  const distSpace = {1..winnow_nelts, 1..winnow_nelts};
  const vectorSpace = {1 .. winnow_nelts};

  var matrix: [randSpace] int;
  var mask: [randSpace] bool;
  var points: [pointSpace] (int, int);
  var dists: [distSpace] real;
  var vector, result: [vectorSpace] real;
}
