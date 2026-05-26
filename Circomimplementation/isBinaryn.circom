template IsBinary(n) {

  // array of n inputs
  signal input in[n];
    
  // n loops: n constraints
  for (var i = 0; i < n; i++) {
    in[i] * (in[i] - 1) === 0;
  }
}

// instantiated w/ 4 inputs & 4 constraints
component main = IsBinary(4);
