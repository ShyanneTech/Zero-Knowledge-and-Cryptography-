pragma circom 2.0.0;

// Constraining 1000 variables to {0,1}
template Constrain1000Example() {
  signal input in[1000];
  
  for (var i = 0; i < 1000; i++) {
    0 === in[i] * (in[i] - 1);
  }
}

component main = Constrain1000Example();
