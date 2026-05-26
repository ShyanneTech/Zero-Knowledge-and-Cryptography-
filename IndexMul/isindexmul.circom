pragma circom 2.0.0;

template IsIndexMultiplied(n) {
  signal input in1[n];
  signal input in2[n];
  
  for (var i = 0; i < n; i++) {
    in1[i] * i === in2[i];
  }
}

component main = IsIndexMultiplied(3);

/* INPUT = {"in1": [0,1,2], "in2": [0,1,4]} */

// accept
// in1[] = [0,1,2]
// in2[] = [0,1,4]

// reject
// in1[] = [0,1,2]
// in2[] = [0,0,2]
