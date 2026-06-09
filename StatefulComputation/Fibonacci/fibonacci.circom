pragma circom 2.1.8;
include "circomlib/circuits/comparators.circom";
include "circomlib/circuits/multiplexer.circom";

template fibonacci(n) {
    assert(n >= 2);

    signal input in;
    signal output out;

    component lessThan = LessThan(252);
    lessThan.in[0] <== in;
    lessThan.in[1] <== n;
    lessThan.out === 1;

    // computing how to get the fibonacci sequence 

    signal fibonacci[n];

    fibonacci[0] <== 1;
    fibonacci[1] <== 1;

    for (var i = 2; i < n; i++) {
        fibonacci[i] <== fibonacci[i-1] + fibonacci[i-2];
        
    }

    component fiboseq = Multiplexer(1,n);
    fiboseq.sel <== in;
    for (var i = 0; i < n; i++) {
        fiboseq.inp[i][0] <== fibonacci[i];

    }

    out <== fiboseq.out[0]; 
  
}

component main = fibonacci(100);