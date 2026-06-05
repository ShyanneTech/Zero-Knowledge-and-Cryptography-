pragma circom 2.0.0;
include "circomlib/circuits/comparators.circom";

template isSorted(n) { // the goal is to confirm that the array is sorted 
    signal input in[n];

    component lt[n-1];
    
    for (var i = 0; i < n-1; i++) {
        lt[i] = LessThan(252);
        lt[i].in[0] <== in[i];
        lt[i].in[1] <== in[i+1];

        lt[i].out === 1;
    }
       
    
}

component main = isSorted(9);
