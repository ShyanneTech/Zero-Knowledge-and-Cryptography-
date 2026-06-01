pragma circom 2.0.0;

include "circomlib/circuits/comparators.circom";

template isEqual() {
    signal input in[2];
    signal output out;

    component isz = IsZero();

    isz.in <== in[0] - in[1];

    isz.out ==> out; 



    
}

component main = isEqual();