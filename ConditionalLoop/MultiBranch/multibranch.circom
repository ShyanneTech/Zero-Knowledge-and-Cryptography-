pragma circom 2.1.8;
include "circomlib/circuits/comparators.circom";

template MultiBranchConditions() {
    signal input x;
    signal output out;

    // CIRCOM SAYS HEHEEEE LEGGOOOO😂

    signal x_eq_5;
    signal x_eq_9;
    signal x_eq_10;
    signal otherwise;

    x_eq_5 <== IsEqual()([x,5]);
    x_eq_9 <== IsEqual()([x,9]);
    x_eq_10 <== IsEqual()([x,10]);
    otherwise <== IsZero()( x_eq_5 +  x_eq_9 +  x_eq_10);

    x_eq_5 +  x_eq_9 +  x_eq_10 + otherwise === 1;
    
    out <==  x_eq_5 * 14 +  x_eq_9 * 22 +  x_eq_10 * 23 + otherwise * 45;
    

}

component main = MultiBranchConditions();