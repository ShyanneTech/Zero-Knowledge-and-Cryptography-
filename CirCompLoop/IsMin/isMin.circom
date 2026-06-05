pragma circom 2.0.0;
include "circomlib/circuits/comparators.circom";

template IsMin(n) {
    signal input in[n];
    signal output out;

    var min = in[0];
    for (var i = 0; i < n; i++) {
        min = min < in[i]? in[i] : min;   
    }

    out <-- min;

    // making sure the supposed min is less than all elements there 

    component lte[n];

    for (var i = 0; i < n; i++) {
        lte[i] = LessEqThan(252);
        lte[i].in[0]<== out;
        lte[i].in[1]<== in[i];
        lte[i].out === 1;
    }


    // confirming it actually belongs to the array of elements 

    component eq[n];
    var sum = 0;

    for (var i = 0; i < n; i++) {
        eq[i] = IsEqual();
        eq[i].in[0] <== out;
        eq[i].in[1] <== in[i];
        sum += eq[i].out  ;
    }

    // if it at least equal to one of the element sum would yield a non zero element and for isz of it would yield zero
    signal isz;
    isz <== IsZero()(sum);
    isz === 0;
    
}

component main = IsMin(5);