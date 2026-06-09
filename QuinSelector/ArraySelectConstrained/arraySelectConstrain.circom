pragma circom 2.1.8;

// I want to constrain the index in such a way that 
include "circomlib/circuits/comparators.circom";

template CalcTotal(n) {
    signal input in[n];
    signal output out;

    signal sums[n];

    sums[0] <== in[0];

    for (var i = 1; i < n; i++) {
        sums[i] <== sums[i-1] + in[i];
    }

    out<== sums[n-1];
    
}

template ArraySelectConstrain(n) {
    signal input in[n];
    signal input index;
    signal output out;

    component lessThan = LessThan(252);
    lessThan.in[0] <== index;
    lessThan.in[1] <== n;

    component neededIndex[n];
    component calcTotal = CalcTotal(n);

    for (var i = 0; i < n; i++) {
        neededIndex[i]= IsEqual();
        neededIndex[i].in[0] <== i ;
        neededIndex[i].in[1] <== index;
        
        calcTotal.in[i] <==  neededIndex[i].out * in[i];
    }


    out <== calcTotal.out;

}

component main = ArraySelectConstrain(7);