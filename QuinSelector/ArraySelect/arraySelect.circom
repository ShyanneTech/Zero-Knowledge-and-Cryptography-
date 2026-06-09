pragma circom 2.1.8;

include "circomlib/circuits/comparators.circom";


// i want the output to be based on the array and index provided

template ArraySelect(n) {
    signal input in[n];
    signal input index;
    signal output out;

    signal prod[n]; // allows the index i don't need to 0 and the one i need to be the index itself

    component neededIndex[n];

    for (var i = 0; i < n; i++) {
        neededIndex[i]= IsEqual();
        neededIndex[i].in[0] <== i ;
        neededIndex[i].in[1] <== index;
        
        prod[i] <==  neededIndex[i].out * in[i];
    }

    var sum = 0;
    for (var i = 0; i < n; i++) {
        sum += prod[i];  
    }

    out <== sum;

}

component main = ArraySelect(5);