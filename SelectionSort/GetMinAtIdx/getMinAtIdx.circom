pragma circom 2.1.8;
include "circomlib/circuits/comparators.circom";
include "circomlib/circuits/multiplexer.circom";
// we build a template that proves we correctly identified the index of the minimum value in a sublist

template GetMinAtIdx(n) {
    signal input in[n];
  
    // compute and constrain min and idx
    // to be the min value in the list
    // and the index of the minimum value
    signal output min;
    signal output idx;

    // compute the minimum and its index
    // outside of the constraints 
    // I want to get the minv and idx 

    var minv = in[0];
    var idxv = 0;

    for (var i = 0; i < n; i++) {
        if (in[i] < minv) {
        minv = in[i];
        idxv = i
          
        }
    }
    min <-- minv ;
    idxv <-- i;

    // constrain that min is ≤ all others, this actually our business 
    //   i want to constrain my min to actually be my min

    component lte[n];
    for (var i = 0; i < n; i++) {
        lte[i] = LessEqThan(252);
        lte[i].in[0] <== min;
        lte[i].in[1] <== in[i];
        lte[i].out === 1;

    }
 
    // assert min is really at in[idx]
    component correctMinIdx = Multiplexer(1,n);
    correctMinIdx.sel <== idx;

    for (var i = 0; i < n; i++) {
        correctMinIdx.inp[i][0] <== in[i];
    }
    correctMinIdx.out[0] === min;
  
}
