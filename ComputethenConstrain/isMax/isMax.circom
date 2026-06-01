pragma circom 2.0.0;
include "circomlib/circuits/comparators.circom";

template isMax(n) {
    signal input in[n];
    signal output out;

    // geeting the assumed max 

    var maxx = in[0];
    for (var i = 0; i < n; i++) {
        if(in[i] > maxx){
            maxx = in[i];
        }
    }

    out <-- maxx;

    // Condition 1
    signal gte[n];

    for (var i = 0; i < n; i++) {
        gte[i] <== GreaterEqThan(252)([out, in[i]]);
        gte[i] === 1;
    }


    // condition 2 

    signal eq[n];
    var sum = 0;

    for (var i = 0; i < n; i++) {
        eq[i] <== IsEqual()([out, in[i]]);
        sum += eq[i];
    }
    signal iz;
    iz <== IsZero()(sum);

    iz === 0;
    
}
component main = isMax(5);