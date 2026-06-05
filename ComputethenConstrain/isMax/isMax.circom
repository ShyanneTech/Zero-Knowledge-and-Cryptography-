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
    component  gte[n]; // when you wan

    for (var i = 0; i < n; i++) {
        // gte[i] <== GreaterEqThan(252)([out, in[i]]);
        // gte[i] === 1; in here i can see i didn't constrain the output and that is bad
        gte[i] <== GreaterEqThan(252);
        gte[i].in[0]<== out;
        gte[i].in[1]<== in[i];
        gte[i].out === 1;


    }


    // condition 2 

    component  eq[n];
    var sum = 0;

    for (var i = 0; i < n; i++) {
        eq[i] <== IsEqual()
        eq[i].in[0]<== out;
        eq[i].in[1]<== in[i];
        sum += eq[i].out;
    }
    signal iz;
    iz <== IsZero()(sum);

    iz === 0;
    
}
component main = isMax(5);