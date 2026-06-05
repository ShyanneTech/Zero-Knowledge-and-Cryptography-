pragma circom 2.1.8;
include "circomlib/circuits/comparators.circom";

template Branching(cond1,cond2,cond3,br1,br2,br3,br4) {
    signal input x;
    signal output out;

    // CIRCOM SAYS HEHEEEE LEGGOOOO 😂
    signal switch1;
    signal switch2;
    signal switch3;
    signal otherwise;

    switch1 <== IsEqual()([x,cond1]);
    switch2 <== IsEqual()([x,cond2]);
    switch3 <== IsEqual()([x,cond3]);
    otherwise <== IsZero()( cond1 + cond2 + cond3 );

    cond1 + cond2 + cond3 + otherwise === 1;
    
    out <==  cond1 * br1 + cond2 * br2 + cond3 * br3 + otherwise * br4 ;
    
}

template MultiBranchConditions() {
    signal input x;
    signal output out;

    component branch = Branching(5,9,10,14,22,23,45);

    branch.x <== x;
    out <== branch.out ;

    
}

component main = MultiBranchConditions();