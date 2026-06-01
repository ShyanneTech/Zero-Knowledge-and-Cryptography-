pragma circom 2.0.0;
include "circomlib/circuits/comparators.circom";
// k is greater than at least 2 of x, y, or z

template AtLeastGreaterThanTwo() {
    signal input k;
    signal input x;
    signal input y;
    signal input z;

    signal GreaterThanX;
    signal GreaterThanY;
    signal GreaterThanZ;

    signal totalGreaterThan;

    GreaterThanX <== GreaterThan(252)([k,x]);
    GreaterThanY <== GreaterThan(252)([k,y]);
    GreaterThanZ <== GreaterThan(252)([k,z]);

    totalGreaterThan <== GreaterThanX + GreaterThanY + GreaterThanZ;

    signal atLeastTwo; 

    atLeastTwo <== GreaterEqThan(252)([totalGreaterThan,2]);

    atLeastTwo === 1;
    
}

component main = AtLeastGreaterThanTwo();