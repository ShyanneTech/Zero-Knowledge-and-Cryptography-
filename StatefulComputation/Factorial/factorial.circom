pragma circom 2.1.8;
include "circomlib/circuits/comparators.circom";
include "circomlib/circuits/multiplexer.circom";


template factorial(n) {
    signal input in;
    signal output out;

    component lessthan = LessThan(252);
    lessthan.in[0]<== in;
    lessthan.in[1]<== n;
    lessthan.out === 1;

    // formula for computing factorial till n 
    signal factorial[n];
    factorial[0] <== 1;
    for (var i = 1; i < n; i++) {
        factorial[i] <== factorial[i-1] * i;  
    }

    component selectInFactorial = Multiplexer(1,n);
    selectInFactorial.sel <== in;
    for (var i = 0; i < n; i++) {
        selectInFactorial.inp[i][0]<== factorial[i];
    }
    

    out <== selectInFactorial.out[0];
    
}

component main = factorial(500);