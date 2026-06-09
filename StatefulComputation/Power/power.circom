pragma circom 2.1.4;
include "circomlib/circuits/comparators.circom";
include "circomlib/circuits/multiplexer.circom";

// Create a circuit which takes an input 'a',(array of length 2 ) , then  implement power modulo 
// and return it using output 'c'.

// HINT: Non Quadratic constraints are not allowed. 

template Pow(n) {
   
   // Your Code here.. 

   signal input a[2]; // this contains base and exponent
   signal output c;

   component lessThan = LessThan(252);
   lessThan.in[0] <== a[1];
   lessThan.in[1] <== n;

    // now proceeds to compute the power equation

    signal Powera[n];
    Powera[0] <== 1;
    for (var i = 1; i < n; i++) {
        Powera[i] <== Powera[i-1] * a[0];
    }

    // proceeds to select only that which i need 

    component exactPowerA = Multiplexer(1,n);
    exactPowerA.sel <== a[1];

    for (var i = 0; i < n; i++) {
        exactPowerA.inp[i][0] <== Powera[i];
    }
    c <== exactPowerA.out[0];

}

component main = Pow(50);
