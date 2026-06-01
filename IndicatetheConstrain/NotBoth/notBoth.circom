pragma circom 2.1.0;
include "circomlib/circuits/comparators.circom";
include "circomlib/circuits/gates.circom";

template NotBoth() {
    signal input x;
    signal input y;
 

    component nand = NAND();
    nand.a <== LessThan(252)([x,100]); 
    nand.b <== LessThan(252)([y,100]) ;

    nand.out === 1;
    
}

component main = NotBoth();