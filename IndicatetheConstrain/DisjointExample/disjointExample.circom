// x is less than 5 or x is greater than 17

pragma circom 2.0.0;
include "circomlib/circuits/comparators.circom";

include "circomlib/circuits/gates.circom";


template DisjointExample() {
    signal input x;

    signal indicator1;
    signal indicator2;

    // using the indicator for x <5 

    indicator1 <== LessThan(252)([x,5]);
    indicator2 <== GreaterThan(252)([x,17]);

    // ORing to see final result 

   component or = OR();
    or.a <== indicator1;
    or.b <== indicator2;
    or.out === 1;


    // Another way of simplifying th code 
    component or2 = OR();
    or2.a <== LessThan(252)([x,5]);
    or2.b <== GreaterThan(252)([x,17]);
    or2.out === 1;
   
}

component main = DisjointExample();