pragma circom 2.0.0;

include "circomlib/circuits/pointbits.circom";
include "circomlib/circuits/comparators.circom";

// Write a Circom function that finds the root of a degree 2 polynomial using the quadratic formula. Remember, everything is done over a finite field, so you need to use the modular square root from the first example.

// Then, write constraints that the two roots (if they exist) satisfy the polynomial. Pass in the polynomial to the Circom template as an array of three coefficients.

template FindQuadRoots() {
    signal input in[3];
    signal output out1;
    signal output out2;



    signal a <== in[2];
    signal b <== in[1];
    signal c <== in[0];



    signal b2;
    b2 <== b * b;

    signal ac;
    ac <== 4 * a * c ;

    signal disc ;

    disc <== b2 - ac;

    signal sideroot1;
    signal sideroot2;

    sideroot1 <-- sqrt(disc);
    sideroot2 <-- sideroot1 * -1;

    sideroot1 * sideroot1 === disc;
    sideroot2 * sideroot2 === disc;

    signal den ;

    den <== 2 *a; 

    signal invden;

    invden <-- den!= 0 ? 1/den : 0;

    
    out1 <-- (-b + sideroot1) * invden;
    out2 <-- (-b + sideroot2) * invden;

    // now constraing it tp ensure i am right 

    // Confirming putting the inputs in 

    signal x1_2;
    signal x2_2;

    x1_2 <== out1 * out1;
    x2_2 <== out2 * out2;

    signal sum1;

    signal ax1_2;
    ax1_2 <== a * x1_2 ;

    sum1<== ax1_2 + b * out1 + c;

    signal sum2;
    signal ax2_2;
    ax2_2 <== a * x2_2 ;

    sum2<== ax2_2 + b * out2 + c;

    signal isz1;
    signal isz2;
     

    isz1 <== IsZero()(sum1);
    isz2 <== IsZero()(sum2);

    isz1 === 1;
    isz2 === 1;




  




    
}

component main = FindQuadRoots();