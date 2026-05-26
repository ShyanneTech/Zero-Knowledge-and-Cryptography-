pragma circom 2.0.0;

// Square template
template Square() {
    signal input in;
    signal output out;

    out <== in * in;

    
}

// Main component 

template Main() {
    signal input a;
    signal input b;
    signal input sumofsquares;

    component a2 = Square();
    component b2 = Square();

    a2.in <== a;
    b2.in <== b;

    a2.out + b2.out === sumofsquares;

    
}

component main = Main();