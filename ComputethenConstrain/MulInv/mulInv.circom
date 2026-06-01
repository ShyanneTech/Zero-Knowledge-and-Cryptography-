pragma circom 2.0.0;

template MulInv() {
    signal input in;
    signal output out;

    var inv = in ** (-2);
    out <-- inv;

    //  out <-- 1 / in;

    out * in === 1;

    
}

component main = MulInv();