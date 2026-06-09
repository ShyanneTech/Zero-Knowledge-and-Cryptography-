pragma circom 2.1.8;

include "circomlib/circuits/comparators.circom";
include "circomlib/circuits/multiplexer.circom";


// The idea is to be able to swap two elemnts based on their positions

template Swap(n) {
    signal input in[n];
    signal input s;
    signal input t;
    signal output out[n];

    // it is possible s = t ... User wahala ;

    signal sEqt;
    sEqt <== IsEqual()([s,t]);

    component positionS = Multiplexer(1,n);
    positionS.sel <== s;
    for (var i = 0; i < n; i++) {
       positionS.inp[i][0] <== in[i]; // output here would give s points to 6
    }

    component positionT = Multiplexer(1,n);
    positionT.sel <== t;
    for (var i = 0; i < n; i++) {
       positionT.inp[i][0] <== in[i]; // output here would give t points to 7
    }

    component ArrayS[n];
    component ArrayT[n];
    component NotST[n];


    for (var i = 0; i < n; i++) {
        ArrayS[i] = IsEqual();
        ArrayS[i].in[0] <== i;
        ArrayS[i].in[1] <== s;
        ArrayT[i] = IsEqual();
        ArrayT[i].in[0] <== i;
        ArrayT[i].in[1] <== t;
        NotST[i] = IsZero();
        NotST[i].in <== ArrayS[i].out + ArrayT[i].out;


    }

    signal BranchS[n];
    signal BranchT[n];
    signal NotBranchST[n];

    for (var i = 0; i < n; i++) {
        BranchS[i] <== ArrayS[i].out * positionT.out[0];
        BranchT[i] <== ArrayT[i].out * positionS.out[0];
        NotBranchST[i] <== NotST[i].out * in[i];

        out[i] <==  (1- sEqt) * BranchS[i]  + BranchT[i]  + NotBranchST[i]; 
        
    }
    
}

component main = Swap(9);