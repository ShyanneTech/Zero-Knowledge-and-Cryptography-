pragma circom 2.0.0;

// Create a circuit that takes an array of four signals
// `in`and a signal s and returns is satisfied if `in`
// is the binary representation of `n`. For example:
// 
// Accept:
// 0,  [0,0,0,0]
// 1,  [1,0,0,0]
// 15, [1,1,1,1]
// 
// Reject:
// 0, [3,0,0,0]
// 
// The circuit is unsatisfiable if n > 15

// n = 8b3 + 4b2+ 2b1+ b0
template FourBitBinary() {
    signal input in[4];
    signal input n;

    assert(n < 16 )

    for ( var i = 0; i< 4;i++){
        0 === in[i] * (in[i] - 1);

    }
    n === 8 * in[3]+4 * in[2]+ 2 * in[1] + 1 * in[0]
;
}

component main{public [n]} = FourBitBinary();