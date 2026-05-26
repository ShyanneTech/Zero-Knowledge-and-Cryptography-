pragma circom 2.0.0;

// Create constraints that enforces all signals
// in `in` are binary, i.e. 0 or 1.

template AllBinary(n) {
    signal input in[n];

    for(var i = 0; i< n; i++){
        in[i] * in[i] === in[i];
    }

}

component main = AllBinary(4);