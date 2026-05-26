pragma circom 2.0.0;

template ValidBinRep(n) {
    signal input in[n];
    signal input k;

    // Limiting in to be just binary 

    for (var i = 0; i < n; i++) {
        in[i] * in[i] === in[i];
    }


    var value;
    for (var i = 0; i < n; i++) {
        value += 2 ** i * in[i];
    }


   value === k; 

    
}

component main = ValidBinRep(5);