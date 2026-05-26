pragma circom 2.0.0;

// Goal is to check that A sum of something equals a said sum 

template  Sum(n) {
    signal input in[n];
    signal input sum;

    var acc;
    for (var i = 0; i < n; i++) {
        acc += in[i];
    }

    acc === sum;

    
}

component main = Sum(5);