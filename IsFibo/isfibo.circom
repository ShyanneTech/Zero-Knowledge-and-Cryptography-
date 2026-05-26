pragma circom 2.0.0;

template isFibo(n){

    assert(n > 1);
    signal input in[n];

    // Creating my variable for fibonacci sequence to test with 

    var CorrectFibo[n];
    CorrectFibo[0] = 0;
    CorrectFibo[1] = 1 ;

    // Creating the logic for the fibonacci sequence 

    for(var i = 2; i < n; i++){
       CorrectFibo[i] = CorrectFibo[i -1 ] + CorrectFibo[i- 2];
    }

    // comparing array in with fibonaaci sequenece 

    for(var i = 0; i < n; i++){
       in[i] === CorrectFibo[i];
    }

}

component main = isFibo(8);