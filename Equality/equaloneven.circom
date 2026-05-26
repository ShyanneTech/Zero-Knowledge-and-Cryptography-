pragma circom 2.0.0;

template EqualOnEven(n){
    signal input in1[n];
    signal input in2[n];

    for(var i = 0; i<n ;i++){
        if(i%2 == 0){
            in1[i] === in2[i];
        }
    }

}

component main = EqualOnEven(4);