pragma circom 2.0.0;

template ModularDiv() {
    signal input in;
    signal output out;

    // compute same as MUlInv by the way 
    out <-- 1 / in;
  
    // then constrain
     out * in === 1;

    
}