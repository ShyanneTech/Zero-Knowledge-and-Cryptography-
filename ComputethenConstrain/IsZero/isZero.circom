pragma circom 2.0.0;

template isZero() {
    signal input in;
    signal output out;

    signal inv;
    inv <-- in != 0 ? 1/in : 0;  // also called advice inputs or non deterministic inputs

    out <== -in*inv + 1;

    in * out === 0;

    
}

component main = isZero();