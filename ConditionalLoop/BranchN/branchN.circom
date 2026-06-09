pragma circom 2.1.0;
include "circomlib/circuits/comparators.circom";
include "circomlib/circuits/multiplexer.circom";

template BranchN(n) {

    assert(n > 1);
    signal input x;
    signal input conds[n-1];
    signal input branches[n];

    signal output out;

    signal switches[n];

    component EqualityChecks[n-1];

    for (var i = 0; i < n-1; i++) {
        EqualityChecks[i] = IsEqual();
        EqualityChecks[i].in[0] <== x;
        EqualityChecks[i].in[1] <== conds[i];
        switches[i]<==  EqualityChecks[i].out;        
    }

  var total = 0;
  for (var i = 0; i < n - 1; i++) {
    total += switches[i];
  }
    switches[n-1]<== IsZero()(total);

    component InnerProduct = EscalarProduct(n);
    for (var i = 0; i < n; i++) {
    InnerProduct.in1[i] <== switches[i];
    InnerProduct.in2[i] <== branches[i];
    }
  
    out <== InnerProduct.out;   
}

template MultiBranchConditions() {
    signal input x;
    signal output out;

    component branchn = BranchN(4);

    var conds[3] = [5, 9, 10];
    var branches[4] = [14, 22, 23, 45];

    for (var i = 0; i < 4; i++) {
     if (i < 3) {
        branchn.conds[i] <== conds[i];
        }
    
        branchn.branches[i] <== branches[i];
    }

     branchn.x <== x;
     branchn.out ==> out;


}

 component main = MultiBranchConditions();