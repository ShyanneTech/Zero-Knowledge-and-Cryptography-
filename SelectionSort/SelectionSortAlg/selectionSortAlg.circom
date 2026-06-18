pragma circom 2.1.8;
include "circomlib/circuits/comparators.circom";
include "circomlib/circuits/Multiplexer.circom";
template Swap(n) {
  signal input in[n];
  signal input s;
  signal input t;
  signal output out[n];

  signal sEqt;
  sEqt <== IsEqual()([s,t]);

  signal valueAtPositionS;
  signal valueAtPositionT;
  

  component getvalueAtPositionS = Multiplexer(1,n);
  component getvalueAtPositionT = Multiplexer(1,n);
  getvalueAtPositionS.sel <== s;
  getvalueAtPositionT.sel <== t;

  for (var i = 0; i < n; i++) {
    getvalueAtPositionS.inp[i][0] <== s; 
    getvalueAtPositionT.inp[i][0] <== t; 
   
  }

  valueAtPositionS <==  getvalueAtPositionS.out[0];
  valueAtPositionT <==  getvalueAtPositionT.out[0];

  component ArrayS[n];
  component ArrayT[n];
  component NotArrayST[n];

  for (var i = 0; i < n; i++) {
    ArrayS[i] = IsEqual();
    ArrayS[i].in[0] <== s;
    ArrayS[i].in[1] <== i; //[0,1,0,0]

    ArrayT[i] = IsEqual();
    ArrayT[i].in[0] <== t;
    ArrayT[i].in[1] <== i; //[0,0,1,0]

    NotArrayST[i] = IsZero();
    NotArrayST[i].in <== ArrayS[i].out + ArrayT[i].out; // [1,0,0,1];  
  }

signal BranchS[n];
signal BranchT[n];
signal NotBranchST[n];

for (var i = 0; i < n; i++) {
    BranchS[i] <== ArrayT[i].out * valueAtPositionS;
    BranchT[i] <== ArrayS[i].out * valueAtPositionT;
    NotBranchST[i] <== NotArrayST[i].out * in[i];

    out[i] <== (1- sEqt) * BranchS[i] + BranchT[i] + NotBranchST[i]; 
 }

}


template Select(n, start) {
  // unsorted list
  signal input in[n];
  
  // index start swapped with the min
  signal output out[n];

  // we will define GetMinAtIdxStartingAt in the next codeblock
  component minIdx0 = GetMinAtIdxStartingAt(n, start);
  for (var i = 0; i < n; i++) {
      minIdx0.in[i] <== in[i];
  }

  component Swap0 = Swap(n);
  Swap0.s <== start; // swap 0 with the min
  Swap0.t <== minIdx0.idx; // with the min (could be idx 0)
  for (var i = 0; i < n; i++) {
    Swap0.in[i] <== in[i];
  } 

  // copy to out
  for (var i = 0; i < n; i++) {
      out[i] <== Swap0.out[i];
  }
}

// formerly GetMinAtIdx
template GetMinAtIdxStartingAt(n, start) {
  signal input in[n];
  signal output min;
  signal output idx;

  // only look for values start and later
  var minv = in[start];
  var idxv = start;
  for (var i = start + 1; i < n; i++) {
    if (in[i] < minv) {
      minv = in[i];
      idxv = i;
    }
  }
  min <-- minv;
  idx <-- idxv;

  // only compare to values start and later
  component lt[n];
  
  // CHANGES HERE: LOOP FROM START TO N-1
  for (var i = start ; i < n; i++) {
    lt[i] = LessEqThan(252);
    lt[i].in[0] <== min;
    lt[i].in[1] <== in[i];
    lt[i].out === 1;
  }

  // assert min is really at in[idx]
    component correctMinIdx = Multiplexer(1,n);
    correctMinIdx.sel <== idx;

    for (var i = 0; i < n; i++) {
        correctMinIdx.inp[i][0] <== in[i];
    }
    correctMinIdx.out[0] === min;
}


template SelectionSort(n) {
  assert(n > 0);

  signal input in[n];
  signal output out[n];

  signal intermediateStates[n][n];

  component SSort[n - 1];
  for (var i = 0; i < n; i++) {
    // copy the input to the first row of
    // intermediateStates. Note that we can do
    // if(i == 0) because i is not a signal
    // and i is known at compile time
    if (i == 0) {
      for (var j = 0; j < n; j++) {
        intermediateStates[0][j] <== in[j];
      }
    }

    else {
      // select sort n items starting at i - 1
      // for i = 1, we compare item at 0 to
      // the rest of the list
      SSort[i - 1] = Select(n, i - 1);
      
      // load in the intermediate state i -1
      for (var j = 0; j < n; j++) {
        SSort[i - 1].in[j] <== intermediateStates[i - 1][j];
      }
      // write the sorted result to row i
      for (var j = 0; j < n; j++) {
        SSort[i - 1].out[j] ==> intermediateStates[i][j];
      }
    }
  }

  // write the final state to the ouput
  for (var i = 0; i < n; i++) {
    intermediateStates[n-1][i] ==> out[i];
  }
}

component main = SelectionSort(9);

/* INPUT = {"in": [3,1,8,2,4,0,1,2,4]} */
