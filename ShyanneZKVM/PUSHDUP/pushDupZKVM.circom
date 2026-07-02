pragma circom 2.1.8;
include "circomlib/circuits/comparators.circom";

// And Template to ensure all conditions are satisfied 
template AND3() {
  signal input in[3];
  signal output out;
  
  signal temp;
  temp <== in[0] * in[1];
  out <== temp * in[2];
}

// I want to deal with the conditions for copy for a column

template ShouldItCopy(j , nbits){
    signal input is_PUSH;
    signal input is_DUP;
    signal input sp;
    signal output out; 

    // sanity checks 
    is_PUSH + is_DUP === 1;
    (1 - is_PUSH ) * is_PUSH === 0;
    (1 - is_DUP ) * is_DUP === 0;

    // you can only copy for a push or dup condition when sp >= 1
    signal spGteOne;
    signal eqZero;
    eqZero <== IsEqual()([sp,0]);
    spGteOne <== 1 - eqZero;

    // now we want to attend the columns to be copied up to sp-1 for both push and dup conditions 

    signal OneBelowSP;
    OneBelowSP <== LessEqThan(nbits)([j,sp-1]);

    component CondA = AND3();
    CondA.in[0] <== (is_PUSH + is_DUP);
    CondA.in[1] <== spGteOne;
    CondA.in[2] <== OneBelowSP;

    out <== CondA.out;
}

// this combines everything together 
template CopyStackInfo(n) {
    signal input is_PUSH;
    signal input is_DUP;
    signal input sp;
    var nbits = 4;
    signal output out[n]; 

    component ShouldCopys[n];
    for (var j = 0; j < n; j++) {
        ShouldCopys[j] = ShouldItCopy(j, nbits);
        ShouldCopys[j].is_PUSH <== is_PUSH;
        ShouldCopys[j].is_DUP <== is_DUP;
        ShouldCopys[j].sp <== sp;

        out[j] <== ShouldCopys[j].out;

    }
    
}

template PushDupZKVM(n) {
    signal input instr[2*n];
    signal output stack[n][n];

    // Instructions I expect 
    var PUSH = 0;
    var DUP = 1;

    signal sp[n + 1];
    signal spBranch[n][1];
    var INC = 0;

    // The table that stored the final of everything and that i can just connect to the main stack output 

    signal metatable[n][3];
    var is_PUSH = 0;
    var is_DUP = 1;
    var ARG = 2;

    // things that are initially set 
    sp[0] <== 0;

    // the first insstruction has to be strictly PUSH which is identified by 0
    
    (PUSH - instr[0]) === 0;

    signal firstIsPUSH;
    firstIsPUSH <== IsEqual()([instr[0], PUSH]);

    stack[0][0]<== firstIsPUSH * instr[1];
    for (var j = 1; j < n; j++) {
        stack[0][j] <== 0;
    }

    sp[1] <== firstIsPUSH;
    spBranch[0][INC] <== firstIsPUSH;

    metatable[0][is_PUSH] <== firstIsPUSH;
    metatable[0][is_DUP] <== 0;
    metatable[0][ARG] <== instr[1];

    // to identify the rest of what exactly is PUSH OR DUP OR EVEN A WRONG OPCODE to update the metatble 

    component EqPUSH[n];
    component EqDUP[n];

    signal previousCellsIfCopy[n][n];
    for (var i = 0; i < n; i++) {
        previousCellsIfCopy[0][i] <== 0;
    }

    component CopyStackInfo[n];

    component eqSp[n][n];
    signal eqSPAndPUSH[n][n];
    signal eqSPAndDUP[n][n];

    for (var i = 0; i < n; i++) {
        eqSPAndPUSH[0][i] <== 0 ;
        eqSPAndDUP[0][i] <== 0 ;
        
    }

    signal eqSPAndDUPWithValue[n][n];
    component eqSPMinusOne [n][n];

    signal eqMinusOneValue[n][n];
    signal duplicatedValue[n][n];

    for (var i = 0; i < n; i++) {
        duplicatedValue[i][0] <== 0;
    }

    for (var i = 1; i < n; i++){
        EqPUSH[i] = IsEqual();
        EqPUSH[i].in[0] <== instr[2 * i];
        EqPUSH[i].in[1] <== PUSH;
        metatable[i][is_PUSH] <== EqPUSH[i].out;

        EqDUP[i] = IsEqual();
        EqDUP[i].in[0] <== instr[2 * i];
        EqDUP[i].in[1] <== DUP;
        metatable[i][is_DUP] <== EqDUP[i].out;
        metatable[i][ARG] <== instr[2 * i + 1];

        // the goal here is to update the spbranch and then update the sp array itself 

        spBranch[i][INC] <== (metatable[i][is_PUSH] + metatable[i][is_DUP]) * (sp[i] + 1);
        sp[i+1] <==  spBranch[i][INC]; 

        CopyStackInfo[i] = CopyStackInfo(n);
        CopyStackInfo[i].sp <== sp[i];
        CopyStackInfo[i].is_PUSH <== is_PUSH;
        CopyStackInfo[i].is_DUP <== is_DUP;


        // now I want to start working on my Stack
        for (var j = 0; j < n; j++) {
            previousCellsIfCopy[i][j] <== stack[i - 1][j] * CopyStackInfo[i].out[j];

            eqSp[i][j]= IsEqual();
            eqSp[i][j].in[0] <== j;
            eqSp[i][j].in[1] <== sp[i];

            eqSPMinusOne[i][j] = IsEqual();
            eqSPMinusOne[i][j].in[0] <== j;
            eqSPMinusOne[i][j].in[1] <== sp[i] - 1;


            eqSPAndPUSH[i][j] <== eqSp[i][j].out * metatable[i][is_PUSH]; 

            eqSPAndDUP[i][j] <== eqSp[i][j].out * metatable[i][is_DUP];


            eqMinusOneValue[i][j] <== eqSPMinusOne[i][j].out * previousCellsIfCopy[i][j];



        } 

         for (var j = 1; j < n; j++) {
            duplicatedValue[i][j] <== eqMinusOneValue[i][j-1]; 
                
        }

        for (var j = 0; j < n; j++) {
            eqSPAndDUPWithValue[i][j] <== eqSPAndDUP[i][j] * duplicatedValue[i][j] ;
            stack[i][j] <==  previousCellsIfCopy[i][j] + (eqSPAndPUSH[i][j] * metatable[i][ARG]) + eqSPAndDUPWithValue[i][j];
        }

    }

    
}
component main  = PushDupZKVM(6);