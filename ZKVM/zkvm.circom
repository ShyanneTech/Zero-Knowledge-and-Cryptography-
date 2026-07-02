pragma circom 2.1.8;

include "circomlib/circuits/comparators.circom";
include "circomlib/circuits/gates.circom";
// now my focus is to work on creating a template that determines if it should copy or not  so it says yes(1) or no (0) 

// 1. if it is push or nop you copy from 0 to 1 , however sp must be from 1 and above 

// 2. if it is pop you copy from 0 to sp - 2 , however sp must be greater than 2 

template AND3() {
  signal input in[3];
  signal output out;
  
  signal temp;
  temp <== in[0] * in[1];
  out <== temp * in[2];
}

template ShouldItCopy(j,bits) {
    signal input sp;
    signal input is_nop;
    signal input is_push;
    signal input is_add;
    signal input is_mul;
    signal output out;

    // sanity checks incase someone wants to play dirty 
    
    is_push + is_mul + is_add + is_nop === 1;
    is_nop * (1 - is_nop) === 0; 
    is_push * (1 - is_push) === 0;
    is_add * (1 - is_add) === 0;
    is_mul * (1 - is_mul) === 0;
  
    // this is focused on the sp value 
    // is sp >= 1 id it isnt equal to zero it will be 

    signal spGteOne;
    signal eqZero;
    eqZero <== IsEqual()([sp,0]);
    spGteOne <== 1 - eqZero;

    // is sp >= 2 id it isnt equal to zero or 1 it will be 

    signal spGteTwo;
    signal eqOne;
    eqOne <== IsEqual()([sp,1]);
    spGteTwo <== 1 - (eqZero + eqOne);


    // to identify each coloumn  
    signal oneBelowSp <== LessEqThan(bits)([j, sp - 1]);
  
    // the current column is 3 or more
    // below the stack pointer
    signal threeBelowSP <== LessEqThan(bits)([j, sp - 3]);

    // Satisfying condidtion 1
    component CondA = AND3(); 
    CondA.in[0] <== spGteOne;
    CondA.in[1] <== oneBelowSp;
    CondA.in[2] <== is_push + is_nop ;

    component CondB = AND3(); 
    CondB.in[0] <== spGteTwo;
    CondB.in[1] <== threeBelowSP;
    CondB.in[2] <== (is_add + is_mul);

    component or = OR();
    or.a <== CondA.out;
    or.b <== CondB.out;

    out <== or.out; 
}


// now we need something that stacks the copy info

template CopyStack(m) {
    var nBits = 4;
    signal input sp;
    signal input is_nop;
    signal input is_push;
    signal input is_add;
    signal input is_mul;
    signal output out[m];

    component ShouldCopyStacks[m];

    for (var j = 0; j < m; j++) {
        ShouldCopyStacks[j] = ShouldItCopy(j,nBits);
        ShouldCopyStacks[j].sp <== sp;
        ShouldCopyStacks[j].is_nop <== is_nop;
        ShouldCopyStacks[j].is_push <== is_push;
        ShouldCopyStacks[j].is_add <== is_add;
        ShouldCopyStacks[j].is_mul <== is_mul;
        out[j] <== ShouldCopyStacks[j].out;      
    }    
}



// my first goal is to be able to have a table that can identify whether i want to push ,pop or nop and determine what my metable is 

// push will be 1 , pop 2 ,nop 0 

template MyZKVM(n) {
    signal input instr[2 * n];
    signal output stack[n][n];
    var NOP = 0;
    var PUSH = 1;
    var ADD = 2;
    var MUL = 3;

    signal sp[n + 1];
    signal spbranch[n][3];

    var SAME  = 0;
    var INC = 1;
    var  DEC = 2;
    

    signal metatable[n][5];
    // the componenets of my metatable 
    var is_NOP = 0;
    var is_PUSH = 1;
    var is_ADD = 2;
    var is_MUL = 3;
    var ARG = 4 ;

    // making initial hardcoded decisions 

    // making sure that the first two instructions are either push or nop

    (PUSH - instr[0]) * (NOP - instr[0]) === 0;

    signal firstispush;
    firstispush <== IsEqual()([instr[0],PUSH]);
    sp[0] <== 0;
    sp[1] <== firstispush;
    
    stack[0][0] <== firstispush * instr[1];
    for (var i = 1; i < n; i++) {
        stack[0][i] <== 0;
    }

    metatable[0][is_PUSH] <== firstispush;
    metatable[0][is_NOP] <== 1 - firstispush;
    metatable[0][is_ADD] <== 0;
    metatable[0][is_MUL] <== 0;
    metatable[0][ARG] <== instr[1];

    spbranch[0][SAME] <== 0 ;
    spbranch[0][INC] <== firstispush;
    spbranch[0][DEC] <== 0;
    // identifies the rest of the instructions and updates metattable from [1]
    component EqPUSH[n];
    component EqNOP[n];
    component EqADD[n];
    component EqMUL[n];

    signal previousCellIfShouldCopy[n][n];
    for (var j = 0; j < n; j++) {
        previousCellIfShouldCopy[0][j] <== 0;
    }

   component  CopyStack[n];
   component eqSP[n][n];
   signal eqSPAndIsPush[n][n];
   for (var i = 0; i < n; i++) {
       eqSPAndIsPush[0][i] <== 0; 
    }

    component eqSPMinus2[n][n];
    signal eqSPMinus2AndIsAdd[n][n];
    signal eqSPMinus2AndIsMul[n][n];  
    for (var i = 0; i < n; i++) {
       eqSPMinus2AndIsAdd[0][i] <== 0; 
       eqSPMinus2AndIsMul[0][i] <== 0; 
    } 
    // (the current column = sp - 2 and is_add) * sum
    signal eqSPMinus2AndIsAddWithValue[n][n];
    signal eqSPMinus2AndIsMulWithValue[n][n];
    signal sum_result[n][n]; // the idea of this is to sum all j and j+1 
    signal mul_result[n][n]; // the idea of this is to multiply all j and j+1 
    for (var i = 0; i < n; i++) {
    eqSPMinus2AndIsAddWithValue[0][i] <== 0;
    eqSPMinus2AndIsMulWithValue[0][i] <== 0;
    sum_result[0][i] <== 0;
    mul_result[0][i] <== 0; 
    }

    

    for (var i = 1; i < n; i++) {
        EqPUSH[i] = IsEqual();
        EqPUSH[i].in[0] <== instr[2*i];
        EqPUSH[i].in[1] <== PUSH;
        metatable[i][is_PUSH] <== EqPUSH[i].out; 

        EqNOP[i] = IsEqual();
        EqNOP[i].in[0] <== instr[2*i];
        EqNOP[i].in[1] <== NOP;
        metatable[i][is_NOP] <== EqNOP[i].out; 

        EqADD[i] = IsEqual();
        EqADD[i].in[0] <== instr[2*i];
        EqADD[i].in[1] <== ADD;
        metatable[i][is_ADD] <== EqADD[i].out; 
        
        EqMUL[i] = IsEqual();
        EqMUL[i].in[0] <== instr[2*i];
        EqMUL[i].in[1] <== MUL;
        metatable[i][is_MUL] <== EqMUL[i].out; 
        
        // get the instruction argument
        metatable[i][ARG] <== instr[2 * i + 1];

        // carry out the sums and muls
        for (var j = 0; j < n - 1; j++) {
            sum_result[i][j] <== stack[i - 1][j] + stack[i - 1][j + 1];
            mul_result[i][j] <== stack[i - 1][j] * stack[i - 1][j + 1];
        }

        // for completeness of array and signals 
        for (var j = n - 1; j < n; j++) {
            sum_result[i][j] <== 0;
            mul_result[i][j] <== 0;
        }


        // let us get the status for the copystack 
        
        CopyStack[i] = CopyStack(n);
        CopyStack[i].sp <== sp[i];
        CopyStack[i].is_push <== metatable[i] [is_PUSH];
        CopyStack[i].is_nop <== metatable[i] [is_NOP];
        CopyStack[i].is_add <== metatable[i] [is_ADD];
        CopyStack[i].is_mul <== metatable[i] [is_MUL];
        
        // using this to update spbranch to get the next sp

        spbranch[i][INC] <==  metatable[i][is_PUSH] * (sp[i] + 1);
        spbranch[i][SAME] <==  metatable[i][is_NOP] * (sp[i]);
        spbranch[i][DEC] <==  (metatable[i][is_ADD] + metatable[i][is_MUL]) * (sp[i] - 1);
        sp[i + 1] <== spbranch[i][INC] + spbranch[i][SAME] + spbranch[i][DEC];

        for (var j = 0; j < n; j++) {
            previousCellIfShouldCopy[i][j] <== CopyStack[i].out[j] * stack[i-1][j];
            
            eqSP[i][j] = IsEqual();
            eqSP[i][j].in[0] <== j;
            eqSP[i][j].in[1] <== sp[i];
            eqSPAndIsPush[i][j] <== eqSP[i][j].out * metatable[i][is_PUSH];

            eqSPMinus2[i][j] = IsEqual();
            eqSPMinus2[i][j].in[0] <== j;
            eqSPMinus2[i][j].in[1] <== sp[i] - 2;
            eqSPMinus2AndIsAdd[i][j] <==  eqSPMinus2[i][j].out * metatable[i][is_ADD];
            eqSPMinus2AndIsMul[i][j] <==  eqSPMinus2[i][j].out * 
            metatable[i][is_MUL];
            eqSPMinus2AndIsAddWithValue[i][j] <== eqSPMinus2AndIsAdd[i][j] * sum_result[i][j];
            eqSPMinus2AndIsMulWithValue[i][j] <== eqSPMinus2AndIsMul[i][j] * mul_result[i][j];



            stack[i][j] <==  (eqSPAndIsPush[i][j] * metatable[i][ARG]) + eqSPMinus2AndIsAddWithValue[i][j] + eqSPMinus2AndIsMulWithValue[i][j] + previousCellIfShouldCopy[i][j];
        }

    }   
}

component main = MyZKVM(5);