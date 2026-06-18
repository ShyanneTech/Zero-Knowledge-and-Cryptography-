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
    signal input is_push;
    signal input is_pop;
    signal input is_nop;
    signal output out;

    // sanity checks incase someone wants to play dirty 
    
    is_push + is_pop + is_nop === 1;
    is_nop * (1 - is_nop) === 0; 
    is_push * (1 - is_push) === 0;
    is_pop * (1 - is_pop) === 0;
  
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
  
    // the current column is 2 or more
    // below the stack pointer
    signal twoBelowSP <== LessEqThan(bits)([j, sp - 2]);

    // Satisfying condidtion 1
    component CondA = AND3(); 
    CondA.in[0] <== spGteOne;
    CondA.in[1] <== oneBelowSp;
    CondA.in[2] <== is_push + is_nop ;

    component CondB = AND3(); 
    CondB.in[0] <== spGteTwo;
    CondB.in[1] <== twoBelowSP;
    CondB.in[2] <== is_pop ;

    component or = OR();
    or.a <== CondA.out;
    or.b <== CondB.out;

    out <== or.out; 
}


// now we need something that stacks the copy info

template CopyStack(m) {
    var nBits;
    signal input sp;
    signal input is_push;
    signal input is_pop;
    signal input is_nop;
    signal output out[m];

    component ShouldCopyStacks[m];

    for (var j = 0; j < m; j++) {
        ShouldCopyStacks[j] = ShouldItCopy(j,nBits);
        ShouldCopyStacks[j].sp <== sp;
        ShouldCopyStacks[j].is_push <== is_push;
        ShouldCopyStacks[j].is_pop <== is_pop;
        ShouldCopyStacks[j].is_nop <== is_nop;
        out[j] <== ShouldCopyStacks[j].out;      
    }    
}



// my first goal is to be able to have a table that can identify whether i want to push ,pop or nop and determine what my metable is 

// push will be 1 , pop 2 ,nop 0 

template MyStackBuilder (n) {
    signal input instr[2 * n];
    signal output stack[n][n];

    var PUSH = 1;
    var POP = 2;
    var NOP = 0;

    signal sp[n + 1];
    signal spbranch[n][3];

    var SAME  = 0;
    var INC = 1;
    var  DEC = 2;
    

    signal metatable[n][4];
    // the componenets of my metatable 
    var is_NOP = 0;
    var is_PUSH = 1;
    var is_POP = 2;
    var ARG = 3 ;

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
    metatable[0][is_POP] <== 0;

    spbranch[0][SAME] <== 0 ;
    spbranch[0][INC] <== firstispush;
    spbranch[0][DEC] <== 0;
    // identifies the rest of the instructions and updates metattable from [1]
    component EqPUSH[n];
    component EqPOP[n];
    component EqNOP[n];

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



    for (var i = 1; i < n; i++) {
        EqPUSH[i] = IsEqual();
        EqPUSH[i].in[0] <== instr[2*i];
        EqPUSH[i].in[1] <== PUSH;
        metatable[i][is_PUSH] <== EqPUSH[i].out; 
        
        EqPOP[i] = IsEqual();
        EqPOP[i].in[0] <== instr[2*i];
        EqPOP[i].in[1] <== POP;
        metatable[i][is_POP] <== EqPOP[i].out; 

        EqNOP[i] = IsEqual();
        EqNOP[i].in[0] <== instr[2*i];
        EqNOP[i].in[1] <== NOP;
        metatable[i][is_NOP] <== EqNOP[i].out; 
        
        // get the instruction argument
        metatable[i][ARG] <== instr[2 * i + 1];

        // let us get the status for the copystack 
        
        CopyStack[i] = CopyStack(n);
        CopyStack[i].sp <== sp[i];
        CopyStack[i].is_push <== metatable[i] [is_PUSH];
        CopyStack[i].is_pop <== metatable[i] [is_POP];
        CopyStack[i].is_nop <== metatable[i] [is_NOP];
        
        // using this to update spbranch to get the next sp

        spbranch[i][INC] <==  metatable[i][is_PUSH] * (sp[i] + 1);
        spbranch[i][DEC] <==  metatable[i][is_POP] * (sp[i] - 1);
        spbranch[i][SAME] <==  metatable[i][is_NOP] * (sp[i]);
        sp[i + 1] <== spbranch[i][INC] + spbranch[i][SAME] + spbranch[i][DEC];

        for (var j = 0; j < n; j++) {
            previousCellIfShouldCopy[i][j] <== CopyStack[i].out[j] * stack[i-1][j];
            
            eqSP[i][j] = IsEqual();
            eqSP[i][j].in[0] <== j;
            eqSP[i][j].in[1] <== sp[i];
            eqSPAndIsPush[i][j] <== eqSP[i][j].out * metatable[i][is_PUSH];

            stack[i][j] <==  eqSPAndIsPush[i][j] * metatable[i][ARG] + previousCellIfShouldCopy[i][j];
        }

    }   
}

component main = MyStackBuilder(3);