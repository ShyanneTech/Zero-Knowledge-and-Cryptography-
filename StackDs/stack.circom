pragma circom 2.1.8;
include "circomlib/circuits/comparators.circom";
include "circomlib/circuits/gates.circom";
// j is the column number
// bits is how many bits we need
// for the LessEqThan component

template AND3() {
  signal input in[3];
  signal output out;
  
  signal temp;
  temp <== in[0] * in[1];
  out <== temp * in[2];
}

template ShouldCopy(j, bits) {
  signal input sp; // this where you confirm if sp >= 1 or sp >= 1 because those are the only conditions for copying  
  signal input is_pop;
  signal input is_push;
  signal input is_nop;
  
  // out = 1 if should copy
  signal output out;
  
  // sanity checks
  is_pop + is_push + is_nop === 1; // only one should be active 
  is_nop * (1 - is_nop) === 0; 
  is_push * (1 - is_push) === 0;
  is_pop * (1 - is_pop) === 0;
  
  // it's cheaper to compute ≠ 0 than > 0 to avoid converting the number to binary
  signal spEqZero;
  signal spGteOne;
  spEqZero <== IsZero()(sp);
  spGteOne <== 1 - spEqZero;
  
  // it's cheaper to compute ≠ 0 and ≠ 1 than ≥ 2
  signal spEqOne;
  signal spGteTwo;
  spEqOne <== IsEqual()([sp, 1]);
  spGteTwo <== (1 - spEqOne) * (1 - spEqZero);
  
  // the current column is 1 or more 
  // below the stack pointer
  signal oneBelowSp <== LessEqThan(bits)([j, sp - 1]);
  
  // the current column is 2 or more
  // below the stack pointer
  signal twoBelowSP <== LessEqThan(bits)([j, sp - 2]);
  
  // condition A
  component a3A = AND3();
  a3A.in[0] <== spGteOne;
  a3A.in[1] <== oneBelowSp;
  a3A.in[2] <== is_push + is_nop;
  
  // condition B
  component a3B = AND3();
  a3B.in[0] <== spGteTwo;
  a3B.in[1] <== twoBelowSP;
  a3B.in[2] <== is_pop;
  
  component or = OR();
  or.a <== a3A.out;
  or.b <== a3B.out;  
  out <== or.out;
}


template CopyStack(m) {
  var nBits = 4;
    signal output out[m];
    signal input sp;
    signal input is_pop;
    signal input is_push;
    signal input is_nop;

    component ShouldCopys[m];
    // signal copy[m];
    
    // loop over the columns
  for (var j = 0; j < m; j++) {
    ShouldCopys[j] = ShouldCopy(j, nBits);
    ShouldCopys[j].sp <== sp;
    ShouldCopys[j].is_pop <== is_pop;
    ShouldCopys[j].is_push <== is_push;
    ShouldCopys[j].is_nop <== is_nop;
    out[j] <== ShouldCopys[j].out;
  }
}



// n is how many instructions we can handle
// since all the instructions might be push,
// our stack needs capacity of up to n
template StackBuilder(n) {

  // to make it easier to identify the opcodes   
  var NOP = 0;
  var PUSH = 1;
  var POP = 2;

  signal input instr[2 * n]; // obviously if you have 3 things to do you would need 6 set of instuction because of the part of opcode and operand [instr,data,instr,data ]
  
  // we add one extra row for sp because
  // our algorithm always writes to the
  // next row and we don't want to conditionally
  // check for an array-out-of-bounds
  signal output sp[n + 1];

  signal output stack[n][n]; // to store based on iteration 
  
  // the purpose of it is to identify and updata metatble  
  var IS_NOP = 0; // the numbering here is to identify the columns yourself
  var IS_PUSH = 1;
  var IS_POP = 2;
  var ARG = 3;
  
  // metaTable is the columns IS_NOP, IS_PUSH, IS_POP, ARG
  signal metaTable[n][4];

  // first instruction must be PUSH or NOP since it starts off being empty...

  (instr[0] - PUSH) * (instr[0] - NOP) === 0; // if this is asserted to be true then it implies that it can be either PUSH or NOP
  // we have to confirm which is which and makes it easy such that if it doesnt satisfy one it automatically satisfies the other 

  signal first_op_is_push;
  first_op_is_push <== IsEqual()([instr[0], PUSH]); 

  // if the first op is NOP, we are forcing the first
  // value to be zero, but this is where the stack
  // pointer is, so it doesn't matter
  stack[0][0] <== first_op_is_push * instr[1]; // this is us hardcoding the 0th row 0th col 

  // initialize the rest of the first stack to be zero for the remaining part of the zeroth row
  for (var i = 1; i < n; i++) {
      stack[0][i] <== 0; // the stack is initially empty [0,0,0,0]
  }

  // we fill out the 0th elements to avoid
  // uninitialzed signals. For a particular
  // execution, we only want one possible witness
  // to correspond to a particular execution

  sp[0] <== 0; // sp is also initially 0 
  sp[1] <== first_op_is_push;
  metaTable[0][IS_PUSH] <== first_op_is_push;
  metaTable[0][IS_POP] <== 0;
  metaTable[0][IS_NOP] <== 1 - first_op_is_push;
  metaTable[0][ARG] <== instr[1];

  // spBranch is what we add to the previous
  // stack pointer based on the opcode.
  // Could be 1, 0, or -1 depending on the
  // opcode. Since the first opcode
  // cannot be POP, -1 is not an option here.
  var SAME = 0; // the numbering here is to identify the columns yourself 
  var INC = 1;
  var DEC = 2;
  signal spBranch[n][3];
  spBranch[0][INC] <== first_op_is_push * 1;
  spBranch[0][SAME] <== (1 - first_op_is_push) * 0;
  spBranch[0][DEC] <== 0;

  // populate the first row of the metaTable
  // and the stack pointer
  // the reason we have it as an array n is that it would check all through the instruction part of the code which will be n in numbers and with that be used to update the metatable 
  component EqPush[n];
  component EqNop[n];
  component EqPop[n];

  component eqSP[n][n];
  signal eqSPAndIsPush[n][n];
  for (var i = 0; i < n; i++) {
      eqSPAndIsPush[0][i] <== 0;
  }

  // signals and components for copying
  component CopyStack[n];
  signal previousCellIfShouldCopy[n][n];
  for (var i = 0; i < n; i++) {
    previousCellIfShouldCopy[0][i] <== 0;
  }

  // These addresses the instructions and which of them is active
  for (var i = 1; i < n; i++) {
    // check which opcode we are executing, the even index part of the instruction input 

    // METATABLE COMPLETELY FILLED
    EqPush[i] = IsEqual();
    EqPush[i].in[0] <== instr[2 * i];
    EqPush[i].in[1] <== PUSH;
    metaTable[i][IS_PUSH] <== EqPush[i].out;

    EqNop[i] = IsEqual();
    EqNop[i].in[0] <== instr[2 * i];
    EqNop[i].in[1] <== NOP;
    metaTable[i][IS_NOP] <== EqNop[i].out;

    EqPop[i] = IsEqual();
    EqPop[i].in[0] <== instr[2 * i];
    EqPop[i].in[1] <== POP;
    metaTable[i][IS_POP] <== EqPop[i].out;

    // get the instruction argument
    metaTable[i][ARG] <== instr[2 * i + 1];

    // This where we handle the copying instructions to fill the stack 
    // if it is a push, write to the stack
    // if it is a copy, write to the stack
    CopyStack[i] = CopyStack(n);
    CopyStack[i].sp <== sp[i];
    CopyStack[i].is_push <== metaTable[i][IS_PUSH];
    CopyStack[i].is_nop <== metaTable[i][IS_NOP];
    CopyStack[i].is_pop <== metaTable[i][IS_POP];
    for (var j = 0; j < n; j++) {
      previousCellIfShouldCopy[i][j] <== CopyStack[i].out[j] * stack[i - 1][j];

      eqSP[i][j] = IsEqual();
      eqSP[i][j].in[0] <== j;
      eqSP[i][j].in[1] <== sp[i];
      eqSPAndIsPush[i][j] <== eqSP[i][j].out * metaTable[i][IS_PUSH];

      // we will either PUSH or COPY or implicilty assign 0
      stack[i][j] <== eqSPAndIsPush[i][j] * metaTable[i][ARG] + previousCellIfShouldCopy[i][j];
    }

    // write to the next row's stack pointer
    spBranch[i][INC] <== metaTable[i][IS_PUSH] * (sp[i] + 1);
    spBranch[i][SAME] <== metaTable[i][IS_NOP] * (sp[i]);
    spBranch[i][DEC] <== metaTable[i][IS_POP] * (sp[i] - 1);
    sp[i + 1] <== spBranch[i][INC] + spBranch[i][SAME] + spBranch[i][DEC];
  }
}

component main = StackBuilder(3);

/* INPUT = {
  "instr": [1, 16, 1, 20, 1, 22]
} */