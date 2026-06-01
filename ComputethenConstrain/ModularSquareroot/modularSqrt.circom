pragma circom 2.0.0;

include "circomlib/circuits/pointbits.circom";


// function sqrt(n) {
//      if (n == 0) {
//         return 0;
//     }

//     // Test that have solution
//     var res = n ** ((-1) >> 1);
// //        if (res!=1) assert(false, "SQRT does not exists");
//     if (res!=1) return 0;

//     var m = 28;
//     var c = 19103219067921713944291392827692070036145651957329286315305642004821462161904;
//     var t = n ** 81540058820840996586704275553141814055101440848469862132140264610111;
//     var r = n ** ((81540058820840996586704275553141814055101440848469862132140264610111+1)>>1);
//     var sq;
//     var i;
//     var b;
//     var j;

//     while ((r != 0)&&(t != 1)) {
//         sq = t*t;
//         i = 1;
//         while (sq!=1) {
//             i++;
//             sq = sq*sq;
//         }

//         // b = c ^ m-i-1
//         b = c;
//         for (j=0; j< m-i-1; j ++) b = b*b;

//         m = i;
//         c = b*b;
//         t = t*c;
//         r = r*b;
//     }

//     if (r < 0 ) {
//         r = -r;
//     }

//     return r;
// }
template ValidSqrt() {
    signal input in;
    signal output out1;
    signal output out2;

    out1 <-- sqrt(in);
    out2 <-- out1 * -1;
    out1 * out1 === in;
    out2 * out2 === in;
  
}

component main = ValidSqrt();
