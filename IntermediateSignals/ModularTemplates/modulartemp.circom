pragma circom 2.0.0;

// separate template
template Mul3() {
  signal input a;
  signal input b;
  signal input c;
  signal input d; // d === a * b * c
  
  // no longer an input
  signal s;
  
  a * b ==> s;
  s * c === d;
}

// main component
template Mul3x2() {
  signal input a;
  signal input b;
  signal input c;
  signal input d; // d === a * b * c
  
  signal input x;
  signal input y;
  signal input z;
  signal input u; // u === x * y * z
  
  component m3_1 = Mul3();
  m3_1.a <== a;
  m3_1.b <== b;
  m3_1.c <== c;
  m3_1.d <== d;
  
  component m3_2 = Mul3();
  m3_2.a <== x;
  m3_2.b <== y;
  m3_2.c <== z;
  m3_2.d <== u;
}
