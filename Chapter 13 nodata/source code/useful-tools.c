#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "array_alloc.h"
#include "useful-tools.h"


// WARNING - inline function was removed from .h because it did not compile well

/* https://gcc.gnu.org/onlinedocs/gcc/Inline.html :
   By declaring a function inline, you can direct GCC to make calls to that function faster. One way GCC can achieve this is to integrate that function's code into the code for its callers. This makes execution faster by eliminating the function-call overhead */
inline int Periodic_int(int a, int A){
  return (a % A + A) % A;
}

inline double Periodic(double x, double X){
  while(x<0)  x+=X;
  while(x>=X) x-=X;
  return x;
}

inline double Shortest(double v, double X_2, double X){
  if(v > X_2)  return v-X;
  if(v < -X_2) return v+X;
  return v;
}

