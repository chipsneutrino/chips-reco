#ifndef WCSIMNVIDIAMATH_HH
#define WCSIMNVIDIAMATH_HH

// We end up calculating a LOT of trig functions
// This namespace contains some good, fast approximate ways
// of evaluating trig and inverse trig functions
//
// Algorithms come from the NVidia CG toolkit
// http://http.developer.nvidia.com/Cg/index_stdlib.html

namespace WCSimNvidiaMath{
  float sin(float radians);
  float cos(float radians); 
};

#endif
