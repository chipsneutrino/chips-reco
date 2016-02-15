#ifndef WCSIMFASTMATH_HH
#define WCSIMFASTMATH_HH

#include <cmath>
// We end up calculating a LOT of trig functions and awkward
// powers, exponentials, error functions, etc.
// This namespace contains some good, fast approximate ways
// of evaluating mathematical expressions to a level of accuracy
// the fitter can live with

namespace WCSimFastMath{
    /*******************************************************************//**
    * \brief Fast cosine function as implemented by Nvidia's CG toolkit
    *
    * Algorithm comes from 
    * http://http.developer.nvidia.com/Cg/index_stdlib.html
    *
    * @param radians Angle to calculate the cosine of, in radians
    * @return cosine of the angle provided, with |error| < 1.8e-7
    */
    inline float cos(float radians)
    {
    	/* C simulation gives a max absolute error of less than 1.8e-7 */
    	float c0[4] = {0.0, 0.5, 1.0, 0.0};
    	float c1[4] = {0.25, -9.0, 0.75, 0.159154943091};
    	float c2[4] = {24.9808039603, -24.9808039603, -60.1458091736,
    			60.1458091736};
    	float c3[4] = {85.4537887573, -85.4537887573, -64.9393539429,
    			64.9393539429};
    	float c4[4] = {19.7392082214, -19.7392082214, -1.0, 1.0};
    
    	/* r0.x = sin(a) */
    	float r0[3] = {0,0,0}, r1[3] = {0,0,0}, r2[3] = {0,0,0};
      r1[0]  = c1[3] * radians;      
      r1[1] = r1[0] - floor(r1[0]);
    
      r2[0] = (float) ( r1[1] < c1[0]);
      r2[1] = (float) ( r1[1] >= c1[1] );
      r2[2] = (float) ( r1[1] >= c1[2] );
    
      r2[1] = (r2[0]*c4[2] + r2[1]*c4[3] + r2[2]*c4[2]);
    
      r0[0] = c0[0] - r1[1];
      r0[1] = c0[1] - r1[1]; 
      r0[2] = c0[2] - r1[1];
    
      r0[0] = r0[0]*r0[0];
      r0[1] = r0[1]*r0[1];
      r0[2] = r0[2]*r0[2];
    
      r1[0] = c2[0] * r0[0] + c2[2];
      r1[1] = c2[1] * r0[1] + c2[3];
      r1[2] = c2[0] * r0[2] + c2[2];
    
      r1[0] = r1[0] * r0[0] + c3[0];
      r1[1] = r1[1] * r0[1] + c3[1]; 
      r1[2] = r1[2] * r0[2] + c3[0];
    
      r1[0] = r1[0] * r0[0] + c3[2];
      r1[1] = r1[1] * r0[1] + c3[3];
      r1[2] = r1[2] * r0[2] + c3[2];
    
      r1[0] = r1[0] * r0[0] + c4[0];
      r1[1] = r1[1] * r0[1] + c4[1];
      r1[2] = r1[2] * r0[2] + c4[0];
    
      r1[0] = r1[0] * r0[0] + c4[2];
      r1[1] = r1[1] * r0[1] + c4[3];
      r1[2] = r1[2] * r0[2] + c4[2]; 
    
      r0[0] = r1[0] * -r2[0] + r1[1] * -r2[1] + r1[2]*-r2[2];
      return r0[0];
    };
    
    /*******************************************************************//**
    * \brief Fast sine function as implemented by Nvidia's CG toolkit
    *
    * Algorithm comes from 
    * http://http.developer.nvidia.com/Cg/index_stdlib.html
    *
    * @param radians Angle to calculate the sine of, in radians
    * @return sine of the angle provided, with |error| < 1.8e-7
    */
    inline float sin(float radians)
    {
    	/* C simulation gives a max absolute error of less than 1.8e-7 */
    	float c0[4] = {0.0, 0.5, 1.0, 0.0};
    	float c1[4] = {0.25, -9.0, 0.75, 0.159154943091};
    	float c2[4] = {24.9808039603, -24.9808039603, -60.1458091736,
    			60.1458091736};
    	float c3[4] = {85.4537887573, -85.4537887573, -64.9393539429,
    			64.9393539429};
    	float c4[4] = {19.7392082214, -19.7392082214, -1.0, 1.0};
    
    	/* r0.x = sin(a) */
    	float r0[3] = {0,0,0}, r1[3] = {0,0,0}, r2[3] = {0,0,0};
      r1[0] = c1[3] * radians - c1[0];
      r1[1] = r1[0] - floor(r1[0]);
    
      r2[0] = (float) ( r1[1] < c1[0]);
      r2[1] = (float) ( r1[1] >= c1[1] );
      r2[2] = (float) ( r1[1] >= c1[2] );
    
      r2[1] = (r2[0]*c4[2] + r2[1]*c4[3] + r2[2]*c4[2]);
    
      r0[0] = c0[0] - r1[1];
      r0[1] = c0[1] - r1[1]; 
      r0[2] = c0[2] - r1[1];
    
      r0[0] = r0[0]*r0[0];
      r0[1] = r0[1]*r0[1];
      r0[2] = r0[2]*r0[2];
    
      r1[0] = c2[0] * r0[0] + c2[2];
      r1[1] = c2[1] * r0[1] + c2[3];
      r1[2] = c2[0] * r0[2] + c2[2];
    
      r1[0] = r1[0] * r0[0] + c3[0];
      r1[1] = r1[1] * r0[1] + c3[1]; 
      r1[2] = r1[2] * r0[2] + c3[0];
    
      r1[0] = r1[0] * r0[0] + c3[2];
      r1[1] = r1[1] * r0[1] + c3[3];
      r1[2] = r1[2] * r0[2] + c3[2];
    
      r1[0] = r1[0] * r0[0] + c4[0];
      r1[1] = r1[1] * r0[1] + c4[1];
      r1[2] = r1[2] * r0[2] + c4[0];
    
      r1[0] = r1[0] * r0[0] + c4[2];
      r1[1] = r1[1] * r0[1] + c4[3];
      r1[2] = r1[2] * r0[2] + c4[2]; 
    
      r0[0] = r1[0] * -r2[0] + r1[1] * -r2[1] + r1[2]*-r2[2];
      return r0[0];
    };
    
    /*******************************************************************//**
    * \brief Fast power function by Martin Ankerl
    *
    * Algorithm comes from 
    * http://martin.ankerl.com/2012/01/25/optimized-approximative-pow-in-c-and-cpp/
    *
    * @param a Number to raise to some power
    * @param b Power to which a should be raised
    * @return a^b 
    */
    inline double pow(double a, double b) {
      // calculate approximation with fraction of the exponent
      if(b < 0) { return pow(a, -b); }

      int e = (int) b;
      union {
        double d;
        int x[2];
      } u = { a };
      u.x[1] = (int)((b - e) * (u.x[1] - 1072632447) + 1072632447);
      u.x[0] = 0;
     
      // exponentiation by squaring with the exponent's integer part
      // double r = u.d makes everything much slower, not sure why
      double r = 1.0;
      while (e) {
        if (e & 1) {
          r *= a;
        }
        a *= a;
        e >>= 1;
      }
     
      return r * u.d;
    };

    /*******************************************************************//**
    * \brief Fast error function by John D. Cook
    *
    * Algorithm comes from 
    *
    * @param a Number to raise to some power
    * @param b Power to which a should be raised
    * @return a^b 
    */
    inline double erf(double x)
    {
        // constants
        double a1 =  0.254829592;
        double a2 = -0.284496736;
        double a3 =  1.421413741;
        double a4 = -1.453152027;
        double a5 =  1.061405429;
        double p  =  0.3275911;

        // Save the sign of x
        int sign = 1;
        if (x < 0)
            sign = -1;
        x = fabs(x);

        // A&S formula 7.1.26
        double t = 1.0/(1.0 + p*x);
        double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

        return sign*y;
    };

};

#endif
