#pragma once

#include <cmath>
#include <cassert>
#include <iostream>
// We end up calculating a LOT of trig functions and awkward
// powers, exponentials, error functions, etc.
// This namespace contains some good, fast approximate ways
// of evaluating mathematical expressions to a level of accuracy
// the fitter can live with

namespace WCSimFastMath
{
/*******************************************************************/ /**
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
	float c2[4] = {24.9808039603, -24.9808039603, -60.1458091736, 60.1458091736};
	float c3[4] = {85.4537887573, -85.4537887573, -64.9393539429, 64.9393539429};
	float c4[4] = {19.7392082214, -19.7392082214, -1.0, 1.0};

	/* r0.x = sin(a) */
	float r0[3] = {0, 0, 0}, r1[3] = {0, 0, 0}, r2[3] = {0, 0, 0};
	r1[0] = c1[3] * radians;
	r1[1] = r1[0] - floor(r1[0]);

	r2[0] = (float)(r1[1] < c1[0]);
	r2[1] = (float)(r1[1] >= c1[1]);
	r2[2] = (float)(r1[1] >= c1[2]);

	r2[1] = (r2[0] * c4[2] + r2[1] * c4[3] + r2[2] * c4[2]);

	r0[0] = c0[0] - r1[1];
	r0[1] = c0[1] - r1[1];
	r0[2] = c0[2] - r1[1];

	r0[0] = r0[0] * r0[0];
	r0[1] = r0[1] * r0[1];
	r0[2] = r0[2] * r0[2];

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

	r0[0] = r1[0] * -r2[0] + r1[1] * -r2[1] + r1[2] * -r2[2];
	return r0[0];
};

/*******************************************************************/ /**
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
	float c2[4] = {24.9808039603, -24.9808039603, -60.1458091736, 60.1458091736};
	float c3[4] = {85.4537887573, -85.4537887573, -64.9393539429, 64.9393539429};
	float c4[4] = {19.7392082214, -19.7392082214, -1.0, 1.0};

	/* r0.x = sin(a) */
	float r0[3] = {0, 0, 0}, r1[3] = {0, 0, 0}, r2[3] = {0, 0, 0};
	r1[0] = c1[3] * radians - c1[0];
	r1[1] = r1[0] - floor(r1[0]);

	r2[0] = (float)(r1[1] < c1[0]);
	r2[1] = (float)(r1[1] >= c1[1]);
	r2[2] = (float)(r1[1] >= c1[2]);

	r2[1] = (r2[0] * c4[2] + r2[1] * c4[3] + r2[2] * c4[2]);

	r0[0] = c0[0] - r1[1];
	r0[1] = c0[1] - r1[1];
	r0[2] = c0[2] - r1[1];

	r0[0] = r0[0] * r0[0];
	r0[1] = r0[1] * r0[1];
	r0[2] = r0[2] * r0[2];

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

	r0[0] = r1[0] * -r2[0] + r1[1] * -r2[1] + r1[2] * -r2[2];
	return r0[0];
};

/*******************************************************************/ /**
 * \brief Fast error function by John D. Cook
 *
 * @param a Number to raise to some power
 * @param b Power to which a should be raised
 * @return a^b
 */
inline double erf(double x)
{
	// constants
	double a1 = 0.254829592;
	double a2 = -0.284496736;
	double a3 = 1.421413741;
	double a4 = -1.453152027;
	double a5 = 1.061405429;
	double p = 0.3275911;

	// Save the sign of x
	int sign = 1;
	if (x < 0)
		sign = -1;
	x = fabs(x);

	// A&S formula 7.1.26
	double t = 1.0 / (1.0 + p * x);
	double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-x * x);

	return sign * y;
};

/**
 * @brief Use a Catmull-Rom spline to smoothly interpolate y = f(x) between the middle two of four evenly-separated x values
 *
 * A Catmull-Rom spline performs a cubic interpolation between the central two points
 * of an array of four (x,y) pairs.  This implementation requires the points be in
 * ascending x order and evenly-separated in x.
 * Catmull-Rom splines are cubic and are guaranteed to be continuous in x and d/dx
 * They are local, requiring only the two points either side of the interpolated value
 * They satisfy f(x1) = y1 and f(x2 = y2)
 * Also, f'(x1) = (y2 - y0)/(x2 - x0) and f'(x2) = (y3 - y1)/(x3 - x1)
 *
 * @param x Array containing four x-coordinates, in ascending order
 * @param y Array containing four corresponding y-coordinates
 * @param xInterp x-coordinate at which to interpolate.  Must lie between x[1] and x[2] inclusive.
 *
 * @return Smoothly interpolated value for y(xInterp)
 */
inline double CatmullRomSpline(double *x, double *y, double xInterp)
{
	// Check for ascending x with even spacing
	double dx1 = x[1] - x[0];
	double dx2 = x[2] - x[1];
	double dx3 = x[3] - x[2];
	if (dx1 < 0)
	{
		assert(dx1 > 0);
	}
	if ((dx3 != dx2) || (dx2 != dx1))
	{
		std::cout << "We want to interpolate at " << xInterp << " and we're using " << x[0] << " " << x[1] << " "
				  << x[2] << " and " << x[3] << std::endl;
		std::cout << "dx1 = " << dx1 << " dx2 = " << dx2 << "dx3 = " << dx3 << std::endl;
		assert((dx3 == dx2) && (dx2 == dx1));
	}

	// Check the value to interpolate at is between x1 and x2
	if (xInterp > x[2] || xInterp < x[1])
	{
		std::cout << "Want to interpolate " << xInterp << " but the bounds we have are " << x[1] << " and " << x[2]
				  << std::endl;
		assert(xInterp <= x[2] && xInterp >= x[1]);
	}

	// Tranform coordinates such that x1->0 and x2->1
	double t = (xInterp - x[1]) / (x[2] - x[1]);

	// Now do the multiplication by the coefficients of the cubic that handles the interpolation
	// Because we transformed coordinates the matrix inversion and multiplication gives nice integer coefficients
	return (0.5 * (2. * y[1] + t * ((y[2] - y[0]) + t * ((2. * y[0] - 5. * y[1] + 4. * y[2] - y[3]) + t * (3. * y[1] - y[0] - 3. * y[2] + y[3])))));
}

}; // namespace WCSimFastMath
