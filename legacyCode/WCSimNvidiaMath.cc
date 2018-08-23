#include "WCSimNvidiaMath.hh"
#include <cmath>
// We end up calculating a LOT of trig functions
// This namespace contains some good, fast approximate ways
// of evaluating trig and inverse trig functions
//
// Algorithms come from the NVidia CG toolkit
// http://http.developer.nvidia.com/Cg/index_stdlib.html

float WCSimNvidiaMath::cos(float radians) {
	/* C simulation gives a max absolute error of less than 1.8e-7 */
	float c0[4] = { 0.0, 0.5, 1.0, 0.0 };
	float c1[4] = { 0.25, -9.0, 0.75, 0.159154943091 };
	float c2[4] = { 24.9808039603, -24.9808039603, -60.1458091736, 60.1458091736 };
	float c3[4] = { 85.4537887573, -85.4537887573, -64.9393539429, 64.9393539429 };
	float c4[4] = { 19.7392082214, -19.7392082214, -1.0, 1.0 };

	/* r0.x = sin(a) */
	float r0[3] = { 0, 0, 0 }, r1[3] = { 0, 0, 0 }, r2[3] = { 0, 0, 0 };
	r1[0] = c1[3] * radians;
	r1[1] = r1[0] - floor(r1[0]);

	r2[0] = (float) (r1[1] < c1[0]);
	r2[1] = (float) (r1[1] >= c1[1]);
	r2[2] = (float) (r1[1] >= c1[2]);

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
}

float WCSimNvidiaMath::sin(float radians) {
	/* C simulation gives a max absolute error of less than 1.8e-7 */
	float c0[4] = { 0.0, 0.5, 1.0, 0.0 };
	float c1[4] = { 0.25, -9.0, 0.75, 0.159154943091 };
	float c2[4] = { 24.9808039603, -24.9808039603, -60.1458091736, 60.1458091736 };
	float c3[4] = { 85.4537887573, -85.4537887573, -64.9393539429, 64.9393539429 };
	float c4[4] = { 19.7392082214, -19.7392082214, -1.0, 1.0 };

	/* r0.x = sin(a) */
	float r0[3] = { 0, 0, 0 }, r1[3] = { 0, 0, 0 }, r2[3] = { 0, 0, 0 };
	r1[0] = c1[3] * radians - c1[0];
	r1[1] = r1[0] - floor(r1[0]);

	r2[0] = (float) (r1[1] < c1[0]);
	r2[1] = (float) (r1[1] >= c1[1]);
	r2[2] = (float) (r1[1] >= c1[2]);

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
}

