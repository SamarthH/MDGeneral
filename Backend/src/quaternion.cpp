#include "quaternion.h"

void qua_add(quaternion q1, quaternion q2, quaternion q)
{
	for (int i = 0; i < 4; ++i)
	{
		q[i] = q1[i] + q2[i];
	}
}

void qua_sub(quaternion q1, quaternion q2, quaternion q)
{
	for (int i = 0; i < 4; ++i)
	{
		q[i] = q1[i] - q2[i];
	}
}

void qua_mult(quaternion q1, quaternion q2, quaternion q)
{
	double out[4];
	out[0] = q1[0]*q2[0] - (q1[1]*q2[1] + q1[2]*q2[2] + q1[3]*q2[3]);
	out[1] = q1[0]*q2[1] + q2[0]*q1[1] + (q1[2]*q2[3] - q1[3]*q2[2]);
	out[2] = q1[0]*q2[2] + q2[0]*q1[2] + (q1[3]*q2[1] - q1[1]*q2[3]);
	out[3] = q1[0]*q2[3] + q2[0]*q1[3] + (q1[1]*q2[2] - q1[2]*q2[1]);
	for (int i = 0; i < 4; ++i)
	{
		q[i] = out[i];
	}
}

void qua_conj(quaternion q1, quaternion q)
{
	q[0] = q1[0];
	for (int i = 1; i < 4; ++i)
	{
		q[i] = -q1[i];
	}
}

double qua_dot(quaternion q1, quaternion q2)
{
	double out = 0;
	for (int i = 0; i < 4; ++i)
	{
		out += q1[i]*q2[i];
	}
	return out;
}

void qua_inv(quaternion q1, quaternion q)
{
	double norm = q1[0]*q1[0] + q1[1]*q1[1] + q1[2]*q1[2] + q1[3]*q1[3];
	for (int i = 0; i < 4; ++i)
	{
		q[i] = q1[i]/norm;
	}
}

void qua_vcross(std::array<double,3> q1, std::array<double,3> q2, quaternion q)
{
	quaternion out;
	out[0] = 0;
	out[1] = q1[1]*q2[2] - q1[2]*q2[1];
	out[2] = q1[2]*q2[0] - q1[0]*q2[2];
	out[3] = q1[0]*q2[1] - q1[1]*q2[0];

	for (int i = 0; i < 4; ++i)
	{
		q[i] = out[i];
	}
}

void qua_normalize(quaternion q1)
{
	double norm = 0;
	for (int i = 0; i < 4; ++i)
	{
		norm += q1[i]*q1[i];
	}
	norm = std::sqrt(norm);
	for (int i = 0; i < 4; ++i)
	{
		q1[i] /= norm;
	}
}