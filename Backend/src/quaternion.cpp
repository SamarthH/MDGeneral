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

void qua_vcross(std::array<double,3> q1, std::array<double,3> q2, std::array<double,3> q)
{
	double out[3];
	out[0] = q1[1]*q2[2] - q1[2]*q2[1];
	out[1] = q1[2]*q2[0] - q1[0]*q2[2];
	out[2] = q1[0]*q2[1] - q1[1]*q2[0];

	for (int i = 0; i < 3; ++i)
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

void matrix_mult_vec(std::array<std::array<double,3>,3> M, std::array<double,3> q, std::array<double,3> out)
{
	double tempout[3] = {0,0,0};
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			tempout[i] += M[i][j]*q[j];
		}
	}
	for (int i = 0; i < 3; ++i)
	{
		out[i] = tempout[i];
	}
}

void generate_rotmat(quaternion q, std::array<std::array<double,3>,3> R)
{
	double prod[4][4];
	for (int i = 0; i < 4; ++i)
	{
		for (int j = 0; j <= i; ++j)
		{
			prod[i][j] = prod[j][i] = q[i]*q[j];
		}
	}

	R[0][0] = 1 - 2*(prod[2][2]+prod[3][3]);
	R[0][1] = 2*(prod[1][2] - prod[3][0]);
	R[0][2] = 2*(prod[1][3] + prod[2][0]);
	R[1][0] = 2*(prod[1][2] + prod[3][0]);
	R[1][1] = 1 - 2*(prod[1][1]+prod[3][3]);
	R[1][2] = 2*(prod[2][3] - prod[1][0]);
	R[2][0] = 2*(prod[1][3] - prod[2][0]);
	R[2][1] = 2*(prod[2][3] + prod[1][0]);
	R[2][2] = 1 - 2*(prod[2][2]+prod[1][1]);
}

void transpose_mat(std::array<std::array<double,3>,3> R, std::array<std::array<double,3>,3> Rt)
{
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			Rt[i][j] = R[j][i];
		}
	}
}

void qua_vec_mult(quaternion q1, std::array<double,3> v, quaternion q)
{
	quaternion out;
	out[0] = 0;
	for (int i = 0; i < 3; ++i)
	{
		out[0] -= q1[i+1]*v[i];
	}

	out[1] = q1[0]*v[0] + (q1[2]*v[2] - q1[3]*v[1]);
	out[2] = q1[0]*v[1] + (q1[3]*v[0] - q1[1]*v[2]);
	out[3] = q1[0]*v[2] + (q1[1]*v[1] - q1[2]*v[0]);

	for (int i = 0; i < 4; ++i)
	{
		q[i] = out[i];
	}
}

double qua_mod(quaternion q)
{
	double mod = 0;
	for (int i = 0; i < 4; ++i)
	{
		mod += q[i]*q[i];
	}
	return mod;
}

void qua_rot_vec_opp(quaternion q1, std::array<double,3> v, std::array<double,3> vout)
{
	quaternion qconj;
	qua_conj(q1,qconj);
	quaternion qout;
	qua_vec_mult(qconj,v,qout);
	qua_mult(qout,q1,qout);
	for (int i = 0; i < 3; ++i)
	{
		vout[i] = qout[i+1];
	}
}