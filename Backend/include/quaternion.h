#ifndef QUATERNION_H
#define QUATERNION_H

#include <cmath>
#include <array>


typedef std::array<double,4> quaternion; ///< defined quaternion

/******************************************************
 * This function adds quaternions q1 and q2, and stores result in q.
*/
void qua_add(quaternion q1, quaternion q2, quaternion q);

/******************************************************
 * This function subtracts quaternions q2 from q1 (performs q1-q2), and stores result in q.
*/
void qua_sub(quaternion q1, quaternion q2, quaternion q);

/******************************************************
 * This function multiplies quaternions q1 and q2 (performs q1q2), and stores result in q.
*/
void qua_mult(quaternion q1, quaternion q2, quaternion q);

/******************************************************
 * This function conjugates quaternion q1, and stores result in q.
*/
void qua_conj(quaternion q1, quaternion q);

/******************************************************
 * This function takes dot product of quaternions q1 and q2 (performs q1 . q2), and returns it
*/
double qua_dot(quaternion q1, quaternion q2);

/******************************************************
 * This function inverts quaternion q1, and stores result in q.
*/
void qua_inv(quaternion q1, quaternion q);

/******************************************************
 * This function takes cross product of vectors q1 and q2 (performs q1.vec \vcross q2.vec), and stores result in q.
*/
void qua_vcross(std::array<double,3> q1, std::array<double,3> q2, std::array<double,3> q);

/******************************************************
 * This function normalizes the quaternion q1
*/
void qua_normalize(quaternion q1);

/*****************************************************
 * This function multiplies a matrix M and a vector double[3] v and stores the result in a quaternion out = Mq
*/

void matrix_mult_vec(std::array<std::array<double,3>,3> M, std::array<double,3> v, std::array<double,3> out);

/*****************************************************
 * This function generates the rotation matrix R corresponding to quaternion q. Only works when ||q|| = 1.
*/

void generate_rotmat(quaternion q, std::array<std::array<double,3>,3> R);

/*****************************************************
 * This function transposes R and stores in Rt
*/

void transpose_mat(std::array<std::array<double,3>,3> R, std::array<std::array<double,3>,3> Rt);

/******************************************************
 * This function multiplies quaternion q1 and vector v (performs q1v), and stores result in q.
*/
void qua_vec_mult(quaternion q1, std::array<double,3> v, quaternion q);

/******************************************************
 * This function returns the modulus of q
*/
double qua_mod(quaternion q);

/******************************************************
 * This function rotates vector v according to vout = q*vq.
*/
void qua_rot_vec_opp(quaternion q1, std::array<double,3> v, std::array<double,3> vout);

#endif