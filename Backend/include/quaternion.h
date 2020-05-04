#ifndef QUATERNION_H
#define QUATERNION_H

#include <cmath>
#include <array>


typedef double quaternion[4]; ///< defined quaternion

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
void qua_vcross(std::array<double,3> q1, std::array<double,3> q2, quaternion q);

/******************************************************
 * This function normalizes the quaternion q1
*/
void qua_normalize(quaternion q1);

#endif