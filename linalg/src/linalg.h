#ifndef LINALG_H
#define LINALG_H

void linalg_transpose(double A[3][3], double A_t[3][3]);
void linalg_sqrmatmul(double A[3][3], double B[3][3], double C[3][3]);
void linalg_vsqrmatmul(double A[3][3], double B[3][1], double C[3][1]);
double linalg_det(double A[3][3]);
void linalg_inverse(double A[3][3], double A_inv[3][3]);
void linalg_pseudoinverse(double A[3][3], double b[3][1], double x[3][1]);

#endif