#include "linalg.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

void linalg_transpose(double A[3][3], double A_t[3][3]){
    /*
    :param A: pointer to the original matrix that will be transposed
    :param A_t: pointer to the to the transposed matrix
    */
   int len = sizeof(A[0]) / sizeof(A[0][0]); //assumes a square matrix
    for(int i=0; i<len; i++){
        for(int j=0; j<len; j++){
            A_t[i][j] = A[j][i]; //transpose
        }
    }   
}

void linalg_sqrmatmul(double A[3][3], double B[3][3], double C[3][3]){
    /*
    :param A: pointer to the first matrix of the multiplication
    :param B: pointer to the first matrix of the multiplication
    :param C: pointer to the result of the multiplication
    */
    int len = sizeof(A[0]) / sizeof(A[0][0]); //assumes a square matrix
    for(int i=0; i<len; i++){
        for(int j=0; j<len; j++){
            for(int k=0; k<len; k++){
                C[i][j] += A[i][k]*B[k][j]; //sum of the products
            }
        }
    }
}

void linalg_vsqrmatmul(double A[3][3], double B[3][1], double C[3][1]){
    /*
    :param A: pointer to the first matrix of the multiplication
    :param B: pointer to the first matrix of the multiplication
    :param C: pointer to the result of the multiplication
    */
    int len = sizeof(A[0]) / sizeof(A[0][0]); //assumes a square matrix
    for(int i=0; i<len; i++){
        for(int j=0; j<len; j++){
            C[i][0] += A[i][j]*B[j][0]; //sum of the products
        }
    }
}

double linalg_det(double A[3][3]){
    /*
    :param A: matrix to compute the determinant
    :return the determinant of the matrix A
    */

    double forward = (A[0][0]*A[1][1]*A[2][2]+A[0][1]*A[1][2]*A[2][0]+A[0][2]*A[1][0]*A[2][1]);
    double backward = (-1*(A[0][1]*A[1][0]*A[2][2])-1*(A[0][0]*A[1][2]*A[2][1])-1*(A[0][2]*A[1][1]*A[2][0]));
    return forward+backward;
}

void linalg_inverse(double A[3][3], double A_inv[3][3]){
    /*
    :param A: pointer to the matrix that will be inversed
    :param A_inv: pointer to the inversed matrix 
    */
    double det = linalg_det(A); //compute the determinant
    
    if(det==0){
        perror("Matrix is singular!");
    }

    double det_inv = 1/det; //inverse of the determinant

    //define the matrix compose by the minor submatrices
    A_inv[0][0] = (A[1][1] * A[2][2] - A[2][1] * A[1][2]) * det_inv;
    A_inv[0][1] = (A[0][2] * A[2][1] - A[0][1] * A[2][2]) * det_inv;
    A_inv[0][2] = (A[0][1] * A[1][2] - A[0][2] * A[1][1]) * det_inv;
    A_inv[1][0] = (A[1][2] * A[2][0] - A[1][0] * A[2][2]) * det_inv;
    A_inv[1][1] = (A[0][0] * A[2][2] - A[0][2] * A[2][0]) * det_inv;
    A_inv[1][2] = (A[1][0] * A[0][2] - A[0][0] * A[1][2]) * det_inv;
    A_inv[2][0] = (A[1][0] * A[2][1] - A[2][0] * A[1][1]) * det_inv;
    A_inv[2][1] = (A[2][0] * A[0][1] - A[0][0] * A[2][1]) * det_inv;
    A_inv[2][2] = (A[0][0] * A[1][1] - A[1][0] * A[0][1]) * det_inv;   
}

void linalg_pseudoinverse(double A[3][3], double b[3][1], double x[3][1]){
    /*
    :param A: square matrix with the coefficients of the quadratic equation
    :param b: array with the function evaluation coefficients
    :param x: array with the function coefficients
    */
    
    //Transposition operation
    double A_t[3][3]; //transposed A
    linalg_transpose(A, A_t);

    //(A.T@A)⁻1
    double AtA[3][3]; 
    linalg_sqrmatmul(A_t, A, AtA);
    double AtA_inv[3][3];
    linalg_inverse(AtA, AtA_inv);

    //(A.T@A)⁻1@A.T
    double A_mult[3][3];
    linalg_sqrmatmul(AtA_inv, A_t, A_mult);

    //(A.T@A)⁻1@A.T@b
    linalg_vsqrmatmul(A_mult, b, x); //pseudo inverse output

}