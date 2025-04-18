#include "./src/linalg.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(){
    double A[3][3] = {{3,5,6},{7,8,9},{12,15,12}};
    double B[3][1] = {{1},{4},{7}};
    double C[3][1];

    linalg_vsqrmatmul(A,B,C);

    //print the result matrix
    for(int i = 0; i<3; i++)
    {
        printf("%f \n", C[i][0]);
    }
}