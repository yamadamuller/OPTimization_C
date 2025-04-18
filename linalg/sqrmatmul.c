#include "./src/linalg.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(){
    double A[3][3] = {{1,2,3},{4,5,6},{7,8,9}};
    double B[3][3] = {{-10, 5.2, 6.7},{-99, 36, 15.5},{62.4, 67.99, 31}};
    double C[3][3];

    linalg_sqrmatmul(A,B,C);

    //print the result matrix
    for(int i = 0; i<3; i++)
    {
        for(int j = 0; j<3; j++)
        {
            printf("%2f ", C[i][j]);
        }
        printf("\n");
    }
}