#include "./src/linalg.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(){
    double A[3][3] = {{0,-3,-2},{1,-4,-2},{-3,4,1}};
    double A_inv[3][3];
    
    linalg_inverse(A, A_inv);

    //print the result matrix
    for(int i = 0; i<3; i++)
    {
        for(int j = 0; j<3; j++)
        {
            printf("%2f ", A_inv[i][j]);
        }
        printf("\n");
    }
}