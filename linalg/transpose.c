#include "./src/linalg.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(){
    double A[3][3] = {{1,1,1},{1,3,9},{1,6,36}};
    double A_t[3][3]; 

    linalg_transpose(A, A_t);

    //print the original matrix
    for(int i = 0; i<3; i++)
    {
        for(int j = 0; j<3; j++)
        {
            printf("%2f ", A[i][j]);
        }
        printf("\n");
    }

    printf("\n");

    //print the transposed matrix
    for(int i = 0; i<3; i++)
    {
        for(int j = 0; j<3; j++)
        {
            printf("%2f ", A_t[i][j]);
        }
        printf("\n");
    }
    

}