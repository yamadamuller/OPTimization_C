#include "./src/linalg.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(){
    double A[3][3] = {{1,2,-1},{2,-3,1},{5,-1,-2}};
    double b[3][1] = {{2},{-1},{-3}};
    double x[3][1];
    
    linalg_pseudoinverse(A, b, x);

    //print the result matrix
    for(int i = 0; i<3; i++)
    {
        printf("%f \n ", x[i][0]);
    }
}