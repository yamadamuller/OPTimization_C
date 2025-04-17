#include "./src/univariate_opt.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double obj_fun(double x){
    return 2*cos(x) -x; //invert the peak
}

int main(){
    double x = 0.5;
    struct OptResults opt_val = gradient_descent_1d(obj_fun, 1., 0.1, 30, 1e-5); 
    printf("%f \n", opt_val.opt_guess);
    printf("%f \n", opt_val.fun_eval);
}