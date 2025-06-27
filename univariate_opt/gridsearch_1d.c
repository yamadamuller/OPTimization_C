#include "./src/univariate_opt.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double obj_fun(double x){
    return 2*cos(x) -x; //invert the peak
}

int main(){
    struct OptResults opt_val = gridsearch_1d(obj_fun, 1., 6., 0.1); 
    printf("%f \n", opt_val.opt_guess);
    printf("%f \n", opt_val.fun_eval);
    printf("%d \n", opt_val.n_iter);
}