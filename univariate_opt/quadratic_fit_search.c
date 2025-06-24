#include "./src/univariate_opt.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double obj_fun(double x){
    return 2*cos(x) -x; //invert the peak
}

int main(){
    struct OptResults opt_val = quadratic_fit_search(obj_fun, 1., 3., 6., 30, 1e-5); //compute the min. value
    printf("%f \n", opt_val.opt_guess);
    printf("%f \n", opt_val.fun_eval);
}