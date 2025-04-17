#include "univariate_opt.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct OptResults gridsearch_1d(double (*opt_fun)(double), double lower_lim, double upper_lim, double step){
    /*
        :param opt_fun: pointer to the objective target function
        :param lower_lim: lower limit of the search interval
        :param upper_lim: upper limit of the search interval
        :param step: step between candidates
    */

   double opt_value = __FLT_MAX__; //variable to store the minimum value
   double opt_guess = lower_lim; //variable to store the optimum guess
   int num_candidates = ceil((upper_lim-lower_lim)/step); //number of candidate values 
   for(int i=0; i<=num_candidates; i++){
        double curr_obj = opt_fun(lower_lim); //compute the objective function for the current guess

        //check if the function has the mininum value in the interval
        if(curr_obj<opt_value){
            opt_value = curr_obj; //comute the opt. value as the current obj. function value
            opt_guess = lower_lim; //update the current optimum value
        }

        lower_lim += step; //update the guess for the next iteration
   }

   struct OptResults opt_output = {opt_value, opt_guess, num_candidates, NAN, NAN}; //store the results in the opt struct

   return opt_output; //return the optimum guess
}

struct OptResults golden_section_search(double (*opt_fun)(double), double lower_lim, double upper_lim, double tol){
    /*
        :param opt_fun: pointer to the objective target function
        :param lower_lim: lower limit of the search interval
        :param upper_lim: upper limit of the search interval
        :param tol: tolerance of the forward error (b-a < tol)
    */

    double c = (-1+sqrt(5))/2; //constant reduction rate
    double c_sub = 1-c; //reduction rate subtracted by 1
    double x1 = lower_lim*c + (upper_lim*c_sub); //first guess of x1
    double x2 = upper_lim*c + (lower_lim*c_sub); //first guess of x2
    double h = upper_lim - lower_lim; //forward error
    int n = (int)(ceil(log(tol/h)/log(c))); //iterations required to convergence

    //gss agorithm
    for(int i=0; i<n; i++){
        double fx1 = opt_fun(x1); //evaluate the obj. function at the x1 guess
        double fx2 = opt_fun(x2); //evaluate the obj. function at the x2 guess
        
        //bracketing
        if(fx1<fx2){
            upper_lim = x2; //update the upper limit
            x2 = x1; //comute the guesses
            x1 = lower_lim*c + upper_lim*c_sub; //update the x1 guess
        }
        else{
            lower_lim = x1; //update the lower limit
            x1 = x2; //comute the guesses
            x2 = upper_lim*c + lower_lim*c_sub; //update the x2 guess
        }
    }

    double opt_guess = (lower_lim+upper_lim)/2; //optimal guess based on the interval
    struct OptResults opt_output = {opt_fun(opt_guess), opt_guess, n, NAN, NAN}; //store the results in the opt struct

    return opt_output; //return the optimum guess
}

double numerical_gradient(double (*opt_fun)(double), double x, double eps){
    /*
        :param opt_fun: pointer to the objective target function
        :param x: the current guess of the minimizer
        :param eps: epsilon
        :return the first order numerical derivative of the objective function
    */
   
    return (opt_fun(x+eps)-opt_fun(x-eps))/(2*eps);
}

double numerical_hessian(double (*opt_fun)(double), double x, double eps){
    /*
        :param opt_fun: pointer to the objective target function
        :param x: the current guess of the minimizer
        :param eps: epsilon
        :return the second order numerical derivative of the objective function
    */

    return (opt_fun(x+eps)-(2*opt_fun(x))+opt_fun(x-eps))/(pow(eps,2));
}

struct OptResults gradient_descent_1d(double (*opt_fun)(double), double x_guess, double gamma, int max_iter, double tol){
    /*
        :param opt_fun: pointer to the objective function
        :param x_guess: the initial guess of the minimizer
        :param gamma: learning rate
        :param max_iter: maximum iteration of the newton-raphson algorithm
        :param tol: tolerance of the minimizer
        :returns a struct with information on the optimization task 
    */

   double curr_grad = 1; //evaluates the evolution of the (1/f'(x)) term
   int n_iter = 0;
   for(int i = 0; i<=max_iter; i++){
        //breaks if the first order derivative is close to zero
        if(abs(curr_grad)<tol){
            break;
        }
        
        curr_grad = numerical_gradient(opt_fun, x_guess, tol); //compute the first order derivative
        x_guess -= gamma*curr_grad; //update the guess with the learning rate
        n_iter = i+1; //update the iteration number
   }

   struct OptResults opt_output = {opt_fun(x_guess), x_guess, n_iter, curr_grad, NAN}; //store the results in the opt struct

   return opt_output; //return the minimum value
}

struct OptResults newton_1d(double (*opt_fun)(double), double x_guess, int max_iter, double tol){
    /*
        :param opt_fun: pointer to the objective function
        :param x_guess: the initial guess of the minimizer
        :param max_iter: maximum iteration of the newton-raphson algorithm
        :param tol: tolerance of the minimizer
        :returns a struct with information on the optimization task 
    */

    double n_term = 1; //evaluates the evolution of the (f'(x)/f''(x)) term
    double curr_grad = 1; //first order derivative
    double curr_hess = 1; //second order derivative
    int n_iter = 0;
    for(int i = 0; i<=max_iter; i++){
        //breaks if the newton term is significantly small
        if (abs(n_term) < tol){
            break; 
        }
        curr_grad = numerical_gradient(opt_fun, x_guess, tol); //compute the first order derivative
        curr_hess = numerical_hessian(opt_fun, x_guess, tol); //compute the second order derivative
        n_term = curr_grad/curr_hess; //compute the newton term
        x_guess -= n_term; //updates the guess
        n_iter = i+1; //update the iteration number
    }

    struct OptResults opt_output = {opt_fun(x_guess), x_guess, n_iter, curr_grad, curr_hess}; //store the results in the opt struct

    return opt_output; //return the guess
}

