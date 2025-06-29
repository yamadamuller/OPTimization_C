#ifndef UNIVARIATE_OPT_H
#define UNIVARIATE_OPT_H

struct OptResults
{
    double fun_eval; //minimum value of the objective function
    double opt_guess; //candidate value that yielded the opt. value
    int n_iter; //number of iterations
    double grad; //gradient at the optimal value
    double hess; //hessian at the optimal value
};

struct OptResults gridsearch_1d(double (*opt_fun)(double), double lower_lim, double upper_lim, double step);
struct OptResults golden_section_search(double (*opt_fun)(double), double lower_lim, double upper_lim, double tol);
double numerical_gradient(double (*opt_fun)(double), double x, double eps);
double numerical_hessian(double (*opt_fun)(double), double x, double eps);
struct OptResults gradient_descent_1d(double (*opt_fun)(double), double x_guess, double gamma, int max_iter, double tol);
struct OptResults newton_1d(double (*opt_fun)(double), double x_guess, int max_iter, double tol);
struct OptResults quadratic_fit_search(double (*opt_fun)(double), double a, double b, double c, int max_inter, double tol);
int evaluate_sign(double x);
struct OptResults bissection_method(double (*opt_fun)(double), double lower_lim, double upper_lim, int max_inter, double tol);

#endif
