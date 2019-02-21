// Compare bisection optimiser from GRATE with GSL Brent's method

#include <iostream>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <functional>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>

// initial guess and tolerance
#define INITIAL_GUESS 0.8
#define TOL 0.00001


int fevals = 0;  // number of function evaluations

// this is the function to be maximised
double func(double x) {
    fevals++;
    return -1.0 * pow(x - 4.8, 2);
//    return sin(x);
}

// gsl does a minimisation, so we multiply the function we want to maximise by -1
double funcmingsl(double x, void* params = 0) {
    (void)(params); // avoid unused parameter warning
    return -1.0 * func(x);
}

// find the maximum using Brent's method from the GSL library
int find_max_brent() {
    // reset function evaluations
    fevals = 0;

    // constants
    const double gold = 1.618034;  // golden ratio

    /* ===== first step to bracket ===== */
    // bracketing adapted from: https://github.com/scipy/scipy/blob/v1.2.1/scipy/optimize/optimize.py#L2263-L2361

    // initial points (two are required)
    double xa = INITIAL_GUESS;
    double xb = xa * 1.1;
    double fa = funcmingsl(xa);
    double fb = funcmingsl(xb);

    // switch xa and xb so we go downhill from a to b
    if (fa < fb) {
        std::swap(xa, xb);
        std::swap(fa, fb);
    }

    // compute new point xc
    double xc = xb + gold * (xb - xa);
    double fc = funcmingsl(xc);
    std::cerr << "starting points: " << xa << ", " << xb << ", " << xc << std::endl;

    // continue until we have bracket the minimum (fa > fb < fc)
    // TODO: try parabolic extrapolation to get new point
    while (fb > fc) {
        // new point
        double w = xc + gold * (xc - xb);
        double fw = funcmingsl(w);

        xa = xb;
        fa = fb;
        xb = xc;
        fb = fc;
        xc = w;
        fc = fw;
    }

    if (xa > xc) {
        std::swap(xa, xc);
    }
    std::cout << "Found bracket: " << xa << ", " << xb << ", " << xc << std::endl;
    std::cout << "             : " << fa << ", " << fb << ", " << fc << std::endl;
    std::cout << "Num function calls after bracketing: " << fevals << std::endl;

    /* ===== second step to converge on minimum ===== */
    // https://www.gnu.org/software/gsl/doc/html/min.html

    // setup the GSL minimiser and function to be minimised
    const gsl_min_fminimizer_type *T = gsl_min_fminimizer_brent;
    gsl_min_fminimizer *s = gsl_min_fminimizer_alloc(T);
    gsl_function F;
    F.function = &funcmingsl;
    F.params = 0;

    // initialise minimiser with values found during bracketing
    gsl_min_fminimizer_set_with_values(s, &F, xb, fb, xa, fa, xc, fc);

    // perform the minimisation
    int status = 0;
    int iter = 0;
    do {
        iter++;
        status = gsl_min_fminimizer_iterate(s);

        xb = gsl_min_fminimizer_x_minimum(s);  // xb is the current minimum point
        xa = gsl_min_fminimizer_x_lower(s);  // xa and xc bracket the minimum
        xc = gsl_min_fminimizer_x_upper(s);

        // check if we have converged
        status = gsl_min_test_interval(xa, xc, TOL, 0.0);
        if (status == GSL_SUCCESS) {
            std::cout << "Converged" << std::endl;
        }
    } while (status == GSL_CONTINUE);

    if (status == GSL_SUCCESS) {
        double opt = xb;
        double optval = func(xb);
        std::cout << "Found maximum at " << opt << ", function value " << optval << std::endl;
        std::cout << "Number of iterations: " << iter << std::endl;
        std::cout << "Total GSL Func evals: " << fevals << std::endl;
    }
    else {
        std::cerr << "Error: failed to converge using GSL Brent method" << std::endl;
    }

    gsl_min_fminimizer_free(s);

    return fevals;
}
    

int find_max_original() {
    double test_plus = 0, test_minus = 0;
    double gradient_1 = 0, gradient_2 = 2;
    double p, p_upper = 0, p_lower = 0;
    double p1, p2;
    double converg;

    // reset function evaluations
    fevals = 0;

    // initial guess
    p = INITIAL_GUESS;

    // calculate gradient numerically
    test_plus = func(p * 1.001);
    test_minus = func(p * 0.999);
    gradient_1 = test_plus - test_minus;
    p1 = p;

    // Now move in the direction of the gradient
    if (gradient_1 > 0)
        p = p + 0.25 * p;  // TODO: what if p is zero???
    else
        p = p - 0.25 * p;

    // calculate gradient numerically
    test_plus = func(p * 1.001);
    test_minus = func(p * 0.999);
    gradient_2 = test_plus - test_minus;
    p2 = p;

    // find the bracket
    while( gradient_1 / gradient_2 > 0 )
    {
        gradient_1 = gradient_2;
        p1 = p;
        if (gradient_1 > 0)
            p = p + 0.25 * p;
        else
            p = p - 0.25 * p;

        test_plus = func(p * 1.001);
        test_minus = func(p * 0.999);
        gradient_2 = test_plus - test_minus;
        p2 = p;
    }

    p_upper = std::max( p1, p2 );
    p_lower = std::min( p1, p2 );
    p = 0.5 * ( p_upper + p_lower );
    converg = ( p_upper - p_lower ) / p;
    std::cout << "Have bracket: (" << p_lower << ", " << p_upper << ")" << std::endl;
    std::cout << "Num function calls after bracketing: " << fevals << std::endl;

    // converge
    while(converg > TOL)
    {
        test_plus = func(p * 1.001);
        test_minus = func(p * 0.999);
        double gradient = test_plus - test_minus;
        if ( gradient > 0 )
            p_lower = p;
        else
            p_upper = p;
        p = 0.5 * ( p_upper + p_lower );
        converg = ( p_upper - p_lower ) / p;
    }

    double opt = p;
    double optval = func(opt);
    std::cout << "Optimum found at p = " << opt << ", f(p) = " << optval << std::endl;
    std::cout << "Total function calls after opt: " << fevals << std::endl;

    return fevals;
}


int main() {
    std::cout << "============= Calling find_max_original =============\n";
    int orig_fevals = find_max_original();

    std::cout << "============= Calling find_max_brent =============\n";
    int brent_fevals = find_max_brent();

    std::cout << "============= Summary =============\n";
    std::cout << "Bisection optimiser used: " << orig_fevals << " force evaluations\n";
    std::cout << "Brent     optimiser used: " << brent_fevals << " force evaluations\n";

    return 0;
}
