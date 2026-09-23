//This file is for various miscellaneous general math functions that we have written for the code and may need to be used elsewhere
///////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "nbody_math_funcs.h"
// ODE solver, orginially written for the king model

// Generic function that numerically solves 2nd order ODEs of the form y'' = f(x, y(x), y'(x)), where y' is dy/dx
// Uses 4th order Runge-Kutta numerical method, function input is of the form of ODE2ndDeriv
// The last parameter returnXWhen0 is a special case boolean where the function will instead return the x value where the otherwise positive y(x) function hits zero/negative
real ODE2ndOrderSolver(real xEval, int stepsPerx, real yInit, real yPrimeInit, ODE2ndDeriv f, Dwarf* params, int returnXWhen0) {
    real nSteps = floor(stepsPerx * xEval);
    real stepRes = stepsPerx*xEval - nSteps;
    real deltax = xEval/nSteps;
    
    real y = yInit; // current value for y(x)
    real z = yPrimeInit; // let z be y'(x)
    real* yCurr = &y;
    real* zCurr = &z;
    real k1, k2, k3, k4, l1, l2, l3, l4;
    real x = deltax;
    // the the last loop needs to be one before the desired number of steps, since the end of a loop gives the values for the next step
    for (int n = 1; n < nSteps + 1; n++) {
        // l factors adjust z, k factors adjust y
        x = n * deltax;

        l1 = deltax * f(x, y, z, params);
        k1 = deltax * z;

        l2 = deltax * f(x + 0.5*deltax, y + 0.5*k1, z + 0.5*l1, params);
        k2 = deltax * (z + 0.5*l1);

        l3 = deltax * f(x + 0.5*deltax, y + 0.5*k2, z + 0.5*l2, params);
        k3 = deltax * (z + 0.5*l2);

        l4 = deltax * f(x + deltax, y + k3, z + l3, params);
        k4 = deltax * (z + l3);

        *yCurr = *yCurr + (1.0/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4);
        *zCurr = *zCurr + (1.0/6.0)*(l1 + 2.0*l2 + 2.0*l3 + l4);

        if ((*yCurr <= 0.0 || isnan(*yCurr) || !isfinite(*yCurr)) && returnXWhen0 == 1) {
            stepRes = 0.0;
            *yCurr = x;
            break;
        }

    }

    // run the iteration one more time if there is a fractional step at the end
    if (stepRes > 0.0) {
        deltax = stepRes/stepsPerx; // smaller than regular deltax
        l1 = deltax * f(x, y, z, params);
        k1 = deltax * z;
        l2 = deltax * f(x + 0.5*deltax, y + 0.5*k1, z + 0.5*l1, params);
        k2 = deltax * (z + 0.5*l1);
        l3 = deltax * f(x + 0.5*deltax, y + 0.5*k2, z + 0.5*l2, params);
        k3 = deltax * (z + 0.5*l2);
        l4 = deltax * f(x + deltax, y + k3, z + l3, params);
        k4 = deltax * (z + l3);
        *yCurr = *yCurr + (1.0/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4);
        *zCurr = *zCurr + (1.0/6.0)*(l1 + 2.0*l2 + 2.0*l3 + l4);
    }
    return *yCurr;
}

//Gamma functions and error functions

// IMPLEMENTATION OF GAMMA FUNCTIONS. COMPLETE AND INCOMPLETE
// These functions are adapted from the Numerical Recipes in C 2nd ed except for gammln.
/* Returns ln(Gamma(z)) via the Lanczos approximation (NR 3rd ed). */
real gammln(const real z) 
{
    /* exact memo: gammln is pure and phase-1 calls it millions of times
     * with only a couple of distinct arguments (delta+2, delta+3). */
    static __thread real gammln_last_z = -1.0e308;
    static __thread real gammln_last_v = 0.0;
    if (z == gammln_last_z) return gammln_last_v;
    //Alogrithm for the calculation of the Lanczos Approx of the complete Gamma function 
    //as implemented in Numerical Recipes 3rd ed, 2007.
    real g = 4.7421875; //g parameter for the gamma function
    real x, tmp, y, A_g;
    
    //these are the cn's
    static const real coeff[14] = {57.1562356658629235,-59.5979603554754912,
                                14.1360979747417471,-0.491913816097620199,.339946499848118887e-4,
                                .465236289270485756e-4,-.983744753048795646e-4,.158088703224912494e-3,
                                -.210264441724104883e-3,.217439618115212643e-3,-.164318106536763890e-3,
                                .844182239838527433e-4,-.261908384015814087e-4,.368991826595316234e-5};
    y = x = z;
    tmp = x + g + 0.5;
    tmp = (x + 0.5) * mw_log(tmp) - tmp;
    A_g = 0.999999999999997092; //this is c0
    
    for (int j = 0; j < 14; j++) 
    {
        A_g += coeff[j] / ++y;
    } //calculates the series approx sum
        
    //sqrt(2 * pi) = 2.5066282746310005
    tmp += mw_log(2.5066282746310005 * A_g / x);//returns the log of the gamma function
    
    return tmp;
}

/*
 * Returns the incomplete gamma function P(a, x), evaluated by its
 * series representation. Also returns ln(Gamma(a)) through gln.
 */
static real gser(const real a, const real x, real* gln)
{
    const int itmax = 100;
    const real eps = (real) 3.0e-7;
    real sum, del, ap;

    *gln = gammln(a);
    if (x <= 0.0)
    {
        return 0.0;
    }

    ap = a;
    del = sum = 1.0 / a;
    for (int n = 1; n <= itmax; ++n)
    {
        ap += 1.0;
        del *= x / ap;
        sum += del;
        if (mw_fabs(del) < mw_fabs(sum) * eps)
        {
            return sum * mw_exp(-x + a * mw_log(x) - (*gln));
        }
    }

    mw_printf("WARNING: gser did not converge (a=%f, x=%f)\n", a, x);
    return sum * mw_exp(-x + a * mw_log(x) - (*gln));
}

/*
 * Returns the incomplete gamma function Q(a, x), evaluated by its
 * continued-fraction representation (modified Lentz method).
 * Also returns ln(Gamma(a)) through gln.
 */
static real gcf(const real a, const real x, real* gln)
{
    const int itmax = 100;
    const real eps = (real) 3.0e-7;
    const real fpmin = (real) 1.0e-30;
    real an, b, c, d, del, h;

    *gln = gammln(a);
    b = x + 1.0 - a;
    if (mw_fabs(b) < fpmin)
    {
        b = fpmin;
    }
    c = 1.0 / fpmin;
    d = 1.0 / b;
    h = d;

    for (int i = 1; i <= itmax; ++i)
    {
        an = -(real) i * ((real) i - a);
        b += 2.0;
        d = an * d + b;
        if (mw_fabs(d) < fpmin)
        {
            d = fpmin;
        }
        c = b + an / c;
        if (mw_fabs(c) < fpmin)
        {
            c = fpmin;
        }
        d = 1.0 / d;
        del = d * c;
        h *= del;
        if (mw_fabs(del - 1.0) < eps)
        {
            return mw_exp(-x + a * mw_log(x) - (*gln)) * h;
        }
    }

    mw_printf("WARNING: gcf did not converge (a=%f, x=%f)\n", a, x);
    return mw_exp(-x + a * mw_log(x) - (*gln)) * h;
}

/*
 * Returns the regularized lower incomplete gamma function:
 * P(a, x) = gamma(a, x) / Gamma(a).
 * Uses series for x < a + 1, otherwise returns 1 - Q(a, x).
 */
real gammp(const real a, const real x)
{
    real gln;

    if (x < 0.0 || a <= 0.0)
    {
        mw_printf("WARNING: Invalid arguments in gammp (a=%f, x=%f)\n", a, x);
        return NAN;
    }

    if (x < (a + 1.0))
    {
        return gser(a, x, &gln);
    }
    return 1.0 - gcf(a, x, &gln);
}

/*
 * Returns the regularized upper incomplete gamma function:
 * Q(a, x) = Gamma(a, x) / Gamma(a) = 1 - P(a, x).
 * Uses series for x < a + 1, otherwise continued fraction.
 */
real gammq(const real a, const real x)
{
    real gln;

    if (!(x >= 0.0 && a > 0.0))
    {
        mw_printf("WARNING: Invalid arguments in gammq (a=%f, x=%f)\n", a, x);
        return NAN;
    }

    if (x < (a + 1.0))
    {
        return 1.0 - gser(a, x, &gln);
    }
    return gcf(a, x, &gln);
}

/* Returns the complete gamma function Gamma(z). */
real GammaFunc(const real z)
{
    return mw_exp(gammln(z));
}

/* Returns the upper incomplete gamma function Gamma(a, x). */
real UpperIncompleteGammaFunc(real a, real x)
{
    return GammaFunc(a) * gammq(a, x);
}

/* Returns the lower incomplete gamma function Gamma(a, x). */
real LowerIncompleteGammaFunc(real a, real x)
{
    return GammaFunc(a) * gammp(a, x);
}

/*
 * Numerical Recipes erff(x):
 * returns erf(x) using the regularized incomplete gamma P(1/2, x^2).
 */
real ErrorFunc(real x)
{
    const real xsq = x * x;
    return (x < 0.0) ? -gammp(0.5, xsq) : gammp(0.5, xsq);
}

/*
 * Numerical Recipes erffc(x):
 * returns erfc(x) using P(1/2, x^2) for x < 0 and Q(1/2, x^2) for x >= 0.
 */
real ComplementaryErrorFunc(real x)
{
    const real xsq = x * x;
    return (x < 0.0) ? 1.0 + gammp(0.5, xsq) : gammq(0.5, xsq);
}

/*
 * Numerical Recipes erfcc(x):
 * fast erfc approximation with fractional error < 1.2e-7.
 */
real ComplementaryErrorFuncApprox(real x)
//Note that this routine is for single precision
//Maybe look into Cody's Rational Chebyshev Approximation if double precision is needed
{
    real z = mw_fabs(x);
    real t = 1.0 / (1.0 + 0.5 * z);
    real ans = t * mw_exp(-z * z - 1.26551223
                          + t * (1.00002368
                          + t * (0.37409196
                          + t * (0.09678418
                          + t * (-0.18628806
                          + t * (0.27886807
                          + t * (-1.13520398
                          + t * (1.48851587
                          + t * (-0.82215223
                          + t * 0.17087277)))))))));
    return (x >= 0.0) ? ans : 2.0 - ans;
}


//First and second derivative functions from mixeddwarf

/*      GENERAL PURPOSE DERIVATIVE, INTEGRATION, MAX FINDING, ROOT FINDING, AND ARRAY SHUFFLER FUNCTIONS        */
real first_derivative(real (*func)(const Dwarf*, real), real x, const Dwarf* comp1)
{
    /* centered 5-point stencil when it stays in-domain; near the origin, switch to a forward stencil so probes never go negative.*/
    const real h = 0.001;

    if (x >= 2.0 * h)
    {
        real p1 =   1.0 * (*func)(comp1, (x - 2.0 * h));
        real p2 = - 8.0 * (*func)(comp1, (x - h) );
        real p3 = - 1.0 * (*func)(comp1, (x + 2.0 * h));
        real p4 =   8.0 * (*func)(comp1, (x + h));
        return (p1 + p2 + p3 + p4) * inv(12.0 * h);
    }

    /* Forward 5-point first derivative, O(h^4) */
    real f0 = (*func)(comp1, x);
    real f1 = (*func)(comp1, x + h);
    real f2 = (*func)(comp1, x + 2.0 * h);
    real f3 = (*func)(comp1, x + 3.0 * h);
    real f4 = (*func)(comp1, x + 4.0 * h);
    return (-25.0 * f0 + 48.0 * f1 - 36.0 * f2 + 16.0 * f3 - 3.0 * f4) * inv(12.0 * h);
}

real second_derivative(real (*func)(const Dwarf*, real), real x, const Dwarf* comp1)
{
    /* Same domain handling as first_derivative: centered away from 0, forward near 0. */
    const real h = 0.001;

    if (x >= 2.0 * h)
    {
        real p1 = - 1.0 * (*func)(comp1, (x + 2.0 * h));
        real p2 =  16.0 * (*func)(comp1, (x + h));
        real p3 = -30.0 * (*func)(comp1, (x));
        real p4 =  16.0 * (*func)(comp1, (x - h));
        real p5 = - 1.0 * (*func)(comp1, (x - 2.0 * h));
        return (p1 + p2 + p3 + p4 + p5) * inv(12.0 * h * h);
    }

    /* Forward 5-point second derivative, O(h^3) */
    real f0 = (*func)(comp1, x);
    real f1 = (*func)(comp1, x + h);
    real f2 = (*func)(comp1, x + 2.0 * h);
    real f3 = (*func)(comp1, x + 3.0 * h);
    real f4 = (*func)(comp1, x + 4.0 * h);
    return (35.0 * f0 - 104.0 * f1 + 114.0 * f2 - 56.0 * f3 + 11.0 * f4) * inv(12.0 * h * h);
}

real gauss_quad(real (*func)(real, const Dwarf*, const Dwarf*, real, mwbool), real lower, real upper, const Dwarf* comp1, const Dwarf* comp2, real energy, mwbool isDark)
{
    /*This is a guassian quadrature routine. It will test to always integrate from the lower to higher of the two limits.
     * If switching the order of the limits was needed to do this then the negative of the integral is returned.
     *
     * Restructured for parallel evaluation: the node positions are
     * enumerated first with EXACTLY the same arithmetic the original
     * sequential loop used, the (pure) integrand is then evaluated at
     * every node in parallel, and finally the partial terms are
     * accumulated in the original order. Bit-identical to the
     * sequential version for any thread count; only wall time changes.
     */
    real intv = 0.0;//initial value of integral
    real a = 0.0, b = 0.0;

    if(lower > upper)
    {
        a = upper;
        b = lower;
    }
    else
    {
        a = lower; 
        b = upper;
    }

    real benchmark = 1.5 * a;
    real Ng = 100.0;//integral resolution
    real hg = (benchmark - a) / (Ng);
    real lowerg = a;
    real upperg = lowerg + hg;


    real coef2 = (lowerg + upperg) / 2.0;//initializes the first coeff to change the function limits
    real coef1 = (upperg - lowerg) / 2.0;//initializes the second coeff to change the function limits
    const real c1 = 5.0 / 9.0;
    const real c2 = 8.0 / 9.0;
    const real c3 = 5.0 / 9.0;
    const real x1 = -sqrt(3.0 / 5.0);
    const real x2 __attribute__((unused)) = 0.00000000000;
    const real x3 = sqrt(3.0 / 5.0);
    real x1n = (coef1 * x1 + coef2);
    /*should be: x2n = (coef1 * x2 + coef2);*/
    real x2n = (coef2);
    real x3n = (coef1 * x3 + coef2);
    int counter = 0;
    /* ---- pass 1: enumerate nodes (same arithmetic, no evaluations) ---- */
    int cap = 1024, nIter = 0;
    real* nx = (real*) malloc((size_t) cap * 4 * sizeof(real));
    if (!nx) return 0.0;
    while (1)
    {
        if (nIter == cap)
        {
            cap *= 2;
            real* nx2 = (real*) realloc(nx, (size_t) cap * 4 * sizeof(real));
            if (!nx2) { free(nx); return 0.0; }
            nx = nx2;
        }
        nx[4*nIter + 0] = x1n;
        nx[4*nIter + 1] = x2n;
        nx[4*nIter + 2] = x3n;
        nx[4*nIter + 3] = coef1;
        nIter++;
        lowerg = upperg;
        upperg = upperg + hg;
        coef2 = (lowerg + upperg) / 2.0;//initializes the first coeff to change the function limits
        coef1 = (upperg - lowerg) / 2.0;

        x1n = ((coef1) * x1 + coef2);
        /*should be: x2n = (coef1 * x2 + coef2);*/
        x2n = (coef2);
        x3n = ((coef1) * x3 + coef2);

        if(lowerg > benchmark)
        {
            Ng = 10.0;//integral resolution
            hg = (b - benchmark) / (Ng);
        }

        if(upper > lower)
        {
            if(lowerg >= upper)//loop termination clause
            {
                break;
            }
        }
        else if(lower > upper)
        {
            if(lowerg >= lower)//loop termination clause
            {
                break;
            }
        }

        if(counter > 100000)
        {
            break;
        }
        else
        {
            counter++;
        }


    }

    /* ---- pass 2: evaluate the integrand at all nodes in parallel ----
     * (*func) is pure; results land in per-node slots, so the outcome
     * is independent of thread count and schedule. */
    real* fv = (real*) malloc((size_t) nIter * 3 * sizeof(real));
    if (!fv) { free(nx); return 0.0; }
    {
        int gi;
      #pragma omp parallel for schedule(static)
        for (gi = 0; gi < 3 * nIter; ++gi)
        {
            fv[gi] = (*func)(nx[4*(gi/3) + (gi%3)], comp1, comp2, energy, isDark);
        }
    }

    /* ---- pass 3: accumulate in the original sequential order ---- */
    for (int it = 0; it < nIter; ++it)
    {
        real cf = nx[4*it + 3];
        intv = intv + c1 * fv[3*it + 0] * cf +
                      c2 * fv[3*it + 1] * cf + 
                      c3 * fv[3*it + 2] * cf;
    }
    free(fv);
    free(nx);

    if(lower > upper)
    {
        intv *= -1.0;
    }

    return intv;
}