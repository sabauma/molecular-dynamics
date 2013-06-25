#ifndef __SOLVER_H__
#define __SOLVER_H__

#include <algorithm>
#include <cassert>
#include <cmath>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>

#define MAX_ORDER 2

using std::max;
using std::min;

inline double MAX(const double a, const double b)
{
    return a > b ? a : b;
}

inline double MIN(const double a, const double b)
{
    return a > b ? b : a;
}

/**
 * The Solver class is templatized by the function it will take. The template
 * parameter is the type of the function object used in place of a function
 * pointer.
 */
template <class Function>
class Solver
{
public:
    // order - number of dependent variables.
    // Starting values must be in the 0th elements of arrays x and y;
    Solver(Function& f, int order, double* x, double* y);

    // void DerivativeFunc(double x, const double* y, double* dydx) is the external
    // routine that computes the derivatives (stored in dydx) from x and y.
    //
    // xMax - final value for x;
    // accuracy - computational acuracy;
    // step - initial guess for the step size (by absolute value);
    // steps - max number of points in arrays x and y.
    // If steps = 0, no intermediate results will be stored and the starting
    // values will be replaced with the final ones;
    // minStep - smallest step size allowed (by absolute value);
    // maxStep - largest step size allowed (by absolute value);
    // xMin - minimum interval at which to record points.

    // double (*DerivativeFunc)(double, const double*, double*)
    int RK4(double xMax, double accuracy,
            double step, int steps, double minStep = 0,
            double maxStep = DBL_MAX, double xMin = 0);

    // The return value for the last RK4 call
    inline int GetSteps() const
    {
        return CurrentStep;
    }

    // Returns suggested stepsize for the next iteration
    inline double GetNextStep() const
    {
        return NextStep;
    }

    // Returns the stepsize used for the last iteration
    inline double GetLastStep() const
    {
        return LastStep;
    }

    // Numerical error types
    enum Error
    {
        STEPSIZE_UNDER_HMIN,
        STEPSIZE_UNDERFLOW
    };

private:
    //derp
    int derp;

    // Number of independent variables
    int Order;

    // Statistics
    int GoodSteps;
    int BadSteps;

    int CurrentStep;
    // double (*Derivatives)(double, const double*, double*);
    Function& Derivatives;

    // These variables are for communication with rkqs
    double x;
    double Accuracy;
    double NextStep;
    double LastStep;
    double Step;

    double yscal[MAX_ORDER];

    // X is the storage for the independent variable's values;
    // space must be allocated by the caller.
    //
    // Y is the storage for the dependent variable;
    // space must be allocated by the caller.
    double* y;
    double* X;
    double* Y;
    double dydx[MAX_ORDER];

    // Output solution
    double* yOut;

    // These variables are used by SingleStep
    double yError[MAX_ORDER];
    double ak2[MAX_ORDER];
    double ak3[MAX_ORDER];
    double ak4[MAX_ORDER];
    double ak5[MAX_ORDER];
    double ak6[MAX_ORDER];
    double ytemp[MAX_ORDER];

    void SingleStep();
    void rkqs();
};

// Constructor
template <class Function>
Solver<Function>::Solver(Function& f, int order, double* x, double* y) :
    Derivatives(f)
{
    assert(order <= MAX_ORDER);
    Order = order;

    X = x;
    Y = y;
    NextStep = 0;
    derp = 0;
}

// Single Runge-Kutta step
template <class Function>
void Solver<Function>::SingleStep()
{
    // Debugging?
    if  (derp == 3000)
    {
        printf("derp = 3000\n");
    }

    // Runge-Kutta constants
    const double a2 = 0.2,
                 a3 = 0.3,
                 a4 = 0.6,
                 a5 = 1.0,
                 a6 = 0.875,

                 b21 = 0.2,
                 b31 = 3.0/40.0,
                 b32 = 9.0/40.0,
                 b41 = 0.3,
                 b42 = -0.9,
                 b43 = 1.2,
                 b51 = -11.0/54.0,
                 b52 = 2.5,
                 b53 = -70.0/27.0,
                 b54 = 35.0/27.0,
                 b61 = 1631.0/55296.0,
                 b62 = 175.0/512.0,
                 b63 = 575.0/13824.0,
                 b64 = 44275.0/110592.0,
                 b65 = 253.0/4096.0,

                 c1 = 37.0/378.0,
                 c3 = 250.0/621.0,
                 c4 = 125.0/594.0,
                 c6 = 512.0/1771.0,

                 dc1 = c1 - 2825.0/27648.0,
                 dc3 = c3 - 18575.0/48384.0,
                 dc4 = c4  -13525.0/55296.0,
                 dc5 = -277.00/14336.0,
                 dc6=c6-0.25;

    double ytemp[MAX_ORDER];

    for ( int i = 0; i < Order; i++ )
        { ytemp[i] = y[i] + b21*LastStep*dydx[i]; }

    Derivatives.Derivatives(x + a2*LastStep, ytemp, ak2);

    for ( int i = 0; i < Order; i++ )
        { ytemp[i] = y[i] + LastStep*(b31*dydx[i] + b32*ak2[i]); }

    Derivatives.Derivatives(x + a3*LastStep, ytemp, ak3);

    for ( int i = 0; i < Order; i++ )
        { ytemp[i] = y[i] + LastStep*(b41*dydx[i] + b42*ak2[i] + b43*ak3[i]); }

    Derivatives.Derivatives(x + a4*LastStep, ytemp, ak4);

    for ( int i = 0; i < Order; i++ )
        { ytemp[i] = y[i] + LastStep*(b51*dydx[i] + b52*ak2[i] + b53*ak3[i] + b54*ak4[i]); }

    Derivatives.Derivatives(x + a5*LastStep, ytemp, ak5);

    for ( int i = 0; i < Order; i++ )
        { ytemp[i] = y[i] + LastStep*(b61*dydx[i] + b62*ak2[i] + b63*ak3[i] + b64*ak4[i] + b65*ak5[i]); }

    Derivatives.Derivatives(x + a6*LastStep, ytemp, ak6);

    // Output solution
    for ( int i = 0; i < Order; i++ )
        { yOut[i] = y[i] + LastStep*(c1*dydx[i] + c3*ak3[i] + c4*ak4[i] + c6*ak6[i]); }

    // Error
    for ( int i = 0; i < Order; i++ )
        { yError[i] = LastStep*(dc1*dydx[i] + dc3*ak3[i] + dc4*ak4[i] + dc5*ak5[i] + dc6*ak6[i]); }
}

template <class Function>
void Solver<Function>::rkqs()
{
    derp++;
    const double SAFETY = 0.9, PGROW = -0.2, PSHRNK = -0.25, ERRCON = 1.89E-4;
    double maxError, tempStep, xnew;

    LastStep = Step;
    for ( ;; )
    {
        SingleStep();
        maxError = 0.0;

        for ( int i = 0; i < Order; i++ )
        {
            maxError = max(maxError, fabs(yError[i] / yscal[i]));
        }

        maxError /= Accuracy;

        if ( maxError <= 1.0 ) { break; }

        tempStep = SAFETY*LastStep*pow(maxError,PSHRNK);
        LastStep = (LastStep >= 0.0 ? max(tempStep, 0.1*LastStep)
                    : min(tempStep, 0.1*LastStep));

        xnew = x + LastStep;

        if (xnew == x)
        {
            throw Error(STEPSIZE_UNDER_HMIN);
        }
    }

    if ( maxError > ERRCON )
    { NextStep = SAFETY*LastStep*pow(maxError, PGROW); }
    else
    { NextStep = 5.0*LastStep; }

    x += LastStep;
}

template <class Function>
int Solver<Function>::RK4(double xMax, double accuracy, double step, int steps,
                          double minStep, double maxStep, double xMin)
{
    const double TINY = 1E-30;

    Accuracy = accuracy;

    double xsav = x = X[0];

    y = Y;
    double xStart = x;

    // Step up or step down?
    Step = (xMax - x) > 0 ? step : -step;

    // Reset step counters
    GoodSteps = BadSteps = 0;

    // Non-zero step?
    if ( steps != 0 )
    {
        CurrentStep = 0;
        // This is where we will write calculated data
        yOut = y + Order;
    }
    else
    {
        // Do nothing if the step size is zero
        CurrentStep = -1;
        yOut = y;
    }

    for ( ;; )
    {
        // Finished integrating because we are outside the interval? Should
        // NEVER happen when we calculate RP because our interval is VERY
        // large
        bool isFinished = xMax > xStart ? (x >= xMax) : (x <= xMax);

        // Debugging?
        if ( derp == 2228 )
        {
            printf("derp = 2228\n");
        }

        if  ( fabs(x - xsav) > xMin || isFinished )
        {
            yOut += Order;
            X[++CurrentStep] = xsav = x;
            if ( CurrentStep == steps || isFinished )
            { return CurrentStep; }
        }

        Derivatives.Derivatives(x, y, dydx);

        for ( int i = 0; i < Order; i++ )
        {
            yscal[i] = fabs(y[i]) + fabs(Step*dydx[i]) + TINY;
        }

        double hleft = xMax - x;
        if ( fabs(Step) > fabs(hleft) ) { Step = hleft; }

        rkqs();
        y = yOut;

        // Keep track of good vs. bad steps
        if ( LastStep == Step )
        { GoodSteps++; }
        else
        { BadSteps++; }

        // Converging too slowly - step too small
        if ( fabs(NextStep) <= minStep )
        { throw Error(STEPSIZE_UNDER_HMIN); }

        // Determine next step
        if ( NextStep > 0 )
        { Step = MIN(NextStep, maxStep); }
        else
        { Step = MAX(NextStep, -maxStep); }
    }
}

#endif /* __SOLVER_H__ */
