#include <iostream>
#include <cmath>
#include <math.h>
#include <iostream>
#include <vector>
#include <array>
#include "APDE.hpp"
#include "PDE1.hpp"

double PDE1::ContinuousAnalyticSolution( double x, double t )
{
     return exp(-4*t)*sin(2*M_PI*x)+x;
}

double PDE1::AnalyticGradientWRTx( double x, double t )
{
    return 2*M_PI*exp(-4*t)*cos(2*M_PI*x)+1;
}
