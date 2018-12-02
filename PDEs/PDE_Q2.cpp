#include <iostream>
#include <cmath>
#include <math.h>
#include <iostream>
#include <vector>
#include <array>
#include "PDE_Q2.hpp"

double PDE_Q2::ContinuousAnalyticSolution( double x, double t )
{
     return exp(-4*t)*sin(2*M_PI*x)+x;
}

double PDE_Q2::AnalyticGradientWRTx( double x, double t )
{
    return 2*M_PI*exp(-4*t)*cos(2*M_PI*x)+1;
}

PDE_Q2::PDE_Q2(  )
{
    a= pow(M_PI, -2), g_0 =0 , g_L = 1, k_0 = pow(10, 300), k_L = pow(10, 300);
}
