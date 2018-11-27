#include <iostream>
#include <cmath>
#include <math.h>
#include <iostream>
#include <vector>
#include <array>
#include "originalPDE.hpp"

double originalPDE::ContinuousAnalyticSolution( double x, double t )
{
     return 6*exp(-t)*sin(M_PI*x);
}

originalPDE::originalPDE( )
{
    a= pow(M_PI, -2), g_0 =0 , g_L = 0, k_0 = pow(10, 300), k_L = pow(10, 300);
}
