#include <iostream>
#include <cmath>
#include <math.h>
#include <iostream>
#include <vector>
#include <array>
#include "APDE.hpp"

//double APDE::ContinuousAnalyticSolution( double x, double t )
//{
//     return 6*exp(-t)*sin(M_PI*x);
//}

double APDE::ContinuousAnalyticSolution( double x, double t )
{
     return exp(-4*t)*sin(2*M_PI*x)+x;
}

