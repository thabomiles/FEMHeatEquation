#include <iostream>
#include <cmath>
#include <math.h>
#include <iostream>
#include <vector>
#include <array>
#include "EllipticPDE2.hpp"

double EllipticPDE2::EllipticalRHSfunction( double x )
{
    return cos(M_PI*x-M_PI*0.5)*pow(M_PI, 2);
}

double EllipticPDE2::ContinuousAnalyticSolution( double x, double t )
{
    return cos(M_PI*x-M_PI*0.5);
}


EllipticPDE2::EllipticPDE2()
{
a= 1, g_0 =0 , g_L = 0, k_0 = pow(10, 300), k_L = pow(10, 300);
}


