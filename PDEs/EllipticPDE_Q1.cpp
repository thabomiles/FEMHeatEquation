#include <iostream>
#include <cmath>
#include <math.h>
#include <iostream>
#include <vector>
#include <array>
#include "EllipticPDE_Q1.hpp"

double EllipticPDE_Q1::EllipticalRHSfunction( double x )
{
    return -sin(M_PI*x)*pow(M_PI, 2);
}

EllipticPDE_Q1::EllipticPDE_Q1()
{
a= 1, g_0 =0 , g_L = 1, k_0 = pow(10, 300), k_L = pow(10, 300);
}

