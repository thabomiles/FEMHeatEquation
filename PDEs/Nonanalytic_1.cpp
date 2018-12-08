#include <iostream>
#include <cmath>
#include <math.h>
#include <iostream>
#include <vector>
#include <array>
#include "Nonanalytic_1.hpp"


double Nonanalytic::FirstBoundary( double t )
{
    return t;
}

double Nonanalytic::SecondBoundary( double t )
{
    return t;
}

double Nonanalytic::InitialCondition ( double x )
{

}

originalPDE::originalPDE( )
{
    a= 1, g_0 =0 , g_L = 1, k_0 = pow(10, 300), k_L = pow(10, 300);
}

