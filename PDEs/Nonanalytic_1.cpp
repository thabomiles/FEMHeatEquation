#include <iostream>
#include <cmath>
#include <math.h>
#include <iostream>
#include <vector>
#include <array>
#include "Nonanalytic_1.hpp"

//double Nonanalytic::ContinuousAnalyticSolution( double x, double t )
//{
//     return exp(-4*t)*sin(2*M_PI*x)+x;
//}
//
//double Nonanalytic::AnalyticGradientWRTx( double x, double t )
//{
//    return 2*M_PI*exp(-4*t)*cos(2*M_PI*x)+1;
//}

double Nonanalytic::FirstBoundary( double t )
{
    return t;
    //return 0;
}

double Nonanalytic::SecondBoundary( double t )
{
    return exp(-t);
}

void Nonanalytic::InitialCondition ( SpaceMesh& a_mesh, std::vector<double>& first_U )
{
    first_U.clear();
    double var;
    for (int i = 0; i<a_mesh.meshsize()+1; i++)
    {
        var = sin(2*M_PI*a_mesh.ReadSpaceNode(i))+a_mesh.ReadSpaceNode(i);
        first_U.push_back(var);
    }
}

Nonanalytic::Nonanalytic( )
{
    a= 0.1, g_0 =0 , g_L = 1, k_0 = pow(10, 300), k_L = pow(10, 300);
}

