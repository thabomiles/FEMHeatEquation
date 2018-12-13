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

void PDE_Q2::InitialCondition ( SpaceMesh& a_mesh, std::vector<double>& first_U )
{
    first_U.clear();
    double var;
    for (int i = 0; i<a_mesh.meshsize()+1; i++)
    {
        var = sin(2*M_PI*a_mesh.ReadSpaceNode(i))+a_mesh.ReadSpaceNode(i);
        first_U.push_back(var);
//        var = ContinuousAnalyticSolution(a_mesh.ReadSpaceNode(i), 0);
//        first_U.push_back(var);
    }
}

double PDE_Q2::FirstBoundary( double t )
{
    return 0;
}

double PDE_Q2::SecondBoundary( double t )
{
    return 1;
}



PDE_Q2::PDE_Q2(  )
{
    a= pow(M_PI, -2), g_0 =0 , g_L = 1, k_0 = pow(10, 300), k_L = pow(10, 300);
}
