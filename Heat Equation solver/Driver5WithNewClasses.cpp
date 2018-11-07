#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <math.h>
#include "heatEquation.hpp"
#include "TriDiagMatrix.hpp"
#include "MassMatrix.hpp"
#include "StiffnessMatrix.hpp"
#include <fstream>
#include <string>
#include <boost/math/quadrature/gauss.hpp>
using namespace std;
using namespace boost::math::quadrature;

void printQuadrature();
double f( double x);
double g (double x);
double h (double x);
const double M_PI = 2*acos(0);


int main(int argc, char* argv[])
{

SpaceMesh smesh;
//smesh.GenerateDefaultSpaceMesh();
//smesh.GloballyBisectSpaceMesh();
smesh.GenerateSpaceMesh({0, 0.25, 0.375, 0.5, 0.75, 1});

TimeMesh tmesh;
tmesh.GenerateUniformTimeMesh(pow(smesh.meshsize(), 2), 1.0);

HeatEquation heat;
heat.SetSpaceTimeMesh( smesh, tmesh, "soultion1.txt");
heat.Solve();

//smesh.PrintSpaceNodes();
heat.PrintSolution();

heat.PrintErrorMesh();

std::cout << heat.GlobalSpaceError();
std::cout << " \n";

}

double f( double x)
{
    double t = 6*sin(M_PI*x)*exp(-5*pow(9, -1));
    double m1 = 4*2.39626*x;
//    return t;
    return pow(t-m1,2);
}

double g (double x)
{
    double t = 6*sin(M_PI*x)*exp(-5*pow(9, -1));
    double m1 = 4*(3.31376-2.39626)*(x-0.25)+2.39626;
    return pow(t-m1,2);
}

double h (double x)
{
    double t = 6*sin(M_PI*x)*exp(-5*pow(9, -1));
    double m1 = 2*(-3.31376)*(x-1);
    return pow(t-m1,2);
}

void printQuadrature()
{
std::cout<<gauss<double, 7>::integrate(f, 0, 0.25) << ", "
        << gauss<double, 7>::integrate(g, 0.25, 0.5) << ", "
        << gauss<double, 7>::integrate(h, 0.5, 1) << ", ";

std::cout << " \n";

}

