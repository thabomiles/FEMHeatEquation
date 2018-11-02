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

double f( double x);
double g (double x);
const double M_PI = 2*acos(0);


int main(int argc, char* argv[])
{

SpaceMesh smesh;
smesh.GenerateDefaultSpaceMesh();

tmesh.GenerateUniformTimeMesh(pow(smesh.meshsize(), 2), 1.0);

HeatEquation heat;
heat.SetSpaceTimeMesh( smesh, tmesh, "soultion1.txt");
heat.Solve();

smesh.PrintSpaceNodes();
heat.PrintSolution();

heat.BuildErrorMesh();
heat.PrintErrorMesh();
std::cout << heat.GlobalSpaceError();
std::cout << " \n";


std::cout << smesh.meshsize()+1;
std::cout << " \n";
//std::cout << sqrt(currentError);
//std::cout << " \n";
//std::cout << " \n";
//
//double M = gauss<double, 7>::integrate(f, 0, 0.25);
//std::cout << M;
//
//std::cout << " \n";
//std::cout << " \n";
//
//double Q = gauss<double, 7>::integrate(g, 0.25, 0.5);
//std::cout << Q;
//
//std::cout << " \n";
//std::cout << " \n";
//
//std::cout << heat.L2ErrorGuass( smesh.ReadSpaceNode(1),smesh.ReadSpaceNode(2))<< ", ";
//
//std::cout << " \n";
//
////std::cout << f (0.2);
//
//std::cout << " \n";
//
////std::cout << heat.ErrorSquared(0.2);
//

}

double f( double x)
{
    double t = 6*sin(M_PI*x)*exp(-1);
    double m1 = 4*1.66329*x;
//    return t;
    return pow(t-m1,2);
}

double g (double x)
{
    double t = 6*sin(M_PI*x)*exp(-1);
    double m1 = 4*(2.29857-1.66329)*(x-0.25)+1.66329;
    return pow(t-m1,2);
}


