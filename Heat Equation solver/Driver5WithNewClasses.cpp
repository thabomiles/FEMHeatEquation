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
smesh.GenerateSpaceMesh({0, 0.15, 0.25, 0.5, 1});

//smesh.GenerateUniformMesh(1, 5);
//smesh.GloballyBisectSpaceMesh();
//smesh.GloballyBisectSpaceMesh();

TimeMesh tmesh;
tmesh.GenerateUniformTimeMesh(pow(smesh.meshsize(), 2), 1.0);

HeatEquation heat;
heat.SetSpaceTimeMesh( smesh, tmesh, "soultion1.txt");
heat.Solve();

smesh.PrintSpaceNodes();
heat.PrintSolution();

//heat.PrintErrorMesh();
//
//std::cout << heat.GlobalSpaceError();
//std::cout << " \n";

}


