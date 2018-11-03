#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <math.h>
#include "AdaptiveHeat.hpp"
#include <fstream>
#include <string>
#include <boost/math/quadrature/gauss.hpp>
using namespace std;
using namespace boost::math::quadrature;

const double M_PI = 2*acos(0);


int main(int argc, char* argv[])
{

SpaceMesh smesh;
smesh.GenerateDefaultSpaceMesh();

TimeMesh tmesh;
tmesh.GenerateUniformTimeMesh(pow(smesh.meshsize(), 2), 1.0);


AdaptiveHeatEquation adaptiveheat;
adaptiveheat.SetSpaceTimeMesh( smesh, tmesh, "soultion1.txt");
adaptiveheat.Solve();

//smesh.PrintSpaceNodes();
adaptiveheat.PrintSolution();

adaptiveheat.BuildErrorMesh();
adaptiveheat.PrintErrorMesh();

std::cout << adaptiveheat.GlobalSpaceError();
std::cout << " \n";

}
