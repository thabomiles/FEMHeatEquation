#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <math.h>
#include "AdaptiveHeat.hpp"
#include "StiffnessMatrix.hpp"
#include <fstream>
#include <string>
#include <boost/math/quadrature/gauss.hpp>
using namespace std;
using namespace boost::math::quadrature;

const double M_PI = 2*acos(0);


int main(int argc, char* argv[])
{

SpaceMesh smesh;
smesh.GenerateSpaceMesh({0, 0.15,  0.25, 0.5, 1.0});
//smesh.GloballyBisectSpaceMesh();

smesh.PrintSpaceNodes();
std::cout<<smesh.meshsize() <<"\n";
std::cout<<"\n";


TimeMesh tmesh;
tmesh.GenerateUniformTimeMesh(pow(smesh.meshsize(), 2), 1.0);

AdaptiveHeatEquation adaptiveheat;
adaptiveheat.SetSpaceTimeMesh( smesh, tmesh, "soultion1.txt");
adaptiveheat.AdaptiveSolver();


//adaptiveheat.Solve();


//adaptiveheat.PrintSolution();

//adaptiveheat.PrintErrorMesh();
adaptiveheat.GlobalSpaceError();

}


