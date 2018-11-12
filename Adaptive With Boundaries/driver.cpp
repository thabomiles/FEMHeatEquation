#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <math.h>
#include "AdaptiveHeat.hpp"
#include "StiffnessMatrix.hpp"
#include "MassMatrix.hpp"
#include "GeneralizedHeat.hpp"
#include <fstream>
#include <string>
#include <boost/math/quadrature/gauss.hpp>
using namespace std;
using namespace boost::math::quadrature;
#include "APDE.hpp"

const double M_PI = 2*acos(0);


int main(int argc, char* argv[])
{

SpaceMesh smesh;
smesh.GenerateSpaceMesh({0, 0.15, 0.25, 0.5, 1.0});
//smesh.GloballyBisectSpaceMesh();
//smesh.GloballyBisectSpaceMesh();
//smesh.GloballyBisectSpaceMesh();


TimeMesh tmesh;
tmesh.GenerateUniformTimeMesh(pow(smesh.meshsize(), 2), 1.0);

APDE apde;

GeneralHeat genheat;
genheat.SetSpaceTimeMesh(smesh, tmesh, apde);

genheat.SolveWithBCs();
//genheat.Solve();

//
//genheat.PrintSolution();
//genheat.PrintErrorMesh();
//genheat.GlobalSpaceError();

}


