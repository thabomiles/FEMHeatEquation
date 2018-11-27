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
#include "PDE_Q2.hpp"
#include "originalPDE.hpp"

const double M_PI = 2*acos(0);


int main(int argc, char* argv[])
{

SpaceMesh smesh;
smesh.GenerateUniformMesh(1, 5);
smesh.GloballyBisectSpaceMesh();
smesh.GloballyBisectSpaceMesh();

TimeMesh tmesh;
tmesh.GenerateUniformTimeMesh(pow(smesh.meshsize(), 2), 1.0);

PDE_Q2 anotherpde;
originalPDE firstpde;

GeneralHeat genheat;
genheat.SetSpaceTimeMesh(smesh, tmesh, firstpde);

genheat.SolveWithBCs();
//
//
//genheat.PrintSolution();
//
//genheat.BuildErrorMesh();
//genheat.PrintErrorMesh();
genheat.GlobalSpaceError();
//
//genheat.UnitTest1();




}


