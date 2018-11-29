#include <iostream>
#include <cmath>
#include <vector>
#include <math.h>
#include "SpaceAdaptiveSolver.hpp"
#include "GeneralizedHeat.hpp"
#include "APDE.hpp"
#include "PDE_Q2.hpp"
#include "originalPDE.hpp"
#include "EllipticPDE_Q1.hpp"
#include "EllipticPDE2.hpp"
#include <string>
#include <boost/math/quadrature/gauss.hpp>
using namespace std;
using namespace boost::math::quadrature;

int main(int argc, char* argv[])
{

SpaceMesh smesh;
smesh.GenerateUniformMesh(1, 5);
//smesh.GloballyBisectSpaceMesh();
//smesh.GloballyBisectSpaceMesh();

TimeMesh tmesh;
tmesh.GenerateUniformTimeMesh(pow(smesh.meshsize(), 2), 1.0);

PDE_Q2 anotherpde;
originalPDE firstpde;
EllipticPDE_Q1 Q1;
EllipticPDE2 elliptic2;

AdaptiveSolver adapt;
adapt.UnitTest();

//GeneralHeat genheat;
//genheat.SetSpaceTimeMesh(smesh, tmesh, elliptic2);
//
//genheat.StationaryHeatEquation();

//genheat.SolveWithBCs();
//
//
//genheat.PrintSolution();
//
//genheat.BuildErrorMesh();
//genheat.PrintErrorMesh();
//genheat.GlobalSpaceError();
//
//genheat.UnitTest1();



}
