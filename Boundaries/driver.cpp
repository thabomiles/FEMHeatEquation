#include <iostream>
#include <cmath>
#include <vector>
#include <math.h>
#include <fstream>
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
//const double M_PI = 2*acos(0);

int main(int argc, char* argv[])
{

SpaceMesh smesh;
smesh.GenerateDefaultSpaceMesh();
//smesh.GloballyBisectSpaceMesh();
//smesh.GloballyBisectSpaceMesh();
//smesh.GloballyBisectSpaceMesh();
//smesh.GloballyBisectSpaceMesh();

TimeMesh tmesh;
tmesh.GenerateUniformTimeMesh(pow(smesh.meshsize(), 2), 1.0);

PDE_Q2 anotherpde;
originalPDE firstpde;
EllipticPDE_Q1 Q1;
EllipticPDE2 elliptic2;

GeneralHeat genheat;
smesh.PrintSpaceNodes();
genheat.SetSpaceTimeMesh(smesh, tmesh, firstpde);

genheat.SolveWithBCs();
genheat.PrintSolution();
//genheat.GlobalSpaceError();

//for (int i=0; i<6; i++)
//{
//    genheat.SetSpaceTimeMesh(smesh, tmesh, firstpde);
//    genheat.SolveWithBCs();
//    //std::cout<< smesh.meshsize() + 1 <<"\n";
//    genheat.GlobalSpaceError();
//    //genheat.H_1Norm();
//    std::cout<< ", " ;
//    smesh.GloballyBisectSpaceMesh();
//    tmesh.GenerateUniformTimeMesh(pow(smesh.meshsize()+1, 2), 1.0);
//}

//genheat.StationaryHeatEquation();

//genheat.SolveWithBCs();

//genheat.EnergyNorm();

//genheat.PrintSolution();
//
//genheat.BuildErrorMesh();
//genheat.PrintErrorMesh();
//genheat.GlobalSpaceError();
//
//genheat.UnitTest1();




}


