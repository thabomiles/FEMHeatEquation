#include <iostream>
#include <fstream>
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

void printTime(SpaceMesh a_smesh, TimeMesh a_tmesh);

int main(int argc, char* argv[])
{

SpaceMesh smesh;
smesh.GenerateUniformMesh(1, 5);
smesh.GloballyBisectSpaceMesh();
smesh.PrintSpaceNodes();


TimeMesh tmesh;
tmesh.GenerateUniformTimeMesh(pow(smesh.meshsize(), 2), 1.0);

PDE_Q2 anotherpde;
originalPDE firstpde;
EllipticPDE_Q1 Q1;
EllipticPDE2 elliptic2;

//printTime(smesh, tmesh);

AdaptiveSolver adapt;
adapt.SetSpaceTimeMesh(smesh, tmesh, firstpde);

adapt.AdaptiveSolve();
adapt.PrintSolution();
}


void printTime(SpaceMesh a_smesh, TimeMesh a_tmesh)
{
    ofstream myfile2;
    myfile2.open ("Y.csv");
    for(int j=1;j<a_tmesh.NumberOfTimeSteps()+1;j++ )
    {
    for(int i = 0; i<a_smesh.meshsize()+1; i++)
    {
        myfile2 << a_tmesh.ReadTimeStep(j) << ", ";
    }
    myfile2 << "\n";
    }
    myfile2.close();
}

