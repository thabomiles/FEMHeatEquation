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

void printTime(SpaceMesh a_smesh, TimeMesh a_tmesh);

int main(int argc, char* argv[])
{

SpaceMesh smesh;
smesh.GenerateDefaultSpaceMesh();
//smesh.GloballyBisectSpaceMesh();
//smesh.GloballyBisectSpaceMesh();
//smesh.GloballyBisectSpaceMesh();
//smesh.GloballyBisectSpaceMesh();

TimeMesh tmesh;
tmesh.GenerateUniformTimeMesh(pow(smesh.meshsize(), 1), 1.0);

PDE_Q2 anotherpde;
originalPDE firstpde;
EllipticPDE_Q1 Q1;
EllipticPDE2 elliptic2;

GeneralHeat genheat;
smesh.PrintSpaceNodes();
genheat.SetSpaceTimeMesh(smesh, tmesh, anotherpde);
//
genheat.SolveWithBCs();
printTime(smesh, tmesh);
genheat.PrintSolution();


//genheat.GlobalSpaceError();

//for (int i=0; i<10; i++)
//{
//    genheat.SetSpaceTimeMesh(smesh, tmesh, anotherpde);
//    genheat.SolveWithBCs();
//    std::cout<< smesh.meshsize() + 1 <<"\n";
//    //genheat.UnitTest1();
//    //genheat.GlobalSpaceError();
//    //genheat.H_1Norm();
//    //genheat.EnergyNorm();
//    std::cout<< ", " ;
//    smesh.GloballyBisectSpaceMesh();
//    tmesh.GenerateUniformTimeMesh(pow(smesh.meshsize()+1, 1), 1.0);
//}


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
