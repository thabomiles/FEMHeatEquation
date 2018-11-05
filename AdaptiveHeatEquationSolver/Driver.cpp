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
smesh.GenerateDefaultSpaceMesh();
smesh.GloballyBisectSpaceMesh();
smesh.GloballyBisectSpaceMesh();
smesh.GloballyBisectSpaceMesh();

smesh.PrintSpaceNodes();


TimeMesh tmesh;
tmesh.GenerateUniformTimeMesh(pow(smesh.meshsize(), 2), 1.0);

TimeMesh tmesh1;




AdaptiveHeatEquation adaptiveheat;
adaptiveheat.SetSpaceTimeMesh( smesh, tmesh, "soultion1.txt");
adaptiveheat.Solve();
//
adaptiveheat.PrintSolution();



adaptiveheat.BuildErrorMesh();
adaptiveheat.PrintErrorMesh();

std::cout << adaptiveheat.GlobalSpaceError();
std::cout << " \n";



}

//std::cout << " \n";
//
//std::vector<double> test = { 2, 2, 2 } ;
//std::cout <<test.size()<< ", "<< smesh.meshsize()+1;
//std::cout << " \n";
//
//std::cout <<adaptiveheat.InterpolantFunction(0.99, test)<< ", ";
//for (int i = 0; i<21; i++)
//{
// //   adaptiveheat.InterpolantFunction(i*0.02, test);
//
//    std::cout <<adaptiveheat.InterpolantFunction(i*0.05, test)<< ", ";
//}


