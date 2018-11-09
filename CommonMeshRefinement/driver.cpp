#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <math.h>
#include "AdaptiveHeat.hpp"
#include "StiffnessMatrix.hpp"
#include "MassMatrix.hpp"
#include <fstream>
#include <string>
#include <boost/math/quadrature/gauss.hpp>
using namespace std;
using namespace boost::math::quadrature;

const double M_PI = 2*acos(0);


int main(int argc, char* argv[])
{

SpaceMesh currentsmesh;
currentsmesh.GenerateUniformMesh(1,6);
//currentsmesh.GenerateSpaceMesh({0,0.125,0.25, 0.5, 1});

SpaceMesh oldsmesh;
oldsmesh.GenerateUniformMesh(1.0, 6);

AdaptiveHeatEquation adaptiveheat;

std::vector<double> u = {2,3,4,2};

std::vector<double> b;
std::vector<double> intervals;
SpaceMesh refinedsmesh;
refinedsmesh.CommonMesh(currentsmesh, oldsmesh);
//refinedsmesh.PrintSpaceNodes();

double integral;
for (int i = 1; i<currentsmesh.meshsize(); i++)
{
    integral=0;
    refinedsmesh.Range(currentsmesh.ReadSpaceNode(i-1), currentsmesh.ReadSpaceNode(i+1), intervals);
    adaptiveheat.PrintVector(intervals);
    for(int j = 0; j<intervals.size()-1; j++)
    {
        integral = integral + adaptiveheat.IntegrateBasisWithU(i, intervals.at(j),
                    intervals.at(j+1), currentsmesh, oldsmesh, u);
                    std::cout<<integral<< "\n";
    }
    b.push_back(integral);
}

adaptiveheat.PrintVector(b);

//std::cout<<adaptiveheat.SolutionTimesBasis(1, 0.25, currentsmesh, oldsmesh, u )<< "\n";


MassMatrix mass;
mass.BuildMassMatrix(currentsmesh);
mass.MatrixVectorMultiplier(u, b);
adaptiveheat.PrintVector(b);

//
//adaptiveheat.PrintVector(b);

}


