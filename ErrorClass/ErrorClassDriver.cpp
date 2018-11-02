#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <math.h>
#include "heatEquation.hpp"
#include "TriDiagMatrix.hpp"
#include "MassMatrix.hpp"
#include "StiffnessMatrix.hpp"
#include <fstream>
#include <string>
#include <boost/math/quadrature/gauss.hpp>
using namespace std;
using namespace boost::math::quadrature;

int main(int argc, char* argv[])
{

double T = 1.0;

    //number of time steps i.e. 1 less than the number of nodes
int m = 32;


double M_PI = 2*acos(0);

double a = pow(M_PI,-2);

std::vector<double> initialSpaceNodes = {0, 0.25, 0.5, 1};
std::vector<double> initialTimeNodes = {0};
std::vector<double> initialSpaceMesh;
std::vector<double> initialTimeMesh;

       //n is the number of elements i.e. 1 less than the number of nodes
int n = initialSpaceNodes.size()-1;


for(int i=1; i<=m; i++)
    {
        initialTimeNodes.push_back( i*(T/m) );
        initialTimeMesh.push_back( initialTimeNodes.at(i)- initialTimeNodes.at(i-1) );
    }



SpaceMesh smesh;
smesh.GenerateSpaceMesh(initialSpaceNodes);
smesh.GloballyBisectSpaceMesh();
smesh.GloballyBisectSpaceMesh();
smesh.GloballyBisectSpaceMesh();
smesh.GloballyBisectSpaceMesh();







}




