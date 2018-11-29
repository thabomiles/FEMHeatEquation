#include <iostream>
#include <cmath>
#include <array>
#include <math.h>
#include "AdaptiveHeat.hpp"
#include "StiffnessMatrix.hpp"
#include "MassMatrix.hpp"
#include <string>
#include <boost/math/quadrature/gauss.hpp>
using namespace std;
using namespace boost::math::quadrature;

const double M_PI = 2*acos(0);


int main(int argc, char* argv[])
{

AdaptiveHeatEquation adaptiveheat;
adaptiveheat.UnitTest();

//SpaceMesh smesh;
//smesh.GenerateUniformMesh(1, 11);
//
//std::array<double, 6> a ={1, 0.98, 0.96, 0.94, 0.92, 0.9};
//
//for (auto i: a)
//    std::cout<< smesh.GeneralTestFunctions(10, i)<< ", " ;


}




