#include <iostream>
#include <cmath>
#include "heatEquation.hpp"
#include <vector>
#include <array>
#include <math.h>

int main(int argc, char* argv[])
{

double T = 4.0;

int n = 3;

int m = 4;

double h = 1/n;

double M_PI = 2*acos(0);

double a = pow(M_PI,-2);

std::vector<double> initialSpaceNodes{0};
std::vector<double> initialTimeNodes{0};


std::vector<double> Diagonal;
std::vector<double> LowerDiag;
std::vector<double> UpperDiag;
std::vector<double> f;
std::vector<double> x;

for(int i=1; i<=n; i++)
    {
        initialSpaceNodes.push_back( i*(1/n) );
    }

 for(int i=1; i<=m; i++)
     {
         initialTimeNodes.push_back( i*(T/m) );
     }
}

