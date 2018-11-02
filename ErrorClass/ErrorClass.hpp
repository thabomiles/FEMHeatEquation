#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <math.h>
#include "SpaceMesh.hpp"
#include "TimeMesh.hpp"
#include <fstream>
#include <string>
#include <boost/math/quadrature/gauss.hpp>
using namespace std;
using namespace boost::math::quadrature;

#ifndef ERRORSHEADERREF
#define ERRORSHEADERREF

class Errors
{
public:

    void GenerateSpaceMesh ( std::vector<double> SpaceNodes );
    void SetSystem (const  )
protected:

    std::vector<double> mpCurrentUh;
    std::vector<double> mpCurrentAnalyticSolution;
    std::vector<double> mpPreviousUh;

};


#endif // ERRORSHEADERREF



