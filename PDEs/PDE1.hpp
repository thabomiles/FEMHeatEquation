#ifndef PDE1HEADERREF
#define PDE1HEADERREF
#include <cmath>
#include <math.h>
#include <vector>
#include <array>
#include "SpaceMesh.hpp"
#include "TimeMesh.hpp"
#include "TriDiagMatrix.hpp"
#include "MassMatrix.hpp"
#include "StiffnessMatrix.hpp"



class PDE1: public APDE
{
    public:
        double ContinuousAnalyticSolution( double x, double t );

    protected:

        double g_0 =0 , g_L = 1, k_0 = pow(10, 300), k_L = pow(10, 300);




};

#endif // PDE1HEADERREF

