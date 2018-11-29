#ifndef ELLIPTICPDE_Q1HEADERREF
#define ELLIPTICPDE_Q1HEADERREF
#include <cmath>
#include <math.h>
#include <vector>
#include <array>
#include "SpaceMesh.hpp"
#include "TimeMesh.hpp"
#include "TriDiagMatrix.hpp"
#include "MassMatrix.hpp"
#include "StiffnessMatrix.hpp"
#include "APDE.hpp"

class EllipticPDE_Q1: public APDE
{
    public:
        friend class GeneralHeat;
        double EllipticalRHSfunction( double x );
        double ContinuousAnalyticSolution( double x, double t );
        EllipticPDE_Q1();

    protected:

};

#endif // ELLIPTICPDE_Q1HEADERREF

