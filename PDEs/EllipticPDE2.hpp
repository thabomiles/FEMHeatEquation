#ifndef ELLIPTICPDE2HEADERREF
#define ELLIPTICPDE2HEADERREF
#include <cmath>
#include <math.h>
#include <vector>
#include <array>
#include "APDE.hpp"

class EllipticPDE2: public APDE
{
    public:
        friend class GeneralHeat;
        double EllipticalRHSfunction( double x );
        double ContinuousAnalyticSolution( double x, double t );

        EllipticPDE2();


    protected:

};

#endif // ELLIPTICPDE2HEADERREF


