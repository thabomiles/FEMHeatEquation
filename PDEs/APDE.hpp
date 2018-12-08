#ifndef APDEHEADERREF
#define APDEHEADERREF
#include <cmath>
#include <math.h>
#include <vector>
#include <array>
#include "SpaceMesh.hpp"
#include "TimeMesh.hpp"
#include "TriDiagMatrix.hpp"
#include "MassMatrix.hpp"
#include "StiffnessMatrix.hpp"

class APDE
{
    public:
        friend class GeneralHeat;
        virtual double ContinuousAnalyticSolution( double x, double t );
        virtual double AnalyticGradientWRTx( double x, double t );
        virtual double EllipticalRHSfunction( double x );
        virtual double FirstBoundary( double t );
        virtual double SecondBoundary( double t );




    protected:
        const double M_PI = 2*acos(0);

        //double a= pow(M_PI, -2), g_0 =0 , g_L = 0, k_0 = pow(10, 300), k_L = pow(10, 300);
        double a, g_0, g_L, k_0, k_L;

};

#endif // APDEHEADERREF
