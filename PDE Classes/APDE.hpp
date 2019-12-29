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
        virtual void InitialCondition ( SpaceMesh& a_mesh, std::vector<double>& first_U );




    protected:
        const double M_PI = 2*acos(0);

        double a, g_0, g_L, k_0, k_L;

};

#endif // APDEHEADERREF
