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
        double ContinuousAnalyticSolution( double x, double t );

    protected:
        const double M_PI = 2*acos(0);

        double g_0 =0 , g_L = 0, k_0 = 0, k_L = 0;




};

#endif // APDEHEADERREF
