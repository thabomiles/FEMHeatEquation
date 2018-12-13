#ifndef PDE_Q2HEADERREF
#define PDE_Q2HEADERREF
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

class PDE_Q2: public APDE
{
    public:
        friend class GeneralHeat;
        double ContinuousAnalyticSolution( double x, double t );
        double AnalyticGradientWRTx( double x, double t );
        void InitialCondition ( SpaceMesh& a_mesh, std::vector<double>& first_U );
        double FirstBoundary( double t );
        double SecondBoundary( double t );
        PDE_Q2( );

    protected:

};

#endif // PDE_Q2HEADERREF
