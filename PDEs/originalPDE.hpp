#ifndef originalPDEHEADERREF
#define originalPDEHEADERREF
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

class originalPDE: public APDE
{
    public:
        friend class GeneralHeat;
        double ContinuousAnalyticSolution( double x, double t );
        originalPDE();

    protected:

};

#endif // originalPDEHEADERREF
