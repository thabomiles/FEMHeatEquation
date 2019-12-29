#ifndef NONANALYTICHEADERREF
#define NONANALYTICHEADERREF
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



class Nonanalytic: public APDE
{
    public:
        friend class GeneralHeat;
        void InitialCondition ( SpaceMesh& a_mesh, std::vector<double>& first_U );
        double FirstBoundary( double t );
        double SecondBoundary( double t );
        Nonanalytic( );

    protected:


};

#endif // NONANALYTICHEADERREF

