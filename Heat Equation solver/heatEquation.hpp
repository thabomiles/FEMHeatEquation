#ifndef HEATEQUATIONHEADERREF
#define HEATEQUATIONHEADERREF
#include <iostream>
#include <vector>
#include <array>
#include "SpaceMesh.hpp"
#include "TimeMesh.hpp"

class HeatEquation
{
public:

    void SolveUniformMesh();

    void SetSpaceTimeMesh( SpaceMesh smesh, TimeMesh tmesh, const std::string outputFileName );

    void SetUniformSystem( double endTime, int numberOfTimeSteps, int numberOfSpaceElements,
                            const std::string outputFileName );

    void Solve();



private:

    int mn, mm;

    double mT;

    std::string moutputFileName;

    SpaceMesh mpSpaceMesh;
    TimeMesh mpTimeMesh;
};


#endif // HEATEQUATIONHEADERREF
