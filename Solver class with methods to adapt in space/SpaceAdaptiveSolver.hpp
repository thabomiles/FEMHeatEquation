#ifndef ADAPTIVESOLVERHEADERREF
#define ADAPTIVESOLVERHEADERREF
#include <iostream>
#include <vector>
#include <array>
#include <functional>
#include "GeneralizedHeat.hpp"
#include "AdaptiveHeat.hpp"
#include "SpaceMesh.hpp"
#include "TimeMesh.hpp"
#include "APDE.hpp"

class AdaptiveSolver: public GeneralHeat
{
public:
    void AdaptiveSolve();
    void BuildRHS();

    double IntegrateBasisWithU( int NodeIndex, double lowerlimit, double upperlimit );
    void RefineMesh();

    void AdaptiveChangingBC();
    void SaveRefinementNodes( std::vector<double>& Nodes_to_insert );
    void SaveIntervalsForRefinement();
    void SaveIntervalsForCoarsening();
    void SaveIntervalsForCoarsening2();
    void MeshRefinement();
    void SetTolerances(const double refineby, const double coarsenby);
    double SolutionTimesBasis( int NodeIndex, double x );
    void UnitTest();

protected:
    double tolerance, coarseningtol;
    std::vector<int> intervalsForRefinement, NodesForRemoval;
    std::vector<double> mpNodeToInsert;
    SpaceMesh oldmesh, refinedsmesh;
};


#endif // ADAPTIVESOLVERHEADERREF
