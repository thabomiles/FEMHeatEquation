#ifndef SPACEMESHHEADERREF
#define SPACEMESHHEADERREF
#include <iostream>
#include <vector>
#include <array>

class SpaceMesh
{
public:
    void Range( double lowerlimit, double upperlimit, std::vector<double>& Nodes);
    void CommonMesh( SpaceMesh& firstmesh, SpaceMesh& secondmesh );
    double TestFunctions(int nodeIndex, double x);
    bool Contained (std::vector<double> SpaceNodes);

    void CopySpaceMesh (const SpaceMesh& oldSpaceMesh);
    void GenerateSpaceMesh ( std::vector<double> SpaceNodes );
    void BisectIntervals (std::vector<int> &intervalsForBisection);
    void CoarsenIntervals (std::vector<int> &intervalsForCoarsening);

    void GenerateUniformMesh(double boundary, int numberofnodes);
    void GenerateDefaultSpaceMesh();
    void GloballyBisectSpaceMesh ();
    void InsertSpaceNode ( double xi );
    void RemoveSpaceNode (int i);
    void RefreshSpaceMesh();
    void PrintSpaceNodes ();
    void PrintSpaceMesh ();
    double ReadSpaceNode (int i);
    double ReadSpaceMesh (int i);
    int meshsize();
    int IndexAbove ( double x );


protected:

    std::vector<double> mpSpaceNodes ;
    std::vector<double> mpSpaceMesh;
    int mpmeshsize;

};


#endif // SPACEMESHHEADERREF


