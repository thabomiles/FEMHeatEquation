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
    double GeneralTestFunctions(int nodeIndex, double x);
    bool Contained (double my_var );

    void CopySpaceMesh (const SpaceMesh& oldSpaceMesh);
    void GenerateSpaceMesh ( std::vector<double> SpaceNodes );
    void BisectIntervals (std::vector<int> &intervalsForBisection);
    void CoarsenIntervals (std::vector<int> &intervalsForCoarsening);
    void ReadLastSpaceNode();

    void GenerateUniformMesh(double boundary, int numberofnodes);
    void GenerateDefaultSpaceMesh();
    void GloballyBisectSpaceMesh ();
    void InsertSpaceNode ( double xi );
    void InsertArray ( std::vector<double> RefinementNodes );
    void RemoveSpaceNode (int i);
    void RefreshSpaceMesh();
    void PrintSpaceNodes ();
    void PrintSpaceMesh ();
    double ReadSpaceNode (int i);
    double ReadSpaceMesh (int i);
    int meshsize();
    int IndexAbove ( double x );

    std::vector<double> mpSpaceNodes;

protected:


    //std::vector<double> mpSpaceMesh;
    int mpmeshsize;

};


#endif // SPACEMESHHEADERREF


