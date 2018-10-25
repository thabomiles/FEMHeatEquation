#ifndef GENERALMESHHEADERREF
#define GENERALMESHHEADERREF
#include <iostream>
#include <vector>
#include <array>

class GeneralMesh
{
public:

    void GeneralSpaceMesh ( std::vector<double> SpaceNodes );
    void GeneralTimeMesh ( std::vector<double> TimeNodes );
    void GloballyBisectSpaceMesh ();
    void GloballyBisectTimeMesh ();
    void InsertSpaceNode ( double xi );
    void InsertTimeNode ( double ti );
    void RefreshSpaceMesh();
    void RefreshTimeMesh();
    void PrintSpaceNodes ();
    void PrintTimeNodes ();
    void PrintSpaceMesh ();
    void PrintTimeMesh ();




protected:

    int mpn;

    std::vector<double> mpSpaceNodes ;
    std::vector<double> mpTimeNodes;
    std::vector<double> mpSpaceMesh;
    std::vector<double> mpTimeMesh;

};


#endif // GENERALMESHHEADERREF


