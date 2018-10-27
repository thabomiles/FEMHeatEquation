#ifndef SPACEMESHHEADERREF
#define SPACEMESHHEADERREF
#include <iostream>
#include <vector>
#include <array>

class SpaceMesh
{
public:

    void GenerateSpaceMesh ( std::vector<double> SpaceNodes );
    void GloballyBisectSpaceMesh ();
    void InsertSpaceNode ( double xi );
    void RefreshSpaceMesh();
    void PrintSpaceNodes ();
    void PrintSpaceMesh ();
    double ReadSpaceNode (int i);
    double ReadSpaceMesh (int i);
    int meshsize();



protected:

    std::vector<double> mpSpaceNodes ;
    std::vector<double> mpSpaceMesh;
    int mpmeshsize;

};


#endif // SPACEMESHHEADERREF


