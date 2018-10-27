#ifndef TIMEMESHHEADERREF
#define TIMEMESHHEADERREF
#include <iostream>
#include <vector>
#include <array>

class TimeMesh
{
public:

    void GenerateTimeMesh ( std::vector<double> TimeNodes );
    void GloballyBisectTimeMesh ();
    void InsertTimeNode ( double ti );
    void RefreshTimeMesh();
    void PrintTimeNodes ();
    void PrintTimeMesh ();
    double ReadTimeStep (int i);
    double ReadTimeMesh (int i);
    int NumberOfTimeSteps();




protected:

    int mpnumberOfTimeSteps;
    std::vector<double> mpTimeNodes;
    std::vector<double> mpTimeMesh;

};


#endif // TIMEMESHHEADERREF


