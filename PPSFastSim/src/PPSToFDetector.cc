#include "FastSimulation/PPSFastSim/interface/PPSToFDetector.h"
#include <math.h>

//=====================================================================================================

PPSToFDetector::PPSToFDetector(double detWidth, double detHeight, double detPosition):
        DetW(detHeight),DetH(detHeight),DetPosition(detPosition){
};

//--------------------------------------------------------------------------------------------------------------//

void PPSToFDetector::AddHit(double x, double y, double tof) {
   if (x>0) return; // The detector is on the negative side
   if (fabs(x)>DetectorWidth+DetectorPosition) return; // hit beyond detector area (W)
   if (fabs(x)<DetectorPosition) return;               // hit falls below detector area
   if (fabs(y)>DetectorHeight) return;                 // hit falls beyond detector area (H)
   X.push_back(x);
   Y.push_back(0.);
   ToF.push_back(tof);
   NHits++;
}

