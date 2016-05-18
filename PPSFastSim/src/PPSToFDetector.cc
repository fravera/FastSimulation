#include "FastSimulation/PPSFastSim/interface/PPSToFDetector.h"
#include <math.h>

//=====================================================================================================

PPSToFDetector::PPSToFDetector(double detWidth, double detHeight, double detPosition):
        DetWidth(detHeight),DetHeight(detHeight),DetPosition(detPosition){
};

//--------------------------------------------------------------------------------------------------------------//

void PPSToFDetector::AddHit(double x, double y, double tof) {
   if (x>0) return; // The detector is on the negative side
   if (fabs(x)>DetWidth+DetPosition) return; // hit beyond detector area (W)
   if (fabs(x)<DetPosition) return;               // hit falls below detector area
   if (fabs(y)>DetHeight) return;                 // hit falls beyond detector area (H)
   X.push_back(x);
   Y.push_back(0.);
   ToF.push_back(tof);
   NHits++;
}

