#include <math.h>
#include <TRandom3.h>
#include "FastSimulation/PPSFastSim/interface/PPSTrkDetector.h"

PPSTrkDetector::PPSTrkDetector(TVector3 detectorPosition, TVector3 detectorRotation):
				fDetectorPosition(detectorPosition),
				fDetectorRotation(detectorRotation)
{
	// Clear();
}


TVector3 PPSTrkDetector::HitSmearing(TVector3 hitVector3)
{
    double x = gRandom->Gaus(hitVector3.X(),fHitSigmaX);
    double y = gRandom->Gaus(hitVector3.Y(),fHitSigmaY);
    double z = gRandom->Gaus(hitVector3.Z(),fHitSigmaZ);
    return TVector3(x,y,z);
}
