#include "FastSimulation/PPSFastSim/interface/PPSPixelDetector.h"

PPSPixelDetector::PPSPixelDetector() 
{
	Clear();
}

void PPSPixelDetector::Clear(){
	fNumberOfHits=0; 
	fHits.clear();

}

void PPSPixelDetector::AddHit(TVector3 hitVector3)
{
// Detector is in the negative side, but DetectorPosition is a positive number
    fHits.push_back(hitVector3);
    fNumberOfHits++;
 	return;

}