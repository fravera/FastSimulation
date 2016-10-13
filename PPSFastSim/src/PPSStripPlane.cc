#include "FastSimulation/PPSFastSim/interface/PPSStripPlane.h"
#include <math.h>
#include <algorithm>
#include <vector>
#include <TRandom3.h>


PPSStripPlane::PPSStripPlane(TVector3 stripPlanePosition,TVector3 stripPlaneRotation, int numberOfStrips, double pitchStrip, double cutSideLength):
			fPlaneId(-1),
			fNumberOfStrips(numberOfStrips),
      		fPitchStrip(pitchStrip),
      		fCutSideLength(cutSideLength){
	Clear();
}


edm::DetSet<TotemRPDigi> PPSStripPlane::FromHitsToDigi(){
	
	std::vector<int> hitStrips;
    edm::DetSet<TotemRPDigi> planeDigi(fPlaneId);

	for(unsigned hIt = 0; hIt<fHits.size(); ++hIt){
		TVector3 localHit = MoveToPlaneSystem(fHits.at(hIt));
		if(!IsInPlane(localHit)) continue;
		localHit.SetY(localHit.Y()+fCutSideLength/2.); //moves the origin in corrispondence of the sensor corner (strip 0)
		localHit.RotateZ(-M_PI/4.); //rotate the sensor to have x axis corresponding to strip 0
		int clusteSize = (int)fClusterSizePlot->GetRandom() + 0.5;
		int centralStrip = localHit.Y()/fPitchStrip;
		int stripMin = -1;
		int stripMax = -1;
		if(clusteSize % 2 == 1){
			stripMin = centralStrip - (clusteSize-1)/2;
			stripMin = stripMin >= 0 ? stripMin : 0;
			stripMax = centralStrip + (clusteSize-1)/2;
			stripMax = stripMax <= fNumberOfStrips  ? stripMax : fNumberOfStrips;
		}
		else{
			stripMin = centralStrip - (clusteSize-2)/2;
			stripMax = centralStrip + (clusteSize-2)/2;
			if(localHit.Y() - centralStrip*fPitchStrip > 0.5) ++stripMax;
			else --stripMin;
			stripMin = stripMin >= 0 ? stripMin : 0;
			stripMax = stripMax <= fNumberOfStrips  ? stripMax : fNumberOfStrips;
		}
		for(int hStr = stripMin; hStr<=stripMax; ++hStr){
			hitStrips.push_back(hStr);
		}

	}

	sort(hitStrips.begin(),hitStrips.end());
	std::vector<int>::iterator it;
	it = unique(hitStrips.begin(),hitStrips.end());
	hitStrips.resize(distance(hitStrips.begin(),it));

	for(unsigned hStr=0; hStr<hitStrips.size(); ++hStr){
		planeDigi.push_back(hitStrips.at(hStr));
	}

	return planeDigi;

}

TVector3 PPSStripPlane::MoveToPlaneSystem(TVector3 globalHit){

	TVector3 localHit = globalHit;
	localHit.RotateZ(-fPlaneRotation.Z()); //takes into account that station near are rotated by 8 degree
	localHit.SetY(localHit.Y()-fPlanePosition.Y()); //moves the origin in corrispondence of the sensor center
	localHit.RotateX(-fPlaneRotation.X()); //takes into account that v station are back flipped
	localHit.SetX(localHit.X()-fPlanePosition.X()); //moves the origin to the strip sensor edge
	localHit.RotateY(-fPlaneRotation.Y()); //rotated the plane, not needed for strips but needed for pixel

	return localHit;
}

bool PPSStripPlane::IsInPlane(TVector3 localHit){
	double xx = localHit.X();
	double yy = localHit.Y();

	if(xx > 0.) return false;
	else if (yy < +xx - fCutSideLength/2.) return false;
	else if (yy > +xx - fCutSideLength/2. + fNumberOfStrips*fPitchStrip*sqrt(2)) return false;
	else if (yy < -xx + fCutSideLength/2. - fNumberOfStrips*fPitchStrip*sqrt(2)) return false;
	else if (yy > -xx + fCutSideLength/2.) return false;
	else return true;
	
}

bool PPSStripPlane::AddHit(TVector3 hit){

	fHits.push_back(hit);
	return IsInPlane(hit);

}