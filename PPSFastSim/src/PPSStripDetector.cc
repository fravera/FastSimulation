#include "FastSimulation/PPSFastSim/interface/PPSStripDetector.h"
#include "DataFormats/Common/interface/DetSet.h"
#include <math.h>
#include <sstream>
#include <string>
#include <set>
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "DataFormats/CTPPSDetId/interface/TotemRPDetId.h"
#include "Math/GenVector/Rotation3D.h"
#include "Math/GenVector/RotationZYX.h"
#include "Math/GenVector/3DConversions.h"

//#include "3DConversions.h"

PPSStripDetector::PPSStripDetector(){
	// BuildDetectorPlanes();
	Clear();
}

void PPSStripDetector::Clear(){
	fNumberOfHits=0; 
	fHits.clear();
	for(std::map<int,PPSStripPlane>::iterator pIt=fStripPlaneMap.begin(); pIt!=fStripPlaneMap.end(); ++pIt){
		pIt->second.Clear();
	}

}

void PPSStripDetector::BuildDetectorPlanes(){
	// double rpLength = 30.; //mm
	// int numberOfPlanes = 10;
	// double interPlanesDistance = (double)rpLength/(numberOfPlanes/2-1);

	// std::cout<<"BuildDetectorPlanes "<<1<<std::endl;
	const std::set<unsigned int> &dets = fRpGeometry->DetsInRP((unsigned int)fDetectorId/10);
	// std::cout<<"BuildDetectorPlanes "<<2<<std::endl;
	std::cout<<"number of detectors in RP"<< dets.size()<<std::endl;

	for (std::set<unsigned int>::iterator dIt = dets.begin(); dIt != dets.end(); ++dIt){
		unsigned int rawId = TotemRPDetId::decToRawId(*dIt);
		std::cout<<*dIt<<" -> "<<rawId<<std::endl;
		DetGeomDesc *stripPlaneDescription = fRpGeometry->GetDetector(rawId);
	 	TVector3 stripPlanePosition(stripPlaneDescription->translation().x(),stripPlaneDescription->translation().y(),stripPlaneDescription->translation().z());
	 	ROOT::Math::Rotation3D rotation3D = stripPlaneDescription->rotation();
	 	ROOT::Math::RotationZYX rotationZYX;
	 	ROOT::Math::gv_detail::convert(rotation3D,rotationZYX);
	 	TVector3 stripPlaneRotation(rotationZYX.Phi(),rotationZYX.Theta(),rotationZYX.Psi());
	 	// PPSStripPlane *stripPlane = new PPSStripPlane(stripPlanePosition,stripPlaneRotation, fNumberOfStrips, fStripPitch, fCutSideLength);
	 	PPSStripPlane *stripPlane = new PPSStripPlane();
	 	stripPlane->SetStripPlanePosition(stripPlanePosition);
        stripPlane->SetStripPlaneRotation(stripPlaneRotation);
        stripPlane->SetNumberOfStrips(fNumberOfStrips);
        stripPlane->SetPitchStrip(fStripPitch);
        stripPlane->SetCutSideLength(fCutSideLength);

	 	stripPlane->SetClusterSizePlot(fClusterSizePlot);
	 	stripPlane->SetPlaneId(rawId);
	 	fStripPlaneMap[rawId] = *stripPlane;

	}
	// for(int i=0; i<numberOfPlanes; i+=2){
	
	// 	stringstream ss;
	// 	ss.str("");
	// 	ss << setw(2) << setfill('0') << i+1;
	// 	string uPlaneNumber = ss.str();
	// 	string uPlaneName = fDetectorName + uPlaneNumber;
	// 	ss.str("");
	// 	ss << setw(2) << setfill('0') << i+2;
	// 	string vPlaneNumber = ss.str();
	// 	string vPlaneName = fDetectorName + vPlaneNumber;
	// 	double zPlanePositon= (double)-rpLength/2.+i/2.*interPlanesDistance;


	// 	TVector3 uStripPlanePosition(fDetectorPosition.X(),fDetectorPosition.Y(),fDetectorPosition.Z()+zPlanePositon);
	// 	TVector3 vStripPlanePosition(fDetectorPosition.X(),fDetectorPosition.Y(),fDetectorPosition.Z()+zPlanePositon);

	// 	TVector3 uStripPlaneRotation(fDetectorRotation.X(),fDetectorRotation.Y(),fDetectorRotation.Z());
	// 	TVector3 vStripPlaneRotation(fDetectorRotation.X()+M_PI,fDetectorRotation.Y(),fDetectorRotation.Z());

	// 	PPSStripPlane *uStripPlane = new PPSStripPlane(uStripPlanePosition,uStripPlaneRotation, fNumberOfStrips, fPitchStrip, fCutSideLength);
	// 	uStripPlane->SetClusterSizePlot(fClusterSizePlot);
	// 	uStripPlane->SetPlaneId(fDetectorId + i);

	// 	PPSStripPlane *vStripPlane = new PPSStripPlane(vStripPlanePosition,vStripPlaneRotation, fNumberOfStrips, fPitchStrip, fCutSideLength);
	// 	vStripPlane->SetClusterSizePlot(fClusterSizePlot);
	// 	vStripPlane->SetPlaneId(fDetectorId + i + 1);

	// 	fStripPlaneMap[uPlaneName] = uStripPlane;
	// 	fStripPlaneMap[vPlaneName] = vStripPlane;

 // 	}

}


std::vector<edm::DetSet<TotemRPDigi> > PPSStripDetector::GetClusters(){

	//  stripZeroPosition = RPTopology::last_strip_to_border_dist_ + (RPTopology::no_of_strips_-1)*RPTopology::pitch_ - RPTopology::y_width_/2.;

	// for (std::set<unsigned int>::iterator dit = dets.begin(); dit != dets.end(); ++dit){
	// 	unsigned int rawId = TotemRPDetId::DecToRawId(*dit);
	// 	DetGeomDesc *ginfo = geom->GetDetector(rawId);
	// 	double zGlobalHit = ginfo->translation().z();
	// 	double xGlobalHit = bx + ax * (z - z0);
	// 	double yGlobalHit = by + ay * (z - z0);

	// 	CLHEP::Hep3Vector globalHit(xGlobalHit,yGlobalHit,zGlobalHit);
	// 	CLHEP::Hep3Vector localHit = geom->GlobalToLocal(rawId, globalHit);
	// 	double uLocalHit = localHit.x();
	// 	double vLocalHit = localHit.y();

	// 	if (!RPTopology::IsHit(uLocalHit, vLocalHit, insensitiveMargin)) continue;
 


	std::vector<edm::DetSet<TotemRPDigi> > rpCluster;
	for(std::map<int,PPSStripPlane>::iterator pIt=fStripPlaneMap.begin(); pIt!=fStripPlaneMap.end(); ++pIt){
		rpCluster.push_back(pIt->second.FromHitsToDigi());
	}

	return rpCluster;
	// poi occorre mergiare i vettori, creare il DetSetVector ed infibe fare lo short
}

void PPSStripDetector::AddHit(TVector3 hitVector3)
{
// Detector is in the negative side, but DetectorPosition is a positive number
    fHits.push_back(hitVector3);
    fNumberOfHits++;
    int numberOfHitPlanes = 0;
    for(std::map<int,PPSStripPlane>::iterator pIt=fStripPlaneMap.begin(); pIt!=fStripPlaneMap.end(); ++pIt){
		bool isHit = pIt->second.AddHit(hitVector3);
		if(isHit) ++numberOfHitPlanes;
	}

	if(numberOfHitPlanes>=6) fSmearedHits.push_back(HitSmearing(hitVector3));

	return;

}