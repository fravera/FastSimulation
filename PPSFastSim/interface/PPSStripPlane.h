#ifndef PPSStripPlane_h
#define PPSStripPlane_h
#include <vector>
#include <TVector3.h>
#include <TH1D.h>
#include "DataFormats/CTPPSDigi/interface/TotemRPDigi.h"
#include "DataFormats/Common/interface/DetSet.h"

class PPSStripPlane{
      public:
            PPSStripPlane(){};
      	PPSStripPlane(TVector3 stripPlanePosition,TVector3 stripPlaneRotation, int numberOfStrips, double pitchStrip, double cutSideLength);
            ~PPSStripPlane() {};
            // virtual ~PPSStripPlane() {};

            void Clear(){fHits.clear();};

            void SetHits(std::vector<TVector3> hits) {fHits=hits;};
            bool AddHit(TVector3 hit);
            void SetClusterSizePlot(TH1D *clusterSizePlot) {fClusterSizePlot=clusterSizePlot;};
            void SetPlaneId(int planeId) {fPlaneId = planeId;};
            void SetStripPlanePosition(TVector3 stripPlanePosition) {fPlanePosition = stripPlanePosition;};
            void SetStripPlaneRotation(TVector3 stripPlaneRotation) {fPlaneRotation = stripPlaneRotation;};
            void SetNumberOfStrips(int numberOfStrips) {fNumberOfStrips=numberOfStrips;};
            void SetPitchStrip(double pitchStrip) {fPitchStrip=pitchStrip;};
            void SetCutSideLength(double cutSideLength) {fCutSideLength=cutSideLength;};

            edm::DetSet<TotemRPDigi> FromHitsToDigi();
            TVector3 MoveToPlaneSystem(TVector3 globalHit);
            bool IsInPlane(TVector3 localHit);

            int GetPlaneId()   {return fPlaneId;};


      private:
            int fPlaneId;
      	int fNumberOfStrips;
      	double fPitchStrip;
      	double fCutSideLength;
      	TVector3 fPlanePosition;
      	TVector3 fPlaneRotation;
            // edm::DetSet<TotemRPDigi> fPlaneDigi;
            std::vector<TVector3> fHits;
            TH1D *fClusterSizePlot;

};

#endif
