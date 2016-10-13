#ifndef PPSStripDetector_h
#define PPSStripDetector_h
#include <vector>
#include <map>
#include <TH1D.h>
#include "FastSimulation/PPSFastSim/interface/PPSTrkDetector.h"
#include "FastSimulation/PPSFastSim/interface/PPSStripPlane.h"
#include "DataFormats/Common/interface/DetSet.h"

class PPSStripDetector : public PPSTrkDetector{
      public:
        PPSStripDetector();
        ~PPSStripDetector() {};
        // virtual ~PPSStripDetector() {};
        void BuildDetectorPlanes();

        std::vector<edm::DetSet<TotemRPDigi> > GetClusters();

        void AddHit(TVector3 hitVector3);
        void Clear();
        
        void SetNumberOfStrips(int numberOfStrips)  { fNumberOfStrips=numberOfStrips;};
        void SetStripPitch(double stripPitch)       { fStripPitch=stripPitch;};
        void SetCutSideLength(double cutSideLength) { fCutSideLength=cutSideLength;};
        void SetClusterSizePlot(TH1D *clusterSizePlot) {fClusterSizePlot=clusterSizePlot;};
        void SetRpPosition(TVector3 rpPosition) {fRpPosition = rpPosition;};

      private:
      	int fNumberOfStrips;
      	double fStripPitch;
      	double fCutSideLength;
      	TH1D *fClusterSizePlot;
        TVector3 fRpPosition;
        std::map<int,PPSStripPlane> fStripPlaneMap;

};

#endif
