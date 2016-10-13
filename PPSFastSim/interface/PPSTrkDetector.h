#ifndef PPSTrkDetector_h
#define PPSTrkDetector_h
#include <vector>
#include <string>
#include <TVector3.h>
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/VeryForwardGeometryBuilder/interface/TotemRPGeometry.h"

class PPSTrkDetector {
      public:
            PPSTrkDetector() {};
            PPSTrkDetector(TVector3 detectorPosition, TVector3 detectorRotation);
            ~PPSTrkDetector() {};
            // virtual ~PPSTrkDetector() {};

            void SetDetectorId(int detectorId) {fDetectorId = detectorId;};
            void SetDetectorName(std::string detectorName) {fDetectorName = detectorName;};
            void SetDetectorResolution(double hitSigmaX, double hitSigmaY, double hitSigmaZ){
                  fHitSigmaX = hitSigmaX;
                  fHitSigmaY = hitSigmaY;
                  fHitSigmaZ = hitSigmaZ; 
            };
            std::string GetDetectorName() {return fDetectorName;};

            // void SetDetectorPosition(TVector3 detectorPosition) {fDetectorPosition=detectorPosition;};
            // void SetDetectorRotation(TVector3 detectorRotation) {fDetectorRotation=detectorRotation;};

            int GetNumberOfHits() {return fNumberOfHits;};
            int GetDetectorId()   {return fDetectorId;};
            std::vector<TVector3> GetHits() {return fHits;};
            std::vector<TVector3> GetSmearedHits() {return fHits;};


            void SetRpGeometry(edm::ESHandle<TotemRPGeometry> rpGeometry) {fRpGeometry=rpGeometry;};

            virtual void Clear() {}; 
            virtual void AddHit(TVector3 hitVector3) {};
            TVector3 HitSmearing(TVector3 hitVector3);

      protected:
            int                 fDetectorId;
            int                 fNumberOfPlanes;
            TVector3            fDetectorPosition;
            TVector3            fDetectorRotation;
            int                 fNumberOfHits;
            std::vector<TVector3>    fHits;
            std::vector<TVector3>    fSmearedHits;
            std::string         fDetectorName;
            edm::ESHandle<TotemRPGeometry> fRpGeometry;
            double fHitSigmaX;
            double fHitSigmaY;
            double fHitSigmaZ; 
};

#endif
