#ifndef PPSPixelDetector_h
#define PPSPixelDetector_h
#include <vector>
#include "FastSimulation/PPSFastSim/interface/PPSTrkDetector.h"

class PPSPixelDetector : public PPSTrkDetector{
      public:
            PPSPixelDetector();
            ~PPSPixelDetector() {};
            // virtual ~PPSPixelDetector() {};
            void AddHit(TVector3 hitVector3);
	        void Clear();
        

};

#endif
