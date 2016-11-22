#ifndef PPSToFDetector_h
#define PPSToFDetector_h
#include <vector>
#include <map>
#include <iterator>
#include <iostream>

class PPSToFDetector {
      public:
            PPSToFDetector(double detWidth, double detHeight, double detPosition);
            ~PPSToFDetector() {};
            // virtual ~PPSToFDetector() {};

            void SetDetectorName(std::string detectorName) {fDetectorName = detectorName;};
            std::string GetDetectorName() {return fDetectorName;};

            double GetHeight()                       {return DetHeight;};
            double GetWidth()                        {return DetWidth;};
            double GetPosition()                     {return DetPosition;};
            
            unsigned get_NHits()                     {return NHits;}; // return the total hit multiplicity (full det)
           

            void AddHit(double x, double y,double tof);
            void clear() {
                  DetId=0;
                  NHits=0;
                  X.clear();
                  Y.clear();
                  ToF.clear();
            };
            
            std::vector<double> X;
            std::vector<double> Y;
            std::vector<double> ToF;

      private:
            double              DetWidth; // detector width
            double              DetHeight; // detector height
            double              DetPosition; // detector position from beam (absolute value)
            std::string         fDetectorName;

            int                 DetId;
            int                 NHits;

};
#endif
