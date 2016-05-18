#ifndef PPSToFDetector_h
#define PPSToFDetector_h
#include <vector>
#include <map>
#include <iterator>
#include <iostream>

class PPSToFDetector {
      public:
            PPSToFDetector(double detWidth, double detHeight, double detPosition);
            virtual ~PPSToFDetector() {};

            double GetHeight()                       {return DetH;};
            double GetWidth()                        {return DetW;};
            double GetPosition()                     {return DetPosition;};
            
            int get_NHits()                          {return NHits;}; // return the total hit multiplicity (full det)
           

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
            double              DetW; // detector width
            double              DetH; // detector height
            double              DetPosition; // detector position from beam (absolute value)

            int                 DetId;
            int                 NHits;

};
#endif
