#ifndef PPSTOFHIT_H
#define PPSTOFHIT_H
#include <vector>
#include "TObject.h"

class PPSToFHit: public TObject {
public:
       PPSToFHit();
       PPSToFHit(double tof, double x,double y) {ToF=tof;X=x;Y=y;};
       virtual ~PPSToFHit() {};
       void set_Hit(double tof,double x, double y){ToF=tof;X=x;Y=y;};

       double X; // x of the cell center
       double Y; // y of the cell center
       double ToF;
ClassDef(PPSToFHit,1);
};
#endif
