#ifndef PPSRECOTRACK_H
#define PPSRECOTRACK_H
#include "FastSimulation/PPSFastObjects/interface/PPSBaseTrack.h"
#include "FastSimulation/PPSFastObjects/interface/PPSTrackerHit.h"
#include "FastSimulation/PPSFastObjects/interface/PPSToFHit.h"
#include "TObject.h"
#include <vector>
class PPSRecoTrack:public PPSBaseTrack {
public:
       PPSRecoTrack(){};
       PPSRecoTrack(const PPSBaseTrack& trk):PPSBaseTrack(trk){};
       PPSRecoTrack(const TLorentzVector& trk,double t, double xi):PPSBaseTrack(trk,t,xi){};

       virtual ~PPSRecoTrack(){};

       void set_HitDet1(double x,double y) {Det1.set_Hit(x,y);};
       void set_HitDet2(double x,double y) {Det2.set_Hit(x,y);};
       void set_HitToF(double tof,double x, double y) {ToF.set_Hit(tof,x,y);};
       void set_ThetaAtIP(double thx,double thy) {thetaX=thx;thetaY=thy;};
       void   set_X0(double x) {X0=x;};
       void   set_Y0(double y) {Y0=y;};
       void   set_Phi(double p) {phi=p;};
       void   set_XTrackChiSquare(double chiSquare) {xTrackChiSquare = chiSquare;};
       void   set_YTrackChiSquare(double chiSquare) {yTrackChiSquare = chiSquare;};
       double set_TimeOfFlight(double time) {timeOfFlight=time;};
       double get_ImpPar() {return sqrt(X0*X0+Y0*Y0);};
       double get_X0() {return X0;}; 
       double get_Y0() {return Y0;};
       double get_XTrackChiSquare() {return xTrackChiSquare;};
       double get_YTrackChiSquare() {return yTrackChiSquare;};
       double get_TimeOfFlight() {return timeOfFlight;};

       PPSTrackerHit Det1;
       PPSTrackerHit Det2;
       PPSToFHit     ToF;
       double        X0;
       double        Y0;
       double        thetaX;   // in urad
       double        thetaY;   // in urad
       double        phi;
       double        xTrackChiSquare;
       double        yTrackChiSquare;
       double        timeOfFlight;

public:
ClassDef(PPSRecoTrack,1);
};
typedef std::vector<PPSRecoTrack> PPSRecoTracks;
#endif
