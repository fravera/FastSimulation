#ifndef PPSSim_h
#define PPSSim_h
// ROOT #includes
#include "TROOT.h"
#include "Rtypes.h"
#include "TH2F.h"
#include "TDirectory.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TGraphErrors.h"
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <math.h> 
#include <string.h>
#include <map>
// Framework includes
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/CTPPSDetId/interface/TotemRPDetId.h"
#include "DataFormats/CTPPSDigi/interface/TotemRPDigi.h"

#include "HepMC/GenEvent.h"
#include "HepMC/GenVertex.h"
#include "HepMC/GenParticle.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "FastSimulation/PPSFastObjects/interface/PPSSpectrometer.h"
#include "FastSimulation/PPSFastSim/interface/PPSTrkDetector.h"
#include "FastSimulation/PPSFastSim/interface/PPSStripDetector.h"
#include "FastSimulation/PPSFastSim/interface/PPSPixelDetector.h"
#include "FastSimulation/PPSFastSim/interface/PPSToFDetector.h"
// Hector #includes
#include "H_BeamLine.h"
#include "H_BeamParticle.h"
#include "H_RecRPObject.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
/*
   The virtuality coming from the reconstruction in Hector (H_RecRPObject) must be summed with
   twice the proton mass squared to recover the quadrimomentum lost of the scattered proton
   */
#include "FastSimulation/PPSFastSim/interface/PPSConstants.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

//=====================================================================================================

typedef ROOT::Math::SMatrix<double,4,4> matrix44Def;
typedef ROOT::Math::SVector<double,4>   vector4Def;
typedef pair<double,vector4Def>         vector5Def;


class PPSSim {
    public:
        PPSSim(bool =false);
        ~PPSSim() {};

        void BeginRun();
        void BeginEvent();
        void EndRun();
        void EndEvent();


        void LoadTwissMatrices();
        std::map<string,std::pair<vector5Def, matrix44Def> > GetMatrixFromTwissFile(string inputTwissFileName, std::vector<string> beamObjects, double &protonMomentumXi);
        map<string,vector5Def> PropagateParticle(int direction, pair<double,vector4Def> protonKinematics);
        int GetXiClassNumber(int direction, double protonXi);

        //void set_Strengths();
        void set_Verbose(bool f)   {fVerbose= f;};
        void set_KickersOFF(bool k){fKickersOFF=k;};
        void set_etaMin(double em) {fEtaMin = em;};
        void set_momentumMin(double pm){fMomentumMin=pm;};
        void set_phiMin(double pm) {fPhiMin = pm;};
        void set_phiMax(double pm) {fPhiMax = pm;};
        void set_xiMin(double xi)  {xi_min  = xi;};
        void set_xiMax(double xi)  {xi_max  = xi;};
        void set_tMin(double t)    {t_min   = t;};
        void set_tMax(double t)    {t_max   = t;};
        void set_CentralMass(double m,double me) {fCentralMass=m;fCentralMassErr=me;};

        //Tracker characteristics
        void set_TrackerEdgeOffset  (double x1F, double x2F, double x1B, double x2B) {
            fTrk1XOffsetF=x1F; fTrk2XOffsetF=x2F; fTrk1XOffsetB=x1B; fTrk2XOffsetB=x2B;
        };
        void set_TrackerLength      (double p)     {
            fTrackerLength=p;
        };
        void set_TrackerSize        (double w,double h) {
            fTrackerWidth=w;fTrackerHeight=h;
        };
        void set_HitSmearing        (bool f=true)    {
            fSmearHit = f;
        };
        void set_TrackerResolution  (double hitSigmaX, double hitSigmaY) {
            fHitSigmaX=hitSigmaX;
            fHitSigmaY=hitSigmaY;
        };
        void set_TrackerGeometry    (string trackerGeometry) {
            fTrackerGeometry = trackerGeometry;
        };
        void set_VerticalShift      (double verticalShift) {
            fVerticalShift = verticalShift;
        };
        void set_TrkDetXOffset      (double trkDetXOffset) {
            fTrkDetXOffset = trkDetXOffset;
        };
        void set_NumberOfStrips     (int numberOfStrips) {
            fNumberOfStrips = numberOfStrips;
        };
        void set_StripPitch         (double stripPitch) {
            fStripPitch = stripPitch;
        };
        void set_CutSideLength      (double cutSideLength) {
            fCutSideLength = cutSideLength;
        };
        void set_ClusterSizePlotFile(string clusterSizePlotFile) {
            fClusterSizePlotFile = clusterSizePlotFile;
        }
        void set_ClusterSizePlotName(string clusterSizePlotName) {
            fClusterSizePlotName = clusterSizePlotName;
        }
        void set_NumberOfRows       (int numberOfRows) {
            fNumberOfRows = numberOfRows;
        };
        void set_NumberOfColumns    (int numberOfColumns) {
            fNumberOfColumns = numberOfColumns;
        };
        void set_DoubleSizeColumn   (vector<int> doubleSizeColumn) {
            fDoubleSizeColumn = doubleSizeColumn;
        };
        void set_DoubleSizeRow      (vector<int> doubleSizeRow) {
            fDoubleSizeRow = doubleSizeRow;
        };
        void set_PixelPitchX        (double pixelPitchX) {
            fPixelPitchX = pixelPitchX;
        };
        void set_PixelPitchY        (double pixelPitchY) {
            fPixelPitchY = pixelPitchY;
        };
        void set_TrackingDetectorIDs(vector<int> trackingDetectorIDs)  {
            fTrackingDetectorIDs = trackingDetectorIDs;
        };
        void set_TrackingDetectorNameMap(std::map<int,string> beam1TrackingDetectorNameMap, std::map<int,string> beam2TrackingDetectorNameMap)  {
            fBeam1TrackingDetectorNameMap = beam1TrackingDetectorNameMap;
            fBeam2TrackingDetectorNameMap = beam2TrackingDetectorNameMap;
        };
        void set_TrackingStationXRotationMap(std::map<int,double> trackingStationXRotationMap) {
            fTrackingStationXRotationMap = trackingStationXRotationMap;
        };
        void set_TrackingStationYRotationMap(std::map<int,double> trackingStationYRotationMap) {
            fTrackingStationYRotationMap = trackingStationYRotationMap;
        };
        void set_TrackingStationZRotationMap(std::map<int,double> trackingStationZRotationMap) {
            fTrackingStationZRotationMap = trackingStationZRotationMap;
        };
        void set_RomanPotGeometry(edm::ESHandle<TotemRPGeometry> rpGeometry) {
            fRpGeometry = rpGeometry;
            // if(fRpGeometry==NULL){
            //     throw cms::Exception("PPSSim.h") << "Very Forward Misaligned Geometry Record not found.";
            // }
        };
       
        //Time of Flight characteristics
        void set_ToFEdgeOffset(double toFXOffsetF, double toFXOffsetB){
            fToFXOffsetF = toFXOffsetF; fToFXOffsetB = toFXOffsetB;};
        void set_ToFSize(double w,double h) {fToFWidth=w;fToFHeight=h;};
        void set_ToFGeometry(std::string tofgeometry)   {fToFGeometry=tofgeometry;};
        void set_ToFDetectorName(string beam1ToFDetectorName, string beam2ToFDetectorName){
            fBeam1ToFDetectorName = beam1ToFDetectorName;
            fBeam2ToFDetectorName = beam2ToFDetectorName;
        };
        void set_ToFResolution(double p, double smearingX, double smearingY)  {fTimeSigma=p; fToFHitSigmaX = smearingX; fToFHitSigmaY = smearingY;};


        void set_FilterHitMap(bool f)        {fFilterHitMap=f;};
        // void set_WindowForTrack(double x, double y,double c)
        //     {fMaxXfromBeam=x;fMaxYfromBeam=y;fDetectorClosestX=c;};
        void set_ApplyFiducialCuts(bool f)   {fApplyFiducialCuts=f;};
        void set_UseToFForTracking(bool f)   {fUseToFForTracking=f;};
        // void set_TCLPosition(const string& tcl,double z1,double z2) {
        //     if (tcl=="TCL4")      {fTCL4Position1=z1;fTCL4Position2=z2;}
        //     else if (tcl=="TCL5") {fTCL5Position1=z1;fTCL5Position2=z2;}
        //     else edm::LogWarning("debug")  <<"WARNING: Unknown Collimator " << tcl ;
        // }

        
        void set_Beam1StationPositions(std::map<string,double> beam1StationPositions){
            fBeam1StationPositions = beam1StationPositions;
        };
        void set_Beam2StationPositions(std::map<string,double> beam2StationPositions){
            fBeam2StationPositions = beam2StationPositions;
        };
        void set_Beam1CenterAtStation(std::map<string,double> beam1CenterAtStation){
            fBeam1CenterAtStation = beam1CenterAtStation;
        };
        void set_Beam1SigmaAtStation(std::map<string,double> beam1SigmaAtStation){
            fBeam1SigmaAtStation = beam1SigmaAtStation;
        };
        void set_Beam1StationSigmaInsertion(std::map<string,double> beam1StationSigmaInsertion){
            fBeam1StationSigmaInsertion = beam1StationSigmaInsertion;
        };
        void set_Beam2CenterAtStation(std::map<string,double> beam2CenterAtStation){
            fBeam2CenterAtStation = beam2CenterAtStation;
        };
        void set_Beam2SigmaAtStation(std::map<string,double> beam2SigmaAtStation){
            fBeam2SigmaAtStation = beam2SigmaAtStation;
        };
        void set_Beam2StationSigmaInsertion(std::map<string,double> beam2StationSigmaInsertion){
            fBeam2StationSigmaInsertion = beam2StationSigmaInsertion;
        };

        void set_TCLNames(vector<string> beam1TCLNames, vector<string> beam2TCLNames){
            fBeam1TCLNames = beam1TCLNames;
            fBeam2TCLNames = beam2TCLNames;
        };

        void set_VertexSmearing(bool f=true) {fSmearVertex = f;};
        // void set_BeamLineFile(std::string b1, std::string b2) {fBeamLine1File=b1;fBeamLine2File=b2;};
        void set_BeamLineFilePath(string beam1FilePath,string beam2FilePath){
            fBeam1FilePath = beam1FilePath;
            fBeam2FilePath = beam2FilePath;
        };
        void set_BeamLineFiles(vector<string> beam1Files,vector<string> beam2Files){
            fBeam1Files = beam1Files;
            fBeam2Files = beam2Files;
        };
        void set_BeamLineFilesXi(vector<double> beam1FilesXi,vector<double> beam2FilesXi){
            fBeam1FilesXi = beam1FilesXi;
            fBeam2FilesXi = beam2FilesXi;
        };
        void set_BeamLineObjects(vector<string> beam1Objects,vector<string> beam2Objects){
            fBeam1Objects = beam1Objects;
            fBeam2Objects = beam2Objects;
        };
        void set_BeamDirection(int b1dir,int b2dir) {fBeam1Direction=b1dir;fBeam2Direction=b2dir;};
        void set_ShowBeamLine()                     {fShowBeamLine=true;};
        void set_GenBeamProfile()                   {fSimBeam=true;};
        void set_CollisionPoint(std::string ip)        {fCollisionPoint=ip;};
        void set_VertexPosition(const double x, const double y, const double z)
        {fVtxMeanX=x;fVtxMeanY=y;fVtxMeanZ=z;};

        int  add_Vertex(const double x, const double y, const double z)
        {fVertex[NVertex++]=TVector3(x,y,z);return NVertex;};
        void add_OutgoingParticle(int ivtx, const TLorentzVector* pF,const TLorentzVector* pB)
        {if (ivtx>=NVertex) {edm::LogWarning("debug")  <<"Error: Adding particle to a non existing vertex.";}
            protonsOut[ivtx]=std::pair<const TLorentzVector*,const TLorentzVector*>(pF,pB);
        };

        void set_CrossingAngleCorrection(bool f=false) {fCrossAngleCorr=f;};
        void set_CrossingAngle(double cr)         {fCrossingAngle=cr;};
        void set_BeamEnergy(double e)             {fBeamEnergy=e;fBeamMomentum = sqrt(fBeamEnergy*fBeamEnergy - ProtonMassSQ);};
        void set_BeamEnergySmearing(bool f=false) {fSmearEnergy=f;};
        void set_BeamEnergyRMS(double rms)        {fBeamEnergyRMS=rms;};
        void set_BeamAngleSmearing(bool f=false)  {fSmearAngle=f;};
        void set_BeamAngleRMS(double rms)         {fBeamAngleRMS=rms;};
        // void set_BeamXSizes(double bsig_ArmF_det1,double bsig_ArmF_det2,double bsig_ArmF_tof,double bsig_ArmB_det1,double bsig_ArmB_det2,double bsig_ArmB_tof) {
        //     fBeamXRMS_ArmB_Trk1=bsig_ArmF_det1;
        //     fBeamXRMS_ArmB_Trk2=bsig_ArmF_det2;
        //     fBeamXRMS_ArmB_ToF=bsig_ArmF_tof;
        //     fBeamXRMS_ArmF_Trk1=bsig_ArmB_det1;
        //     fBeamXRMS_ArmF_Trk2=bsig_ArmB_det2;
        //     fBeamXRMS_ArmF_ToF=bsig_ArmB_tof;
        // };
        // void set_BeamXRMSTracker(vector<double> beamXRMS_ArmF_Trk,vector<double> beamXRMS_ArmB_Trk){
        //     fBeamXRMS_ArmF_Trk = beamXRMS_ArmF_Trk;
        //     fBeamXRMS_ArmB_Trk = beamXRMS_ArmB_Trk;
        // };

        void set_TrackImpactParameterCut(double rip){fTrackImpactParameterCut=rip;}
        void set_ThetaXRangeatDet1(double thx_min,double thx_max){fMinThetaXatDet1=thx_min;fMaxThetaXatDet1=thx_max;}
        void set_ThetaYRangeatDet1(double thy_min,double thy_max){fMinThetaYatDet1=thy_min;fMaxThetaYatDet1=thy_max;}

        void ReadGenEvent(const std::vector<TrackingParticle>*);

        void ReadGenEvent(const HepMC::GenEvent* ); // read data from HepMC::GenEvent
        void ReadGenEvent(const std::vector<reco::GenParticle>* ); // read data from reco::GenParticleCollection
        void set_GenData();   // to be used when the generation is done by external generator
        std::map<string,TH2F*> GenBeamProfile(int  direction);
        void Generation();
        void Simulation();
        void Reconstruction();
        bool SearchTrack(TGraphErrors *XLineProjection , TGraphErrors *YLineProjection ,int Direction,double& xi,double& t,double& partP,double& pt,double& thx,double& thy,double& x0,double& y0, double &xChiSquare, double &yChiSquare);
        void TrackReco(int Direction,H_RecRPObject* station,PPSBaseData* arm);
        void VertexReco();
        void Digitization();
        void TrackerDigi();
        void ToFDigi(int Direction, const PPSBaseData*,PPSToFDetector*);
        void ReconstructArm(H_RecRPObject* pps_station, double x1,double y1,double x2,double y2, double& tx, double& ty,double& eloss);
        void Get_t_and_xi(const TLorentzVector* p,double& t, double& xi);

        // void Propagate(H_BeamParticle* p1,int Direction) ; // returns true if particle has stopped before detector position
        void GenSingleParticle(double& , double& ,double&);
        void GenCentralMass(double& , double& ,double&,double&,double&,double&);

        void Run();
        void SmearVertexPosition(double& ,double& ,double&);
        void CrossingAngleCorrection(H_BeamParticle& p_out, const int Direction);
        void CrossingAngleCorrection(TLorentzVector& );
        void ApplyBeamSmearing(TLorentzVector&);
        void LorentzBoost(TLorentzVector& p_out, const string& frame);
        TLorentzVector shoot(const double& t, const double& xi,const double& phi,const int);
        void LoadParameters();
        void PrintParameters();
        // void HitSmearing(double& x, double& y, double& z);
        void ToFSmearing(double& t) {if (fSmearHit) t = gRandom3->Gaus(t,fTimeSigma);};
        void set_XTrackChiSquareCut(double xTrackChiSquareCut) {fXTrackChiSquareCut=xTrackChiSquareCut;};
        void set_YTrackChiSquareCut(double yTrackChiSquareCut) {fYTrackChiSquareCut=yTrackChiSquareCut;};
        double Minimum_t(const double& xi);

        bool isPhysical(const double& xi)    { return (Minimum_t(xi) < t_max)&&xi<=1.&&xi>=0; }
        PPSSpectrometer<Gen> * get_GenDataHolder() {return fGen;};
        PPSSpectrometer<Sim> * get_SimDataHolder() {return fSim;};
        PPSSpectrometer<Reco> * get_RecoDataHolder(){return fReco;};
        edm::DetSetVector<TotemRPDigi> * get_TotemDigi() {return fTotemDigi;};

        double get_BeamMomentum()            {return fBeamMomentum;};
        double get_BeamEnergy()              {return fBeamEnergy;};
        TH2F* get_Beam1Profile()             {return beam1profile;};
        TH2F* get_Beam2Profile()             {return beam2profile;};
        
    private:
        bool         fExternalGenerator;
        bool         fVerbose;
        int          NEvent;
        std::string  fGenMode;
        TRandom3 *gRandom3;
        TH2F*            beam1profile;
        TH2F*            beam2profile;
        PPSSpectrometer<Gen>* fGen;
        PPSSpectrometer<Sim>* fSim;
        PPSSpectrometer<Reco>* fReco;
        // Hector objects
        H_BeamLine*    beamlineF;
        H_BeamLine*    beamlineB;
        H_RecRPObject* pps_stationF;
        H_RecRPObject* pps_stationB;
        // LHC and det parameters
        // std::string         fBeamLine1File;
        // std::string         fBeamLine2File;
        string          fBeam1FilePath;
        string          fBeam2FilePath;
        vector<string>  fBeam1Files;
        vector<string>  fBeam2Files;
        vector<double>  fBeam1FilesXi;
        vector<double>  fBeam2FilesXi;
        vector<string>  fBeam1Objects;
        vector<string>  fBeam2Objects;

        // beam map : < xiClassNumber, map < stationName, pair <beamPositionAtElement, TwissMatrix> > >
        // using beamPositionAtElement = (S(m), X(mm), thetaX(mrad), Y(mm), thetaY(mrad)) from matrix files
        std::map< int, std::map<string,std::pair<vector5Def, matrix44Def> > > fBeam1ParameterCollection;
        std::map< int, std::map<string,std::pair<vector5Def, matrix44Def> > > fBeam2ParameterCollection;
        std::map< int, pair<double, double> > fBeam1XiClasses;
        std::map< int, pair<double, double> > fBeam2XiClasses;

        int                 fBeam1Direction;
        int                 fBeam2Direction;
        bool                fShowBeamLine;
        std::string         fCollisionPoint;
        float          fBeamLineLength;
        double         fBeamEnergy;
        double         fBeamMomentum;
        double         fCrossingAngle; // in micro radians
        bool           fCrossAngleCorr;
        bool           fKickersOFF;
        double         fDetectorClosestX;
        double         fMaxXfromBeam;
        double         fMaxYfromBeam;

        //Tracker characteristics
        string         fTrackerGeometry;
        double         fTrackerZPosition;
        double         fTrackerLength;
        double         fTrackerWidth;
        double         fTrackerHeight;
        double         fVerticalShift;
        double         fTrkDetXOffset;
        int            fNumberOfStrips;
        double         fStripPitch;
        double         fCutSideLength;
        string         fClusterSizePlotFile;
        string         fClusterSizePlot;
        int            fNumberOfRows;
        int            fNumberOfColumns;
        vector<int>    fDoubleSizeColumn;
        vector<int>    fDoubleSizeRow;
        double         fPixelPitchX;
        double         fPixelPitchY;
        vector<int>     fTrackingDetectorIDs;
        std::map<int,string> fBeam1TrackingDetectorNameMap;
        std::map<int,string> fBeam2TrackingDetectorNameMap;
        std::map<int,double> fTrackingStationXRotationMap;
        std::map<int,double> fTrackingStationYRotationMap;
        std::map<int,double> fTrackingStationZRotationMap;
        std::string fClusterSizePlotName;
        edm::ESHandle<TotemRPGeometry> fRpGeometry;


        //Time of Flight characteristics
        double         fToFWidth;
        double         fToFHeight;
        string         fToFGeometry;
        double         fToFZPosition;
        string         fBeam1ToFDetectorName;
        string         fBeam2ToFDetectorName;
        
        std::map<string,double> fBeam1StationPositions;
        std::map<string,double> fBeam2StationPositions;
        std::map<string,double> fBeam1CenterAtStation;
        std::map<string,double> fBeam2CenterAtStation;
        std::map<string,double> fBeam1SigmaAtStation;
        std::map<string,double> fBeam2SigmaAtStation;
        std::map<string,double> fBeam1StationSigmaInsertion;
        std::map<string,double> fBeam2StationSigmaInsertion;
        vector<string> fBeam1TCLNames;
        vector<string> fBeam2TCLNames;

        double         fTCL4Position1;
        double         fTCL4Position2;
        double         fTCL5Position1;
        double         fTCL5Position2;
       
        std::map<int,PPSTrkDetector*> fTrackingStationForward;
        std::map<int,PPSTrkDetector*> fTrackingStationBackward;
        // PPSTrkStation* TrkStation_F; // auxiliary object with the tracker geometry
        // PPSTrkStation* TrkStation_B; // auxiliary object with the tracker geometry
        PPSToFDetector* ToFDet_F;  // idem for the ToF detector
        PPSToFDetector* ToFDet_B;  // idem for the ToF detector
        // Parameters for vertex smearing
        bool   fSmearVertex;
        double fVtxMeanX;
        double fVtxMeanY;
        double fVtxMeanZ;
        double fVtxSigmaX;
        double fVtxSigmaY;
        double fVtxSigmaZ;

        // Parameters for hit smearing
        double fSmearHit;
        double fHitSigmaX;
        double fHitSigmaY;
        double fHitSigmaZ;
        double fToFHitSigmaX;
        double fToFHitSigmaY;
        double fTimeSigma;

        // Parameters for the detector missalignment
        double fTrk1XOffsetF; // distance from RP window to traker 1 Arm F sensor sensitive area
        double fTrk2XOffsetF; // distance from RP window to traker 2 Arm F sensor sensitive area
        double fTrk1XOffsetB; // distance from RP window to traker 1 Arm B sensor sensitive area
        double fTrk2XOffsetB; // distance from RP window to traker 2 Arm B sensor sensitive area
        double fToFXOffsetF ; // distance from RP window to ToF Arm F sensor sensitive area
        double fToFXOffsetB ; // distance from RP window to ToF Arm B sensor sensitive area

        // Parameter for time smearing

        // Parameter for angular smearing
        bool   fSmearAngle;
        double fBeamAngleRMS; // in micro radians

        // Parameter for energy smearing
        bool   fSmearEnergy;
        double fBeamEnergyRMS;

        // Parameter for the Reconstruction
        bool   fFilterHitMap;
        bool   fApplyFiducialCuts;
        double fTrackImpactParameterCut; // maximum impact parameter at IP in cm
        double fMinThetaXatDet1;         // minimum thetaX at first tracker detector in urad
        double fMaxThetaXatDet1;         // maximum thetaX at first tracker detector in urad
        double fMinThetaYatDet1;         // minimum thetaY at first tracker detector in urad
        double fMaxThetaYatDet1;         // maximum thetaY at first tracker detector in urad
        double fXTrackChiSquareCut;
        double fYTrackChiSquareCut;
        bool   fUseToFForTracking;

        // Parameters for the simulation
        double xi_min;
        double xi_max;
        double t_min;
        double t_max;
        double fPhiMin;
        double fPhiMax;
        double fEtaMin;
        double fMomentumMin;
        double fCentralMass;
        double fCentralMassErr;
        std::vector<double> CheckPoints;
        // Generated proton
        int NVertex;
        std::map<int,TVector3> fVertex;
        std::map<int,pair<const TLorentzVector*,const TLorentzVector*> > protonsOut;
        std::map<int,pair<bool,bool> >                       fHasStopped;
        bool fHasStoppedF;
        bool fHasStoppedB;
        const TLorentzVector* protonF;
        const TLorentzVector* protonB;
        //
        TLorentzVector Beam1;
        TLorentzVector Beam2;
        // Simulated hits
        bool fSimBeam;
        edm::DetSetVector<TotemRPDigi> * fTotemDigi;
        // private methods
};
#endif