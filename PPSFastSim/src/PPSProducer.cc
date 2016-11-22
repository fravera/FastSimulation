// -*- C++ -*-
//
// Package:    PPSProducer
// Class:      PPSProducer
// 
/**\class PPSProducer PPSProducer.cc FastSimulation/PPSFastSim/src/PPSProducer.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Luiz Martins Mundim Filho,22 1-020,+41227677686,
//         Created:  Thu Apr 25 17:37:14 CEST 2013
// $Id$
// 
//


// system include files
#include <memory>

// user include files

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "HepMC/GenEvent.h"
#include "HepMC/GenVertex.h"
#include "HepMC/GenParticle.h"
//
#include "FastSimulation/PPSFastObjects/interface/PPSSpectrometer.h"
#include "FastSimulation/PPSFastSim/interface/PPSSim.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Geometry/Records/interface/VeryForwardMisalignedGeometryRecord.h"
#include "Geometry/Records/interface/VeryForwardRealGeometryRecord.h"
#include "DataFormats/CTPPSDigi/interface/TotemRPDigi.h"


#include <map>
//
// class declaration
//

class PPSProducer : public edm::EDProducer {
    public:
        explicit PPSProducer(const edm::ParameterSet&);
        ~PPSProducer();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    private:
        virtual void beginJob() ;
        virtual void beginRun();
        // virtual void beginRun(edm::Run&, edm::EventSetup const&);
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endRun(edm::Run&, edm::EventSetup const&);
        virtual void endJob() ;

        virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
        virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

        // ----------member data ---------------------------
        bool    fVerbose;
        bool   fUseHepMCProducer;
        bool fBeginRun;
        PPSSim* pps;
        //Including GetToken 
        edm::InputTag gensrc;
        edm::EDGetTokenT<edm::HepMCProduct> gensrcTokenHepMC;

        edm::InputTag gensrcPart;
        edm::EDGetTokenT<std::vector<reco::GenParticle>> gensrcPartToken;		
};

//--------------------------------------------------------------------------------------------------------------//

PPSProducer::PPSProducer(const edm::ParameterSet& iConfig):fVerbose(false)
{
    //register your products
    produces<PPSSpectrometer<Sim> >("PPSSim");
    produces<PPSSpectrometer<Gen> >("PPSGen");
    produces<PPSSpectrometer<Reco> >("PPSReco");
    produces<edm::DetSetVector<TotemRPDigi> >("TotemRPDigi");
    
    //now do what ever other initialization is needed
    pps = NULL;
    // edm::FileInPath beam1filename_  = iConfig.getParameter<edm::FileInPath> ("Beam1File");
    // edm::FileInPath beam2filename_  = iConfig.getParameter<edm::FileInPath> ("Beam2File");

    string fBeam1FilePath             = iConfig.getParameter<string>("Beam1FilePath");
    std::vector<string> fBeam1Files   = iConfig.getParameter<vector<string> >("Beam1Files");
    std::vector<double> fBeam1FilesXi = iConfig.getParameter<vector<double> >("Beam1FilesXi");
    std::vector<string> fBeam1Objects = iConfig.getParameter<vector<string> >("Beam1Objects");
    string fBeam2FilePath             = iConfig.getParameter<string>("Beam2FilePath");
    std::vector<string> fBeam2Files   = iConfig.getParameter<vector<string> >("Beam2Files");
    std::vector<double> fBeam2FilesXi = iConfig.getParameter<vector<double> >("Beam2FilesXi");
    std::vector<string> fBeam2Objects = iConfig.getParameter<vector<string> >("Beam2Objects");

    std::vector<double> fBeam1CenterAtStation       = iConfig.getParameter<vector<double> >("Beam1CenterAtStation");
    std::vector<double> fBeam1SigmaAtStation        = iConfig.getParameter<vector<double> >("Beam1SigmaAtStation");
    std::vector<double> fBeam1StationPositions      = iConfig.getParameter<vector<double> >("Beam1StationPositions");
    std::vector<double> fBeam1StationSigmaInsertion = iConfig.getParameter<vector<double> >("Beam1StationSigmaInsertion");
    std::vector<double> fBeam2CenterAtStation       = iConfig.getParameter<vector<double> >("Beam2CenterAtStation");
    std::vector<double> fBeam2SigmaAtStation        = iConfig.getParameter<vector<double> >("Beam2SigmaAtStation");
    std::vector<double> fBeam2StationPositions      = iConfig.getParameter<vector<double> >("Beam2StationPositions");
    std::vector<double> fBeam2StationSigmaInsertion = iConfig.getParameter<vector<double> >("Beam2StationSigmaInsertion");


    int    beam1dir            = iConfig.getParameter<int>("Beam1Direction");
    int    beam2dir            = iConfig.getParameter<int>("Beam2Direction");
    bool   showbeam            = iConfig.getUntrackedParameter<bool>("ShowBeamLine",false);
    bool   simbeam             = iConfig.getUntrackedParameter<bool>("SimBeamProfile",false);
    double fVtxMeanX           = iConfig.getParameter<double>("VtxMeanX");
    double fVtxMeanY           = iConfig.getParameter<double>("VtxMeanY");
    double fVtxMeanZ           = iConfig.getParameter<double>("VtxMeanZ");
    string ip                  = iConfig.getParameter<string>("CollisionPoint");
    fVerbose                   = iConfig.getUntrackedParameter<int>("Verbosity",0);
    // double fTrackerLength      = iConfig.getParameter<double>("TrackerLength");
    // double fTrackerWidth       = iConfig.getParameter<double>("TrackerWidth"); // tracker width in mm
    // double fTrackerHeight      = iConfig.getParameter<double>("TrackerHeight"); // tracker height in mm
    double fToFWidth           = iConfig.getParameter<double>("ToFWidth"); // Tof width in mm
    double fToFHeight          = iConfig.getParameter<double>("ToFHeight"); // Tof height in mm
    string tofgeometry         = iConfig.getParameter<string>("ToFGeometry"); 
    // double fTrk1XOffset        = iConfig.getParameter<double>("TrkDet1XOffset");
    // double fTrk2XOffset        = iConfig.getParameter<double>("TrkDet2XOffset");
    //double fTrkXOffset         = iConfig.getParameter<double>("TrkDetXOffset");
    double fToFXOffset         = iConfig.getParameter<double>("ToFDetXOffset");
    string fBeam1ToFDetectorName    = iConfig.getParameter<string>("Beam1ToFDetectorName");
    string fBeam2ToFDetectorName    = iConfig.getParameter<string>("Beam2ToFDetectorName");
    // vector<double> fBeam1TCLPositions  = iConfig.getUntrackedParameter<vector<double> >("Beam1TCLPositions");
    // vector<double> fBeam2TCLPositions  = iConfig.getUntrackedParameter<vector<double> >("Beam2TCLPositions");
    // vector<double> fBeam1CenterAtTCL   = iConfig.getUntrackedParameter<vector<double> >("fBeam1CenterAtTCL");
    // vector<double> fBeam1SigmaAtTCL   = iConfig.getUntrackedParameter<vector<double> >("fBeam1SigmaAtTCL");
    // vector<double> fBeam1TCLSigmaInsertion   = iConfig.getUntrackedParameter<vector<double> >("Beam1TCLSigmaInsertion");
    // vector<double> fBeam2CenterAtTCL   = iConfig.getUntrackedParameter<vector<double> >("fBeam2CenterAtTCL");
    // vector<double> fBeam2SigmaAtTCL   = iConfig.getUntrackedParameter<vector<double> >("fBeam2SigmaAtTCL");
    // vector<double> fBeam2TCLSigmaInsertion   = iConfig.getUntrackedParameter<vector<double> >("Beam2TCLSigmaInsertion");
    
    vector<string> fBeam1TCLNames           = iConfig.getParameter<vector<string> >("Beam1TCLNames");
    vector<string> fBeam2TCLNames           = iConfig.getParameter<vector<string> >("Beam2TCLNames");

    // double fTCL4Position       = iConfig.getUntrackedParameter<double>("TCL4Position",0.);
    // double fTCL5Position       = iConfig.getUntrackedParameter<double>("TCL5Position",0.);
    bool   fSmearHit           = iConfig.getParameter<bool>("SmearHit");
    double fHitSigmaX          = iConfig.getParameter<double>("HitSigmaX");
    double fHitSigmaY          = iConfig.getParameter<double>("HitSigmaY");
    double fTimeSigma          = iConfig.getParameter<double>("TimeSigma");
    double fToFHitSigmaX       = iConfig.getParameter<double>("ToFHitSigmaX");
    double fToFHitSigmaY       = iConfig.getParameter<double>("ToFHitSigmaY");
    bool   fSmearAngle         = iConfig.getParameter<bool>("SmearAngle");
    double fBeamAngleRMS       = iConfig.getParameter<double>("BeamAngleRMS"); // beam angular dispersion in urad
    bool   fSmearEnergy        = iConfig.getParameter<bool>("SmearEnergy");
    double fBeamEnergy         = iConfig.getParameter<double>("BeamEnergy");    // beam energy in GeV
    double fBeamEnergyRMS      = iConfig.getParameter<double>("BeamEnergyRMS");    // beam energy dispersion in GeV
    // vector<double> fBeam1CenterAtStation = iConfig.getParameter<vector<double> >("Beam1CenterAtStation"); // beam sigma(X) at Arm Forward tracker stations in mm
    // vector<double> fBeam1SigmaAtStation  = iConfig.getParameter<vector<double> >("Beam1SigmaAtStation"); // beam sigma(X) at Arm Forward tracker stations in mm
    // vector<double> fBeam2CenterAtStation = iConfig.getParameter<vector<double> >("Beam2CenterAtStation"); // beam sigma(X) at Arm Forward tracker stations in mm
    // vector<double> fBeam2SigmaAtStation  = iConfig.getParameter<vector<double> >("Beam2SigmaAtStation"); // beam sigma(X) at Arm Forward tracker stations in mm
    
    double fPhiMin             = iConfig.getParameter<double>("PhiMin");
    double fPhiMax             = iConfig.getParameter<double>("PhiMax");
    double fCentralMass        = iConfig.getParameter<double>("CentralMass");
    double fCentralMassErr     = iConfig.getParameter<double>("CentralMassErr");
    bool   fKickersOFF         = iConfig.getUntrackedParameter<bool>("KickersOFF",true);
    bool   fCrossAngleCorr     = iConfig.getParameter<bool>("CrossAngleCorr");
    double fCrossingAngle      = iConfig.getParameter<double>("CrossingAngle");
    double fEtaMin             = iConfig.getParameter<double>("EtaMin");
    double fMomentumMin        = iConfig.getParameter<double>("MomentumMin");
    double fImpParcut          = iConfig.getParameter<double>("TrackImpactParameterCut"); // exclude hit combination that lead to high imp. par. reco tracks (in cm)
    double fMinthx             = iConfig.getParameter<double>("MinThetaXatDet1"); // minimum thetaX at first tracker detector (in urad)
    double fMaxthx             = iConfig.getParameter<double>("MaxThetaXatDet1"); // maximum thetaX at first tracker detector (in urad)
    double fMinthy             = iConfig.getParameter<double>("MinThetaYatDet1"); // minimum thetaY at first tracker detector (in urad)
    double fMaxthy             = iConfig.getParameter<double>("MaxThetaYatDet1"); // maximum thetaY at first tracker detector (in urad)
    // double fMaxXfromBeam       = iConfig.getParameter<double>("MaxXfromBeam");    // maximum distance (X) from beam a hit is accepted (in mm, negative)
    // double fMaxYfromBeam       = iConfig.getParameter<double>("MaxYfromBeam");    // maximum distance (Y) from beam a hit is accepted (in mm, positive, simetric)
    // double fDetectorClosestX   = iConfig.getParameter<double>("DetectorClosestX");// minimum distance (X) from beam a hit is accepted (in mm, negative)
    // bool   fFilterHitMap       = iConfig.getParameter<bool>("FilterHitMap");       // apply geometrical cuts in the hit position (RP window+distance from beam)
    double   fXTrackChiSquareCut = iConfig.getParameter<double>("XTrackChiSquareCut");
    double   fYTrackChiSquareCut = iConfig.getParameter<double>("YTrackChiSquareCut");
    bool   fApplyFiducialCuts  = iConfig.getParameter<bool>("ApplyFiducialCuts");  // apply geometrical cuts in the hit position (Detector size)
    bool   fUseToFForTracking  = iConfig.getParameter<bool>("UseToFForTracking"); 

    vector<string> fBeam1TrkDetNames       = iConfig.getParameter<vector<string> >("Beam1TrkDetNames");
    vector<string> fBeam2TrkDetNames       = iConfig.getParameter<vector<string> >("Beam2TrkDetNames");
    vector<int>    fTrkDetIds         = iConfig.getParameter<vector<int> >("TrkDetIds");
    vector<double> fTrkDetXRotation   = iConfig.getParameter<vector<double> >("TrkDetXRotation");
    vector<double> fTrkDetYRotation   = iConfig.getParameter<vector<double> >("TrkDetYRotation");
    vector<double> fTrkDetZRotation   = iConfig.getParameter<vector<double> >("TrkDetZRotation");
    
    string fTrackerGeometry = iConfig.getParameter<string>("TrackerGeometryType");
    double fVerticalShift   = iConfig.getParameter<double>("VerticalShift");
    double fTrkDetXOffset   = iConfig.getParameter<double>("TrkDetXOffset");

    int fNumberOfStrips     = -1;
    double fStripPitch      = -1.;
    double fCutSideLength   = -1.;
    string fClusterSizePlotFile = "";
    string fClusterSizePlotName = "";
    int fNumberOfRows       = -1;
    int fNumberOfColumns    = -1;
    vector<int> fDoubleSizeColumn;
    vector<int> fDoubleSizeRow;
    double fPixelPitchX     = -1.;  // pixel pitch in X in microns
    double fPixelPitchY     = -1.;  // pixel pitch in Y in microns
    
    if(fTrackerGeometry=="TOTEMStrip"){
        fNumberOfStrips = iConfig.getParameter<int>("NumberOfStrips");
        fStripPitch     = iConfig.getParameter<double>("StripPitch");
        fCutSideLength  = iConfig.getParameter<double>("CutSideLength");
        fClusterSizePlotFile  = iConfig.getParameter<string>("ClusterSizePlotFile");
        fClusterSizePlotName  = iConfig.getParameter<string>("ClusterSizePlotName");
    }

    if(fTrackerGeometry=="Pixel"){
        fNumberOfRows             = iConfig.getParameter<int>("NumberOfRows");
        fNumberOfColumns          = iConfig.getParameter<int>("NumberOfColumns");
        fDoubleSizeRow    = iConfig.getParameter<vector<int> >("DoubleSizeRow");
        fDoubleSizeColumn = iConfig.getParameter<vector<int> >("DoubleSizeColumn");
        fPixelPitchX           = iConfig.getParameter<double>("PixelPitchX");
        fPixelPitchY           = iConfig.getParameter<double>("PixelPitchY");
    }

    gensrc              = iConfig.getParameter<edm::InputTag>("genSource");
    gensrcTokenHepMC = consumes<edm::HepMCProduct>(gensrc);

    gensrcPart           = iConfig.getParameter<edm::InputTag>("genSource");
    gensrcPartToken = consumes<std::vector<reco::GenParticle>>(gensrcPart);

    // New Flag to use the HepMCProducer (True)  or GenParticle + PU Particles (false)
    fUseHepMCProducer         = iConfig.getUntrackedParameter<bool>("UseHepMCProducer",false);

    // To get Path for Optics files
    // std::string beam1filename = beam1filename_.fullPath();
    // std::string beam2filename = beam2filename_.fullPath();
    // if (beam1filename.find(".") == 0 && beam2filename.find(".") ==0  ){
    //     beam1filename.erase(0, 2);
    //     beam2filename.erase(0, 2);		}

    pps = new PPSSim(true); // instanciate PPSSim with External Generator
    pps->set_KickersOFF(fKickersOFF);
    // pps->set_BeamLineFile(beam1filename,beam2filename);
    pps->set_BeamLineFilePath(fBeam1FilePath,fBeam2FilePath);
    pps->set_BeamLineFiles(fBeam1Files,fBeam2Files);
    pps->set_BeamLineFilesXi(fBeam1FilesXi,fBeam2FilesXi);
    pps->set_BeamLineObjects(fBeam1Objects,fBeam2Objects);

    pps->set_BeamDirection(beam1dir,beam2dir);
    pps->set_BeamEnergySmearing(fSmearEnergy);
    pps->set_BeamEnergy(fBeamEnergy);
    pps->set_BeamEnergyRMS(fBeamEnergyRMS);
    pps->set_BeamAngleSmearing(fSmearAngle);
    pps->set_BeamAngleRMS(fBeamAngleRMS);

    map<string, double> fBeam1StationPositionMap;
    map<string, double> fBeam2StationPositionMap;
    map<string, double> fBeam1CenterAtStationMap;
    map<string, double> fBeam2CenterAtStationMap;
    map<string, double> fBeam1SigmaAtStationMap;
    map<string, double> fBeam2SigmaAtStationMap;
    map<string, double> fBeam1StationSigmaInsertionMap;
    map<string, double> fBeam2StationSigmaInsertionMap;

    for(unsigned i=0; i<fBeam1Objects.size(); ++i){
        if (fBeam1StationPositionMap.find(fBeam1Objects.at(i)) != fBeam1StationPositionMap.end()){
            throw cms::Exception("PPSProducer.cc") << "Station name "<<fBeam1Objects.at(i)<<" is repeeted";
        }
        fBeam1StationPositionMap[fBeam1Objects.at(i)]    = fBeam1StationPositions.at(i);
        fBeam1CenterAtStationMap[fBeam1Objects.at(i)]    = fBeam1CenterAtStation.at(i);
        fBeam1SigmaAtStationMap[fBeam1Objects.at(i)]     = fBeam1SigmaAtStation.at(i);
        fBeam1StationSigmaInsertionMap[fBeam1Objects.at(i)] = fBeam1StationSigmaInsertion.at(i);
    }

    for(unsigned i=0; i<fBeam2Objects.size(); ++i){
        if (fBeam2StationPositionMap.find(fBeam2Objects.at(i)) != fBeam2StationPositionMap.end()){
            throw cms::Exception("PPSProducer.cc") << "Station name "<<fBeam2Objects.at(i)<<" is repeeted";
        }
        fBeam2StationPositionMap[fBeam2Objects.at(i)]    = fBeam2StationPositions.at(i);
        fBeam2CenterAtStationMap[fBeam2Objects.at(i)]    = fBeam2CenterAtStation.at(i);
        fBeam2SigmaAtStationMap[fBeam2Objects.at(i)]     = fBeam2SigmaAtStation.at(i);
        fBeam2StationSigmaInsertionMap[fBeam2Objects.at(i)] = fBeam2StationSigmaInsertion.at(i);
    }

    pps->set_Beam1StationPositions(fBeam1StationPositionMap);
    pps->set_Beam2StationPositions(fBeam2StationPositionMap);
    pps->set_Beam1CenterAtStation(fBeam1CenterAtStationMap);
    pps->set_Beam1SigmaAtStation(fBeam1SigmaAtStationMap);
    pps->set_Beam1StationSigmaInsertion(fBeam1StationSigmaInsertionMap);
    pps->set_Beam2CenterAtStation(fBeam2CenterAtStationMap);
    pps->set_Beam2SigmaAtStation(fBeam2SigmaAtStationMap);
    pps->set_Beam2StationSigmaInsertion(fBeam2StationSigmaInsertionMap);

    pps->set_TCLNames(fBeam1TCLNames,fBeam2TCLNames);

    if (showbeam) pps->set_ShowBeamLine();
    if (simbeam)  pps->set_GenBeamProfile();
    pps->set_VertexPosition(fVtxMeanX,fVtxMeanY,fVtxMeanZ);
    pps->set_CollisionPoint(ip);

    //Tracker characteristics
    if(fBeam1TrkDetNames.size() != fTrkDetIds.size()){
        edm::LogError("debug") << "Number of Tracking detector names of beam 1 and tracking detector Z positions are not equal";
        exit(1);
    }
    if(fBeam2TrkDetNames.size() != fTrkDetIds.size()){
        edm::LogError("debug") << "Number of Tracking detector names of beam 2 and tracking detector Z positions are not equal";
        exit(1);
    }

    map<int, string> fBeam1TrackingDetectorNamesMap;
    map<int, string> fBeam2TrackingDetectorNamesMap;
    map<int, double> fTrackingDetectorXRotationMap;
    map<int, double> fTrackingDetectorYRotationMap;
    map<int, double> fTrackingDetectorZRotationMap;
    map<int,double>  fBeamXRMSArmFTrackerMap;
    map<int,double>  fBeamXRMSArmBTrackerMap;
    map<int,double>  fBeamXCenterArmFTrackerMap;
    map<int,double>  fBeamXCenterArmBTrackerMap;

    for(unsigned i=0; i<fTrkDetIds.size(); ++i){
        if (fBeam1TrackingDetectorNamesMap.find(fTrkDetIds.at(i)) != fBeam1TrackingDetectorNamesMap.end()){
            edm::LogError("debug") << "Tracking station name "<<fTrkDetIds.at(i)<<" is repeeted";
            edm::LogError("debug") <<"Two station must not have identical names";
            exit(1);
        }
        fBeam1TrackingDetectorNamesMap[fTrkDetIds.at(i)] = fBeam1TrkDetNames.at(i);
        fBeam2TrackingDetectorNamesMap[fTrkDetIds.at(i)] = fBeam2TrkDetNames.at(i);
        fTrackingDetectorXRotationMap[fTrkDetIds.at(i)] = fTrkDetXRotation.at(i);
        fTrackingDetectorYRotationMap[fTrkDetIds.at(i)] = fTrkDetYRotation.at(i);
        fTrackingDetectorZRotationMap[fTrkDetIds.at(i)] = fTrkDetZRotation.at(i);
    }

    pps->set_TrackingDetectorIDs(fTrkDetIds);
    pps->set_TrackingDetectorNameMap(fBeam1TrackingDetectorNamesMap,fBeam2TrackingDetectorNamesMap);
    pps->set_TrackingStationXRotationMap(fTrackingDetectorXRotationMap);
    pps->set_TrackingStationYRotationMap(fTrackingDetectorYRotationMap);
    pps->set_TrackingStationZRotationMap(fTrackingDetectorZRotationMap);
    
    pps->set_TrackerGeometry(fTrackerGeometry);
    pps->set_TrackerEdgeOffset(fTrkDetXOffset,fTrkDetXOffset,fTrkDetXOffset,fTrkDetXOffset); // use the same offset for the whole tracker
    pps->set_VerticalShift(fVerticalShift);
    pps->set_TrkDetXOffset(fTrkDetXOffset);
    pps->set_NumberOfStrips(fNumberOfStrips);
    pps->set_StripPitch(fStripPitch);
    pps->set_CutSideLength(fCutSideLength);
    pps->set_ClusterSizePlotFile(fClusterSizePlotFile);
    pps->set_ClusterSizePlotName(fClusterSizePlotName);
    pps->set_NumberOfRows(fNumberOfRows);
    pps->set_NumberOfColumns(fNumberOfColumns);
    pps->set_DoubleSizeColumn(fDoubleSizeColumn);
    pps->set_DoubleSizeRow(fDoubleSizeRow);
    pps->set_PixelPitchX(fPixelPitchX);
    pps->set_PixelPitchY(fPixelPitchY);
  

    //Time of Flight characteristics
    pps->set_ToFGeometry(tofgeometry);
    pps->set_ToFSize(fToFWidth,fToFHeight);
    pps->set_ToFDetectorName(fBeam1ToFDetectorName,fBeam2ToFDetectorName);
    pps->set_ToFResolution(fTimeSigma,fToFHitSigmaX,fToFHitSigmaY);
    pps->set_ToFEdgeOffset    (fToFXOffset,fToFXOffset); // use the same offset for the forward and backward arm
    
    pps->set_HitSmearing(fSmearHit);
    pps->set_VertexSmearing(false); // when using cmssw, vertex smearing is done somewhere else
    pps->set_phiMin(fPhiMin);
    pps->set_phiMax(fPhiMax);
    pps->set_etaMin(fEtaMin);
    pps->set_momentumMin(fMomentumMin);
    pps->set_CentralMass(fCentralMass,fCentralMassErr);
    // pps->set_HitSmearing(fSmearHit);
    pps->set_TrackerResolution(fHitSigmaX,fHitSigmaY);
    pps->set_Verbose(fVerbose);
    pps->set_CrossingAngleCorrection(fCrossAngleCorr);
    pps->set_CrossingAngle(fCrossingAngle);
    pps->set_TrackImpactParameterCut(fImpParcut);
    pps->set_ThetaXRangeatDet1(fMinthx,fMaxthx);
    pps->set_ThetaYRangeatDet1(fMinthy,fMaxthy);
    // pps->set_WindowForTrack(fMaxXfromBeam,fMaxYfromBeam,fDetectorClosestX);
    // pps->set_FilterHitMap(fFilterHitMap);
    pps->set_XTrackChiSquareCut(fXTrackChiSquareCut);
    pps->set_YTrackChiSquareCut(fYTrackChiSquareCut);
    pps->set_ApplyFiducialCuts(fApplyFiducialCuts);
    pps->set_UseToFForTracking(fUseToFForTracking);

}


//--------------------------------------------------------------------------------------------------------------//

//
// member functions
//

// ------------ method called to produce the data  ------------
void PPSProducer::produce(edm::Event& iEvent, const edm::EventSetup& eventSetup)
{

    std::cout<<"!!!----- PPSProducer::produce -----!!!"<<std::endl;

    using namespace edm;
    using namespace reco;
    using namespace HepMC;
    using namespace CLHEP;

    if (pps && fBeginRun){
    edm::ESHandle<TotemRPGeometry> rpGeometry;
    eventSetup.get<VeryForwardMisalignedGeometryRecord>().get(rpGeometry);
    //std::cout<<5.1<<std::endl;
    //rpGeometry->GetDetector(1020);
    // std::cout<<"Reading RP Geometry"<<std::endl;
        pps->set_RomanPotGeometry(rpGeometry);

        pps->PrintParameters();
        pps->BeginRun();
        pps->PrintParameters();
        // fBeginRun = false;
    }
   
    pps->BeginEvent();
    if(fUseHepMCProducer){
        try{
            Handle<HepMCProduct> genHandle;
            iEvent.getByToken(gensrcTokenHepMC,genHandle);
            const HepMC::GenEvent* evt = genHandle->GetEvent();
            pps->ReadGenEvent(evt);
        }
        catch(const Exception&){
            edm::LogWarning("debug") <<"PPSProducer::produce: No MC information";
            return;
        }
    } else {
        try{
            // Try first genParticles with PileUp events
            Handle<std::vector<reco::GenParticle> > genParticles;
            iEvent.getByToken(gensrcPartToken, genParticles);
            pps->ReadGenEvent(&(*genParticles));
        }
        catch(const Exception&){
            edm::LogWarning("debug") <<"PPSProducer::produce: No GenParticle information. Trying HepMC::Event.";
            return;
        }
    }
    pps->Run();
    pps->EndEvent();
    std::unique_ptr<PPSSpectrometer<Gen> > fGen(new PPSSpectrometer<Gen>(*(pps->get_GenDataHolder())));
    std::unique_ptr<PPSSpectrometer<Sim> > fSim(new PPSSpectrometer<Sim>(*(pps->get_SimDataHolder())));
    std::unique_ptr<PPSSpectrometer<Reco> > fReco(new PPSSpectrometer<Reco>(*(pps->get_RecoDataHolder())));
    std::unique_ptr<edm::DetSetVector<TotemRPDigi> > fTotemDigi(new edm::DetSetVector<TotemRPDigi>(*(pps->get_TotemDigi())));
    // iEvent.put(fGen);
    // iEvent.put(fSim);
    // iEvent.put(fReco);
    iEvent.put(std::move(fGen),"PPSGen");
    iEvent.put(std::move(fSim),"PPSSim");
    iEvent.put(std::move(fReco),"PPSReco");
    iEvent.put(std::move(fTotemDigi),"TotemRPDigi");
}

//--------------------------------------------------------------------------------------------------------------//

PPSProducer::~PPSProducer()
{

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)

}

//--------------------------------------------------------------------------------------------------------------//

// ------------ method called once each job just before starting event loop  ------------
void PPSProducer::beginJob()
// void PPSProducer::beginJob(edm::EventSetup const& eventSetup)
{
    std::cout<<"!!!----- PPSProducer::beginJob -----!!!"<<std::endl;
    fBeginRun = true;
    pps->LoadTwissMatrices();
    // edm::ESHandle<TotemRPGeometry> rpGeometry;
    // eventSetup.get<VeryForwardMisalignedGeometryRecord>().get(rpGeometry);
    // //std::cout<<5.1<<std::endl;
    // //rpGeometry->GetDetector(1020);
    // std::cout<<"Cout reading RP Geometry"<<std::endl;
    // pps->set_RomanPotGeometry(rpGeometry);

    // if (pps){
    //     pps->PrintParameters();
    //     pps->BeginRun();
    //     pps->PrintParameters();
    // }
   
}

//--------------------------------------------------------------------------------------------------------------//

// ------------ method called once each job just after ending the event loop  ------------
void PPSProducer::endJob() {

    std::cout<<"!!!----- PPSProducer::endJob -----!!!"<<std::endl;
    if(pps) pps->EndRun();
}

//--------------------------------------------------------------------------------------------------------------//

// ------------ method called when starting to processes a run  ------------
// void PPSProducer::beginRun(edm::Run&, edm::EventSetup const& eventSetup)
void PPSProducer::beginRun()
{
    std::cout<<"!!!----- PPSProducer::beginRun -----!!!"<<std::endl;


}

//--------------------------------------------------------------------------------------------------------------//

// ------------ method called when ending the processing of a run  ------------
void PPSProducer::endRun(edm::Run&, edm::EventSetup const&)
{
    std::cout<<"!!!----- PPSProducer::endRun -----!!!"<<std::endl;

}

//--------------------------------------------------------------------------------------------------------------//

// ------------ method called when starting to processes a luminosity block  ------------
void PPSProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
    std::cout<<"!!!----- PPSProducer::beginLuminosityBlock -----!!!"<<std::endl;

}

//--------------------------------------------------------------------------------------------------------------//

// ------------ method called when ending the processing of a luminosity block  ------------
void PPSProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
    std::cout<<"!!!----- PPSProducer::endLuminosityBlock -----!!!"<<std::endl;

}

//--------------------------------------------------------------------------------------------------------------//

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void PPSProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//--------------------------------------------------------------------------------------------------------------//

//define this as a plug-in
DEFINE_FWK_MODULE(PPSProducer);
