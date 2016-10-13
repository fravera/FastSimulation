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
    //now do what ever other initialization is needed
    pps = NULL;
    edm::FileInPath beam1filename_  = iConfig.getParameter<edm::FileInPath> ("Beam1File");
    edm::FileInPath beam2filename_  = iConfig.getParameter<edm::FileInPath> ("Beam2File");
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
    double fToFInsertion       = iConfig.getParameter<double>("ToFInsertion");     // Number of sigms (X) from the beam for the tof
    string tofgeometry         = iConfig.getParameter<string>("ToFGeometry"); 
    // double fTrk1XOffset        = iConfig.getParameter<double>("TrkDet1XOffset");
    // double fTrk2XOffset        = iConfig.getParameter<double>("TrkDet2XOffset");
    //double fTrkXOffset         = iConfig.getParameter<double>("TrkDetXOffset");
    double fToFXOffset         = iConfig.getParameter<double>("ToFDetXOffset");
    double fToFZPosition       = iConfig.getParameter<double>("ToFZPosition");
    double fTCL4Position       = iConfig.getUntrackedParameter<double>("TCL4Position",0.);
    double fTCL5Position       = iConfig.getUntrackedParameter<double>("TCL5Position",0.);
    bool   fSmearHit           = iConfig.getParameter<bool>("SmearHit");
    // double fHitSigmaX          = iConfig.getParameter<double>("HitSigmaX");
    // double fHitSigmaY          = iConfig.getParameter<double>("HitSigmaY");
    double fTimeSigma          = iConfig.getParameter<double>("TimeSigma");
    double fToFHitSigmaX       = iConfig.getParameter<double>("ToFHitSigmaX");
    double fToFHitSigmaY       = iConfig.getParameter<double>("ToFHitSigmaY");
    bool   fSmearAngle         = iConfig.getParameter<bool>("SmearAngle");
    double fBeamAngleRMS       = iConfig.getParameter<double>("BeamAngleRMS"); // beam angular dispersion in urad
    bool   fSmearEnergy        = iConfig.getParameter<bool>("SmearEnergy");
    double fBeamEnergy         = iConfig.getParameter<double>("BeamEnergy");    // beam energy in GeV
    double fBeamEnergyRMS      = iConfig.getParameter<double>("BeamEnergyRMS");    // beam energy dispersion in GeV
    vector<double> fBeamSizes_ArmF_Trk     = iConfig.getParameter<vector<double> >("BeamSizes_ArmF_Trk"); // beam sigma(X) at Arm Forward tracker stations in mm
    double fBeamSize_ArmF_ToF      = iConfig.getParameter<double>("BeamSize_ArmF_ToF" ); // beam sigma(X) at Arm Forward timing station in mm
    vector<double> fBeamSizes_ArmB_Trk     = iConfig.getParameter<vector<double> >("BeamSizes_ArmB_Trk"); // beam sigma(X) at Arm Backward first tracker station in mm
    double fBeamSize_ArmB_ToF      = iConfig.getParameter<double>("BeamSize_ArmB_ToF" ); // beam sigma(X) at Arm Backward timing station in mm
    vector<double> fBeamCenters_ArmF_Trk     = iConfig.getParameter<vector<double> >("BeamCenters_ArmF_Trk"); // beam center(X) at Arm Forward tracker stations in mm
    double fBeamCenter_ArmF_ToF      = iConfig.getParameter<double>("BeamCenter_ArmF_ToF" ); // beam center(X) at Arm Forward timing station in mm
    vector<double> fBeamCenters_ArmB_Trk     = iConfig.getParameter<vector<double> >("BeamCenters_ArmB_Trk"); // beam center(X) at Arm Backward first tracker station in mm
    double fBeamCenter_ArmB_ToF      = iConfig.getParameter<double>("BeamCenter_ArmB_ToF" ); // beam center(X) at Arm Backward timing station in mm
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


    double fTrackerInsertion   = iConfig.getParameter<double>("TrackerInsertion"); // Number of sigms (X) from the beam for the tracker
    vector<string> fTrkDetNames       = iConfig.getParameter<vector<string> >("TrkDetNames");
    vector<int>    fTrkDetIds         = iConfig.getParameter<vector<int> >("TrkDetIds");
    vector<double> fTrkDetZPosition   = iConfig.getParameter<vector<double> >("TrkDetZPosition");
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
    std::string beam1filename = beam1filename_.fullPath();
    std::string beam2filename = beam2filename_.fullPath();
    if (beam1filename.find(".") == 0 && beam2filename.find(".") ==0  ){
        beam1filename.erase(0, 2);
        beam2filename.erase(0, 2);		}

    pps = new PPSSim(true); // instanciate PPSSim with External Generator
    pps->set_KickersOFF(fKickersOFF);
    pps->set_BeamLineFile(beam1filename,beam2filename);
    pps->set_BeamDirection(beam1dir,beam2dir);
    pps->set_BeamEnergySmearing(fSmearEnergy);
    pps->set_BeamEnergy(fBeamEnergy);
    pps->set_BeamEnergyRMS(fBeamEnergyRMS);
    pps->set_BeamAngleSmearing(fSmearAngle);
    pps->set_BeamAngleRMS(fBeamAngleRMS);

    pps->set_BeamXRMSToF(fBeamSize_ArmF_ToF,fBeamSize_ArmB_ToF);
    pps->set_BeamXCenterToF(fBeamCenter_ArmF_ToF,fBeamCenter_ArmB_ToF);
    // pps->set_BeamXSizes(fBeamSize_ArmF_Trk1,fBeamSize_ArmF_Trk2,fBeamSize_ArmF_ToF,fBeamSize_ArmB_Trk1,fBeamSize_ArmB_Trk2,fBeamSize_ArmB_ToF);
    pps->set_TCLPosition("TCL4",fTCL4Position,fTCL4Position);
    pps->set_TCLPosition("TCL5",fTCL5Position,fTCL5Position);
    if (showbeam) pps->set_ShowBeamLine();
    if (simbeam)  pps->set_GenBeamProfile();
    pps->set_VertexPosition(fVtxMeanX,fVtxMeanY,fVtxMeanZ);
    pps->set_CollisionPoint(ip);

    //Tracker characteristics
    if(fTrkDetZPosition.size() != fTrkDetIds.size()){
        edm::LogError("debug") << "Number of Tracking detector names and tracking detector Z positions are not equal";
        exit(1);
    }
    if(fTrkDetNames.size() != fTrkDetIds.size()){
        edm::LogError("debug") << "Number of Tracking detector names and tracking detector Z positions are not equal";
        exit(1);
    }

    map<int, string> fTrackingDetectorNamesMap;
    map<int, double> fTrackingDetectorZPositionMap;
    map<int, double> fTrackingDetectorXRotationMap;
    map<int, double> fTrackingDetectorYRotationMap;
    map<int, double> fTrackingDetectorZRotationMap;
    map<int,double>  fBeamXRMSArmFTrackerMap;
    map<int,double>  fBeamXRMSArmBTrackerMap;
    map<int,double>  fBeamXCenterArmFTrackerMap;
    map<int,double>  fBeamXCenterArmBTrackerMap;

    for(unsigned i=0; i<fTrkDetIds.size(); ++i){
        if (fTrackingDetectorZPositionMap.find(fTrkDetIds.at(i)) != fTrackingDetectorZPositionMap.end()){
            edm::LogError("debug") << "Tracking station name "<<fTrkDetIds.at(i)<<" is repeeted";
            edm::LogError("debug") <<"Two station must not have identical names";
            exit(1);
        }
        fTrackingDetectorNamesMap[fTrkDetIds.at(i)] = fTrkDetNames.at(i);
        fTrackingDetectorZPositionMap[fTrkDetIds.at(i)] = fTrkDetZPosition.at(i);
        fTrackingDetectorXRotationMap[fTrkDetIds.at(i)] = fTrkDetXRotation.at(i);
        fTrackingDetectorYRotationMap[fTrkDetIds.at(i)] = fTrkDetYRotation.at(i);
        fTrackingDetectorZRotationMap[fTrkDetIds.at(i)] = fTrkDetZRotation.at(i);
        fBeamXRMSArmFTrackerMap[fTrkDetIds.at(i)] = fBeamSizes_ArmF_Trk.at(i);
        fBeamXRMSArmBTrackerMap[fTrkDetIds.at(i)] = fBeamSizes_ArmB_Trk.at(i);
        fBeamXCenterArmFTrackerMap[fTrkDetIds.at(i)] = fBeamCenters_ArmF_Trk.at(i);
        fBeamXCenterArmBTrackerMap[fTrkDetIds.at(i)] = fBeamCenters_ArmB_Trk.at(i);

    }

    pps->set_TrackingDetectorIDs(fTrkDetIds);
    pps->set_TrackingDetectorNameMap(fTrackingDetectorNamesMap);
    pps->set_TrackingStationZPositionMap(fTrackingDetectorZPositionMap);
    pps->set_TrackingStationXRotationMap(fTrackingDetectorXRotationMap);
    pps->set_TrackingStationYRotationMap(fTrackingDetectorYRotationMap);
    pps->set_TrackingStationZRotationMap(fTrackingDetectorZRotationMap);
    pps->set_BeamXSizeTracker(fBeamXRMSArmFTrackerMap,fBeamXRMSArmBTrackerMap);
    pps->set_BeamXCenterTracker(fBeamXCenterArmFTrackerMap,fBeamXCenterArmBTrackerMap);

    pps->set_TrackerGeometry(fTrackerGeometry);
    // pps->set_TrackerZPosition(fTrackerZPosition);
    pps->set_TrackerInsertion(fTrackerInsertion);
    // pps->set_TrackerLength(fTrackerLength);
    // pps->set_TrackerSize(fTrackerWidth,fTrackerHeight);
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
    pps->set_ToFInsertion(fToFInsertion);
    pps->set_ToFSize(fToFWidth,fToFHeight);
    pps->set_ToFZPosition(fToFZPosition);
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
    // pps->set_TrackerResolution((fHitSigmaX+fHitSigmaY)/2.*um_to_mm);
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

    edm::ESHandle<TotemRPGeometry> rpGeometry;
    eventSetup.get<VeryForwardMisalignedGeometryRecord>().get(rpGeometry);
    //std::cout<<5.1<<std::endl;
    //rpGeometry->GetDetector(1020);
    std::cout<<"Reading RP Geometry"<<std::endl;
    pps->set_RomanPotGeometry(rpGeometry);

    if (pps){
        pps->PrintParameters();
        pps->BeginRun();
        pps->PrintParameters();
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
    // iEvent.put(fGen);
    // iEvent.put(fSim);
    // iEvent.put(fReco);
    iEvent.put(std::move(fGen),"PPSGen");
    iEvent.put(std::move(fSim),"PPSSim");
    iEvent.put(std::move(fReco),"PPSReco");
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
