// ROOT #includes
#include "FastSimulation/PPSFastSim/interface/PPSSim.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <TMatrixD.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>

//=====================================================================================================

PPSSim::PPSSim(bool ext_gen): fExternalGenerator(ext_gen),
    fVerbose(false),NEvent(0),fGenMode(""),
    fBeam1FilePath(""),fBeam2FilePath(""),fBeam1Files(),fBeam2Files(),fBeam1FilesXi(),fBeam2FilesXi(),fBeam1Objects(),fBeam2Objects(),
    fBeam1Direction(1),fBeam2Direction(1),fShowBeamLine(false),
    fCollisionPoint(""),fBeamLineLength(500),fBeamEnergy(0),fBeamMomentum(0),
    fCrossingAngle(0.),fCrossAngleCorr(false),fKickersOFF(false),
    fDetectorClosestX(-2.),fMaxXfromBeam(-25),fMaxYfromBeam(10), fTrackerGeometry(""),
    fTrackerLength(0.),fTrackerWidth(0.),fTrackerHeight(0.),
    fVerticalShift(0.),fTrkDetXOffset(0.),fNumberOfStrips(-1),
    fStripPitch(0.),fCutSideLength(0.),fNumberOfRows(0),fNumberOfColumns(0),
    fDoubleSizeColumn(),fDoubleSizeRow(),fPixelPitchX(0.),fPixelPitchY(0.),
    fToFWidth(0.),fToFHeight(0.),fToFGeometry(""),
    fSmearVertex(false),fVtxMeanX(0.),fVtxMeanY(0.),fVtxMeanZ(0.),fVtxSigmaX(0.),fVtxSigmaY(0.),fVtxSigmaZ(0.),
    fSmearHit(1.),fHitSigmaX(0.),fHitSigmaY(0.),fHitSigmaZ(0.),fTimeSigma(0.),
    fTrk1XOffsetF(0.),fTrk2XOffsetF(0.),fTrk1XOffsetB(0.),fTrk2XOffsetB(0.),
    fSmearAngle(false),fBeamAngleRMS(0.),fSmearEnergy(false),fBeamEnergyRMS(0.),
    fFilterHitMap(true),fApplyFiducialCuts(true),
    fTrackImpactParameterCut(0.),fMinThetaXatDet1(-200),fMaxThetaXatDet1(200),
    fMinThetaYatDet1(-200),fMaxThetaYatDet1(200),
    fXTrackChiSquareCut(0.),fYTrackChiSquareCut(0.),fUseToFForTracking(false),
    xi_min(0.),xi_max(0.),t_min(0.),t_max(0.),fPhiMin(-TMath::Pi()),fPhiMax(TMath::Pi()),
    fEtaMin(7.),fMomentumMin(3000.),fCentralMass(0.),fCentralMassErr(0.),
    CheckPoints(),fSimBeam(false)
{
    beam1profile=NULL;
    beam2profile=NULL;
    gRandom3  = new TRandom3(0);
}

//--------------------------------------------------------------------------------------------------------------//

void PPSSim::BeginRun()
{    
    // if (fVerbose) edm::LogWarning("debug") << "fBeamLine1File: " << fBeamLine1File ;
    // if (fVerbose) edm::LogWarning("debug") << "fBeamLine2File: " << fBeamLine2File ;    
    extern int kickers_on;
    kickers_on = (fKickersOFF)?0:1;
    beamlineF = new H_BeamLine(-1,fBeamLineLength);
    beamlineB = new H_BeamLine( 1,fBeamLineLength);
    beamlineF->fill(fBeam2FilePath + "/" + fBeam2Files[0],fBeam2Direction,fCollisionPoint);
    beamlineB->fill(fBeam1FilePath + "/" + fBeam1Files[0],fBeam1Direction,fCollisionPoint);
    //
    //set_Strengths();
    beamlineF->offsetElements( 120, 0.097);
    beamlineB->offsetElements( 120, 0.097);
    beamlineF->calcMatrix(); beamlineB->calcMatrix();
    if (fShowBeamLine) {
        beamlineF->showElements();
        beamlineB->showElements();
    }
    // Create a particle to get the beam energy from the beam file
    //
    
    pps_stationF = new H_RecRPObject(fBeam2StationPositions[fBeam2TrackingDetectorNameMap[fTrackingDetectorIDs.at(0)]],fBeam2StationPositions[fBeam2TrackingDetectorNameMap[fTrackingDetectorIDs.at(1)]],*beamlineF);
    pps_stationB = new H_RecRPObject(fBeam1StationPositions[fBeam1TrackingDetectorNameMap[fTrackingDetectorIDs.at(0)]],fBeam1StationPositions[fBeam1TrackingDetectorNameMap[fTrackingDetectorIDs.at(1)]],*beamlineB);
    //
    // check the kinematic limits in case it is requested to generate the higgs mass in the central system
    if (fGenMode=="M_X") { // check the kinematic limits
        if (xi_min*(2.*fBeamEnergy)>fCentralMass+fCentralMassErr||xi_max*(2.*fBeamEnergy)<fCentralMass-fCentralMassErr) {
            edm::LogWarning("debug") << "xi limits outside kinematic limits for the given Central mass. Stopping..." ;
            exit(1);
        }
    }
    if (fSimBeam) {

        TFile beamProfilesFile("BeamProfiles.root","RECREATE");

        std::map<string,TH2F*> beam1Profile = GenBeamProfile(-1);
        std::map<string,TH2F*> beam2Profile = GenBeamProfile( 1);

        for(vector<string>::iterator vIt=fBeam1Objects.begin(); vIt!=fBeam1Objects.end(); ++vIt){
           beam1Profile[*vIt]->Write();
        }

        for(vector<string>::iterator vIt=fBeam2Objects.begin(); vIt!=fBeam2Objects.end(); ++vIt){
            beam2Profile[*vIt]->Write();
        }

        beamProfilesFile.Close();
        fSimBeam = false;
    }
    // 
    fGen = new PPSSpectrometer<Gen>();
    fSim = new PPSSpectrometer<Sim>();
    fReco= new PPSSpectrometer<Reco>();
    //

    fTrackerZPosition = fBeam1StationPositions[fBeam1TrackingDetectorNameMap[fTrackingDetectorIDs.at(0)]];
    fTrackerLength = fBeam1StationPositions[fBeam1TrackingDetectorNameMap[fTrackingDetectorIDs.at(1)]] - fBeam1StationPositions[fBeam1TrackingDetectorNameMap[fTrackingDetectorIDs.at(0)]];;

    for(unsigned i=0; i< fTrackingDetectorIDs.size(); ++i){
        int detectorId = fTrackingDetectorIDs.at(i);
        edm::LogWarning("debug") << "Building detector "<< detectorId;
        string beam1DetectorName = fBeam1TrackingDetectorNameMap[i];
        string beam2DetectorName = fBeam2TrackingDetectorNameMap[i];
        PPSTrkDetector *trackingDetector_ArmForward = NULL;
        PPSTrkDetector *trackingDetector_ArmBackward = NULL;

        if(fTrackerGeometry=="TOTEMStrip"){
            TFile *clusterSizePlotFile = new TFile(fClusterSizePlotFile.data());
            if(clusterSizePlotFile == NULL){
                edm::LogError("debug") << "File for strip cluster with name " << fClusterSizePlotFile << " not found" ;
                exit(1);
            }
            TH1D *clusterSizePlot = (TH1D*)clusterSizePlotFile->Get(fClusterSizePlotName.data());
            if(clusterSizePlot == NULL){
                edm::LogError("debug") << "Histogram for strip cluster with name " << fClusterSizePlotName << " not found in file " << fClusterSizePlotFile;
                exit(1);
            }
            // std::cout<<1<<std::endl;
            clusterSizePlot->SetDirectory(0);
            // std::cout<<2<<std::endl;

            PPSStripDetector *stripDetector_ArmForward = new PPSStripDetector();
            // std::cout<<3<<std::endl;
            stripDetector_ArmForward->SetDetectorId(0+detectorId);
            // std::cout<<4<<std::endl;
            stripDetector_ArmForward->SetDetectorName(beam2DetectorName);
            // std::cout<<5<<std::endl;
            // unsigned int rawId = TotemRPDetId::decToRawId(1000+detectorId);
            // std::cout<<5.1<<std::endl;
            // fRpGeometry->GetDetector(rawId);
            // if(fRpGeometry==NULL){
            //     throw cms::Exception("PPSSim.h") << "Very Forward Misaligned Geometry Record not found.";
            // }
            // std::cout<<5.2<<std::endl;
            // DetGeomDesc *stripPlaneDescription_ArmForward = fRpGeometry->GetDetector(TotemRPDetId::decToRawId(0+detectorId));
            // // std::cout<<6<<std::endl;
            // if(stripPlaneDescription_ArmForward == NULL){
            //     edm::LogError("debug") << "Strip detector with ID " << 0+detectorId << " not found in file ";
            //     exit(1);
            // }
            // std::cout<<7<<std::endl;
            TVector3 rpPosition_ArmForward(
                -fTrkDetXOffset+fBeam2CenterAtStation[beam2DetectorName]-fBeam2SigmaAtStation[beam2DetectorName]*fBeam2StationSigmaInsertion[beam2DetectorName],
                +fVerticalShift,
                -fBeam2StationPositions[beam2DetectorName]*m_to_mm);
            // std::cout<<"Forward Pot Position: x = "<<rpPosition_ArmForward.X()<<" y = "<<rpPosition_ArmForward.Y()<<" z = "<<rpPosition_ArmForward.Z()<<std::endl;
            // std::cout<<8<<std::endl;
            stripDetector_ArmForward->SetRpPosition(rpPosition_ArmForward);
            // std::cout<<9<<std::endl;
            stripDetector_ArmForward->SetDetectorResolution(fHitSigmaX, fHitSigmaY, fHitSigmaZ);
            // std::cout<<10<<std::endl;
            stripDetector_ArmForward->SetNumberOfStrips(fNumberOfStrips);
            // std::cout<<11<<std::endl;
            stripDetector_ArmForward->SetStripPitch(fStripPitch);
            // std::cout<<12<<std::endl;
            stripDetector_ArmForward->SetCutSideLength(fCutSideLength);
            // std::cout<<13<<std::endl;
            stripDetector_ArmForward->SetClusterSizePlot(clusterSizePlot);
            // std::cout<<13.1<<std::endl;
            stripDetector_ArmForward->SetRpGeometry(fRpGeometry);
            // std::cout<<13.2<<std::endl;
            stripDetector_ArmForward->BuildDetectorPlanes();
            // std::cout<<14<<std::endl;
            // std::cout<<15<<std::endl;
            // stripDetector_ArmForward->AddHit(TVector3(0.,0.,0.));
            // std::cout<<16<<std::endl;
            trackingDetector_ArmForward = stripDetector_ArmForward;
            // std::cout<<17<<std::endl;
            // trackingDetector_ArmForward->AddHit(TVector3(0.,0.,0.));
            // std::cout<<18<<std::endl;

            PPSStripDetector *stripDetector_ArmBackward = new PPSStripDetector();
            // std::cout<<15<<std::endl;
            stripDetector_ArmBackward->SetDetectorId(1000+detectorId);
            // std::cout<<16<<std::endl;
            stripDetector_ArmBackward->SetDetectorName(beam1DetectorName);
            // std::cout<<17<<std::endl;
            // DetGeomDesc *stripPlaneDescription_ArmBackward = fRpGeometry->GetDetector(TotemRPDetId::decToRawId(1000+detectorId));
            // // std::cout<<18<<std::endl;
            // if(stripPlaneDescription_ArmBackward == NULL){
            //     edm::LogError("debug") << "Strip detector with ID " << 0+detectorId << " not found in file ";
            //     exit(1);
            // }
            // std::cout<<19<<std::endl;
            TVector3 rpPosition_ArmBackward(
                -fTrkDetXOffset+fBeam1CenterAtStation[beam1DetectorName]-fBeam1SigmaAtStation[beam1DetectorName]*fBeam1StationSigmaInsertion[beam1DetectorName],
                +fVerticalShift,
                -fBeam1StationPositions[beam1DetectorName]*m_to_mm);
            // std::cout<<"Backward Pot Position: x = "<<rpPosition_ArmBackward.X()<<" y = "<<rpPosition_ArmBackward.Y()<<" z = "<<rpPosition_ArmBackward.Z()<<std::endl;
            // std::cout<<20<<std::endl;
            stripDetector_ArmBackward->SetRpPosition(rpPosition_ArmBackward);
            // std::cout<<21<<std::endl;
            stripDetector_ArmBackward->SetDetectorResolution(fHitSigmaX, fHitSigmaY, fHitSigmaZ);
            // std::cout<<22<<std::endl;
            stripDetector_ArmBackward->SetNumberOfStrips(fNumberOfStrips);
            // std::cout<<23<<std::endl;
            stripDetector_ArmBackward->SetStripPitch(fStripPitch);
            // std::cout<<24<<std::endl;
            stripDetector_ArmBackward->SetCutSideLength(fCutSideLength);
            // std::cout<<25<<std::endl;
            stripDetector_ArmBackward->SetClusterSizePlot(clusterSizePlot);
            // std::cout<<25.1<<std::endl;
            stripDetector_ArmBackward->SetRpGeometry(fRpGeometry);
            // std::cout<<25.2<<std::endl;
            stripDetector_ArmBackward->BuildDetectorPlanes();
            // std::cout<<26<<std::endl;
            trackingDetector_ArmBackward = stripDetector_ArmBackward;
            // std::cout<<27<<std::endl;

            clusterSizePlotFile->Close();
            // std::cout<<28<<std::endl;
            delete clusterSizePlotFile;
            // std::cout<<29<<std::endl;
        }

        if(fTrackerGeometry=="Pixel"){
            PPSPixelDetector *pixelDetector_ArmForward = new PPSPixelDetector();
            pixelDetector_ArmForward->SetDetectorId(1000+detectorId);
            pixelDetector_ArmForward->SetDetectorName(beam2DetectorName);
            trackingDetector_ArmForward = pixelDetector_ArmForward;

            PPSPixelDetector *pixelDetector_ArmBackward = new PPSPixelDetector();
            pixelDetector_ArmForward->SetDetectorId(0+detectorId);
            pixelDetector_ArmBackward->SetDetectorName(beam1DetectorName);
            trackingDetector_ArmBackward = pixelDetector_ArmBackward;
        }

        fTrackingStationForward[detectorId]=trackingDetector_ArmForward;
        fTrackingStationBackward[detectorId]=trackingDetector_ArmBackward;
        // PPSStripDetector *stripDetectorF = (stripDetectorF*)fTrackingStationForward[detectorId];
        // std::cout<<"Direction = F hit x = "<<0.<<" hit y = "<<0.<<" hit z = "<<0.<<std::endl;
        // // stripDetectorF->AddHit(TVector3(0.,0.,0.));
        // fTrackingStationForward[detectorId]->AddHit(TVector3(0.,0.,0.));
        // PPSStripDetector *stripDetectorB = (PPSStripDetector*)&fTrackingStationBackward[detectorId];
        // std::cout<<"Direction = B hit x = "<<0.<<" hit y = "<<0.<<" hit z = "<<0.<<std::endl;
        // stripDetectorB->AddHit(TVector3(0.,0.,0.));
        
    }

    if(fToFGeometry=="diamond") {
        ToFDet_F  = new PPSToFDetector(fToFWidth,fToFHeight,fBeam2StationSigmaInsertion[fBeam2ToFDetectorName]*fBeam2SigmaAtStation[fBeam2ToFDetectorName]+fToFXOffsetF+fBeam2CenterAtStation[fBeam2ToFDetectorName]);
        ToFDet_F->SetDetectorName(fBeam2ToFDetectorName);
        ToFDet_B  = new PPSToFDetector(fToFWidth,fToFHeight,fBeam1StationSigmaInsertion[fBeam1ToFDetectorName]*fBeam1SigmaAtStation[fBeam1ToFDetectorName]+fToFXOffsetF+fBeam1CenterAtStation[fBeam1ToFDetectorName]);
        ToFDet_B->SetDetectorName(fBeam1ToFDetectorName);
    }
         
    //   Check the overall kinematic limits
    if (fGenMode=="M_X") {
        if (xi_min*(2.*fBeamEnergy)>fCentralMass+fCentralMassErr||xi_max*(2.*fBeamEnergy)<fCentralMass-fCentralMassErr) {
            edm::LogWarning("debug") << "xi limits outside kinematic limits for the given Central mass. Stopping..." ;
            exit(1);
        }
    }
}

//--------------------------------------------------------------------------------------------------------------//

void PPSSim::EndRun()
{
}

//--------------------------------------------------------------------------------------------------------------//

void PPSSim::BeginEvent()
{
    std::cout<< " !!! ------ PPSSim::BeginEvent ------ !!!\n";
    fGen->clear();
    fSim->clear();
    fReco->clear();
    protonF = NULL;
    protonB = NULL;
    fHasStoppedF = false;
    fHasStoppedB = false;
    NVertex=0;
    fVertex.clear();
    protonsOut.clear();
    fHasStopped.clear();

    for(unsigned i=0; i< fTrackingDetectorIDs.size(); ++i){
        int detectorId = fTrackingDetectorIDs.at(i);
        fTrackingStationForward[detectorId]->Clear();

    }
    
    ToFDet_F->clear();
    ToFDet_B->clear();

}

//--------------------------------------------------------------------------------------------------------------//

void PPSSim::EndEvent()
{
    std::cout<< " !!! ------ PPSSim::EndEvent ------ !!!\n";

    for(int i=0;i<NVertex;i++) {
        protonF = (protonsOut[i].first); protonB = (protonsOut[i].second);
        if (protonF) delete protonF;
        if (protonB) delete protonB;
    }
    protonsOut.clear();
    protonF=NULL;
    protonB=NULL;
}

//--------------------------------------------------------------------------------------------------------------//

void PPSSim::Run() 
{
    std::cout<< " !!! ------ PPSSim::Run ------ !!!\n";

    if (!fExternalGenerator) Generation();
    if (fVerbose) edm::LogWarning("debug") << "PPSSim: Starting Simulation step."; 
    Simulation();
    if (fVerbose) edm::LogWarning("debug") << "PPSSim: Starting Digitization step.";
    Digitization(); // a fake digitization procedure
    if (fVerbose) edm::LogWarning("debug") << "PPSSim: Starting Reconstruction step.";
    Reconstruction();
}

//--------------------------------------------------------------------------------------------------------------//

void PPSSim::Generation()
{
    // Uses the CMS units in the vertex and mm for the PPS parameters
    // sorts a proton in the forward direction
    double t1,xi1,phi1;
    double t2,xi2,phi2;
    if (fGenMode=="M_X") {
        if (fCentralMass==0) {
            edm::LogWarning("debug") << "PPSSim::Generation: Central mass not defined. Exiting...";
            exit(1);
        }
        GenCentralMass(t1,t2,xi1,xi2,phi1,phi2);
    }
    else {
        GenSingleParticle(t1,xi1,phi1);
        GenSingleParticle(t2,xi2,phi2);
    }
    int nvtx=add_Vertex(fVtxMeanX,fVtxMeanY,fVtxMeanZ);
    protonF=new TLorentzVector(shoot(t1,xi1,phi1,1)); //  1 = positive/forward direction
    protonB=new TLorentzVector(shoot(t2,xi2,phi2,-1));
    add_OutgoingParticle(nvtx-1,protonF,protonB);
    fHasStoppedF=false;fHasStoppedB=false;
    set_GenData();
}

//--------------------------------------------------------------------------------------------------------------//

//
// fill the tree structures 
void PPSSim::ReadGenEvent(const std::vector<TrackingParticle>* gentrackingP)
{
}

//--------------------------------------------------------------------------------------------------------------//

//Using reco gen particle
void PPSSim::ReadGenEvent(const std::vector<reco::GenParticle>* genP)
{
    if (!genP) return;
    int ivtx = -1;
    int colId = -1;
    double momB = 0.;
    double momF = 0.;
    TLorentzVector* pF = NULL;
    TLorentzVector* pB = NULL;
    double vtxX=0.;
    double vtxY=0.;
    double vtxZ=0.;
    for(size_t i=0;i<genP->size();++i) {
        const reco::GenParticle& p = (*genP)[i];
        if (p.pdgId()!=2212) continue;
        if (p.status()  !=1) continue;
        //double pz = p.pt()*sinh (p.eta());
        double px = p.px();double py=p.py();double pz = p.pz();
        if (ivtx<0) {
            ivtx=0;
            vtxX=p.vx();vtxY=p.vy();vtxZ=p.vz();// Contrary to HepMC, reco::genParticle already uses cm, so no convertion is needed
            colId = p.collisionId();
        }

        if (colId!=p.collisionId()) {
            int nvtx=add_Vertex(vtxX,vtxY,vtxZ);
            if (ivtx!=nvtx-1) {edm::LogWarning("debug")<< "WARNING: unexpected vertex number.";}
            add_OutgoingParticle(nvtx-1,pF,pB);
            colId = p.collisionId();
            vtxX=p.vx();vtxY=p.vy();vtxZ=p.vz();
            ivtx++;
            momF=0.; momB=0.;
            pF=NULL;pB=NULL;
        } else {
            // verify the vertex consistency
            if (vtxX!=p.vx()||vtxY!=p.vy()||vtxZ!=p.vz()) {
                edm::LogWarning("debug") << "WARNING: unexpected new vertex position";
            }
        }
        if (p.eta()>0&&momF<pz) {
            momF=pz; pF = new TLorentzVector(px,py,pz,sqrt(px*px+py*py+pz*pz+ProtonMassSQ));
        } else if (p.eta()<0&&momB<fabs(pz)) {
            momB=fabs(pz);pB = new TLorentzVector(px,py,pz,sqrt(px*px+py*py+pz*pz+ProtonMassSQ));
        }
        // this is the  last particle, add it anyway..
        if (i==genP->size()-1) {
            int nvtx=add_Vertex(vtxX,vtxY,vtxZ);
            if (ivtx!=nvtx-1) {edm::LogWarning("debug") << "WARNING: unexpected vertex number.";}
            if (fVerbose) {if(pF) pF->Print();if (pB) pB->Print();}
            add_OutgoingParticle(nvtx-1,pF,pB);
        }
    }
    set_GenData();
}

//--------------------------------------------------------------------------------------------------------------//

void PPSSim::ReadGenEvent(const HepMC::GenEvent* evt)
{
    using namespace CLHEP;
    TLorentzVector* pF = NULL;
    TLorentzVector* pB = NULL;
    if (!evt) return;
    int nvtx =0;
    double vtxX = 0.;
    double vtxY = 0.;
    double vtxZ = 0.;
    for(HepMC::GenEvent::vertex_const_iterator ivtx = evt->vertices_begin();ivtx!=evt->vertices_end();ivtx++) {
        if (fVerbose) (*ivtx)->print();
        if ((*ivtx)->id()!=0) continue;
        vtxX = (*ivtx)->position().x()/cm; // CMS uses cm but HepMC uses mm
        vtxY = (*ivtx)->position().y()/cm;
        vtxZ = (*ivtx)->position().z()/cm;

        // choose the highest momentum particle on each side to be propagated
        double momF = 0.; double momB = 0.;

        for(HepMC::GenVertex::particles_out_const_iterator pItr = (*ivtx)->particles_out_const_begin();
                pItr!= (*ivtx)->particles_out_const_end();pItr++) {
            if (fVerbose) (*pItr)->print();
            if ((*pItr)->status() != 1) continue; // this is not a final state particle
            if ((*pItr)->pdg_id()!=2212) continue; // only protons to be processed
            if (fabs((*pItr)->momentum().eta()) < fEtaMin) continue; 
            if ((*pItr)->momentum().e()<fMomentumMin) continue; 
            if ((*pItr)->momentum().pz()>0&&(*pItr)->momentum().pz()>momF) {
                momF= (*pItr)->momentum().pz();
                pF  = new TLorentzVector((*pItr)->momentum().px(),(*pItr)->momentum().py(),(*pItr)->momentum().pz(),
                        sqrt(ProtonMassSQ+pow((*pItr)->momentum().px(),2)+pow((*pItr)->momentum().py(),2)+pow((*pItr)->momentum().pz(),2)));
            } else if ((*pItr)->momentum().pz()<0&&fabs((*pItr)->momentum().pz())>momB){
                momB = fabs((*pItr)->momentum().pz());
                pB   = new TLorentzVector((*pItr)->momentum().px(),(*pItr)->momentum().py(),(*pItr)->momentum().pz(),
                        sqrt(ProtonMassSQ+pow((*pItr)->momentum().px(),2)+pow((*pItr)->momentum().py(),2)+pow((*pItr)->momentum().pz(),2)));
            }
        }
    }
    if (!pF&&!pB) return;
    nvtx=add_Vertex(vtxX,vtxY,vtxZ);
    add_OutgoingParticle(nvtx-1,pF,pB);
    if (fVerbose) {
        if (pF) {pF->Print();}
        if (pB) {pB->Print();}
    }
    set_GenData();
}

//--------------------------------------------------------------------------------------------------------------//

void PPSSim::set_GenData()
{
    for(int i=0;i<NVertex;i++) {
        protonF = (protonsOut[i].first); protonB = (protonsOut[i].second);
        int tF=-1;if (protonF) tF=i;
        int tB=-1;if (protonB) tB=i;
        fGen->Vertices->Add(fVertex[i].x(),fVertex[i].y(),fVertex[i].z(),tF,tB);
        double t,xi;
        if (protonF) {Get_t_and_xi(protonF,t,xi); fGen->ArmF.addParticle(*protonF,t,xi);}
        if (protonB) {Get_t_and_xi(protonB,t,xi); fGen->ArmB.addParticle(*protonB,t,xi);}
        if (protonF){
            ApplyBeamSmearing(const_cast<TLorentzVector&>(*protonF));
        }
        if (protonB) {
            ApplyBeamSmearing(const_cast<TLorentzVector&>(*protonB));
        }
    }
}

//--------------------------------------------------------------------------------------------------------------//

void PPSSim::Get_t_and_xi(const TLorentzVector* proton,double& t,double& xi) 
{
    t = 0.;
    xi = -1.;
    if (!proton) return;
    double mom    = proton->P();
    if (mom>fBeamMomentum) mom=fBeamMomentum;
    double energy = proton->E();
    double theta  = (proton->Pz()>0)?proton->Theta():TMath::Pi()-proton->Theta();
    t      = -2.*(ProtonMassSQ-fBeamEnergy*energy+fBeamMomentum*mom*cos(theta));
    xi     = (1.0-energy/fBeamEnergy);
}

//--------------------------------------------------------------------------------------------------------------//

void PPSSim::Simulation()
{
    for(int i=0;i<NVertex;i++) {
        double vtxX=fVertex[i].x();
        double vtxY=fVertex[i].y();
        double vtxZ=fVertex[i].z();
        protonF = (protonsOut[i].first); 
        protonB = (protonsOut[i].second);
        
        fHasStoppedF=false;
        fHasStoppedB=false;
        
        // At this point, one should be using the CMS units (cm)
        int tF=-1; 
        if (protonF) tF=i;
        int tB=-1; 
        if (protonB) tB=i;

        fSim->Vertices->Add(vtxX,vtxY,vtxZ,tF,tB);
        // FIRST, propagate to the positive(forward) direction, then to the other side
        //
        // Propagate until PPS, filling pps_station accordingly for the reconstruction if needed
        // Remember: 
        //         HECTOR uses um for X and Y coordinates and m for Z
        //         is backward for LHC, which means when propagating to the CMS forward direction, one needs to rotate it
        //         by means of doing x_LHC = -x_CMS and z = -z 

        H_BeamParticle *part = NULL;
        if (protonF) {
            int Direction = 1;
            double t,xi;
            Get_t_and_xi(protonF,t,xi);
            fSim->ArmF.AddTrack(*protonF,t,xi);
            //
            if (fCrossAngleCorr) LorentzBoost(const_cast<TLorentzVector&>(*protonF),"LAB");
            //
            part = new H_BeamParticle(ProtonMass,1);
            part->setPosition(-(vtxX-fVtxMeanX)*cm_to_um,(vtxY-fVtxMeanY)*cm_to_um,0.,0.,-(vtxZ)*cm_to_m);
            part->set4Momentum(-protonF->Px(),protonF->Py(),-protonF->Pz(),protonF->E());
            part->computePath(beamlineF);
            // Propagate(part,Direction);
            if (part) {delete part;part = NULL;}
        } 
        //
        //  Propagate to the negative/backward direction
        //
        if (protonB) {
            int Direction = -1;
            double t,xi;
            Get_t_and_xi(protonB,t,xi);
            fSim->ArmB.AddTrack(*protonB,t,xi);
            //
            if (fCrossAngleCorr) LorentzBoost(const_cast<TLorentzVector&>(*protonB),"LAB");
            //
            part = new H_BeamParticle(ProtonMass,1);
            part->setPosition(-(vtxX-fVtxMeanX)*cm_to_um,(vtxY-fVtxMeanY)*cm_to_um,0.,0.,-(vtxZ)*cm_to_m);
            part->set4Momentum(-protonB->Px(),protonB->Py(),-protonB->Pz(),protonB->E()); // HECTOR uses always positive z momentum
            part->computePath(beamlineB);
            // Propagate(part,Direction);
            if (part) {delete part;part = NULL;}
        } 
    }
}

//--------------------------------------------------------------------------------------------------------------//

void PPSSim::Reconstruction()
{
    int Direction;
    Direction=1;
    TrackReco(Direction,pps_stationF,&(fReco->ArmF));
    Direction=-1;
    TrackReco(Direction,pps_stationB,&(fReco->ArmB));
    VertexReco();
}

//--------------------------------------------------------------------------------------------------------------//

bool PPSSim::SearchTrack( TGraphErrors *xLineProjection, TGraphErrors *yLineProjection ,int Direction,double& xi,double& t,double& partP,double& pt,double& thx,double& thy,double& x0, double& y0, double &xChiSquare, double &yChiSquare)
{
    std::cout<< " !!! ------ PPSSim::SearchTrack ------ !!!\n";
    double theta=0.;
    xi = 0; t=0; partP=0; pt=0; x0=0.;y0=0.;
    H_RecRPObject*  station = NULL;
    if (Direction>0) {
        station = pps_stationF;
    } else {
        station = pps_stationB;
    }

    double x1 = 0.; double y1 = 0.;
    double x2 = 0.; double y2 = 0.;

    TF1 *xLine = new TF1("xLine","pol1",fTrackerZPosition,fToFZPosition);
    xLineProjection->Fit(xLine,"Q0");
    if(xLineProjection->GetN()>2) xChiSquare = xLine->GetChisquare();
    else xChiSquare = 0.;
    if (fVerbose) edm::LogWarning("debug") << "Direction "<<Direction<<" Line x chiSquare "<<xChiSquare<<" number of points "<<xLineProjection->GetN();
    if(xLineProjection->GetN()>2 && xChiSquare>fXTrackChiSquareCut) return false;

    TF1 *yLine = new TF1("yLine","pol1",fTrackerZPosition,fToFZPosition);
    yLineProjection->Fit(yLine,"Q0");
    if(yLineProjection->GetN()>2) yChiSquare = yLine->GetChisquare();
    else yChiSquare = 0.;
    if (fVerbose) edm::LogWarning("debug") <<"Direction "<<Direction<<" Line y chiSquare "<<yChiSquare<<" number of points "<<yLineProjection->GetN();
    if(yLineProjection->GetN()>2 && yChiSquare>fYTrackChiSquareCut) return false;

    x1 = xLine->Eval(fTrackerZPosition);
    x2 = xLine->Eval(fTrackerZPosition+fTrackerLength);
    if (fVerbose) edm::LogWarning("debug") << "Direction "<<Direction<<" xTracker1 "<<x1<<" xTracker2 "<<x2;

    y1 = yLine->Eval(fTrackerZPosition);
    y2 = yLine->Eval(fTrackerZPosition+fTrackerLength);
    if (fVerbose) edm::LogWarning("debug") << "Direction "<<Direction<<" yTracker1 "<<y1<<" yTracker2 "<<y2;

    double eloss;
    ReconstructArm(station, x1,y1,x2,y2,thx,thy,eloss);
    // Protect for unphysical results
    if (std::isnan(eloss)||std::isinf(eloss)||
            std::isnan(thx)  || std::isinf(thx) ||
            std::isnan(thy)  || std::isinf(thy)) return false;
    //
    if (-thx<-100||-thx>300) return false;
    if (thy<-200||thy>200) return false;
    //
    x0 = -station->getX0()*um_to_cm;
    y0 = station->getY0()*um_to_cm;
    double ImpPar=sqrt(x0*x0+y0*y0);
    if (fTrackImpactParameterCut>0.) {
        if (ImpPar>fTrackImpactParameterCut) return false;
    }
    if (eloss<0||eloss>fBeamEnergy) return false;
    theta = sqrt(thx*thx+thy*thy)*urad;
    xi    = eloss/fBeamEnergy;
    double energy= fBeamEnergy*(1.-xi);
    partP = sqrt(energy*energy-ProtonMassSQ);
    t     = -2.*(ProtonMassSQ - fBeamEnergy*energy+fBeamMomentum*partP*cos(theta));
    pt    = sqrt(pow(partP*thx*urad,2)+pow(partP*thy*urad,2));
    if (xi<0.||xi>1.||t<0.||t>10.||pt<=0.) {
        xi = 0; t=0; partP=0; pt=0; theta=0; x0=0.;y0=0.;
        return false; // unphysical values 
    }
    return true;
}

//--------------------------------------------------------------------------------------------------------------//

void PPSSim::TrackReco(int Direction,H_RecRPObject* station,PPSBaseData* arm_base)
{
    std::cout<< " !!! ------ PPSSim::TrackReco ------ !!!\n";
    //
    PPSRecoData* arm = dynamic_cast<PPSRecoData*>(arm_base);
    double xi,t,partP,pt,phi,theta,x0,y0,thx,thy,xChiSquare, yChiSquare;
    std::map<int,PPSTrkDetector*> trakingDetectorMap;
    PPSToFDetector* ToF  = NULL;
    if (Direction>0) {
        trakingDetectorMap = fTrackingStationForward;
        ToF=ToFDet_F;
    } else {
        trakingDetectorMap = fTrackingStationBackward;
        ToF=ToFDet_B;
    }

    std::vector<TVector3> hitTracker1 = trakingDetectorMap[fTrackingDetectorIDs.at(0)]->GetSmearedHits(); 
    std::vector<TVector3> hitTracker2 = trakingDetectorMap[fTrackingDetectorIDs.at(1)]->GetSmearedHits(); 
    
    for(unsigned i=0;i<hitTracker1.size();i++){
        arm->AddHitTrk1(hitTracker1.at(i).X(),hitTracker1.at(i).Y());
        // std::cout<<"Direction = "<<Direction<<" Recorded hit x = "<<hitTracker1.at(i).X()<<" hit y = "<<hitTracker1.at(i).Y()<<" hit z = "<<hitTracker1.at(i).Z()<<std::endl;
    }
    for(unsigned i=0;i<hitTracker2.size();i++){
        arm->AddHitTrk2(hitTracker2.at(i).X(),hitTracker2.at(i).Y());
        // std::cout<<"Direction = "<<Direction<<" Recorded hit x = "<<hitTracker2.at(i).X()<<" hit y = "<<hitTracker2.at(i).Y()<<" hit z = "<<hitTracker2.at(i).Z()<<std::endl;
    }  
    unsigned numberOfTrk1Points = (hitTracker1.size()==0) ? 1 : hitTracker1.size();
    unsigned numberOfTrk2Points = (hitTracker2.size()==0) ? 1 : hitTracker2.size();
    unsigned numberOfToFPoints = (ToF->get_NHits()==0 || !fUseToFForTracking) ? 1 : ToF->get_NHits();

    for(unsigned i=0;i<numberOfTrk1Points;i++) {
        for(unsigned j=0;j<numberOfTrk2Points;j++){
            for(unsigned k=0; k<numberOfToFPoints; k++){
                vector<double> xPoints;
                vector<double> xPointsError;
                vector<double> yPoints;
                vector<double> yPointsError;
                vector<double> zPoints;
                if(i<hitTracker1.size()){
                    xPoints.push_back(hitTracker1.at(i).X());
                    yPoints.push_back(hitTracker1.at(i).Y());
                    zPoints.push_back(fTrackerZPosition);
                    xPointsError.push_back(fHitSigmaX);
                    yPointsError.push_back(fHitSigmaY);
                    if (fVerbose) edm::LogWarning("debug") <<"Direction "<<Direction<<" added point for Tracker 1: x = "<< xPoints.back() <<" +o- "<<xPointsError.back()<<" y = "<< yPoints.back() <<" +o- "<<yPointsError.back();
                }
                if(j<hitTracker2.size()){
                    xPoints.push_back(hitTracker2.at(i).X());
                    yPoints.push_back(hitTracker2.at(j).Y());
                    zPoints.push_back(fTrackerZPosition+fTrackerLength);
                    xPointsError.push_back(fHitSigmaX);
                    yPointsError.push_back(fHitSigmaY);
                    if (fVerbose) edm::LogWarning("debug") <<"Direction "<<Direction<<" added point for Tracker 2: x = "<< xPoints.back() <<" +o- "<<xPointsError.back()<<" y = "<< yPoints.back() <<" +o- "<<yPointsError.back();
                }
                if(k<ToF->get_NHits() && fUseToFForTracking){
                    // std::cout<<1<<std::endl;
                    xPoints.push_back(ToF->X.at(k));
                    // std::cout<<2<<std::endl;
                    yPoints.push_back(ToF->Y.at(k));
                    // std::cout<<3<<std::endl;
                    zPoints.push_back(fToFZPosition);
                    // std::cout<<4<<std::endl;
                    xPointsError.push_back(fToFHitSigmaX);
                    // std::cout<<5<<std::endl;
                    yPointsError.push_back(fToFHitSigmaY);
                    // std::cout<<6<<std::endl;
                    if (fVerbose) edm::LogWarning("debug") <<"Direction "<<Direction<<" added point for ToF      : x = "<< xPoints.back() <<" +o- "<<xPointsError.back()<<" y = "<< yPoints.back() <<" +o- "<<yPointsError.back();
                }

                if(xPoints.size()<2) continue;
                // std::cout<<xPoints.size()<<std::endl;

                TGraphErrors *xLineProjection = new TGraphErrors(xPoints.size(), &(zPoints[0]), &(xPoints[0]), 0, &(xPointsError[0]));
                TGraphErrors *yLineProjection = new TGraphErrors(xPoints.size(), &(zPoints[0]), &(yPoints[0]), 0, &(yPointsError[0]));

                if (SearchTrack(xLineProjection,yLineProjection,Direction,xi,t,partP,pt,thx,thy,x0,y0,xChiSquare,yChiSquare)) {
                    if (fVerbose) edm::LogWarning("debug") <<"Direction "<<Direction<<" added track!!!";
                    theta = sqrt(thx*thx+thy*thy)*urad;
                    phi   = (Direction>0)?-atan2(thy,-thx):atan2(thy,thx); // defined according to the positive direction
                    if (Direction<0) { theta=TMath::Pi()-theta; }
                    double px = partP*sin(theta)*cos(phi);
                    double py = partP*sin(theta)*sin(phi);
                    double pz = partP*cos(theta);
                    double  e = sqrt(partP*partP+ProtonMassSQ);
                    TLorentzVector p(px,py,pz,e);
                    if (fCrossAngleCorr) LorentzBoost(p,"MC");
                    Get_t_and_xi(const_cast<TLorentzVector*>(&p),t,xi);
                    arm->AddTrack(p,t,xi);

                    if(i<hitTracker1.size()) arm->get_Track().set_HitDet1(hitTracker1.at(i).X(),hitTracker1.at(i).Y());
                    if(j<hitTracker2.size()) arm->get_Track().set_HitDet2(hitTracker2.at(i).X(),hitTracker2.at(j).Y());
                    if(k<ToF->get_NHits()){
                        arm->get_Track().set_HitToF (ToF->ToF.at(k),ToF->X.at(k),ToF->X.at(k));
                        arm->get_Track().set_TimeOfFlight(ToF->ToF.at(k));
                    }
                    else{
                        arm->get_Track().set_TimeOfFlight(-1.);
                    }
                    arm->get_Track().set_X0(x0);
                    arm->get_Track().set_Y0(y0);
                    arm->get_Track().set_Phi(phi);
                    arm->get_Track().set_ThetaAtIP(thx,thy); // thx is given in CMS coordinates
                    arm->get_Track().set_XTrackChiSquare(xChiSquare);
                    arm->get_Track().set_YTrackChiSquare(yChiSquare);
                       
                }
                else if (fVerbose) edm::LogWarning("debug") <<"Direction "<<Direction<<" no track found!!!";

                delete xLineProjection;
                delete yLineProjection;

            }
        } 
    }
    if (fVerbose) edm::LogWarning("debug") <<"Direction "<<Direction<<" track number: "<<arm->Tracks.size();

}

//--------------------------------------------------------------------------------------------------------------//

void PPSSim::VertexReco()
{
    if (fVerbose) edm::LogWarning("debug") <<"Starting Vertex Reco";
    PPSRecoTracks* tracksF=&(fReco->ArmF.Tracks);
    if (fVerbose) edm::LogWarning("debug") <<"ArmF track number: "<<tracksF->size();
    if (tracksF->size()==0) return;
    PPSRecoTracks* tracksB=&(fReco->ArmB.Tracks);
    if (fVerbose) edm::LogWarning("debug") <<"ArmB track number: "<<tracksB->size();
    if (tracksB->size()==0) return;
    PPSRecoVertex* vtxs = (PPSRecoVertex*)fReco->Vertices;
    vtxs->clear();

    double vtxX,vtxY,vtxZ;
    double tofF,tofB,ToFtot,d_ToFtot;
    // double xtF,ytF,xtB,ytB;
    // int cellidF=0,cellidB=0;

    d_ToFtot = sqrt(2.)*fTimeSigma; // uncertainty on ToFtot due to detector resolution
    int Nsigma = 3.0;              // # of sigmas (CL for vertex reconstruction)

    for(int i=0;i<(int)tracksF->size();i++){
        tofF=tracksF->at(i).get_TimeOfFlight();
        if(tofF==-1.) continue;
        for(int j=0;j<(int)tracksB->size();j++) {
            tofB=tracksB->at(j).get_TimeOfFlight();
            if(tofB==-1.) continue;
            // if (ToFDet_B->GetADC(cellidB,l)==0) edm::LogWarning("debug") << "WARNING: no ADC found";
            ToFtot = tofF+tofB;
            cout<<fabs(ToFtot-2*fToFZPosition/c_light_ns)<<"  "<<Nsigma*d_ToFtot<<endl;
            if (fabs(ToFtot-2*fToFZPosition/c_light_ns)>Nsigma*d_ToFtot) continue;
            vtxZ=-c_light_ns*(tofF-tofB)/2.0*m_to_cm;
            vtxX=(tracksF->at(i).get_X0()+tracksB->at(j).get_X0())/2.; // this is not very meaningful, there is not enough precision
            vtxY=(tracksF->at(i).get_Y0()+tracksB->at(j).get_Y0())/2.; // idem
            vtxs->Add(vtxX,vtxY,vtxZ,i,j);
        } 
    } 
}

//--------------------------------------------------------------------------------------------------------------//

void PPSSim::ReconstructArm(H_RecRPObject* pps_station, double x1, double y1, double x2, double y2, double& tx, double& ty, double& eloss)
{
    tx=0.;
    ty=0.;
    eloss=0.;
    if (!pps_station) return;
    // Change the orientation and units according to Hector
    x1*=-mm_to_um;
    x2*=-mm_to_um;
    y1*= mm_to_um;
    y2*= mm_to_um;
    pps_station->setPositions(x1,y1,x2,y2);
    double energy = pps_station->getE(AM); // dummy call needed to calculate some Hector internal parameter
    if (std::isnan(energy)||std::isinf(energy)) return;
    tx =  -pps_station->getTXIP();  // change orientation to CMS
    ty =  pps_station->getTYIP();
    eloss = pps_station->getE();
}

//--------------------------------------------------------------------------------------------------------------//

void PPSSim::Digitization()
{
    //    Fake method to mimic a digitization procedure
    //    Just copy the information from the fSim branch and smear the hit according to a given
    //    detector resolution;
    TrackerDigi();

    int Direction=1;
    ToFDigi(Direction,&(fSim->ArmF),ToFDet_F);
    Direction=-1;
    ToFDigi(Direction,&(fSim->ArmB),ToFDet_B);

}

//--------------------------------------------------------------------------------------------------------------//

void PPSSim::TrackerDigi()
{

    std::vector<edm::DetSet<TotemRPDigi> > rpStripClusters;

    for(map<int,PPSTrkDetector*>::iterator dIt=fTrackingStationForward.begin(); dIt!=fTrackingStationForward.end(); ++dIt){
        if(fTrackerGeometry == "TOTEMStrip"){
            PPSStripDetector *trackingStation = (PPSStripDetector*)dIt->second;
            std::vector<edm::DetSet<TotemRPDigi> > stationDigi = trackingStation->GetClusters();
            rpStripClusters.insert(rpStripClusters.end(), stationDigi.begin(), stationDigi.end());
        }
        // sum clusters
    }

    for(map<int,PPSTrkDetector*>::iterator dIt=fTrackingStationBackward.begin(); dIt!=fTrackingStationBackward.end(); ++dIt){
        if(fTrackerGeometry == "TOTEMStrip"){
            PPSStripDetector *trackingStation = (PPSStripDetector*)dIt->second;
            std::vector<edm::DetSet<TotemRPDigi> > stationDigi = trackingStation->GetClusters();
            rpStripClusters.insert(rpStripClusters.end(), stationDigi.begin(), stationDigi.end());
        }
        // sum clusters
    }

    // std::cout<<"Number of hit strips = "<<rpStripClusters.size()<<std::endl;
    fTotemDigi = new edm::DetSetVector<TotemRPDigi>(rpStripClusters);

}

//--------------------------------------------------------------------------------------------------------------//

void PPSSim::ToFDigi(int Direction, const PPSBaseData* arm_sim,PPSToFDetector* ToFDet)
{
    
    double toFOffset=0.;
    double beamXRMS_ToF=0.;
    double toFInsertion=0.;
    if (Direction>0) {
        toFOffset = fToFXOffsetF;
        beamXRMS_ToF = fBeam2SigmaAtStation[fBeam2ToFDetectorName];
        toFInsertion = fBeam2StationSigmaInsertion[fBeam2ToFDetectorName];
    } else {
        toFOffset = fToFXOffsetB;
        beamXRMS_ToF = fBeam1SigmaAtStation[fBeam1ToFDetectorName];
        toFInsertion = fBeam1StationSigmaInsertion[fBeam1ToFDetectorName];
    }

    if(!arm_sim||!ToFDet) return;
    // what direction?
    PPSRecoData* arm_reco=NULL;
    ToFDet->clear();
    if (ToFDet==ToFDet_F) arm_reco = &(fReco->ArmF);
    if (ToFDet==ToFDet_B) arm_reco = &(fReco->ArmB);
    //
    for(int i=0;i<const_cast<PPSBaseData*>(arm_sim)->ToFDet.NHits();i++){
        double x = arm_sim->ToFDet.at(i).X;
        double y = arm_sim->ToFDet.at(i).Y;
        // if (fFilterHitMap&&(x>fDetectorClosestX||x<fMaxXfromBeam||fabs(y)>fabs(fMaxYfromBeam))) {
        //     continue;
        // }
        if (fApplyFiducialCuts) {
            double xmin = toFInsertion*beamXRMS_ToF + toFOffset;
            double xmax = xmin+fToFWidth;
            if (fabs(x)<xmin||fabs(x)>xmax||fabs(y)>fabs(fToFHeight/2)) { // use ABS because the detector are on the negative X side
                continue;
            }
        }
        double t = arm_sim->ToFDet.at(i).ToF;
        if (t>0) ToFSmearing(t);
        ToFDet->AddHit(x,y,t);
        //cout<<"ToF Digi x = "<<x<<endl;
        if (!arm_reco) continue;
        //int cellid = ToFDet->findCellId(x,y);
        //if (cellid==0) continue;
        // find x,y of the center of the cell
        x = gRandom3->Gaus(x,fToFHitSigmaX);
        // double xc=0;
        // double yc=0;
        // if (ToFDet->get_CellCenter(cellid,xc,yc)) arm_reco->AddHitToF(cellid,t,xc,yc);
        // else arm_reco->AddHitToF(cellid,t,0.,0.);
        arm_reco->AddHitToF(t,x,0.);
        //cout<<"ToF Digi smeared x = "<<x<<endl;
    }
}

//--------------------------------------------------------------------------------------------------------------//

void PPSSim::GenSingleParticle(double& t, double& xi, double& phi)
{
    phi = gRandom3->Uniform(fPhiMin,fPhiMax);
    if (fGenMode=="linear"||fGenMode=="uniform") {
        xi = gRandom3->Uniform(xi_min,xi_max);
        t  = gRandom3->Uniform(t_min,t_max);
    }
    else if (fGenMode=="log"){
        if (t_min==0) t_min = 1e-6; // limit t to 1 MeV
        xi = pow(10,gRandom3->Uniform(log10(xi_min),log10(xi_max)));
        t  = pow(10,gRandom3->Uniform(log10(t_min),log10(t_max)));
    }
    double min_t = Minimum_t(xi);
    if (t<min_t) t = min_t;
}

//--------------------------------------------------------------------------------------------------------------//

void PPSSim::GenCentralMass(double& t1, double& t2, double& xi1, double& xi2, double& phi1, double& phi2)
{
    if (fCentralMassErr>0) {
        double m_h = gRandom3->Gaus(fCentralMass,fCentralMassErr);
        while(1) {
            xi1 = gRandom3->Uniform(xi_min,xi_max);
            xi2 = gRandom3->Uniform(xi_min,xi_max);
            double mh_2 = sqrt(xi1*xi2)*2.*fBeamEnergy;
            if ((fabs(m_h-mh_2)<fCentralMassErr) &&
                    (isPhysical(xi1)&&isPhysical(xi2))) break;// check validity of kinematic region
        }
    } else {
        while(1) {
            xi1 = gRandom3->Uniform(xi_min,xi_max);
            xi2 = pow(0.5*fCentralMass/fBeamEnergy,2)/xi1;
            if (isPhysical(xi1)&&isPhysical(xi2)) break;
        }
    }

    phi1 = gRandom3->Uniform(fPhiMin,fPhiMax);
    phi2 = gRandom3->Uniform(fPhiMin,fPhiMax);
    t1   = gRandom3->Uniform(Minimum_t(xi1),t_max);
    t2   = gRandom3->Uniform(Minimum_t(xi2),t_max);
}

//--------------------------------------------------------------------------------------------------------------//

void PPSSim::LorentzBoost(TLorentzVector& p_out, const string& frame)
{
    // Use a matrix
    double microrad = 1.e-6;
    TMatrixD tmpboost(4,4);
    double alpha_ = 0.;
    double phi_  = fCrossingAngle*microrad;
    if (p_out.Pz()<0) phi_*=-1;
    tmpboost(0,0) = 1./cos(phi_);
    tmpboost(0,1) = - cos(alpha_)*sin(phi_);
    tmpboost(0,2) = - tan(phi_)*sin(phi_);
    tmpboost(0,3) = - sin(alpha_)*sin(phi_);
    tmpboost(1,0) = - cos(alpha_)*tan(phi_);
    tmpboost(1,1) = 1.;
    tmpboost(1,2) = cos(alpha_)*tan(phi_);
    tmpboost(1,3) = 0.;
    tmpboost(2,0) = 0.;
    tmpboost(2,1) = - cos(alpha_)*sin(phi_);
    tmpboost(2,2) = cos(phi_);
    tmpboost(2,3) = - sin(alpha_)*sin(phi_);
    tmpboost(3,0) = - sin(alpha_)*tan(phi_);
    tmpboost(3,1) = 0.;
    tmpboost(3,2) = sin(alpha_)*tan(phi_);
    tmpboost(3,3) = 1.;

    if(frame=="LAB") tmpboost.Invert();

    TMatrixD p4(4,1);
    p4(0,0) = p_out.E();
    p4(1,0) = p_out.Px();
    p4(2,0) = p_out.Py();
    p4(3,0) = p_out.Pz();
    TMatrixD p4lab(4,1);
    p4lab = tmpboost * p4;
    p_out.SetPxPyPzE(p4lab(1,0),p4lab(2,0),p4lab(3,0),p4lab(0,0));
}

//--------------------------------------------------------------------------------------------------------------//

void PPSSim::ApplyBeamSmearing(TLorentzVector& p_out)
{
    double microrad = 1.e-6;
    double theta = p_out.Theta(); if (p_out.Pz()<0) theta=TMath::Pi()-theta;
    double dtheta_x = (double)(fSmearAngle)?gRandom3->Gaus(0.,fBeamAngleRMS):0;
    double dtheta_y = (double)(fSmearAngle)?gRandom3->Gaus(0.,fBeamAngleRMS):0;
    double denergy  = (double)(fSmearEnergy)?gRandom3->Gaus(0.,fBeamEnergyRMS):0.;

    double px = p_out.P()*sin(theta+dtheta_x*microrad)*cos(p_out.Phi());
    double py = p_out.P()*sin(theta+dtheta_y*microrad)*sin(p_out.Phi());
    double pz = p_out.P()*(cos(theta)+denergy);

    if (p_out.Pz()<0) pz*=-1;

    double e  = sqrt(px*px+py*py+pz*pz+ProtonMassSQ);
    p_out.SetPxPyPzE(px,py,pz,e);
}

//--------------------------------------------------------------------------------------------------------------//

void PPSSim::CrossingAngleCorrection(TLorentzVector& p_out)
{
    double microrad = 1.e-6;
    double theta = p_out.Theta(); if (p_out.Pz()<0) theta=TMath::Pi()-theta;
    double dtheta_x = (double)((fSmearAngle)?gRandom3->Gaus(0.,fBeamAngleRMS):0+
            (p_out.Pz()>0)?fCrossingAngle:-fCrossingAngle);
    double dtheta_y = (double)(fSmearAngle)?gRandom3->Gaus(0.,fBeamAngleRMS):0;
    double denergy  = (double)(fSmearEnergy)?gRandom3->Gaus(0.,fBeamEnergyRMS):0.;

    double px = p_out.P()*(theta*cos(p_out.Phi())+dtheta_x*microrad);
    double py = p_out.P()*(theta*sin(p_out.Phi())+dtheta_y*microrad);
    double pz = p_out.P()*(cos(theta)+denergy);

    if (p_out.Pz()<0) pz*=-1;

    double e  = sqrt(px*px+py*py+pz*pz+ProtonMassSQ);
    p_out.SetPxPyPzE(px,py,pz,e);
}

//--------------------------------------------------------------------------------------------------------------//

void PPSSim::CrossingAngleCorrection(H_BeamParticle& p_out, const int Direction)
{
    // 
    // Remember: Hector  used X,Z inverted in ref. to CMS, but pz is always positive
    double partP = sqrt(pow(p_out.getE(),2)-ProtonMassSQ);
    double px = -Direction*partP*p_out.getTX()*urad;
    double py = partP*p_out.getTY()*urad;
    double pz = Direction*partP*cos(sqrt(pow(p_out.getTX(),2)+pow(p_out.getTY(),2))*urad);
    TLorentzVector p(px,py,pz,p_out.getE());
    CrossingAngleCorrection(p);
    p_out.set4Momentum(-Direction*p.Px(),p.Py(),Direction*p.Pz(),p.E());
    return;
}

//--------------------------------------------------------------------------------------------------------------//

TLorentzVector PPSSim::shoot(const double& t,const double& xi, const double& phi,const int Direction)
{
    long double energy=fBeamEnergy*(1.-xi);
    long double partP = sqrt((long double)(energy*energy-ProtonMassSQ));
    long double theta = acos((-t/2. - ProtonMassSQ + fBeamEnergy*energy)/(fBeamMomentum*partP)); // this angle is the scattering one
    long double px = partP*sin(theta)*cos((long double)phi)*Direction;
    long double py = partP*sin(theta)*sin((long double)phi);
    long double pz = partP*cos(theta)*Direction;
    return TLorentzVector((double)px,(double)py,(double)pz,(double)energy);
}

//--------------------------------------------------------------------------------------------------------------//

// void PPSSim::Propagate(H_BeamParticle* pbeam,int Direction) {
//     PPSSimData* arm = NULL;
//     H_BeamLine* beamline = NULL;
//     double startZ = -pbeam->getS(); // in the CMS ref. frame, in meters
//     double tcl4pos = 0;
//     double tcl5pos = 0;

//     std::map<int,PPSTrkDetector*> trakingDetectorMap;

//     if (Direction>0) {
//         arm = &(fSim->ArmF);
//         beamline=beamlineF;
//         tcl4pos=fTCL4Position2;
//         tcl5pos=fTCL5Position2;
//         trakingDetectorMap = fTrackingStationForward;
//     }
//     if (Direction<0) {
//         arm = &(fSim->ArmB);
//         beamline=beamlineB;
//         tcl4pos=fTCL4Position1;
//         tcl5pos=fTCL5Position1;
//         trakingDetectorMap = fTrackingStationBackward;
//     }
//     // Propagate until TCL4 and 5
//     if (tcl4pos>0) {
//         double beampos = (Direction<0)?fBeam1PosAtTCL4.first:fBeam2PosAtTCL4.first;
//         double beamrms = (Direction<0)?fBeam1RMSAtTCL4.first:fBeam2RMSAtTCL4.first;
//         pbeam->propagate(tcl4pos);
//         double xpos = -pbeam->getX()*um_to_mm;
//         arm->get_Track().set_XatTCL4(fabs(xpos-beampos)/beamrms);
//     }
//     if (tcl5pos>0) {
//         double beampos = (Direction<0)?fBeam1PosAtTCL5.first:fBeam2PosAtTCL5.first;
//         double beamrms = (Direction<0)?fBeam1RMSAtTCL5.first:fBeam2RMSAtTCL5.first;
//         pbeam->propagate(tcl5pos);double xpos=-pbeam->getX()*um_to_mm;
//         arm->get_Track().set_XatTCL5(fabs(xpos-beampos)/beamrms);
//     }

//     for(vector<int>::iterator dIt=fTrackingDetectorIDs.begin(); dIt!=fTrackingDetectorIDs.end(); ++dIt){
//         double hitZ = fTrackingStationZPositionMap[*dIt];
//         pbeam->propagate(abs(hitZ));

//         int stopped = (pbeam->stopped(beamline) && pbeam->getStoppingElement()->getS()<abs(hitZ))?1:0;
//         if (stopped) continue;
//         // uses mm for X,Y and m for Z in the PPS station
//         double hitX = -pbeam->getX()*um_to_mm;
//         double hitY = pbeam->getY()*um_to_mm;
//         //
//         // PPSStripDetector *stripDetector = (PPSStripDetector*)&trakingDetectorMap[*dIt];
//         // if(stripDetector==NULL) std::cout<<"NULL\n";
//         // std::cout<<"Direction = "<<Direction<<" hit x = "<<hitX<<" hit y = "<<hitY<<" hit z = "<<hitZ<<std::endl;
//         trakingDetectorMap[*dIt]->AddHit(TVector3(hitX,hitY,hitZ*m_to_mm));
//         // stripDetector->AddHit(TVector3(hitX,hitY,hitZ));
//         // ((PPSStripDetector)trakingDetectorMap[*dIt]).AddHit(TVector3(hitX,hitY,hitZ));
//     }

//     pbeam->propagate(fTrackerZPosition);
//     // std::cout<<"fTrackerZPosition "<<fTrackerZPosition<<std::endl;

//     int stopped = (pbeam->stopped(beamline) && pbeam->getStoppingElement()->getS()<fTrackerZPosition)?1:0;
//     if (!stopped){

//         // uses mm for X,Y and m for Z in the PPS station
//         double x1 = -pbeam->getX()*um_to_mm;
//         double y1 = pbeam->getY()*um_to_mm;
//         //
//         // std::cout<<"Old Direction = "<<Direction<<" hit x = "<<x1<<" hit y = "<<y1<<" hit z = "<<fTrackerZPosition<<std::endl;
//         arm->get_Track().set_HitDet1(x1,y1);
//         arm->AddHitTrk1(x1,y1);
//     }

//     pbeam->propagate(fTrackerZPosition+fTrackerLength);

//     stopped=(pbeam->stopped(beamline) && pbeam->getStoppingElement()->getS()<fTrackerZPosition+fTrackerLength)?1:0;
//     if (!stopped){
//         double x2 = -pbeam->getX()*um_to_mm;
//         double y2 = pbeam->getY()*um_to_mm;
        
//         // std::cout<<"Old Direction = "<<Direction<<" hit x = "<<x2<<" hit y = "<<y2<<" hit z = "<<fTrackerZPosition+fTrackerLength<<std::endl;
//         arm->get_Track().set_HitDet2(x2,y2);
//         arm->AddHitTrk2(x2,y2);
//     }

//     // Propagate until Time detector
//     pbeam->propagate(fToFZPosition);

//     double xt = -pbeam->getX()*um_to_mm;
//     double yt = pbeam->getY()*um_to_mm;
//     stopped=(pbeam->stopped(beamline) && pbeam->getStoppingElement()->getS()<fToFZPosition)?1:0;
//     if (!stopped){
//         //
//         double tof = (fToFZPosition-Direction*startZ)/c_light_ns;
//         arm->get_Track().set_HitToF(tof,xt,yt);
//         arm->AddHitToF(tof,xt,yt);
//     }
// }

//--------------------------------------------------------------------------------------------------------------//

void PPSSim::SmearVertexPosition(double& vtxX,double& vtxY, double& vtxZ)
{
    vtxX = fVtxMeanX;
    vtxY = fVtxMeanY;
    vtxZ = fVtxMeanZ;
    if (fSmearVertex) {
        vtxX=gRandom3->Gaus(fVtxMeanX,fVtxSigmaX); // in cm
        vtxY=gRandom3->Gaus(fVtxMeanY,fVtxSigmaY); // in cm
        vtxZ=gRandom3->Gaus(fVtxMeanZ,fVtxSigmaZ); // in cm
    }
}

//--------------------------------------------------------------------------------------------------------------//

//--------------------------------------------------------------------------------------------------------------//

double PPSSim::Minimum_t(const double& xi)
{
    double partE = fBeamEnergy*(1.- xi);
    double partP = sqrt(partE*partE-ProtonMassSQ);
    return -2.*(fBeamMomentum*partP-fBeamEnergy*partE+ProtonMassSQ);
}

//--------------------------------------------------------------------------------------------------------------//

void PPSSim::PrintParameters()
{
    edm::LogWarning("debug") << "Running with:\n"
        << "TrackerPosition     = " <<  fTrackerZPosition << "\n"
        << "TrackerLength       = " <<  fTrackerLength << "\n"
        << "TrackerLength       = " <<  fTrackerLength << "\n"
        << "BeamLineLength      = " <<  fBeamLineLength << "\n"
        << "SmearVertex         = " <<  fSmearVertex << "\n"
        << "VtxMeanX            = " <<  fVtxMeanX << "\n"
        << "VtxMeanY            = " <<  fVtxMeanY << "\n"
        << "VtxMeanZ            = " <<  fVtxMeanZ << "\n"
        << "VtxSigmaX           = " <<  fVtxSigmaX << "\n"
        << "VtxSigmaY           = " <<  fVtxSigmaY << "\n"
        << "VtxSigmaZ           = " <<  fVtxSigmaZ << "\n"
        << "VtxMeanZ            = " <<  fVtxMeanZ << "\n"
        << "VtxSigmaX           = " <<  fVtxSigmaX << "\n"
        << "VtxSigmaY           = " <<  fVtxSigmaY << "\n"
        << "VtxSigmaZ           = " <<  fVtxSigmaZ << "\n"
        << "SmearHit            = " <<  fSmearHit << "\n"
        << "HitSigmaX           = " <<  fHitSigmaX << "\n"
        << "HitSigmaY           = " <<  fHitSigmaY << "\n"
        << "HitSigmaZ           = " <<  fHitSigmaZ << "\n"
        << "TimeSigma           = " <<  fTimeSigma << "\n"
        << "SimBeam             = " <<  fSimBeam   << "\n"
        << "PhiMin              = " <<  fPhiMin    << "\n"
        << "PhiMax              = " <<  fPhiMax    << "\n"
        << "EtaMin              = " <<  fEtaMin    << "\n"
        << "MomentumMin         = " <<  fMomentumMin << "\n"
        << "CrossAngleCorr      = " <<  fCrossAngleCorr << "\n"
        << "KickersOFF          = " <<  fKickersOFF << "\n"
        << "Central Mass        = " <<  fCentralMass << " +- " << fCentralMassErr << "\n"
        << "TrackImpactParameterCut  = " << fTrackImpactParameterCut << "\n"
        << "MinThetaXatDet1     = " <<fMinThetaXatDet1 << "\n"
        << "MaxThetaXatDet1     = " <<fMaxThetaXatDet1 << "\n"
        << "MinThetaYatDet1     = " <<fMinThetaYatDet1 << "\n"
        << "MaxThetaYatDet1     = " <<fMaxThetaYatDet1 << "\n";

    edm::LogWarning("debug") << "\nBeam1 station parameters\n";
    edm::LogWarning("debug") << "Name\tZ\tXCenter\tXSigma\tInsertionSigma  \n";
    for(vector<string>::iterator vIt=fBeam1Objects.begin(); vIt!=fBeam1Objects.end(); ++vIt){
        edm::LogWarning("debug") << *vIt << "\t" << fBeam1StationPositions[*vIt] << "\t" << fBeam1CenterAtStation[*vIt] << "\t" << fBeam1SigmaAtStation[*vIt] << "\t" << fBeam1StationSigmaInsertion[*vIt];
    }

    edm::LogWarning("debug") << "\nBeam2 station parameters\n";
    edm::LogWarning("debug") << "Name\tZ\tXCenter\tXSigma\tInsertionSigma  \n";
    for(vector<string>::iterator vIt=fBeam2Objects.begin(); vIt!=fBeam2Objects.end(); ++vIt){
        edm::LogWarning("debug") << *vIt << "\t" << fBeam2StationPositions[*vIt] << "\t" << fBeam2CenterAtStation[*vIt] << "\t" << fBeam2SigmaAtStation[*vIt] << "\t" << fBeam2StationSigmaInsertion[*vIt];
    }

}

//--------------------------------------------------------------------------------------------------------------//

std::map<string,TH2F*> PPSSim::GenBeamProfile(int  direction)
{
    float beamp_w = 20.0;//beam pipe width
    vector<string> beamObjects;
    vector<string> beamTCLNames;


    if(direction < 0){
        beamObjects = fBeam1Objects;
        beamTCLNames=fBeam1TCLNames;
    }
    if(direction > 0){
        beamObjects = fBeam2Objects;
        beamTCLNames=fBeam2TCLNames;
    }

    //double protonBeamMomentum = TMath::Sqrt(ProtonMassSQ-fBeamEnergy*fBeamEnergy);
    int   nbins = 500;
    TH2F* beamprofile = (TH2F*)gDirectory->FindObject("BeamProfile");
    if (beamprofile) delete beamprofile;
    std::map<string,TH2F*> beamProfilesHistograms;


    for(vector<string>::iterator vIt=beamObjects.begin(); vIt!=beamObjects.end(); ++vIt){
        beamProfilesHistograms[*vIt] = new TH2F(Form("Beam%iProfile_%s",direction,vIt->data()),Form("Beam %i Profile at %s",direction,vIt->data()),nbins,-beamp_w,beamp_w,nbins,-beamp_w,beamp_w);
    }

    for(vector<string>::iterator vIt=beamTCLNames.begin(); vIt!=beamTCLNames.end(); ++vIt){
        for(int n=0;n<10000;n++) {
        // for(int n=0;n<100000;n++) {
            vector5Def protonKinematicAndEnergy;
            TVector3 protonBeamVertex(gRandom->Gaus(0.,fVtxMeanX*mm_to_m),gRandom->Gaus(0.,fVtxMeanY*mm_to_m),0.);
            // TVector3 protonBeamVertex(0.,0.,0.);
            TVector3 protonBeamAngle(gRandom->Gaus(0.,fBeamAngleRMS*urad_to_rad),gRandom->Gaus(0.,fBeamAngleRMS*urad_to_rad),0.);
            // TVector3 protonBeamAngle(0.,0.,0.);
            vector4Def protonKinematic(protonBeamVertex.X(),protonBeamAngle.X(),protonBeamVertex.Y(),protonBeamAngle.Y());
            protonKinematicAndEnergy = std::pair<double,vector4Def>(0.,protonKinematic);
            map<string,vector5Def> protonKinematicsAtStations = PPSSim::PropagateParticle(direction, protonKinematicAndEnergy);
            if(protonKinematicsAtStations[*vIt].first < 0.) continue;
            vector4Def protonKinematicsAtStation = protonKinematicsAtStations[*vIt].second;
            beamProfilesHistograms[*vIt]->Fill(protonKinematicsAtStation(0)*m_to_mm,protonKinematicsAtStation(2)*m_to_mm);
        }
        
        if(direction < 0){
            fBeam1CenterAtStation[*vIt] = beamProfilesHistograms[*vIt]->GetMean(1);
            fBeam1SigmaAtStation[*vIt] = beamProfilesHistograms[*vIt]->GetRMS(1);
        }
        if(direction > 0){
            fBeam2CenterAtStation[*vIt] = beamProfilesHistograms[*vIt]->GetMean(1);
            fBeam2SigmaAtStation[*vIt] = beamProfilesHistograms[*vIt]->GetRMS(1);
        }

    }

    for(int n=0;n<10000;n++) {
    // for(int n=0;n<100000;n++) {
        vector5Def protonKinematicAndEnergy;
        TVector3 protonBeamVertex(gRandom->Gaus(0.,fVtxMeanX*mm_to_m),gRandom->Gaus(0.,fVtxMeanY*mm_to_m),0.);
        // TVector3 protonBeamVertex(0.,0.,0.);
        TVector3 protonBeamAngle(gRandom->Gaus(0.,fBeamAngleRMS*urad_to_rad),gRandom->Gaus(0.,fBeamAngleRMS*urad_to_rad),0.);
        // TVector3 protonBeamAngle(0.,0.,0.);
        vector4Def protonKinematic(protonBeamVertex.X(),protonBeamAngle.X(),protonBeamVertex.Y(),protonBeamAngle.Y());
        protonKinematicAndEnergy = std::pair<double,vector4Def>(0.,protonKinematic);
        map<string,vector5Def> protonKinematicsAtStations = PPSSim::PropagateParticle(direction, protonKinematicAndEnergy);
        for(map<string,vector5Def>::iterator mIt=protonKinematicsAtStations.begin(); mIt!=protonKinematicsAtStations.end(); ++mIt){
            if(find(beamTCLNames.begin(),beamTCLNames.end(),mIt->first) != beamTCLNames.end()) continue;
            if(mIt->second.first < 0.) continue;
            vector4Def protonKinematicsAtStation = mIt->second.second;
            beamProfilesHistograms[mIt->first]->Fill(protonKinematicsAtStation(0)*m_to_mm,protonKinematicsAtStation(2)*m_to_mm);
        }
    }


    for(vector<string>::iterator vIt=beamObjects.begin(); vIt!=beamObjects.end(); ++vIt){
        if(find(beamTCLNames.begin(),beamTCLNames.end(),*vIt) != beamTCLNames.end()) continue;
        if(direction < 0){
            fBeam1CenterAtStation[*vIt] = beamProfilesHistograms[*vIt]->GetMean(1);
            fBeam1SigmaAtStation[*vIt] = beamProfilesHistograms[*vIt]->GetRMS(1);
        }
        if(direction > 0){
            fBeam2CenterAtStation[*vIt] = beamProfilesHistograms[*vIt]->GetMean(1);
            fBeam2SigmaAtStation[*vIt] = beamProfilesHistograms[*vIt]->GetRMS(1);
        }
    }

    return beamProfilesHistograms;
}

//--------------------------------------------------------------------------------------------------------------//

map<string,vector5Def> PPSSim::PropagateParticle(int direction, vector5Def protonKinematic){
    
    vector<string> beamObjects;
    map<string,double> beamCenterAtTCLMap;
    map<string,double> beamSigmaAtTCLMap;
    vector<string> beamTCLNames;
    std::map< int, std::map<string,std::pair<vector5Def, matrix44Def> > > beamParameterCollection;
    std::map<string,double> beamCenterAtStation;
    std::map<string,double> beamSigmaAtStation;
    std::map<string,double> beamStationSigmaInsertion;

    if(direction < 0){
        beamObjects = fBeam1Objects;
        beamTCLNames=fBeam1TCLNames;
        beamParameterCollection = fBeam1ParameterCollection;
        beamCenterAtStation = fBeam1CenterAtStation;
        beamSigmaAtStation = fBeam1SigmaAtStation;
        beamStationSigmaInsertion = fBeam1StationSigmaInsertion;
    }
    if(direction > 0){
        beamObjects = fBeam2Objects;
        beamTCLNames=fBeam2TCLNames;
        beamParameterCollection =  fBeam2ParameterCollection;
        beamCenterAtStation = fBeam2CenterAtStation;
        beamSigmaAtStation = fBeam2SigmaAtStation;
        beamStationSigmaInsertion = fBeam2StationSigmaInsertion;
    }

    map<string,vector5Def> protonKinematicsAtStations;
    for(vector<string>::iterator vIt=beamObjects.begin(); vIt!=beamObjects.end(); ++vIt){
        vector5Def defaultKinematics = std::pair<double,vector4Def>(-1.,vector4Def(0.,0.,0.,0.));
        protonKinematicsAtStations[*vIt] = defaultKinematics;
    }

    double protonXi = protonKinematic.first;

    for(vector<string>::iterator vIt=beamObjects.begin(); vIt!=beamObjects.end(); ++vIt){
        std::pair<vector5Def, matrix44Def> opticForStationProton = beamParameterCollection[GetXiClassNumber(direction, protonXi)][*vIt];
        // std::pair<vector5Def, matrix44Def> opticForStationBeam = beamParameterCollection[GetXiClassNumber(direction, 0.)][*vIt];
        vector4Def protonKinematicCalculated =  opticForStationProton.first.second + opticForStationProton.second * protonKinematic.second;
        // if(*vIt=="XRPH.E6R5.B1") std::cout<< opticForStationProton.first.second(0)<<"\t"<<protonKinematic.second(0)<<"\n";
        // vector4Def protonKinematicCalculated = - opticForStationBeam.first.second + opticForStationProton.first.second + opticForStationProton.second * protonKinematic.second;
        // if(*vIt==beamObjects.back() && direction==-1) std::cout<<(opticForStationBeam.first.second)(0)<<std::endl;
        protonKinematicsAtStations[*vIt] = std::pair<double,vector4Def>(protonXi,protonKinematicCalculated);
        // std::cout<<*vIt<<std::endl;
        if(find(beamTCLNames.begin(),beamTCLNames.end(),*vIt) != beamTCLNames.end() && 
            (protonKinematicCalculated(0)*m_to_mm<beamCenterAtStation[*vIt]-beamSigmaAtStation[*vIt]*beamStationSigmaInsertion[*vIt] ||
            protonKinematicCalculated(0)*m_to_mm>beamCenterAtStation[*vIt]+beamSigmaAtStation[*vIt]*beamStationSigmaInsertion[*vIt] ))  return protonKinematicsAtStations;
    }

    return protonKinematicsAtStations;

}

//--------------------------------------------------------------------------------------------------------------//

std::map<string,std::pair<vector5Def, matrix44Def> > PPSSim::GetMatrixFromTwissFile(string inputTwissFileName, std::vector<string> beamObjects, double &protonMomentumXi){

    std::map<string,std::pair<vector5Def, matrix44Def> > matrixMap;

    FILE *inputTwissFile;
    inputTwissFile = fopen(inputTwissFileName.data(),"r");

    if(inputTwissFile == NULL){
        throw cms::Exception("PPSSim.h") << "Twiss file " << inputTwissFileName << " does not exist";
    }

    map<string, int> nameVsColumn;

    double beamMomentum=-1.;
    double protonMomentumLoss=-1.;

    // get positions in the columns
    while(!feof(inputTwissFile)){
        char cInput[1000];
        fgets(cInput,1000,inputTwissFile);
        if(cInput[0] =='@'){
            if(strncmp(cInput,"@ PC",4)==0){
                string headerString(cInput);
                stringstream headerStream;
                string dummy;
                headerStream<<headerString;
                headerStream>>dummy>>dummy>>dummy>>beamMomentum;

            }
            if(strncmp(cInput,"@ DELTAP'",8)==0){
                string headerString(cInput);
                stringstream headerStream;
                string dummy;
                headerStream<<headerString;
                headerStream>>dummy>>dummy>>dummy>>protonMomentumXi;

            }
            continue;
        }
        string columnNamesString(cInput);
        stringstream columnNameStream;
        columnNameStream<<columnNamesString;
        string columnName;
        int counter = -2; //starting with -1 because there is al "*" as a first column name
        while((columnNameStream >> columnName)){
            // std::cout<<columnName<<"   "<<counter<<std::endl;
            nameVsColumn[columnName]=counter++;
        }
        fgets(cInput,1000,inputTwissFile);
        break;
    }    

    // protonMomentumXi=(beamMomentum-protonMomentumLoss)/beamMomentum;
    protonMomentumXi = -protonMomentumXi;
    double protonPz = beamMomentum*(1-protonMomentumXi);

    while(!feof(inputTwissFile)){
        char cInput[1000];
        fgets(cInput,1000,inputTwissFile);
        if(feof(inputTwissFile)) break;
        string columnNamesString(cInput);
        stringstream columnNameStream;
        columnNameStream<<columnNamesString;
        string valueString;
        columnNameStream >> valueString;
        bool foundNeededStation=false;
        for(std::vector<string>::iterator vIt=beamObjects.begin(); vIt!=beamObjects.end(); ++vIt){
            if(valueString.find(*vIt)==1){
                foundNeededStation=true;
                break;
            }
        }
        if(!foundNeededStation) continue;
        double columnValue;
        std::vector<double> columnValues;
        while((columnNameStream >> columnValue)) {
            columnValues.push_back(columnValue);
        }

        vector4Def beamPositionAndAngle;
        beamPositionAndAngle(0) = columnValues[nameVsColumn["X"]];
        beamPositionAndAngle(1) = TMath::ATan(columnValues[nameVsColumn["PX"]]/beamMomentum);
        beamPositionAndAngle(2) = columnValues[nameVsColumn["Y"]];
        beamPositionAndAngle(3) = TMath::ATan(columnValues[nameVsColumn["PY"]]/beamMomentum);
        vector5Def beamPositionAtElement = vector5Def(columnValues[nameVsColumn["S"]],beamPositionAndAngle);
        // for(unsigned int i=0; i<columnValues.size(); ++i){
        //     std::cout<< i << "  " << columnValues[i]<<std::endl;
        // }
        // std::cout<<"S"<< "  " <<nameVsColumn["S"] << "  " << columnValues[nameVsColumn["S"]]<<std::endl;
        // std::cout<<"X"<< "  " <<nameVsColumn["X"] << "  " << columnValues[nameVsColumn["X"]]<<std::endl;
        // std::cout<<"Y"<< "  " <<nameVsColumn["Y"] << "  " << columnValues[nameVsColumn["Y"]]<<std::endl;
        matrix44Def transportMatrix;
        for(int r=1; r<=4; ++r){
            for(int c=1; c<=4; ++c){
                transportMatrix(r-1,c-1) = columnValues[nameVsColumn[Form("RE%i%i",r,c)]];
                //std::cout<<nameVsColumn[Form("R%i%i",r,c)]<<std::endl;
            }
        }

        std::pair<vector5Def, matrix44Def> elementParameter(beamPositionAtElement,transportMatrix);
        matrixMap[valueString.substr(1,valueString.length()-2)] = elementParameter;
        // if("XRPH.E6R5.B1"==valueString.substr(1,valueString.length()-2)) std::cout<<inputTwissFileName<<"\t"<<std::setprecision(10)<<columnValues[nameVsColumn["X"]]<<std::endl;
    }

    return matrixMap;
 }

 //--------------------------------------------------------------------------------------------------------------//

void PPSSim::LoadTwissMatrices(){

    fBeam1ParameterCollection.clear();
    fBeam1XiClasses.clear();

    std::vector<double> twissFileBeam1Xis;
    for (unsigned int i = 0; i < fBeam1Files.size(); ++i){
        double protonMomentumXi;
        fBeam1ParameterCollection[i] = GetMatrixFromTwissFile(fBeam1FilePath + "/" + fBeam1Files.at(i),fBeam1Objects, protonMomentumXi);
        // if(i==0) protonMomentumXi = 0;
        twissFileBeam1Xis.push_back(protonMomentumXi);
        // std::cout<<"Class number " << i << "\n";
        // for(int r=1; r<=4; ++r){
        //     for(int c=1; c<=4; ++c){
        //         std::cout<<fBeam1ParameterCollection[i]["TCL.4R5.B1"].second(r-1,c-1)<<"\t";
        //     }
        //     std::cout<<"\n";
        // }
        // std::cout<< std::setprecision(10) <<protonMomentumXi<<std::endl;
    }

    //std::cout<< fBeam1ParameterCollection.size() <<std::endl;

    // std::cout<<"beam1:\n";
    for (unsigned int i = 0; i < twissFileBeam1Xis.size(); ++i){
        //double protonMomentumXi = twissFileBeam1Xis.at(i);
        double xiMin = -1;
        double xiMax = -1;

        if(i==0) xiMin=-1.;
        else xiMin = (twissFileBeam1Xis.at(i-1) + twissFileBeam1Xis.at(i))/2.;
        if(i==fBeam1Files.size()-1) xiMax = 1.;
        else xiMax = (twissFileBeam1Xis.at(i) + twissFileBeam1Xis.at(i+1))/2.;
        // std::cout<<"class "<<i<<" -> "<<std::setprecision(7)<<xiMin<<" : "<<xiMax<<std::endl;

        fBeam1XiClasses[i] = pair<double,double> (xiMin,xiMax);
    }

    //std::cout<< fBeam1XiClasses.size() <<std::endl;

    std::vector<double> twissFileBeam2Xis;
    for (unsigned int i = 0; i < fBeam2Files.size(); ++i){
        double protonMomentumXi;
        fBeam2ParameterCollection[i] = GetMatrixFromTwissFile(fBeam2FilePath + "/" + fBeam2Files.at(i),fBeam2Objects, protonMomentumXi);
        // if(i==0) protonMomentumXi = 0;
        twissFileBeam2Xis.push_back(protonMomentumXi);
        // std::cout<< std::setprecision(10) <<protonMomentumXi<<std::endl;
    }

    //std::cout<< fBeam2ParameterCollection.size() <<std::endl;

    // std::cout<<"beam2:\n";
    for (unsigned int i = 0; i < twissFileBeam2Xis.size(); ++i){
        //double protonMomentumXi = twissFileBeam2Xis.at(i);
        double xiMin = -1;
        double xiMax = -1;

        if(i==0) xiMin=-1.;
        else xiMin = (twissFileBeam2Xis.at(i-1) + twissFileBeam2Xis.at(i))/2.;
        if(i==fBeam2Files.size()-1) xiMax = 1.;
        else xiMax = (twissFileBeam2Xis.at(i) + twissFileBeam2Xis.at(i+1))/2.;
        // std::cout<<"class "<<i<<" -> "<<std::setprecision(7)<<xiMin<<" : "<<xiMax<<std::endl;

        fBeam2XiClasses[i] = pair<double,double> (xiMin,xiMax);
    }

    //std::cout<< fBeam2XiClasses.size() <<std::endl;

    return;

}

int PPSSim::GetXiClassNumber(int direction, double protonXi){

    std::map< int, pair<double, double> > beamXiClasses;

    if(direction<0) beamXiClasses = fBeam1XiClasses;
    if(direction>0) beamXiClasses = fBeam2XiClasses;

    int xiClassIterator=-1;
    for(std::map< int, pair<double, double> >::iterator mIt=beamXiClasses.begin(); mIt!=beamXiClasses.end(); ++mIt){
        xiClassIterator = mIt->first;
        if(mIt->second.second >= protonXi) break;
    }
    // std::cout<<xiClassIterator<<" ";
    return xiClassIterator;
 }