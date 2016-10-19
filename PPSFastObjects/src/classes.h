#include "FastSimulation/PPSFastObjects/interface/PPSSpectrometer.h"
#include "FastSimulation/PPSFastObjects/interface/PPSGenData.h"
#include "FastSimulation/PPSFastObjects/interface/PPSRecoData.h"
#include "FastSimulation/PPSFastObjects/interface/PPSGenParticles.h"
#include "FastSimulation/PPSFastObjects/interface/PPSSimTracks.h"
#include "FastSimulation/PPSFastObjects/interface/PPSRecoTracks.h"
#include "FastSimulation/PPSFastObjects/interface/PPSRecoVertex.h"
#include "FastSimulation/PPSFastObjects/interface/PPSGenVertex.h"
#include "FastSimulation/PPSFastObjects/interface/PPSVertex.h"
#include "DataFormats/Common/interface/Wrapper.h"

#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/CTPPSDigi/interface/TotemRPDigi.h"

namespace { struct dictionaryPPS {
    PPSVertex                    vtx;
    edm::Wrapper<PPSVertex>      vtxs;
    PPSGenVertex                 genvtx;
    edm::Wrapper<PPSGenVertex>   genWvtx;
    PPSRecoVertex                recovtx;
    edm::Wrapper<PPSRecoVertex>  recoWvtx;
    PPSTrackerHit                trkhit;
    PPSTrackerHits               trkhits;
    PPSToFHit                    tofhit;
    PPSToFHits                   tofhits;
    PPSBaseTrack                 basetrk;
    PPSSimTrack                  simtrk;
    PPSSimTracks                 simtrks;
    PPSRecoTrack                 recotrk;
    PPSRecoTracks                recotrks;
    PPSBaseData                  basedata;
    PPSGenData                   genData;
    PPSSimData                   PpsSimData;
    PPSRecoData                  PpsRecoData;
    PPSSpectrometer<PPSGenData>  genPpsSpectrometer;
    PPSSpectrometer<PPSSimData>  simPpsSpectrometer;
    PPSSpectrometer<PPSRecoData> recoPpsSpectrometer;
    edm::Wrapper<PPSSpectrometer<PPSGenData> > gdummy;
    edm::Wrapper<PPSSpectrometer<PPSSimData> > sdummy;
    edm::Wrapper<PPSSpectrometer<PPSRecoData> > rdummy;

    // TotemRPDigi rp_str_dig;
    // edm::DetSet<TotemRPDigi> ds_rp_str_dig;
    // std::vector<TotemRPDigi> vec_rp_str_dig;
    // edm::DetSetVector<TotemRPDigi> dsv_rp_str_dig;
    // std::vector<edm::DetSet<TotemRPDigi> > vec_ds_rp_str_dig;
    // edm::Wrapper<edm::DetSet<TotemRPDigi> > wds_rp_str_dig;
    // edm::Wrapper<edm::DetSetVector<TotemRPDigi> > wdsv_rp_str_dig;
};
}
