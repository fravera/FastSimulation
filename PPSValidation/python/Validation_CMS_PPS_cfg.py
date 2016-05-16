import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'file:GluGluTo2Jets_M_300_2000_13TeV_exhume_cff_py_GEN_SIM_RECOBEFMIX_DIGI_RECO_NoPU.root'
        #'file:GluGluTo2Jets_M_300_2000_13TeV_exhume_cff_py_GEN_SIM_RECOBEFMIX_DIGI_RECO_PU.root'
        #'file:TTbar_13TeV_TuneCUETP8M1_cfi_GEN_SIM_RECOBEFMIX_DIGI_RECO_PU.root'
        #'file:GluGluTo2Jets_M_300_2000_13TeV_exhume_GEN_SIM_RECOBEFMIX_DIGI_RECO.root'
        #'file:GluGluTo2Jets_M_300_2000_13TeV_exhume_GEN_SIM_RECOBEFMIX_DIGI_RECO_NoPU.root'
        #'file:GluGluTo2Jets_M_300_2000_13TeV_exhume_GEN_SIM_RECOBEFMIX_DIGI_RECO_PU.root'
        #'file:GluGluTo2Jets_M_100_7TeV_exhume_cff_py_GEN_SIM_RECOBEFMIX_DIGI_RECO.root',
        'file:GluGluTo2Jets_M_300_2000_13TeV_exhume_RECO_PU.root'
    )
)

process.demo = cms.EDAnalyzer('Validation'
		, jets          = cms.InputTag('ak4PFJets')
		, vertices      = cms.InputTag('offlinePrimaryVertices')
		, muons         = cms.InputTag('muons')
		, electrons     = cms.InputTag('gedGsfElectrons')
		, ppsGen        = cms.InputTag('ppssim:PPSGen')
		, ppsSim        = cms.InputTag('ppssim:PPSSim')
		, ppsReco       = cms.InputTag('ppssim:PPSReco')
        , ToFCellWidth  = cms.untracked.vdouble(0.81, 0.91, 1.02, 1.16, 1.75, 2.35, 4.2, 4.2) # tof cell width in mm #move to vector - diamond geometry
        , BeamSizeAtToF = cms.double(0.113) # beam sigma (X) at timing station in mm
        , ToFInsertion  = cms.double(15), # Number of sigmas (X) from the beam for the tof
)

process.TFileService = cms.Service("TFileService",
                #fileName = cms.string('histo_validation_FastSim_CTPPS_NoPU.root')
                #fileName = cms.string('histo_validation_FastSim_CTPPS_PU.root')
                #fileName = cms.string('histo_validation_FastSim_CTPPS_ttbar.root')
                #fileName = cms.string('histo_validation_GluGlu_NoPU.root')
                #fileName = cms.string('histo_validation_GluGlu_PU.root')
                fileName = cms.string('histo_validation.root')
        )

process.p = cms.Path(process.demo)
