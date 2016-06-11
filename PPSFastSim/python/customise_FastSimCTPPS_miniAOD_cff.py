import FWCore.ParameterSet.Config as cms

def customise_pu_protons_ctpps(process):
	
	# process=customise(process)

	# process.mix.mixObjects.mixHepMC.makeCrossingFrame = cms.untracked.bool(True)

	# process.ppssim.UseHepMCProducer = cms.untracked.bool(False)
	# PU gen particle       
	# process.genParticlesPU = cms.EDProducer("GenParticleProducer",
	# saveBarCodes = cms.untracked.bool(True),
	# mix = cms.string("mix"),
	# abortOnUnknownPDGCode = cms.untracked.bool(False),
	# useCrossingFrame = cms.untracked.bool(True)
	# )
	process.genProtonsPU = cms.EDFilter("GenParticleSelector",filter = cms.bool(False),src = cms.InputTag("genParticlesPU"),cut = cms.string(''))
	process.MINIAODSIMoutput.outputCommands.extend(cms.untracked.vstring('keep PPSGenDataPPSSpectrometer_*_*_*','keep PPSSimDataPPSSpectrometer_*_*_*','keep PPSRecoDataPPSSpectrometer_*_*_*'))## add CTPPS eventContent
	# process.genProtonsPU.cut = 'status = 1 & pdgId == 2212 & abs(pz) >= %f' % ( 0.5*13000./2.0)
	outputCommandsPU = [ 'keep *_genParticlesPU_*_*', 'keep *_genProtonsPU_*_*']
	
	# process.ppssim.genSource = cms.InputTag("genProtonsPU") # for Pile-up events

	#process.digitisation_step.replace(process.pdigi_valid, process.pdigi_valid * process.genParticlesPU * process.genProtonsPU)
	# process.digitisation_step.replace(process.ppssim, process.genParticlesPU * process.genProtonsPU * process.ppssim)

	return (process)
