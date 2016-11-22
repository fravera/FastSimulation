import FWCore.ParameterSet.Config as cms
import math

#common stuff here
det1 = 203.827   #position of first tracker detector
det2 = 212.550
tof  = 215.7   #position of time of flight detector
trklen = det2-det1
hit_smear= True
Ang_smear= True
E_smear = True
vtx_smear = False
run_with_CR = True
useCR = True
ecms = 13000.
#det1xoffset = 0.25 + 0. #distance from sensitive edge + RP security margin
#det2xoffset = 0.25 + 0. #distance from sensitive edge + RP security margin
tofxoffset  = 0.6  + 0. #distance from sensitive edge + RP security margin
phi_min = -math.pi
phi_max =  math.pi

def customise_ctpps(process):
######## CTPPS
    process.load('PhysicsTools.HepMCCandAlgos.genParticles_cfi')
    process.load("FWCore.MessageLogger.MessageLogger_cfi")
    process.MessageLogger = cms.Service("MessageLogger")

    VertexX = cms.double(0.)
    VertexY = cms.double(0.)
    VertexZ = cms.double(0.)

    if hasattr(process,"VtxSmeared") and hasattr(process.VtxSmeared,"X0"):
        VertexX = process.VtxSmeared.X0
        VertexY = process.VtxSmeared.Y0
        VertexZ = process.VtxSmeared.Z0

    if hasattr(process,"VtxSmeared") and hasattr(process.VtxSmeared,"MeanX"):
        VertexX = process.VtxSmeared.MeanX
        VertexY = process.VtxSmeared.MeanY
        VertexZ = process.VtxSmeared.MeanZ
    

    print 'Setting CT-PPS FastSim'
    ppssim_beam_options = cms.PSet(
                        Verbosity = cms.untracked.int32(0),
               	        # Beam1File = cms.FileInPath("FastSimulation/PPSFastSim/data/LHCB1_Beta0.40_6.5TeV_CR185_2016.tfs"),
                        # Beam2File = cms.FileInPath("FastSimulation/PPSFastSim/data/LHCB2_Beta0.40_6.5TeV_CR185_2016.tfs"),
                        Beam1FilePath  = cms.string("data/"),
                        Beam2FilePath  = cms.string("data/"),
                        Beam1Files     = cms.vstring(
                            "twiss_lhcb1.tfs0000","twiss_lhcb1.tfs0001","twiss_lhcb1.tfs0002","twiss_lhcb1.tfs0003","twiss_lhcb1.tfs0004",
                            "twiss_lhcb1.tfs0005","twiss_lhcb1.tfs0006","twiss_lhcb1.tfs0007","twiss_lhcb1.tfs0008","twiss_lhcb1.tfs0009","twiss_lhcb1.tfs0010","twiss_lhcb1.tfs0011","twiss_lhcb1.tfs0012","twiss_lhcb1.tfs0013","twiss_lhcb1.tfs0014",
                            "twiss_lhcb1.tfs0015","twiss_lhcb1.tfs0016","twiss_lhcb1.tfs0017","twiss_lhcb1.tfs0018","twiss_lhcb1.tfs0019","twiss_lhcb1.tfs0020","twiss_lhcb1.tfs0021","twiss_lhcb1.tfs0022","twiss_lhcb1.tfs0023","twiss_lhcb1.tfs0024",
                            "twiss_lhcb1.tfs0025","twiss_lhcb1.tfs0026","twiss_lhcb1.tfs0027","twiss_lhcb1.tfs0028","twiss_lhcb1.tfs0029","twiss_lhcb1.tfs0030","twiss_lhcb1.tfs0031","twiss_lhcb1.tfs0032","twiss_lhcb1.tfs0033","twiss_lhcb1.tfs0034",
                            "twiss_lhcb1.tfs0035","twiss_lhcb1.tfs0036","twiss_lhcb1.tfs0037","twiss_lhcb1.tfs0038","twiss_lhcb1.tfs0039","twiss_lhcb1.tfs0040","twiss_lhcb1.tfs0041","twiss_lhcb1.tfs0042","twiss_lhcb1.tfs0043","twiss_lhcb1.tfs0044",
                            "twiss_lhcb1.tfs0045","twiss_lhcb1.tfs0046","twiss_lhcb1.tfs0047","twiss_lhcb1.tfs0048","twiss_lhcb1.tfs0049","twiss_lhcb1.tfs0050","twiss_lhcb1.tfs0051","twiss_lhcb1.tfs0052","twiss_lhcb1.tfs0053","twiss_lhcb1.tfs0054",
                            "twiss_lhcb1.tfs0055","twiss_lhcb1.tfs0056","twiss_lhcb1.tfs0057","twiss_lhcb1.tfs0058","twiss_lhcb1.tfs0059","twiss_lhcb1.tfs0060","twiss_lhcb1.tfs0061","twiss_lhcb1.tfs0062","twiss_lhcb1.tfs0063","twiss_lhcb1.tfs0064",
                            "twiss_lhcb1.tfs0065","twiss_lhcb1.tfs0066","twiss_lhcb1.tfs0067","twiss_lhcb1.tfs0068","twiss_lhcb1.tfs0069","twiss_lhcb1.tfs0070","twiss_lhcb1.tfs0071","twiss_lhcb1.tfs0072","twiss_lhcb1.tfs0073","twiss_lhcb1.tfs0074",
                            "twiss_lhcb1.tfs0075","twiss_lhcb1.tfs0076","twiss_lhcb1.tfs0077","twiss_lhcb1.tfs0078","twiss_lhcb1.tfs0079","twiss_lhcb1.tfs0080","twiss_lhcb1.tfs0081","twiss_lhcb1.tfs0082","twiss_lhcb1.tfs0083","twiss_lhcb1.tfs0084",
                            "twiss_lhcb1.tfs0085","twiss_lhcb1.tfs0086","twiss_lhcb1.tfs0087","twiss_lhcb1.tfs0088","twiss_lhcb1.tfs0089","twiss_lhcb1.tfs0090","twiss_lhcb1.tfs0091","twiss_lhcb1.tfs0092","twiss_lhcb1.tfs0093","twiss_lhcb1.tfs0094",
                            "twiss_lhcb1.tfs0095","twiss_lhcb1.tfs0096","twiss_lhcb1.tfs0097","twiss_lhcb1.tfs0098","twiss_lhcb1.tfs0099","twiss_lhcb1.tfs0100","twiss_lhcb1.tfs0101","twiss_lhcb1.tfs0102","twiss_lhcb1.tfs0103","twiss_lhcb1.tfs0104",
                            "twiss_lhcb1.tfs0105","twiss_lhcb1.tfs0106","twiss_lhcb1.tfs0107","twiss_lhcb1.tfs0108","twiss_lhcb1.tfs0109","twiss_lhcb1.tfs0110","twiss_lhcb1.tfs0111","twiss_lhcb1.tfs0112","twiss_lhcb1.tfs0113","twiss_lhcb1.tfs0114",
                            "twiss_lhcb1.tfs0115","twiss_lhcb1.tfs0116","twiss_lhcb1.tfs0117","twiss_lhcb1.tfs0118","twiss_lhcb1.tfs0119","twiss_lhcb1.tfs0120","twiss_lhcb1.tfs0121","twiss_lhcb1.tfs0122","twiss_lhcb1.tfs0123","twiss_lhcb1.tfs0124",
                            "twiss_lhcb1.tfs0125","twiss_lhcb1.tfs0126","twiss_lhcb1.tfs0127","twiss_lhcb1.tfs0128","twiss_lhcb1.tfs0129","twiss_lhcb1.tfs0130","twiss_lhcb1.tfs0131","twiss_lhcb1.tfs0132","twiss_lhcb1.tfs0133","twiss_lhcb1.tfs0134",
                            "twiss_lhcb1.tfs0135","twiss_lhcb1.tfs0136","twiss_lhcb1.tfs0137","twiss_lhcb1.tfs0138","twiss_lhcb1.tfs0139","twiss_lhcb1.tfs0140","twiss_lhcb1.tfs0141","twiss_lhcb1.tfs0142","twiss_lhcb1.tfs0143","twiss_lhcb1.tfs0144",
                            "twiss_lhcb1.tfs0145","twiss_lhcb1.tfs0146","twiss_lhcb1.tfs0147","twiss_lhcb1.tfs0148","twiss_lhcb1.tfs0149","twiss_lhcb1.tfs0150","twiss_lhcb1.tfs0151","twiss_lhcb1.tfs0152","twiss_lhcb1.tfs0153","twiss_lhcb1.tfs0154",
                            "twiss_lhcb1.tfs0155","twiss_lhcb1.tfs0156","twiss_lhcb1.tfs0157","twiss_lhcb1.tfs0158","twiss_lhcb1.tfs0159","twiss_lhcb1.tfs0160","twiss_lhcb1.tfs0161","twiss_lhcb1.tfs0162","twiss_lhcb1.tfs0163","twiss_lhcb1.tfs0164",
                            "twiss_lhcb1.tfs0165","twiss_lhcb1.tfs0166","twiss_lhcb1.tfs0167","twiss_lhcb1.tfs0168","twiss_lhcb1.tfs0169","twiss_lhcb1.tfs0170","twiss_lhcb1.tfs0171","twiss_lhcb1.tfs0172","twiss_lhcb1.tfs0173","twiss_lhcb1.tfs0174",
                            "twiss_lhcb1.tfs0175","twiss_lhcb1.tfs0176","twiss_lhcb1.tfs0177","twiss_lhcb1.tfs0178","twiss_lhcb1.tfs0179","twiss_lhcb1.tfs0180","twiss_lhcb1.tfs0181","twiss_lhcb1.tfs0182","twiss_lhcb1.tfs0183","twiss_lhcb1.tfs0184",
                            "twiss_lhcb1.tfs0185","twiss_lhcb1.tfs0186","twiss_lhcb1.tfs0187","twiss_lhcb1.tfs0188","twiss_lhcb1.tfs0189","twiss_lhcb1.tfs0190","twiss_lhcb1.tfs0191","twiss_lhcb1.tfs0192","twiss_lhcb1.tfs0193","twiss_lhcb1.tfs0194",
                            "twiss_lhcb1.tfs0195","twiss_lhcb1.tfs0196","twiss_lhcb1.tfs0197","twiss_lhcb1.tfs0198","twiss_lhcb1.tfs0199","twiss_lhcb1.tfs0200"
                            ),
                        Beam1FilesXi   = cms.vdouble(0.,0.001),
                        Beam1Objects   = cms.vstring("TCL.4R5.B1","TCL.5R5.B1","XRPH.C6R5.B1","XRPH.D6R5.B1","XRPH.E6R5.B1","XRPH.B6R5.B1"),
                        Beam2Files     = cms.vstring(
                            "twiss_lhcb2.tfs0000",
                            "twiss_lhcb2.tfs0001","twiss_lhcb2.tfs0002","twiss_lhcb2.tfs0003","twiss_lhcb2.tfs0004",
                            "twiss_lhcb2.tfs0005","twiss_lhcb2.tfs0006","twiss_lhcb2.tfs0007","twiss_lhcb2.tfs0008","twiss_lhcb2.tfs0009","twiss_lhcb2.tfs0010","twiss_lhcb2.tfs0011","twiss_lhcb2.tfs0012","twiss_lhcb2.tfs0013","twiss_lhcb2.tfs0014",
                            "twiss_lhcb2.tfs0015","twiss_lhcb2.tfs0016","twiss_lhcb2.tfs0017","twiss_lhcb2.tfs0018","twiss_lhcb2.tfs0019","twiss_lhcb2.tfs0020","twiss_lhcb2.tfs0021","twiss_lhcb2.tfs0022","twiss_lhcb2.tfs0023","twiss_lhcb2.tfs0024",
                            "twiss_lhcb2.tfs0025","twiss_lhcb2.tfs0026","twiss_lhcb2.tfs0027","twiss_lhcb2.tfs0028","twiss_lhcb2.tfs0029","twiss_lhcb2.tfs0030","twiss_lhcb2.tfs0031","twiss_lhcb2.tfs0032","twiss_lhcb2.tfs0033","twiss_lhcb2.tfs0034",
                            "twiss_lhcb2.tfs0035","twiss_lhcb2.tfs0036","twiss_lhcb2.tfs0037","twiss_lhcb2.tfs0038","twiss_lhcb2.tfs0039","twiss_lhcb2.tfs0040","twiss_lhcb2.tfs0041","twiss_lhcb2.tfs0042","twiss_lhcb2.tfs0043","twiss_lhcb2.tfs0044",
                            "twiss_lhcb2.tfs0045","twiss_lhcb2.tfs0046","twiss_lhcb2.tfs0047","twiss_lhcb2.tfs0048","twiss_lhcb2.tfs0049","twiss_lhcb2.tfs0050","twiss_lhcb2.tfs0051","twiss_lhcb2.tfs0052","twiss_lhcb2.tfs0053","twiss_lhcb2.tfs0054",
                            "twiss_lhcb2.tfs0055","twiss_lhcb2.tfs0056","twiss_lhcb2.tfs0057","twiss_lhcb2.tfs0058","twiss_lhcb2.tfs0059","twiss_lhcb2.tfs0060","twiss_lhcb2.tfs0061","twiss_lhcb2.tfs0062","twiss_lhcb2.tfs0063","twiss_lhcb2.tfs0064",
                            "twiss_lhcb2.tfs0065","twiss_lhcb2.tfs0066","twiss_lhcb2.tfs0067","twiss_lhcb2.tfs0068","twiss_lhcb2.tfs0069","twiss_lhcb2.tfs0070","twiss_lhcb2.tfs0071","twiss_lhcb2.tfs0072","twiss_lhcb2.tfs0073","twiss_lhcb2.tfs0074",
                            "twiss_lhcb2.tfs0075","twiss_lhcb2.tfs0076","twiss_lhcb2.tfs0077","twiss_lhcb2.tfs0078","twiss_lhcb2.tfs0079","twiss_lhcb2.tfs0080","twiss_lhcb2.tfs0081","twiss_lhcb2.tfs0082","twiss_lhcb2.tfs0083","twiss_lhcb2.tfs0084",
                            "twiss_lhcb2.tfs0085","twiss_lhcb2.tfs0086","twiss_lhcb2.tfs0087","twiss_lhcb2.tfs0088","twiss_lhcb2.tfs0089","twiss_lhcb2.tfs0090","twiss_lhcb2.tfs0091","twiss_lhcb2.tfs0092","twiss_lhcb2.tfs0093","twiss_lhcb2.tfs0094",
                            "twiss_lhcb2.tfs0095","twiss_lhcb2.tfs0096","twiss_lhcb2.tfs0097","twiss_lhcb2.tfs0098","twiss_lhcb2.tfs0099","twiss_lhcb2.tfs0100","twiss_lhcb2.tfs0101","twiss_lhcb2.tfs0102","twiss_lhcb2.tfs0103","twiss_lhcb2.tfs0104",
                            "twiss_lhcb2.tfs0105","twiss_lhcb2.tfs0106","twiss_lhcb2.tfs0107","twiss_lhcb2.tfs0108","twiss_lhcb2.tfs0109","twiss_lhcb2.tfs0110","twiss_lhcb2.tfs0111","twiss_lhcb2.tfs0112","twiss_lhcb2.tfs0113","twiss_lhcb2.tfs0114",
                            "twiss_lhcb2.tfs0115","twiss_lhcb2.tfs0116","twiss_lhcb2.tfs0117","twiss_lhcb2.tfs0118","twiss_lhcb2.tfs0119","twiss_lhcb2.tfs0120","twiss_lhcb2.tfs0121","twiss_lhcb2.tfs0122","twiss_lhcb2.tfs0123","twiss_lhcb2.tfs0124",
                            "twiss_lhcb2.tfs0125","twiss_lhcb2.tfs0126","twiss_lhcb2.tfs0127","twiss_lhcb2.tfs0128","twiss_lhcb2.tfs0129","twiss_lhcb2.tfs0130","twiss_lhcb2.tfs0131","twiss_lhcb2.tfs0132","twiss_lhcb2.tfs0133","twiss_lhcb2.tfs0134",
                            "twiss_lhcb2.tfs0135","twiss_lhcb2.tfs0136","twiss_lhcb2.tfs0137","twiss_lhcb2.tfs0138","twiss_lhcb2.tfs0139","twiss_lhcb2.tfs0140","twiss_lhcb2.tfs0141","twiss_lhcb2.tfs0142","twiss_lhcb2.tfs0143","twiss_lhcb2.tfs0144",
                            "twiss_lhcb2.tfs0145","twiss_lhcb2.tfs0146","twiss_lhcb2.tfs0147","twiss_lhcb2.tfs0148","twiss_lhcb2.tfs0149","twiss_lhcb2.tfs0150","twiss_lhcb2.tfs0151","twiss_lhcb2.tfs0152","twiss_lhcb2.tfs0153","twiss_lhcb2.tfs0154",
                            "twiss_lhcb2.tfs0155","twiss_lhcb2.tfs0156","twiss_lhcb2.tfs0157","twiss_lhcb2.tfs0158","twiss_lhcb2.tfs0159","twiss_lhcb2.tfs0160","twiss_lhcb2.tfs0161","twiss_lhcb2.tfs0162","twiss_lhcb2.tfs0163","twiss_lhcb2.tfs0164",
                            "twiss_lhcb2.tfs0165","twiss_lhcb2.tfs0166","twiss_lhcb2.tfs0167","twiss_lhcb2.tfs0168","twiss_lhcb2.tfs0169","twiss_lhcb2.tfs0170","twiss_lhcb2.tfs0171","twiss_lhcb2.tfs0172","twiss_lhcb2.tfs0173","twiss_lhcb2.tfs0174",
                            "twiss_lhcb2.tfs0175","twiss_lhcb2.tfs0176","twiss_lhcb2.tfs0177","twiss_lhcb2.tfs0178","twiss_lhcb2.tfs0179","twiss_lhcb2.tfs0180","twiss_lhcb2.tfs0181","twiss_lhcb2.tfs0182","twiss_lhcb2.tfs0183","twiss_lhcb2.tfs0184",
                            "twiss_lhcb2.tfs0185","twiss_lhcb2.tfs0186","twiss_lhcb2.tfs0187","twiss_lhcb2.tfs0188","twiss_lhcb2.tfs0189","twiss_lhcb2.tfs0190","twiss_lhcb2.tfs0191","twiss_lhcb2.tfs0192","twiss_lhcb2.tfs0193","twiss_lhcb2.tfs0194",
                            "twiss_lhcb2.tfs0195","twiss_lhcb2.tfs0196","twiss_lhcb2.tfs0197","twiss_lhcb2.tfs0198","twiss_lhcb2.tfs0199","twiss_lhcb2.tfs0200"
                            ),
                        Beam2FilesXi   = cms.vdouble(0.,0.001),
                        Beam2Objects   = cms.vstring("TCL.4L5.B2","TCL.5L5.B2","XRPH.C6L5.B2","XRPH.D6L5.B2","XRPH.E6L5.B2","XRPH.B6L5.B2"),
                        Beam1Direction = cms.int32(1),
                        Beam2Direction = cms.int32(1),
                        SmearEnergy    = cms.bool(E_smear),
                        SmearAngle     = cms.bool(Ang_smear),
                        BeamEnergy     = cms.double(ecms/2.0),
                        BeamEnergyRMS  = cms.double(1.11e-4),
                        BeamAngleRMS   = cms.double(35.355), # in urad
                        Beam1CenterAtStation = cms.vdouble(0.,0.,0.,0.,0.,0.), # beam sigma (X) at Arm Forward tracker stations in mm
                        Beam1SigmaAtStation  = cms.vdouble(0.,0.,0.266876,0.336219,0.361528,0.361528), # beam sigma (X) at Arm Forward tracker stations in mm
                        Beam1StationPositions = cms.vdouble(143.0,183.8,203.827,212.550,215.7,219.551), #Z position in m of the stations
                        Beam1StationSigmaInsertion = cms.vdouble(15.,15.,15.,15.,15.,15.),
                        Beam2CenterAtStation = cms.vdouble(0.,0.,0.,0.,0.,0.), # beam sigma (X) at Arm Forward tracker stations in mm
                        Beam2SigmaAtStation  = cms.vdouble(0.,0.,0.149449,0.216493,0.241087,0.361528), # beam sigma (X) at Arm Forward tracker stations in mm
                        Beam2StationPositions = cms.vdouble(143.0,183.8,203.827,212.550,215.7,219.551), #Z position in m of the stations
                        Beam2StationSigmaInsertion = cms.vdouble(15.,15.,15.,15.,15.,15.),
                        ShowBeamLine   = cms.untracked.bool(False),
                        SimBeamProfile = cms.untracked.bool(True),
                        CrossAngleCorr = cms.bool(useCR),
                        CrossingAngle  = cms.double(205.0) #in mrad
                        )
    ppssim_tofDiamond = cms.PSet(
                        ToFGeometry       = cms.string("diamond"),
                        ToFWidth          = cms.double(4.2*4.),
                        ToFHeight         = cms.double(4.2), # tof cell height in mm
                        )
    ppssim_tofQuartz = cms.PSet(
                        ToFGeometry       = cms.string("quartz"),
                        ToFCellWidth      = cms.untracked.vdouble(3.0), # tof cell width in mm #move to vector - diamond geometry
                        ToFCellHeight     = cms.double(3.0), # tof cell height in mm
                        ToFNCellX         = cms.int32(5),      # number of cells in X
                        ToFNCellY         = cms.int32(4),      # number of cells in Y
                        )
    ppssim_detector_options = cms.PSet(
                        TrackerInsertion    = cms.double(15), # Number of sigmas (X) from the beam for the tracker
                        ToFInsertion        = cms.double(15), # Number of sigmas (X) from the beam for the tof
                        # TrackerZPosition    = cms.double(det1),
                        # TrackerLength       = cms.double(trklen),
                        Beam1TrkDetNames    = cms.vstring("XRPH.C6R5.B1","XRPH.D6R5.B1","XRPH.B6R5.B1"),
                        Beam2TrkDetNames    = cms.vstring("XRPH.C6L5.B2","XRPH.D6L5.B2","XRPH.B6L5.B2"),
                        # TrkBeamObjectsName  = cms.vstring("XRPV.C6R5","XRPH.D6R5"),
                        TrkDetIds           = cms.vint32(20,30,120),
                        ToFDetXOffset       = cms.double(tofxoffset), # ToF sensitive area distance from RP edge
                        ToFZPosition        = cms.double(tof),
                        Beam1ToFDetectorName= cms.string("XRPH.E6R5.B1"),
                        Beam2ToFDetectorName= cms.string("XRPH.E6L5.B2"),
                        SmearHit            = cms.bool(hit_smear),
                        HitSigmaX           = cms.double(0.01),
                        HitSigmaY           = cms.double(0.01),
                        HitSigmaZ           = cms.double(0),
                        ToFHitSigmaX        = cms.double(0.15),
                        ToFHitSigmaY        = cms.double(4.2/math.sqrt(12)),
                        TimeSigma           = cms.double(0.08), #in ns
                        MinThetaXatDet1     = cms.double(-500.), #min. theta x at first tracker in urad
                        MaxThetaXatDet1     = cms.double(500.),   #max. theta x at first tracker in urad
                        MinThetaYatDet1     = cms.double(-500.), #min. theta y at first tracker in urad
                        MaxThetaYatDet1     = cms.double(500.), #max. theta y at first tracker in urad
                        XTrackChiSquareCut  = cms.double(2.),     #X track ChiSquare cut 
                        YTrackChiSquareCut  = cms.double(2.),     #X track ChiSquare cut 
                        ApplyFiducialCuts   = cms.bool(True),     #apply geometrical (X,Y) in the hits
                        UseToFForTracking   = cms.bool(False) #not working!! problem with added point for Tracker 1
                        )

    ppssim_tracker_strip = cms.PSet(
                        TrackerGeometryType = cms.string("TOTEMStrip"),
                        NumberOfStrips      = cms.int32(512),
                        StripPitch          = cms.double(0.066), #mm
                        CutSideLength       = cms.double(2.), #mm length of the cut side of the strip sensor
                        VerticalShift       = cms.double(0.), #vertical shift of the detector with respect to the origin
                        TrkDetXOffset       = cms.double(0.2+0.05), # tracker sensitive area distance from RP edge in mm
                        TrkDetXRotation     = cms.vdouble(0.,0.,0.), #Rotation in degree along X axis of the tracking station in the same orded of TrkDetName
                        TrkDetYRotation     = cms.vdouble(0.,0.,0.), #Rotation in degree along Y axis of the tracking station in the same orded of TrkDetName
                        TrkDetZRotation     = cms.vdouble(8.,0.,0.), #Rotation in degree along Z axis of the tracking station in the same orded of TrkDetName
                        ClusterSizePlotFile = cms.string("FastSimulation/PPSFastSim/data/StripClusterSizeDistribution.root"),
                        ClusterSizePlotName = cms.string("hStripClusterSize")
                        )

    #         / \
    #       /     \
    #     /         \
    #   /             |
    #                 |  o 
    #   \             |
    #     \         /
    #       \     / 
    #         \ /
    #
    #       y|  /z
    #        | /
    #        |/____x
    
    ppssim_tracker_pixel = cms.PSet(
                        TrackerGeometryType = cms.string("Pixel"),
                        NumberOfRows        = cms.int32(160),
                        NumberOfColumns     = cms.int32(156),
                        DoubleSizeColumn    = cms.vint32(51,52,103,104),
                        DoubleSizeRow       = cms.vint32(79,80),
                        PixelPitchX         = cms.double(0.150),  # pitch between cells in X in mm
                        PixelPitchY         = cms.double(0.100),  # pitch between cells in Y in mm
                        VerticalShift       = cms.double(0.), #vertical shift of the detector with respect to the origin
                        TrkDetXOffset       = cms.double(0.200+0.200), # tracker sensitive area distance from RP edge 
                        TrkDetXRotation     = cms.vdouble(0.,0.,0.), #Rotation in degree along X axis of the tracking station in the same orded of TrkDetName
                        TrkDetYRotation     = cms.vdouble(0.,0.,0.), #Rotation in degree along Y axis of the tracking station in the same orded of TrkDetName
                        TrkDetZRotation     = cms.vdouble(8.,0.,0.) #Rotation in degree along Z axis of the tracking station in the same orded of TrkDetName
                        )

    #        column <--
    #    ______________ row
    #   |              | | 
    #   |              | V o
    #   |______________| 
    #               
    #
    #       y|  /z
    #        | /
    #        |/____x
    
    ppssim_general_options = cms.PSet(
                        UseHepMCProducer        = cms.untracked.bool(True), 
                        VtxMeanX                = VertexX,
                        VtxMeanY                = VertexY,
                        VtxMeanZ                = VertexZ,
                        CollisionPoint          = cms.string("IP5"),
                        Beam1TCLSigmaInsertion  = cms.untracked.vdouble(15.,15.),
                        Beam2TCLSigmaInsertion  = cms.untracked.vdouble(15.,15.),
                        Beam1TCLNames           = cms.vstring("TCL.4R5.B1","TCL.5R5.B1"),
                        Beam2TCLNames           = cms.vstring("TCL.4L5.B2","TCL.5L5.B2"),
                        Beam1CenterAtTCL        = cms.vdouble(0.,0.), # beam center (X) at Arm Forward tracker stations in mm
                        Beam1SizeAtTCL          = cms.vdouble(0.,0.), # beam sigma (X) at Arm Forward tracker stations in mm
                        Beam2CenterAtTCL        = cms.vdouble(0.,0.), # beam center (X) at Arm Forward tracker stations in mm
                        Beam2SizeAtTCL          = cms.vdouble(0.,0.), # beam sigma (X) at Arm Forward tracker stations in mm
                        PhiMin                  = cms.double(phi_min),
                        PhiMax                  = cms.double(phi_max),
                        CentralMass             = cms.double(125.7),
                        CentralMassErr          = cms.double(0.4),
                        EtaMin                  = cms.double(7.0),  # min eta to be tracked by HECTOR
                        MomentumMin             = cms.double(3.000), # min mom. to be tracked by HECTOR
                        TrackImpactParameterCut = cms.double(0.5) # max. imp. par. for reco tracks
                        )


    process.ppssim = cms.EDProducer('PPSProducer',
                        ppssim_beam_options,
                        ppssim_tofDiamond,
                        #ppssim_tofQuartz,
                        ppssim_tracker_strip,
                        ppssim_detector_options,
                        ppssim_general_options,
                        genSource = cms.InputTag("generatorSmeared") # for HepMC event -> no pileup events
                        )

    # Adding CT-PPS 
    process.AODSIMoutput.outputCommands.extend(cms.untracked.vstring('keep PPSGenDataPPSSpectrometer_*_*_*','keep PPSSimDataPPSSpectrometer_*_*_*','keep PPSRecoDataPPSSpectrometer_*_*_*','keep TotemRPDigiedmDetSetVector_*_*_*'))## add CTPPS eventContent
    # Path and EndPath definitions
    #process.ppssim_step = cms.Path(process.ppssim)#CTPPS
    #process.generation_step.replace(process.generator,process.generator * process.VtxSmeared)#CTPPS
    # process.simulation_step.replace(process.psim, process.psim * process.ppssim)
    process.digitisation_step.replace(process.pdigi_valid, process.pdigi_valid * process.ppssim)

    # CMSSW_7_6_0_pre1 or greather 
    #process.schedule.insert(6,process.ppssim_step)	

    return (process)

def customise_pu_protons_ctpps(process):

    process=customise(process)
        
    process.mix.mixObjects.mixHepMC.makeCrossingFrame = cms.untracked.bool(True)

    process.ppssim.UseHepMCProducer = cms.untracked.bool(False)
    # PU gen particle   		
    process.genParticlesPU = cms.EDProducer("GenParticleProducer",
    saveBarCodes = cms.untracked.bool(True),
    mix = cms.string("mix"),
    abortOnUnknownPDGCode = cms.untracked.bool(False),
    useCrossingFrame = cms.untracked.bool(True)
    )
    process.genProtonsPU = cms.EDFilter("GenParticleSelector",
        filter = cms.bool(False),
        src = cms.InputTag("genParticlesPU"),
        cut = cms.string('')
    )

    process.genProtonsPU.cut = 'status = 1 & pdgId == 2212 & abs(pz) >= %f' % ( 0.5*13000./2.0)
    outputCommandsPU = [ 'keep *_genParticlesPU_*_*', 'keep *_genProtonsPU_*_*']

    process.ppssim.genSource = cms.InputTag("genProtonsPU") # for Pile-up events

    #process.digitisation_step.replace(process.pdigi_valid, process.pdigi_valid * process.genParticlesPU * process.genProtonsPU)
    process.digitisation_step.replace(process.ppssim, process.genParticlesPU * process.genProtonsPU * process.ppssim)

    return (process)
