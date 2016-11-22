
def customise(process):
    #X0 = cms.double(0.0322)
    #Y0 = cms.double(0.0)
    #Z0 = cms.double(0.0)
    process.ppssim.VtxMeanX = 0.0142
    process.ppssim.VtxMeanY = 0.0142
    process.ppssim.VtxMeanZ = 0.
    return process
