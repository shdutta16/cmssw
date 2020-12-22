import os

import FWCore.ParameterSet.Config as cms

processName = "Demo"

from Configuration.StandardSequences.Eras import eras
process = cms.Process(processName, eras.Phase2C9)

# import of standard configurations
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.RecoSim_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, "auto:phase2_realistic_T15", "")
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D49_cff')


############################## Parse arguments ##############################

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing("analysis")

options.register("debugFile",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Create debug file" # Description
)

options.parseArguments()


process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))

sourceFile = "list_of_input_files.txt"

fNames = []

if (len(options.inputFiles)) :
    
    fNames = options.inputFiles

else :
    
    with open(sourceFile) as f:
        
        fNames = f.readlines()


for iFile, fName in enumerate(fNames) :
    
    if (
        "file:" not in fName and
        "root:" not in fName
    ) :
        
        fNames[iFile] = "file:%s" %(fName)


outFile = "ntupleTree.root"


sourceFileNames = cms.untracked.vstring(fNames)

process.source = cms.Source("PoolSource",
    fileNames = sourceFileNames,
)


# HGCal variables
from MyTools.EDProducers.producers_cfi import *

process.HGCalElectronClusIso = HGCalElectronClusIsoProducer.clone()
process.HGCalElectronHoverE = HGCalElectronHoverEProducer.clone()
process.HGCalElectronRvar = HGCalElectronRvarProducer.clone()
process.HGCalElectronPCA = HGCalElectronPCAProducer.clone()

process.HGCalVar_seq = cms.Sequence(
    process.HGCalElectronClusIso *
    process.HGCalElectronHoverE *
    process.HGCalElectronRvar *
    process.HGCalElectronPCA
)

# Output
process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string(outFile)
)


process.schedule = cms.Schedule()


# Aging
from SLHCUpgradeSimulations.Configuration.aging import customise_aging_1000
customise_aging_1000(process)


process.reco_seq = cms.Sequence(
    process.RawToDigi *
    process.L1Reco *
    process.reconstruction
)


###### PixelCPE issue
process.TrackProducer.TTRHBuilder = "WithTrackAngle"
process.PixelCPEGenericESProducer.UseErrorsFromTemplates = False
process.PixelCPEGenericESProducer.LoadTemplatesFromDB = False
process.PixelCPEGenericESProducer.TruncatePixelCharge = False
process.PixelCPEGenericESProducer.IrradiationBiasCorrection = False
process.PixelCPEGenericESProducer.DoCosmics = False
process.PixelCPEGenericESProducer.Upgrade = cms.bool(True) 
######


process.p = cms.Path(
    process.reco_seq *
    process.HGCalVar_seq
)


process.schedule = cms.Schedule()
process.schedule.insert(0, process.p)


# Debug
if (options.debugFile) :
    
    process.out = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string("debug.root")
    )
    
    process.output_step = cms.EndPath(process.out)
    process.schedule.extend([process.output_step])
