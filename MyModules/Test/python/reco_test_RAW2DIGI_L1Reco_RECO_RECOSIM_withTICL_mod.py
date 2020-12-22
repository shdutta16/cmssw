# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: reco --filein file:Py8PtGun_cfi_py_GEN_SIM_DIGI_L1_L1TrackTrigger_DIGI2RAW_HLT.root --conditions auto:phase2_realistic -n 10 --era Phase2C8_timing_layer_bar --eventcontent FEVTDEBUGHLT --runUnscheduled -s RAW2DIGI,L1Reco,RECO,RECOSIM --datatier GEN-SIM-RECO --geometry Extended2023D41 --no_exec
import os

import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

#process = cms.Process('RECO',eras.Phase2C8_timing_layer_bar)
process = cms.Process('RECO',eras.Phase2C9)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
#process.load('Configuration.Geometry.GeometryExtended2026D41Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.RecoSim_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


############################## Parse arguments ##############################
import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing("analysis")

#options.outputFile = "output_GEN-SIM-RECO.root"

options.register("sourceFile",
    "", # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.string, # string, int, or float
    "File containing list of input files" # Description
)

options.register("outputDir",
    "", # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.string, # string, int, or float
    "Output directory" # Description
)

options.register("outFileNumber",
    -1, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "File number (will be added to the filename if >= 0)" # Description
)

options.register("isCRAB",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "CRAB run or not" # Description
)

options.register("eventRange",
    "", # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.string, # string, int, or float
    "Syntax: Run1:Event1-Run2:Event2 (includes both)" # Description
)

options.register("debugFile",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Create debug file" # Description
)

options.register("trace",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Trace modules" # Description
)

options.register("memoryCheck",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Check memory usage" # Description
)

options.register("printTime",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Print timing information" # Description
)


options.parseArguments()


if (options.isCRAB) :
    
    options.debugFile = 0

fNames = []


if (not options.isCRAB) :
    
    #sourceFile = "sourceFiles/RelValElectronGunPt2To100_CMSSW_10_6_0_pre4-106X_upgrade2023_realistic_v2_2023D41noPU-v1_GEN-SIM-DIGI-RAW/RelValElectronGunPt2To100_CMSSW_10_6_0_pre4-106X_upgrade2023_realistic_v2_2023D41noPU-v1_GEN-SIM-DIGI-RAW.txt"
    sourceFile = "sourceFiles/SingleElectronFlatPtGun_pT-15_eta-1p5-3p0_sobhatta-crab_SingleElectronFlatPtGun_pT-15_eta-1p5-3p0_GEN-SIM-DIGI-RAW-284165e958c955242461cd9651f5a03c_USER/SingleElectronFlatPtGun_pT-15_eta-1p5-3p0_sobhatta-crab_SingleElectronFlatPtGun_pT-15_eta-1p5-3p0_GEN-SIM-DIGI-RAW-284165e958c955242461cd9651f5a03c_USER_mod.txt"
    
    if (not(len(options.inputFiles))) :

        if (len(options.sourceFile)) :
    
            sourceFile = options.sourceFile
    
        with open(sourceFile) as f:
    
            fNames = f.readlines()

                
        options.inputFiles = fNames
        #print "Error: Enter a valid inputFiles."
        #exit(1)
    
    #if (not len(options.outputFile)) :
    #    
    #    #print "Error: Enter a valid outputFile."
    #    #exit(1)
    
    
    for iFile in range(0, len(options.inputFiles)) :
        
        fileName = options.inputFiles[iFile]
        
        if (fileName.find("file:") < 0 and fileName.find("root:") < 0) :
            
            options.inputFiles[iFile] = "file:%s" %(fileName)


outFileSuffix = ""

if (options.outFileNumber >= 0) :
    
    outFileSuffix = "%s_%d" %(outFileSuffix, options.outFileNumber)


outFile = "output_GEN-SIM-RECO%s.root" %(outFileSuffix)

if (len(options.outputDir)) :
    
    os.system("mkdir -p %s" %(options.outputDir))
    
    outFile = "%s/%s" %(options.outputDir, outFile)


#print options.maxEvents
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)


# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
    secondaryFileNames = cms.untracked.vstring()
)

if (len(options.eventRange)) :
    
    process.source.eventsToProcess = cms.untracked.VEventRange(options.eventRange)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('reco nevts:10'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output file name modification
if (outFile.find("/eos/cms") ==  0) :
    
    outFile = outFile.replace("/eos/cms", "root://eoscms.cern.ch//eos/cms")

print "Output file: %s" %(outFile)


# Output definition
process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        #dataTier = cms.untracked.string('AODSIM'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string(outFile),
    
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    #outputCommands = process.RECOSIMEventContent.outputCommands,
    #outputCommands = process.AODSIMEventContent.outputCommands,
    
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T15', '')


from Configuration.StandardSequences.Reconstruction_cff import *


## Local reco
#print "*"*50
#print "localreco    :", process.localreco
#localreco_mod = process.localreco.copyAndExclude([
#    muonlocalreco,
#    castorreco,
#])
#print ""
#print "localreco_mod:", localreco_mod
#print ""
#
#
## Global reco
##from RecoTracker.IterativeTracking.MuonSeededStep_cff import *
#from RecoTracker.IterativeTracking.iterativeTk_cff import *
#from RecoVertex.Configuration.RecoVertex_cff import *
#
#print "*"*50
#print "globalreco_tracking    :", process.globalreco_tracking
#globalreco_tracking_mod = process.globalreco_tracking.copyAndExclude([
#    standalonemuontracking,
#    muonSeededStepTask,
#    
#    unsortedOfflinePrimaryVertices4Dfastsim,
#    trackWithVertexRefSelectorBeforeSorting4Dfastsim,
#    trackRefsForJetsBeforeSorting4Dfastsim,
#    offlinePrimaryVertices4Dfastsim,
#    offlinePrimaryVertices4DfastsimWithBS,
#    unsortedOfflinePrimaryVertices4Dfastsim,
#    trackWithVertexRefSelectorBeforeSorting4Dfastsim,
#    trackRefsForJetsBeforeSorting4Dfastsim,
#    offlinePrimaryVertices4Dfastsim,
#    offlinePrimaryVertices4DfastsimWithBS,
#])
#print ""
#print "globalreco_tracking_mod:", globalreco_tracking_mod
#print ""
#
#print "*"*50
#print "globalreco    :", process.globalreco
#globalreco_mod = process.globalreco.copyAndExclude([
#    jetGlobalReco,
#    muonGlobalReco,
#    muoncosmicreco,
#    CastorFullReco,
#])
#
#globalreco_mod.replace(globalreco_tracking, globalreco_tracking_mod)
#print ""
#print "globalreco_mod:", globalreco_mod
#print ""
#
#
## High level reco
#print "*"*50
#print "highlevelreco    :", process.highlevelreco
#highlevelreco_mod = process.highlevelreco.copyAndExclude([
#    particleFlowReco,
#    muoncosmichighlevelreco,
#    muonshighlevelreco,
#    jetHighLevelReco,
#    metrecoPlusHCALNoise,
#    btagging,
#    recoPFMET,
#    PFTau,
#    cosmicDCTracksSeq,
#])
#print ""
#print "highlevelreco_mod:", highlevelreco_mod
#print ""
#
#
## Reco
#print "*"*50
#print "reconstruction    :", process.reconstruction
#
#from RecoTracker.IterativeTracking.MuonSeededStep_cff import *
#
#reconstruction_mod = process.reconstruction.copy()
#
#reconstruction_mod.replace(localreco, localreco_mod)
#reconstruction_mod.replace(globalreco, globalreco_mod)
##reconstruction_mod.replace(highlevelreco, highlevelreco_mod)
#reconstruction_mod.replace(highlevelreco, cms.Sequence(process.egammaHighLevelRecoPrePF))
#print ""
#print "reconstruction_mod:", reconstruction_mod
#print ""
#
##reconstruction_mod = cms.Sequence(
##    localreco_mod *
##    globalreco_mod *
##    highlevelreco_mod *
##    
##    process.logErrorHarvester
##)
#
#
# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
#process.reconstruction_step = cms.Path(reconstruction_mod)
process.recosim_step = cms.Path(process.recosim)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)


# Aging
from SLHCUpgradeSimulations.Configuration.aging import customise_aging_1000
customise_aging_1000(process)


# Import TICL
#from ticl_iterations import TICL_iterations
from RecoHGCal.TICL.ticl_iterations import TICL_iterations_withReco

# Attach a sequence to process: process.TICL
#process = TICL_iterations(process)
#TICL_iterations_withReco(process)


process.FEVTDEBUGHLTEventContent.outputCommands.extend([
#process.RECOSIMEventContent.outputCommands.extend([
#process.AODSIMEventContent.outputCommands.extend([
    "keep *_*particleFlowRecHitHGC*_*_*",
    "keep *_hgcalDigis_*_*",
    "keep *_particleFlowClusterECAL_*_*",
    "keep *_siPhase2Clusters_*_*",
    "keep *_siStripDigis_*_*",
    "keep *_siPixelClusters_*_*",
    "keep *_siPixelClustersCache_*_*", 
    "keep *_siPixelClusterShapeCache_*_*",
    "keep *_siPixelRecHits*_*_*",
    "keep *_trackerDrivenElectronSeeds_*_*",
    "keep *_pfTrackElec_*_*",
    
    "keep *_tracksters*_*_*",
    "keep *_ecalDrivenGsfElectrons*_*_*",
    "keep *_*FromTICL*_*_*",
    
    #"drop *",
    
    #"drop FEDRawDataCollection_*_*_*",
    #
    #"drop *_g4SimHits_*_*",
    #"keep *_g4SimHits_HGCHits*_*",
    #
    #"drop *_*ecalDigi*_*_*",
    #"drop *_*ecal*RecHit*_*_*",
    #
    #"drop *Hcal*_*_*_*",
    #"drop *_*Hcal*_*_*",
    #
    #"drop *_*horeco*_*_*",
    #
    #"drop *Muon*_*_*_*",
    #"drop *_*Muon*_*_*",
    #
    #"drop *Tau*_*_*_*",
    #"drop *_*Tau*_*_*",
    #
    #"drop *Vector*_*Jet*_*_*",
    #
    #"drop *_*mtdRecHit*_*_*",
    #
    #"drop *_*particleFlowRecHitECAL*_*_*",
    #"drop *_*particleFlowRecHitHBHE*_*_*",
    #"drop *_*particleFlowRecHitHF*_*_*",
    #
    #"drop *_*ctpps*_*_*",
    #"drop *_*totem*_*_*",
])


#from MyModules.Test.ecalDrivenGsfElectronsFromTICL_cff import ecalDrivenGsfElectronsFromTICL_customizeProcess
#process = ecalDrivenGsfElectronsFromTICL_customizeProcess(process)


########## Filters ##########

from EDFilters.MyFilters.GenParticleFilter_cfi import *

process.GenParticleFilter_ele = GenParticleFilter.clone()
process.GenParticleFilter_ele.atLeastN = cms.int32(2)
process.GenParticleFilter_ele.pdgId = cms.int32(11)
process.GenParticleFilter_ele.minEta = cms.double(1.479)
process.GenParticleFilter_ele.maxEta = cms.double(3.0)
#process.GenParticleFilter_ele.debug = cms.bool(True)


process.filter_seq = cms.Sequence(
    process.GenParticleFilter_ele
)

process.filter_path = cms.Path(
    process.filter_seq
)

process.FEVTDEBUGHLToutput.SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring("filter_path")
)

#process.reco_seq = cms.Sequence(
#    process.RawToDigi *
#    process.L1Reco *
#    process.reconstruction *
#    #process.recosim *
#    process.endOfProcess
#)
#
#process.p = cms.Path(
#    process.filter_seq *
#    process.reco_seq
#)


###### PixelCPE issue
process.TrackProducer.TTRHBuilder = "WithTrackAngle"
process.PixelCPEGenericESProducer.UseErrorsFromTemplates = False
process.PixelCPEGenericESProducer.LoadTemplatesFromDB = False
process.PixelCPEGenericESProducer.TruncatePixelCharge = False
process.PixelCPEGenericESProducer.IrradiationBiasCorrection = False
process.PixelCPEGenericESProducer.DoCosmics = False
process.PixelCPEGenericESProducer.Upgrade = cms.bool(True) 
######


# Schedule definition
process.schedule = cms.Schedule(
    process.filter_path,
    process.raw2digi_step,
    process.L1Reco_step,
    process.reconstruction_step,
    process.recosim_step,
    #process.ecalDrivenGsfElectronsFromTICL_step,
    
    #process.p,
    
    process.endjob_step,
    process.FEVTDEBUGHLToutput_step
)

#TICL_iterations_withReco(process)


# Tracer
if (options.trace) :
    
    process.Tracer = cms.Service("Tracer")


#Timing
if (options.printTime) :
    
    process.Timing = cms.Service("Timing",
        summaryOnly = cms.untracked.bool(False),
        useJobReport = cms.untracked.bool(True)
    )


# Debug
if (options.debugFile) :
    
    process.out = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string("debug.root")
    )
    
    process.output_step = cms.EndPath(process.out)
    process.schedule.extend([process.output_step])

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)


# Customisation from command line

#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion

