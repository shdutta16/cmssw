import os

import FWCore.ParameterSet.Config as cms

processName = "Demo"

from Configuration.StandardSequences.Eras import eras
#process = cms.Process("RECO", eras.Run2_2017)
#process = cms.Process("RECO", eras.Run2_2018)

#process = cms.Process(processName, eras.phase2_hgcal)
#process = cms.Process(processName, eras.Phase2C8_timing_layer_bar)
process = cms.Process(processName, eras.Phase2C9)

#process = cms.Process("RECO", eras.Phase2C8_timing_layer_bar)

#process.load("FWCore.MessageService.MessageLogger_cfi")
#MessageLogger = cms.Service("MessageLogger")


# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
##process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
##process.load('SimGeneral.MixingModule.mixNoPU_cfi')
##process.load('Configuration.Geometry.GeometryExtended2023D41Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.RecoSim_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
##process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
#process.GlobalTag.globaltag = "94X_mc2017_realistic_v10"
#process.GlobalTag.globaltag = "100X_upgrade2018_realistic_v10"

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, "auto:phase2_realistic", "")
process.GlobalTag = GlobalTag(process.GlobalTag, "auto:phase2_realistic_T15", "")
#process.GlobalTag = GlobalTag(process.GlobalTag, "94X_mc2017_realistic_v10")
#process.GlobalTag = GlobalTag(process.GlobalTag, "94X_mc2017_realistic_v11")


#process.load('Configuration.Geometry.GeometryExtended2023D41Reco_cff')
#process.load('Configuration.Geometry.GeometryExtended2026D41Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')

#process.load('Configuration.Geometry.GeometryExtended2023D41_cff')
#process.load('Configuration.Geometry.GeometryExtended2026D41_cff')
process.load('Configuration.Geometry.GeometryExtended2026D49_cff')


############################## Parse arguments ##############################

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing("analysis")


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

options.register("rerunTICL",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Whether to rerun TICL" # Description
)

options.register("modTICLele",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Whether to use modified TICL-electron sequence" # Description
)

options.register("onRaw",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Running on RAW" # Description
)

options.register("modTICLeleWithRerunTICL",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Whether to use the rerun TICL for the modified TICL-electron sequence" # Description
)

options.register("storeSimHit",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Store sim-hits" # Description
)

options.register("storeRecHit",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Store rec-hits" # Description
)

options.register("storeHGCALlayerClus",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Store HGCal layer clusters" # Description
)

options.register("storeSuperClusTICLclus",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Store info about TICL-electron SC, SC seed, and TICL-cluster matches" # Description
)

options.register("TICLeleGenMatchDR",
    99999, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.float, # string, int, or float
    "DeltaR to use for TICL-electron gen-matching (will store only the gen-matched ones)" # Description
)

options.register("TICLphoGenMatchDR",
    99999, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.float, # string, int, or float
    "DeltaR to use for TICL-photon gen-matching (will store only the gen-matched ones)" # Description
)

options.register("isGunSample",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Is it a particle gun sample" # Description
)

options.register("genEleFilter",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Apply gen-electron filter" # Description
)

options.register("genPhoFilter",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Apply gen-photon filter" # Description
)

options.register("genPartonFilter",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Apply gen-parton filter" # Description
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

options.register("depGraph",
    0, # Default value
    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    VarParsing.VarParsing.varType.int, # string, int, or float
    "Produce dependency graph only" # Description
)

options.parseArguments()


#maxEvents = -1
#options.maxEvents = 15
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))


#sourceFile = "sourceFiles/SingleElectronFlatPtGun_pT-15_rajdeep/SingleElectronFlatPtGun_pT-15_rajdeep.txt"
#sourceFile = "sourceFiles/SingleElectronFlatPtGun_pT-15_eta-1p5-3p0_sobhatta-crab_SingleElectronFlatPtGun_pT-15_eta-1p5-3p0_GEN-SIM-RECO-ffc2278112c688bef3890fc698a39794_USER/SingleElectronFlatPtGun_pT-15_eta-1p5-3p0_sobhatta-crab_SingleElectronFlatPtGun_pT-15_eta-1p5-3p0_GEN-SIM-RECO-ffc2278112c688bef3890fc698a39794_USER.txt"
#sourceFile = "sourceFiles/SingleElectronFlatPtGun_pT-15_eta-1p5-3p0_withTICLfractions/SingleElectronFlatPtGun_pT-15_eta-1p5-3p0_withTICLfractions.txt"

#sourceFile = "sourceFiles/SinglePi0FlatPtGun_pT-15_eta-1p5-3p0_GEN-SIM-RECO/SinglePi0FlatPtGun_pT-15_eta-1p5-3p0_GEN-SIM-RECO.txt"
#sourceFile = "sourceFiles/SingleZprimeToEEflatPtGun_m-5_pT-35_eta-1p5-3p0_GEN-SIM-RECO/SingleZprimeToEEflatPtGun_m-5_pT-35_eta-1p5-3p0_GEN-SIM-RECO.txt"

#sourceFile = "sourceFiles/RelValElectronGunPt2To100_CMSSW_10_6_0_pre4-106X_upgrade2023_realistic_v2_2023D41noPU-v1_GEN-SIM-RECO/RelValElectronGunPt2To100_CMSSW_10_6_0_pre4-106X_upgrade2023_realistic_v2_2023D41noPU-v1_GEN-SIM-RECO.txt"
#sourceFile = "sourceFiles/RelValElectronGunPt2To100_CMSSW_10_6_0_patch2-106X_upgrade2023_realistic_v3_2023D41noPU-v1_GEN-SIM-RECO/RelValElectronGunPt2To100_CMSSW_10_6_0_patch2-106X_upgrade2023_realistic_v3_2023D41noPU-v1_GEN-SIM-RECO.txt"

sourceFile = "sourceFiles/SingleElectronFlatPtGun_fpantale_pT-0-200_eta-1p5-3p0_GEN-SIM-RECO/SingleElectronFlatPtGun_fpantale_pT-0-200_eta-1p5-3p0_GEN-SIM-RECO_mod.txt"

if (len(options.sourceFile)) :
    
    sourceFile = options.sourceFile


fNames = []

if (len(options.inputFiles)) :
    
    fNames = options.inputFiles

else :
    
    with open(sourceFile) as f:
        
        fNames = f.readlines()

#fNames = ["file:/afs/cern.ch/work/s/sobhatta/private/HGCal_ele-reco/CMSSW_11_0_0_pre4/src/output_GEN-SIM-RECO_numEvent20.root"]
#fNames = ["file:/afs/cern.ch/work/s/sobhatta/private/HGCal_ele-reco/CMSSW_11_0_0_pre4/src/output_GEN-SIM-RECO_numEvent100.root"]
#fNames = ["file:/afs/cern.ch/work/s/sobhatta/private/HGCal_ele-reco/CMSSW_11_0_0_pre4/src/output_GEN-SIM-RECO_numEvent2000.root"]
#fNames = ["file:/eos/cms/store/group/phys_egamma/fpantale/output_GEN-SIM-RECO.root"]


for iFile, fName in enumerate(fNames) :
    
    if (
        "file:" not in fName and
        "root:" not in fName
    ) :
        
        fNames[iFile] = "file:%s" %(fName)


outFileSuffix = ""

if (options.rerunTICL) :
    
    outFileSuffix = "%s_rerunTICL" %(outFileSuffix)


if (options.modTICLele) :
    
    if (options.rerunTICL and options.modTICLeleWithRerunTICL) :
        
        outFileSuffix = "%s_modTICLeleWithRerunTICL" %(outFileSuffix)
        
    else :
        
        outFileSuffix = "%s_modTICLele" %(outFileSuffix)


if (options.onRaw) :
    
    outFileSuffix = "%s_onRaw" %(outFileSuffix)


if (options.outFileNumber >= 0) :
    
    outFileSuffix = "%s_%d" %(outFileSuffix, options.outFileNumber)


outFile = "ntupleTree%s.root" %(outFileSuffix)

if (len(options.outputDir)) :
    
    os.system("mkdir -p %s" %(options.outputDir))
    
    outFile = "%s/%s" %(options.outputDir, outFile)


sourceFileNames = cms.untracked.vstring(fNames)
#print sourceFileNames

process.source = cms.Source("PoolSource",
    fileNames = sourceFileNames,
    
    # Run1:Event1 to Run2:Event2
    #eventsToProcess = cms.untracked.VEventRange("1:78722-1:78722"),
    
    #duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
)


if (len(options.eventRange)) :
    
    process.source.eventsToProcess = cms.untracked.VEventRange(options.eventRange)


#process.options = cms.untracked.PSet(
#    #SkipEvent = cms.untracked.vstring("ProductNotFound"),
#    
#    #printDependencies = cms.untracked.bool(True),
#)


if (options.depGraph) :
    
    process.DependencyGraph = cms.Service("DependencyGraph")
    process.source = cms.Source("EmptySource")
    process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(0))


if (options.onRaw) :
    
    #from Configuration.StandardSequences.Reconstruction_cff import *
    #
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
    #    
    #    #egammaGlobalReco,
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
    
    
    # Reco
    #print "*"*50
    #print "reconstruction    :", process.reconstruction
    
    #from RecoTracker.IterativeTracking.MuonSeededStep_cff import *
    
    process.reconstruction_mod = process.reconstruction.copy()
    
    #process.reconstruction_mod.replace(localreco, localreco_mod)
    #process.reconstruction_mod.replace(globalreco, globalreco_mod)
    ###process.reconstruction_mod.replace(highlevelreco, highlevelreco_mod)
    ###process.reconstruction_mod.replace(highlevelreco, cms.Sequence(process.egammaHighLevelRecoPrePF))
    #process.reconstruction_mod.replace(highlevelreco, cms.Sequence())
    
    #print ""
    #print "process.reconstruction_mod:", process.reconstruction_mod
    #print ""
    
    
    options.rerunTICL = 1
    options.modTICLele = 1
    options.modTICLeleWithRerunTICL = 1




#from EDProducers.EnergySharedTICLmultiClusterProducer.EnergySharedTICLmultiCluster_cfi import *

#energySharingAlgo = "Mean"
energySharingAlgo = "Expo"
#energySharingAlgo = "Gaus"

#distanceType = "2Ddist"
distanceType = "3Ddist"

#process.EnergySharedTICLmultiClusters = EnergySharedTICLmultiClusters.clone()
#process.EnergySharedTICLmultiClusters.algoTypeStr = energySharingAlgo
#process.EnergySharedTICLmultiClusters.distTypeStr = distanceType


# TICL-electron ambiguity resolver
#import EDProducers.HGCalElectronAmbiguityResolver.HGCalElectronAmbiguityResolver_cfi

#process.ambiguityResolvedTICLelectrons = HGCalElectronAmbiguityResolver_cfi.AmbiguityResolvedElectrons


# Rerun TICL
#label_TICLtrackster = cms.InputTag("trackstersEM", "", "RECO")
#label_TICLmultiCluster = cms.InputTag("multiClustersFromTrackstersEM", "MultiClustersFromTracksterByCA", "RECO")

label_TICLtrackster = cms.InputTag("ticlTrackstersMerge", "", "RECO")
label_TICLmultiCluster = cms.InputTag("ticlMultiClustersFromTrackstersMerge", "", "RECO")

if (options.rerunTICL) :
    
    label_TICLtrackster = cms.InputTag("ticlTrackstersMerge", "", processName)
    label_TICLmultiCluster = cms.InputTag("ticlMultiClustersFromTrackstersMerge", "", processName)
    #label_TICLmultiCluster = cms.InputTag("hgcalMultiClusters", "", processName)


# Mod TICL-ele
#label_gsfEleFromTICL = cms.InputTag("ecalDrivenGsfElectronsFromTICL", "", "RECO")
label_gsfEleFromTICL = cms.InputTag("ecalDrivenGsfElectronsFromMultiCl", "", "RECO")

if (options.modTICLele) :
    
    label_gsfEleFromTICL = cms.InputTag("ecalDrivenGsfElectronsFromMultiCl", "", processName)
    
    label_phoFromTICL = cms.InputTag("photonsFromMultiCl", "", "Demo")


# TICL-ele variables
from MyTools.EDProducers.producers_cfi import *

l_var_TICLele = []
l_mapProdVars_TICLele = []
TICLeleVar_task = cms.Task()

l_var_TICLpho = []
l_mapProdVars_TICLpho = []
TICLphoVar_task = cms.Task()

l_rad = [1.8, 2.0, 2.8, 3.5]
#l_rad = [2.8]
l_hitEnCut = [0.0, 0.02, 0.5, 1.0, 1.5]

l_isoCone = [0.12, 0.15, 0.2, 0.25]
l_clusEnCut = [0.0, 0.02, 0.5, 1.0, 1.5]


for iCone, cone in enumerate(l_isoCone) :
    
    for iEnCut, enCut in enumerate(l_clusEnCut) :
        
        radStr = "R%s" %(str(cone).replace(".", "p"))
        enCutStr = "En%s" %(str(enCut).replace(".", "p"))
        
        
        # Ele H/E
        prodLabel = "TICLeleHoverEProducer%s%s" %(radStr, enCutStr)
        
        setattr(
            process,
            prodLabel,
            HGCalElectronHoverEProducer.clone(
                electrons = label_gsfEleFromTICL,
                coneDR = cone,
                minClusE = enCut,
                minPt = 10.0,
                #debug = cms.bool(True),
            )
        )
        
        TICLeleVar_task.add(getattr(process, prodLabel))
        
        l_mapProdVars_TICLele.append(cms.InputTag(prodLabel, getattr(process, prodLabel).instanceName.value(), processName))
        l_var_TICLele.append("%s_%s" %(prodLabel, getattr(process, prodLabel).instanceName.value()))
        
        
        # Pho H/E
        prodLabel = "TICLphoHoverEProducer%s%s" %(radStr, enCutStr)
        
        setattr(
            process,
            prodLabel,
            HGCalPhotonHoverEProducer.clone(
                photons = label_phoFromTICL,
                coneDR = cone,
                minClusE = enCut,
                minPt = 10.0,
                #debug = cms.bool(True),
            )
        )
        
        TICLphoVar_task.add(getattr(process, prodLabel))
        
        l_mapProdVars_TICLpho.append(cms.InputTag(prodLabel, getattr(process, prodLabel).instanceName.value(), processName))
        l_var_TICLpho.append("%s_%s" %(prodLabel, getattr(process, prodLabel).instanceName.value()))
    
    
    ## Track-iso
    #prodLabel = "TICLeleTrackIsoProducer%s" %(radStr)
    #
    #setattr(
    #    process,
    #    prodLabel,
    #    HGCalElectronTrackIsoProducer.clone(
    #        electrons = label_gsfEleFromTICL,
    #        isoConeDR = cone,
    #        minTrackPt = 0.0,
    #        maxTrackEleDz = 99999.0,
    #        #debug = cms.bool(True),
    #    )
    #)
    #
    #TICLeleVar_task.add(getattr(process, prodLabel))
    #
    #l_mapProdVars_TICLele.append(
    #    cms.InputTag(prodLabel, getattr(process, prodLabel).instanceName.value(), processName)
    #)
    #
    #
    #prodLabel = "TICLeleTrackIsoProducer%sTrkDz0p15TrkPt1" %(radStr)
    #
    #setattr(
    #    process,
    #    prodLabel,
    #    HGCalElectronTrackIsoProducer.clone(
    #        electrons = label_gsfEleFromTICL,
    #        isoConeDR = cone,
    #        #debug = cms.bool(True),
    #    )
    #)
    #
    #TICLeleVar_task.add(getattr(process, prodLabel))
    #
    #l_mapProdVars_TICLele.append(
    #    cms.InputTag(prodLabel, getattr(process, prodLabel).instanceName.value(), processName)
    #)


for iRad, rad in enumerate(l_rad) :
    
    for iEnCut, enCut in enumerate(l_hitEnCut) :
        
        radStr = "R%s" %(str(rad).replace(".", "p"))
        enCutStr = "En%s" %(str(enCut).replace(".", "p"))
        
        
        # Ele Rvar
        prodLabel = "TICLeleRvarProducer%s%s" %(radStr, enCutStr)
        
        setattr(
            process,
            prodLabel,
            HGCalElectronRvarProducer.clone(
                electrons = label_gsfEleFromTICL,
                cylinderR = rad,
                minHitE = enCut,
                minPt = 10.0,
                #debug = cms.bool(True),
            )
        )
        
        TICLeleVar_task.add(getattr(process, prodLabel))
        
        l_mapProdVars_TICLele.append(cms.InputTag(prodLabel, getattr(process, prodLabel).instanceName.value(), processName))
        l_var_TICLele.append("%s_%s" %(prodLabel, getattr(process, prodLabel).instanceName.value()))
        
        
        # Pho Rvar
        prodLabel = "TICLphoRvarProducer%s%s" %(radStr, enCutStr)
        
        setattr(
            process,
            prodLabel,
            HGCalPhotonRvarProducer.clone(
                photons = label_phoFromTICL,
                cylinderR = rad,
                minHitE = enCut,
                minPt = 10.0,
                #debug = cms.bool(True),
            )
        )
        
        TICLphoVar_task.add(getattr(process, prodLabel))
        
        l_mapProdVars_TICLpho.append(cms.InputTag(prodLabel, getattr(process, prodLabel).instanceName.value(), processName))
        l_var_TICLpho.append("%s_%s" %(prodLabel, getattr(process, prodLabel).instanceName.value()))
        
        
        # Ele PCA
        prodLabel = "TICLelePCAProducer%s%s" %(radStr, enCutStr)
        
        setattr(
            process,
            prodLabel,
            HGCalElectronPCAProducer.clone(
                electrons = label_gsfEleFromTICL,
                cylinderR = rad,
                minHitE = enCut,
                minPt = 10.0,
                #debug = cms.bool(True),
            )
        )
        
        TICLeleVar_task.add(getattr(process, prodLabel))
        
        l_mapProdVars_TICLele.append(cms.InputTag(prodLabel, getattr(process, prodLabel).instanceName.value()+"Sigma2UU", processName))
        l_mapProdVars_TICLele.append(cms.InputTag(prodLabel, getattr(process, prodLabel).instanceName.value()+"Sigma2VV", processName))
        l_mapProdVars_TICLele.append(cms.InputTag(prodLabel, getattr(process, prodLabel).instanceName.value()+"Sigma2WW", processName))
        
        l_var_TICLele.append("%s_%s" %(prodLabel, getattr(process, prodLabel).instanceName.value()+"Sigma2UU"))
        l_var_TICLele.append("%s_%s" %(prodLabel, getattr(process, prodLabel).instanceName.value()+"Sigma2VV"))
        l_var_TICLele.append("%s_%s" %(prodLabel, getattr(process, prodLabel).instanceName.value()+"Sigma2WW"))
        
        
        # Pho PCA
        prodLabel = "TICLphoPCAProducer%s%s" %(radStr, enCutStr)
        
        setattr(
            process,
            prodLabel,
            HGCalPhotonPCAProducer.clone(
                photons = label_phoFromTICL,
                cylinderR = rad,
                minHitE = enCut,
                minPt = 10.0,
                #debug = cms.bool(True),
            )
        )
        
        TICLphoVar_task.add(getattr(process, prodLabel))
        
        l_mapProdVars_TICLpho.append(cms.InputTag(prodLabel, getattr(process, prodLabel).instanceName.value()+"Sigma2UU", processName))
        l_mapProdVars_TICLpho.append(cms.InputTag(prodLabel, getattr(process, prodLabel).instanceName.value()+"Sigma2VV", processName))
        l_mapProdVars_TICLpho.append(cms.InputTag(prodLabel, getattr(process, prodLabel).instanceName.value()+"Sigma2WW", processName))
        
        l_var_TICLpho.append("%s_%s" %(prodLabel, getattr(process, prodLabel).instanceName.value()+"Sigma2UU"))
        l_var_TICLpho.append("%s_%s" %(prodLabel, getattr(process, prodLabel).instanceName.value()+"Sigma2VV"))
        l_var_TICLpho.append("%s_%s" %(prodLabel, getattr(process, prodLabel).instanceName.value()+"Sigma2WW"))


# TICL ele
process.HGCalElectronVarMap = mapProducer.clone(
    collections = cms.VInputTag(l_mapProdVars_TICLele),
    
    debug = cms.bool(False),
)

# TICL pho
process.HGCalPhotonVarMap = mapProducer.clone(
    collections = cms.VInputTag(l_mapProdVars_TICLpho),
    
    debug = cms.bool(False),
)


process.HGCalVar_seq = cms.Sequence(
    process.HGCalElectronVarMap *
    process.HGCalPhotonVarMap
)

process.HGCalVar_seq.associate(TICLeleVar_task)
process.HGCalVar_seq.associate(TICLphoVar_task)


process.treeMaker = cms.EDAnalyzer(
    "TreeMaker",
    
    ############################## My stuff ##############################
    debug = cms.bool(False),
    
    isGunSample = cms.bool(bool(options.isGunSample)),
    
    storeSimHit = cms.bool(bool(options.storeSimHit)),
    storeRecHit = cms.bool(bool(options.storeRecHit)),
    storeHGCALlayerClus = cms.bool(bool(options.storeHGCALlayerClus)),
    storeSuperClusTICLclus = cms.bool(bool(options.storeSuperClusTICLclus)),
    
    TICLeleGenMatchDR = cms.double(options.TICLeleGenMatchDR),
    TICLphoGenMatchDR = cms.double(options.TICLphoGenMatchDR),
    
    
    ############################## GEN ##############################
    
    label_generator = cms.InputTag("generator"),
    label_genParticle = cms.InputTag("genParticles"),
    
    
    ############################## RECO ##############################
    
    label_pileup = cms.InputTag("addPileupInfo"),
    label_rho = cms.InputTag("fixedGridRhoFastjetAll"),
    
    label_HGCEESimHit = cms.InputTag("g4SimHits", "HGCHitsEE"),
    #label_HGCEESimHit = cms.InputTag("g4SimHits", "HGCHitsEE", "SIM"),
    
    label_HGCEERecHit = cms.InputTag("HGCalRecHit", "HGCEERecHits"),
    label_HGCHEFRecHit = cms.InputTag("HGCalRecHit", "HGCHEFRecHits"),
    label_HGCHEBRecHit = cms.InputTag("HGCalRecHit", "HGCHEBRecHits"),
    
    label_HGCALlayerCluster = cms.InputTag("hgcalLayerClusters"),
    
    #label_TICLmultiCluster = cms.InputTag("MultiClustersFromTracksters", "MultiClustersFromTracksterByCA", "RECO"),
    #label_TICLmultiCluster = cms.InputTag("EnergySharedTICLmultiClusters", "EnergySharedTICLmultiClusters%s%s" %(energySharingAlgo, distanceType)),
    
    label_TICLtrackster = label_TICLtrackster,
    label_TICLmultiCluster = label_TICLmultiCluster,
    
    #label_TICLmultiClusterMIP = cms.InputTag("MultiClustersFromTrackstersMIP", "MIPMultiClustersFromTracksterByCA"),
    
    label_PFRecHit = cms.InputTag("particleFlowRecHitHGC", "Cleaned"),
    
    label_caloParticle = cms.InputTag("mix", "MergedCaloTruth"),
    
    label_gsfEleFromMultiClus = cms.InputTag("ecalDrivenGsfElectronsFromMultiCl", ""),
    
    label_gsfEleFromTICL = label_gsfEleFromTICL,
    label_gsfEleFromTICLvarList = cms.vstring(l_var_TICLele),
    label_gsfEleFromTICLvarMap = cms.InputTag("HGCalElectronVarMap", process.HGCalElectronVarMap.instanceName.value(), processName),
    
    label_generalTrack = cms.InputTag("generalTracks"),
    
    label_phoFromMultiClus = cms.InputTag("photonsFromMultiCl", "", "RECO"),
    
    label_phoFromTICL = label_phoFromTICL,
    label_phoFromTICLvarList = cms.vstring(l_var_TICLpho),
    label_phoFromTICLvarMap = cms.InputTag("HGCalPhotonVarMap", process.HGCalPhotonVarMap.instanceName.value(), processName),
    
    
    ########## AK4 jet ##########
    
    #label_ak4PFjet = cms.InputTag("ak4PFJets"),
    
    
)



#if (options.modTICLele) :
#    
#    ###from MyModules.Test.ecalDrivenGsfElectronsFromTICL_cff_orig import ecalDrivenGsfElectronsFromTICL_customizeProcess
#    from MyModules.Test.ecalDrivenGsfElectronsFromTICL_cff import ecalDrivenGsfElectronsFromTICL_customizeProcess
#    
#    process = ecalDrivenGsfElectronsFromTICL_customizeProcess(process, onReco = (not options.onRaw))
#    
#    if (options.rerunTICL and options.modTICLeleWithRerunTICL) :
#        
#        #print label_TICLmultiCluster
#        #print label_TICLmultiCluster.__dict__
#        
#        #process.particleFlowClusterHGCalFromTICL.initialClusteringStep.clusterSrc = label_TICLmultiCluster
#        
#        process.particleFlowClusterHGCalFromTICL.initialClusteringStep.clusterSrc = cms.InputTag(
#            label_TICLmultiCluster.__dict__["_InputTag__moduleLabel"],
#            label_TICLmultiCluster.__dict__["_InputTag__productInstance"],
#            label_TICLmultiCluster.__dict__["_InputTag__processName"],
#        )


########## Filters ##########

from EDFilters.MyFilters.GenParticleFilter_cfi import *

# Gen-ele filter
process.GenParticleFilter_ele = GenParticleFilter.clone()
process.GenParticleFilter_ele.atLeastN = cms.int32(1)
process.GenParticleFilter_ele.pdgIds = cms.vint32(11)
process.GenParticleFilter_ele.minPt = cms.double(10)
process.GenParticleFilter_ele.minEta = cms.double(1.479)
process.GenParticleFilter_ele.maxEta = cms.double(3.1)
process.GenParticleFilter_ele.isGunSample = cms.bool(bool(options.isGunSample))
#process.GenParticleFilter_ele.debug = cms.bool(True)

process.filter_seq_genEle = cms.Sequence()

if (options.genEleFilter) :

    process.filter_seq_genEle = cms.Sequence(process.GenParticleFilter_ele)


# Gen-pho filter
process.GenParticleFilter_pho = GenParticleFilter.clone()
process.GenParticleFilter_pho.atLeastN = cms.int32(1)
process.GenParticleFilter_pho.pdgIds = cms.vint32(22)
process.GenParticleFilter_pho.minPt = cms.double(10)
process.GenParticleFilter_pho.minEta = cms.double(1.479)
process.GenParticleFilter_pho.maxEta = cms.double(3.1)
process.GenParticleFilter_pho.isGunSample = cms.bool(bool(options.isGunSample))
#process.GenParticleFilter_pho.debug = cms.bool(True)

process.filter_seq_genPho = cms.Sequence()

if (options.genPhoFilter) :

    process.filter_seq_genPho = cms.Sequence(process.GenParticleFilter_pho)


# Gen-parton filter
process.GenParticleFilter_part = GenParticleFilter.clone()
process.GenParticleFilter_part.atLeastN = cms.int32(1)
process.GenParticleFilter_part.pdgIds = cms.vint32(1, 2, 3, 4, 5, 21)
process.GenParticleFilter_part.minPt = cms.double(10)
process.GenParticleFilter_part.minEta = cms.double(1.479)
process.GenParticleFilter_part.maxEta = cms.double(3.1)
process.GenParticleFilter_part.isGunSample = cms.bool(bool(options.isGunSample))
#process.GenParticleFilter_part.debug = cms.bool(True)

process.filter_seq_genPart = cms.Sequence()

if (options.genPartonFilter) :

    process.filter_seq_genPart = cms.Sequence(process.GenParticleFilter_part)



#process.filter_path = cms.Path(
#    process.filter_seq_genEle *
#    process.filter_seq_genPart
#)


print "Deleting existing output file."
os.system("rm %s" %(outFile))


# Output file name modification
if (outFile.find("/eos/cms") ==  0) :
    
    outFile = outFile.replace("/eos/cms", "root://eoscms.cern.ch//eos/cms")


# Output
process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string(outFile)
)


process.schedule = cms.Schedule()


# Aging
from SLHCUpgradeSimulations.Configuration.aging import customise_aging_1000
customise_aging_1000(process)


process.reco_seq = cms.Sequence()

if (options.onRaw) :
    
    process.reco_seq = cms.Sequence(
        process.RawToDigi *
        process.L1Reco *
        process.reconstruction_mod
    )


## TICL
#if (options.rerunTICL) :
#    
#    #from RecoHGCal.TICL.ticl_iterations import TICL_iterations
#    #TICL_iterations(process)
#    
#    from RecoHGCal.TICL.ticl_iterations import TICL_iterations_withReco
#    TICL_iterations_withReco(process)
#
#
#process.TICLele_seq = cms.Sequence()
#
#if (options.modTICLele) :
#    
#    process.TICLele_seq = cms.Sequence(process.ecalDrivenGsfElectronsFromTICL_step)


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
    #process.filter_seq *
    
    process.filter_seq_genEle *
    process.filter_seq_genPho *
    process.filter_seq_genPart *
    
    process.reco_seq *
    #process.TICLele_seq *
    process.HGCalVar_seq *
    process.treeMaker
)

process.schedule.insert(0, process.p)

print "\n"
print "*"*50
print "process.schedule:", process.schedule
print "*"*50
#print "process.schedule.__dict__:", process.schedule.__dict__
#print "*"*50
print "\n"


# Tracer
if (options.trace) :
    
    process.Tracer = cms.Service("Tracer")


if (options.memoryCheck) :
    
    process.SimpleMemoryCheck = cms.Service(
        "SimpleMemoryCheck",
        moduleMemorySummary = cms.untracked.bool(True),
    )


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


process.MessageLogger = cms.Service(
    "MessageLogger",
    destinations = cms.untracked.vstring(
        "cerr",
    ),
    cerr = cms.untracked.PSet(
        #threshold  = cms.untracked.string("ERROR"),
        DEBUG = cms.untracked.PSet(limit = cms.untracked.int32(0)),
        WARNING = cms.untracked.PSet(limit = cms.untracked.int32(0)),
        ERROR = cms.untracked.PSet(limit = cms.untracked.int32(0)),
    )
)


from FWCore.ParameterSet.Utilities import convertToUnscheduled
process = convertToUnscheduled(process)


# Add early deletion of temporary data products to reduce peak memory need
#from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
#process = customiseEarlyDelete(process)
# End adding early deletion
