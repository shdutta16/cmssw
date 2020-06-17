import FWCore.ParameterSet.Config as cms


AmbiguityResolvedElectrons = cms.EDProducer(
    "EnergySharedTICLmultiClusterProducer",
    
    instanceName = cms.string("HGCalElectronAmbiguityResolver"),
    
    label_electron = cms.untracked.InputTag("ecalDrivenGsfElectronsFromTICL"),
    
    label_HGCEERecHit = cms.untracked.InputTag("HGCalRecHit" , "HGCEERecHits" ),
    label_HGCHEFRecHit = cms.untracked.InputTag("HGCalRecHit", "HGCHEFRecHits"),
    label_HGCHEBRecHit = cms.untracked.InputTag("HGCalRecHit", "HGCHEBRecHits"),
    
    
    debug = cms.bool(False),
)
