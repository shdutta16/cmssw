import FWCore.ParameterSet.Config as cms


HGCalElectronClusIsoProducer = cms.EDProducer(
    "HGCalElectronClusterIsoProducer",
    
    instanceName = cms.string("HGCalElectronClusIso"),
    
    electrons = cms.InputTag("ecalDrivenGsfElectronsFromMultiCl"),
    layerClusters = cms.InputTag("hgcalLayerClusters"),
    
    coneDRmin = cms.double(0.02),
    coneDRmax = cms.double(0.2),
    
    minEmClusE = cms.double(0.02),
    minEmClusET = cms.double(0.0),
    
    minHadClusE = cms.double(0.07),
    minHadClusET = cms.double(0.0),
    
    minPt = cms.double(0.0),
    
    debug = cms.bool(False),
)


HGCalPhotonClusIsoProducer = cms.EDProducer(
    "HGCalPhotonClusterIsoProducer",
    
    instanceName = cms.string("HGCalPhotonClusIso"),
    
    photons = cms.InputTag("photonsFromMultiCl"),
    layerClusters = cms.InputTag("hgcalLayerClusters"),
    
    coneDRmin = cms.double(0.02),
    coneDRmax = cms.double(0.2),
    
    minEmClusE = cms.double(0.02),
    minEmClusET = cms.double(0.0),
    
    minHadClusE = cms.double(0.07),
    minHadClusET = cms.double(0.0),
    
    minPt = cms.double(0.0),
    
    debug = cms.bool(False),
)


HGCalElectronHoverEProducer = cms.EDProducer(
    "HGCalElectronHoverEProducer",
    
    instanceName = cms.string("HGCalElectronHoverE"),
    
    electrons = cms.InputTag("ecalDrivenGsfElectronsFromMultiCl"),
    layerClusters = cms.InputTag("hgcalLayerClusters"),
    
    coneDR = cms.double(0.15),
    
    minClusE = cms.double(0.0),
    minClusET = cms.double(0.0),
    
    minPt = cms.double(0.0),
    
    debug = cms.bool(False),
)


HGCalPhotonHoverEProducer = cms.EDProducer(
    "HGCalPhotonHoverEProducer",
    
    instanceName = cms.string("HGCalPhotonHoverE"),
    
    photons = cms.InputTag("photonsFromMultiCl"),
    layerClusters = cms.InputTag("hgcalLayerClusters"),
    
    coneDR = cms.double(0.15),
    
    minClusE = cms.double(0.0),
    minClusET = cms.double(0.0),
    
    minPt = cms.double(0.0),
    
    debug = cms.bool(False),
)


HGCalElectronTrackIsoProducer = cms.EDProducer(
    "HGCalElectronTrackIsoProducer",
    
    instanceName = cms.string("HGCalElectronTrackIso"),
    
    electrons = cms.InputTag("ecalDrivenGsfElectronsFromMultiCl"),
    tracks = cms.InputTag("generalTracks"),
    
    isoConeDR = cms.double(0.3),
    vetoConeDR = cms.double(0.01),
    vetoPhiStripDeta = cms.double(0.03),
    
    minTrackPt = cms.double(1.0),
    maxTrackEleDz = cms.double(0.15),
    
    debug = cms.bool(False),
)


HGCalElectronRvarProducer = cms.EDProducer(
    "HGCalElectronRvarProducer",
    
    instanceName = cms.string("HGCalElectronRvar"),
    
    electrons = cms.InputTag("ecalDrivenGsfElectronsFromMultiCl"),
    
    PFRecHits = cms.InputTag("particleFlowRecHitHGC"),
    HGCEERecHits = cms.untracked.InputTag("HGCalRecHit", "HGCEERecHits"),
    HGCHEFRecHits = cms.untracked.InputTag("HGCalRecHit", "HGCHEFRecHits"),
    HGCHEBRecHits = cms.untracked.InputTag("HGCalRecHit", "HGCHEBRecHits"),
    
    cylinderR = cms.double(2.8),
    nLayer = cms.int32(28),
    
    minHitE = cms.double(0.0),
    minHitET = cms.double(0.0),
    
    minPt = cms.double(0.0),
    
    debug = cms.bool(False),
)


HGCalPhotonRvarProducer = cms.EDProducer(
    "HGCalPhotonRvarProducer",
    
    instanceName = cms.string("HGCalPhotonRvar"),
    
    photons = cms.InputTag("photonsFromMultiCl"),
    
    PFRecHits = cms.InputTag("particleFlowRecHitHGC"),
    HGCEERecHits = cms.untracked.InputTag("HGCalRecHit", "HGCEERecHits"),
    HGCHEFRecHits = cms.untracked.InputTag("HGCalRecHit", "HGCHEFRecHits"),
    HGCHEBRecHits = cms.untracked.InputTag("HGCalRecHit", "HGCHEBRecHits"),
    
    cylinderR = cms.double(2.8),
    nLayer = cms.int32(28),
    
    minHitE = cms.double(0.0),
    minHitET = cms.double(0.0),
    
    minPt = cms.double(0.0),
    
    debug = cms.bool(False),
)


HGCalElectronPCAProducer = cms.EDProducer(
    "HGCalElectronPCAProducer",
    
    instanceName = cms.string("HGCalElectronPCA"),
    
    electrons = cms.InputTag("ecalDrivenGsfElectronsFromMultiCl"),
    
    PFRecHits = cms.InputTag("particleFlowRecHitHGC"),
    HGCEERecHits = cms.untracked.InputTag("HGCalRecHit", "HGCEERecHits"),
    HGCHEFRecHits = cms.untracked.InputTag("HGCalRecHit", "HGCHEFRecHits"),
    HGCHEBRecHits = cms.untracked.InputTag("HGCalRecHit", "HGCHEBRecHits"),
    
    cylinderR = cms.double(2.8),
    nLayer = cms.int32(28),
    
    minHitE = cms.double(0.0),
    minHitET = cms.double(0.0),
    
    minPt = cms.double(0.0),
    
    debug = cms.bool(False),
)


HGCalPhotonPCAProducer = cms.EDProducer(
    "HGCalPhotonPCAProducer",
    
    instanceName = cms.string("HGCalPhotonPCA"),
    
    photons = cms.InputTag("photonsFromMultiCl"),
    
    PFRecHits = cms.InputTag("particleFlowRecHitHGC"),
    HGCEERecHits = cms.untracked.InputTag("HGCalRecHit", "HGCEERecHits"),
    HGCHEFRecHits = cms.untracked.InputTag("HGCalRecHit", "HGCHEFRecHits"),
    HGCHEBRecHits = cms.untracked.InputTag("HGCalRecHit", "HGCHEBRecHits"),
    
    cylinderR = cms.double(2.8),
    nLayer = cms.int32(28),
    
    minHitE = cms.double(0.0),
    minHitET = cms.double(0.0),
    
    minPt = cms.double(0.0),
    
    debug = cms.bool(False),
)


mapProducer = cms.EDProducer(
    "MapProducer",
    
    instanceName = cms.string("MapProducer"),
    useProcessName = cms.bool(False),
    
    collections = cms.VInputTag(),
    
    debug = cms.bool(False),
)
