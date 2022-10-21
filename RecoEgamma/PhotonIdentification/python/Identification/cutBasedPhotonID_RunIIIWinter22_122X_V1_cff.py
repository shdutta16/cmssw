from PhysicsTools.SelectorUtils.centralIDRegistry import central_id_registry

import FWCore.ParameterSet.Config as cms

# Common functions and classes for ID definition are imported here:
from RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_tools \
    import ( WorkingPoint_V3,
             IsolationCutInputs,
             ClusterIsolationCutInputs,
             HoverECutInputs,
             configureVIDCutBasedPhoID_V6 )


# DOCUMENTATION UPDATE WILL BE DONE SOON


#
# First, define cut values
#


# Loose working point Barrel and Endcap
idName = "cutBasedPhotonID-RunIIIWinter22-122X-V1-loose"
WP_Loose_EB = WorkingPoint_V3(
    idName     ,  # idName
    0.011145   ,  # full5x5_SigmaIEtaIEtaCut
    0.12999    ,  # hOverEWithEACut
# Isolation cuts are generally absIso < C1 + pt*C2, except for HCalClus is < C1 + pt*C2 + pt*pt*C3
    1.8852     ,  # absPFChgHadIsoWithEACut_C1
    0.0        ,  # absPFChgHadIsoWithEACut_C2
    0.70379    ,  # absPFECalClusIsoWithEACut_C1
    0.00065204 ,  # absPFECalClusIsoWithEACut_C2
    6.3440     ,  # absPFHCalClusIsoWithEACut_C1
    0.010055   ,  # absPFHCalClusIsoWithEACut_C2
    0.00005783    # absPFHCalClusIsoWithEACut_C3
    )
WP_Loose_EE = WorkingPoint_V3(
    idName     ,  # idName
    0.027674   ,  # full5x5_SigmaIEtaIEtaCut
    0.15343    ,  # hOverEWithEACut
# Isolation cuts are generally absIso < C1 + pt*C2, except for HCalClus is < C1 + pt*C2 + pt*pt*C3
    1.6540     ,  # absPFChgHadIsoWithEACut_C1
    0.0        ,  # absPFChgHadIsoWithEACut_C2
    6.61585    ,  # absPFECalClusIsoWithEACut_C1
    0.00019549 ,  # absPFECalClusIsoWithEACut_C2
    1.8588     ,  # absPFHCalClusIsoWithEACut_C1
    0.01170    ,  # absPFHCalClusIsoWithEACut_C2
    0.00007476    # absPFHCalClusIsoWithEACut_C3
    )



# Medium working point Barrel and Endcap
idName = "cutBasedPhotonID-RunIIIWinter22-122X-V1-medium"
WP_Medium_EB = WorkingPoint_V3(
    idName     ,  # idName
    0.0100     ,  # full5x5_SigmaIEtaIEtaCut
    0.058305   ,  # hOverEWithEACut
# Isolation cuts are generally absIso < C1 + pt*C2, except for HCalClus is < C1 + pt*C2 + pt*pt*C3
    0.93929    ,  # absPFChgHadIsoWithEACut_C1
    0.0        ,  # absPFChgHadIsoWithEACut_C2
    0.22770    ,  # absPFECalClusIsoWithEACut_C1
    0.0065204  ,  # absPFECalClusIsoWithEACut_C2
    2.1890     ,  # absPFHCalClusIsoWithEACut_C1
    0.010055   ,  # absPFHCalClusIsoWithEACut_C2
    0.00005783    # absPFHCalClusIsoWithEACut_C3
    )

WP_Medium_EE = WorkingPoint_V3(
    idName     ,  #idName
    0.02687    ,  # full5x5_SigmaIEtaIEtaCut
    0.005181   ,  # hOverECutWithEA
# Isolation cuts are generally absIso < C1 + pt*C2, except for HCalClus is < C1 + pt*C2 + pt*pt*C3
    0.97029    ,  # absPFChgHadIsoWithEACut_C1
    0.0        ,  # absPFChaHadIsoWithEACut_C2
    1.124      ,  # absPFECalClusIsoWithEACut_C1
    0.00019549 ,  # absPFECalClusIsoWithEACut_C2
    0.033670   ,  # absPFHCalClusIsowithEACut_C1
    0.01170    ,  # absPFHCalClusIsoWithEACut_C2
    0.00007476    # absPFHCalClusIsoWithEACut_C3
    )



# Tight working point Barrel and Endcap
idName = "cutBasedPhotonID-RunIIIWinter22-122X-V1-tight"
WP_Tight_EB = WorkingPoint_V3(
    idName     ,  # idName
    0.009993   ,  # full5x5_SigmaIEtaIEtaCut
    0.0417588  ,  # hOverECutWithEA
# Isolation cuts are generally absIso < C1 + pt*C2, except for HCalClus is < C1 + pt*C2 + pt*pt*C3
    0.3163     ,  # absPFChgHadIsoWithEACut_C1
    0.0        ,  # absPFChgHadIsoWithEACut_C2
    0.14189    ,  # absPFECalClusIsoWithEACut_C1
    0.00065204 ,  # absPFECalClusIsoWithEACut_C2
    0.010055   ,  # absPFHCalClusIsoWithEACut_C1
    2.044      ,  # absPFHCalClusIsoWithEACut_C2
    0.00005783    # absPFHCalClusIsoWithEACut_C3
    )

WP_Tight_EE = WorkingPoint_V3(
    idName     ,  # idName
    0.02687    ,  # full5x5_SigmaIEtaIEtaCut
    0.0025426  ,  # hOverECutWithEA
# Isolation cuts are generally absIso < C1 + pt*C2, except for HCalClus is < C1 + pt*C2 + pt*pt*C3
    0.29266    ,  # absPFChgHadIsoWithEACut_C1
    0.0        ,  # absPFChgHadIsoWithEACut_C2
    1.04269    ,  # absPFECalClusIsoWithEACut_C1
    0.00019549 ,  # absPFECalClusIsoWithEACut_C2
    0.029262   ,  # absPFHCalClusIsowithEACut_C1
    0.01170    ,  # absPFHCalClusIsoWithEACut_C2
    0.00007476    # absPFHCalClusIsoWithEACut_C3
    )


# Second, define where to find the precomputed isolations and what effective
# areas to use for pile-up correction
isoInputs = IsolationCutInputs(
    # chHadIsolationMapName
    'photonIDValueMapProducer:phoChargedIsolation' ,
    # chHadIsolationEffAreas
    "RecoEgamma/PhotonIdentification/data/RunIII_Winter22/effectiveArea_ChgHadronIso_95percentBased.txt",
    # neuHadIsolationMapName
    'photonIDValueMapProducer:phoNeutralHadronIsolation' ,
    # neuHadIsolationEffAreas
    "RecoEgamma/PhotonIdentification/data/RunIII_Winter22/effectiveArea_NeuHadronIso_95percentBased.txt" ,
    # phoIsolationMapName
    'photonIDValueMapProducer:phoPhotonIsolation' ,
    # phoIsolationEffAreas
    "RecoEgamma/PhotonIdentification/data/RunIII_Winter22/effectiveArea_PhotonIso_95percentBased.txt"
)

clusterIsoInputs = ClusterIsolationCutInputs(
    # trkIsolationMapName
    'photonIDValueMapProducer:phoTrkIsolation' ,
    # trkIsolationEffAreas
    "RecoEgamma/PhotonIdentification/data/RunIII_Winter22/effectiveArea_TrackerIso_95percentBased.txt",
    # ecalClusIsolationMapName
    'photonIDValueMapProducer:phoEcalPFClIsolation' ,
    # ecalClusIsolationEffAreas
    "RecoEgamma/PhotonIdentification/data/RunIII_Winter22/effectiveArea_ECalClusterIso_95percentBased.txt",
    # hcalClusIsolationMapName
    'photonIDValueMapProducer:phoHcalPFClIsolation' ,
    # hcalClusIsolationEffAreas
    "RecoEgamma/PhotonIdentification/data/RunIII_Winter22/effectiveArea_HCalClusterIso_95percentBased.txt"
)

hOverEInputs = HoverECutInputs(
    # hOverEEffAreas
    "RecoEgamma/PhotonIdentification/data/RunIII_Winter22/effectiveArea_coneBasedHoverE_95percentBased.txt"
)

#
# Finally, set up VID configuration for all cuts
#
cutBasedPhotonID_RunIIIWinter22_122X_V1_loose  = configureVIDCutBasedPhoID_V6 ( WP_Loose_EB, WP_Loose_EE, isoInputs, clusterIsoInputs, hOverEInputs)
cutBasedPhotonID_RunIIIWinter22_122X_V1_medium = configureVIDCutBasedPhoID_V6 ( WP_Medium_EB, WP_Medium_EE, isoInputs, clusterIsoInputs, hOverEInputs)
cutBasedPhotonID_RunIIIWinter22_122X_V1_tight  = configureVIDCutBasedPhoID_V6 ( WP_Tight_EB, WP_Tight_EE, isoInputs, clusterIsoInputs, hOverEInputs)

## The MD5 sum numbers below reflect the exact set of cut variables
# and values above. If anything changes, one has to
# 1) comment out the lines below about the registry,
# 2) run "calculateIdMD5 <this file name> <one of the VID config names just above>
# 3) update the MD5 sum strings below and uncomment the lines again.
#

central_id_registry.register(cutBasedPhotonID_RunIIIWinter22_122X_V1_loose.idName,
                             '7c76eb5b7876b6c42e0a8dfb5bd68b9048853401')
central_id_registry.register(cutBasedPhotonID_RunIIIWinter22_122X_V1_medium.idName,
                             '7477f41fac6e7b3d35b48cdfffa675fd3a2ad1ad')
central_id_registry.register(cutBasedPhotonID_RunIIIWinter22_122X_V1_tight.idName,
                             '53416e10005a7a274591f90c1bcccb58bffaa650')

cutBasedPhotonID_RunIIIWinter22_122X_V1_loose.isPOGApproved = cms.untracked.bool(True)
cutBasedPhotonID_RunIIIWinter22_122X_V1_medium.isPOGApproved = cms.untracked.bool(True)
cutBasedPhotonID_RunIIIWinter22_122X_V1_tight.isPOGApproved = cms.untracked.bool(True)

