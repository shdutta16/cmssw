import FWCore.ParameterSet.Config as cms

from RecoHGCal.TICL.TICLSeedingRegions_cff import ticlSeedingGlobal, ticlSeedingGlobalHFNose
from RecoHGCal.TICL.trackstersProducer_cfi import trackstersProducer as _trackstersProducer
from RecoHGCal.TICL.filteredLayerClustersProducer_cfi import filteredLayerClustersProducer as _filteredLayerClustersProducer
from RecoHGCal.TICL.multiClustersFromTrackstersProducer_cfi import multiClustersFromTrackstersProducer as _multiClustersFromTrackstersProducer

# CLUSTER FILTERING/MASKING

filteredLayerClustersEM = _filteredLayerClustersProducer.clone(
    clusterFilter = "ClusterFilterByAlgoAndSizeAndLayerRange",
    min_cluster_size = 3, # inclusive
    max_layerId = 30, # inclusive
    algo_number = 8,
    LayerClustersInputMask = 'ticlTrackstersTrkEM',
    iteration_label = "EM"
)

filteredLayerClustersEMnoTrk = filteredLayerClustersEM.clone(
    LayerClustersInputMask = cms.InputTag("hgcalLayerClusters","InitialLayerClustersMask")
)

# CA - PATTERN RECOGNITION

ticlTrackstersEM = _trackstersProducer.clone(
    filtered_mask = "filteredLayerClustersEM:EM",
    original_mask = 'ticlTrackstersTrkEM',
    seeding_regions = "ticlSeedingGlobal",
    filter_on_categories = [0, 1],
    pid_threshold = 0.5,
    energy_em_over_total_threshold = 0.9,
    max_longitudinal_sigmaPCA = 10,
    shower_start_max_layer = 5, #inclusive
    max_out_in_hops = 1,
    skip_layers = 2,
    max_missing_layers_in_trackster = 1,
    min_layers_per_trackster = 10,
    min_cos_theta = 0.97,  # ~14 degrees
    min_cos_pointing = 0.9, # ~25 degrees
    max_delta_time = 3.,
    itername = "EM",
    algo_verbosity = 0,
)

ticlTrackstersEMnoTrk = ticlTrackstersEM.clone(
    filtered_mask = "filteredLayerClustersEMnoTrk:EM",
    original_mask = cms.InputTag("hgcalLayerClusters","InitialLayerClustersMask"),
)

# MULTICLUSTERS

ticlMultiClustersFromTrackstersEM = _multiClustersFromTrackstersProducer.clone(
    Tracksters = "ticlTrackstersEM"
)

ticlMultiClustersFromTrackstersEMnoTrk = ticlMultiClustersFromTrackstersEM.clone(Tracksters = "ticlTrackstersEMnoTrk")

ticlEMStepTask = cms.Task(ticlSeedingGlobal
    ,filteredLayerClustersEM
    ,ticlTrackstersEM
    ,ticlMultiClustersFromTrackstersEM
    ,filteredLayerClustersEMnoTrk
    ,ticlTrackstersEMnoTrk
    ,ticlMultiClustersFromTrackstersEMnoTrk
)

filteredLayerClustersHFNoseEM = filteredLayerClustersEM.clone(
    LayerClusters = 'hgcalLayerClustersHFNose',
    LayerClustersInputMask = "hgcalLayerClustersHFNose:InitialLayerClustersMask",
    iteration_label = "EMn",
    algo_number = 9
#no tracking mask for EM for now
)

ticlTrackstersHFNoseEM = ticlTrackstersEM.clone(
    detector = "HFNose",
    layer_clusters = "hgcalLayerClustersHFNose",
    layer_clusters_hfnose_tiles = "ticlLayerTileHFNose",
    original_mask = "hgcalLayerClustersHFNose:InitialLayerClustersMask",
    filtered_mask = "filteredLayerClustersHFNoseEM:EMn",
    seeding_regions = "ticlSeedingGlobalHFNose",
    time_layerclusters = "hgcalLayerClustersHFNose:timeLayerCluster",
    min_layers_per_trackster = 6
)

ticlHFNoseEMStepTask = cms.Task(ticlSeedingGlobalHFNose
                              ,filteredLayerClustersHFNoseEM
                              ,ticlTrackstersHFNoseEM
)
