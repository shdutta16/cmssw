# ifndef HGCalAlgoSuperClusterHoverE_H
# define HGCalAlgoSuperClusterHoverE_H


// system include files
#include <memory>

// user include files
# include "DataFormats/CaloRecHit/interface/CaloCluster.h"
# include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
# include "DataFormats/FWLite/interface/ESHandle.h"
# include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
# include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
# include "DataFormats/Math/interface/deltaR.h"
# include "DataFormats/Math/interface/LorentzVector.h"
# include "FWCore/Framework/interface/Event.h"
# include "FWCore/Framework/interface/Frameworkfwd.h"
# include "FWCore/Framework/interface/MakerMacros.h"
# include "FWCore/Framework/interface/stream/EDProducer.h"
# include "FWCore/ParameterSet/interface/ParameterSet.h"
# include "FWCore/Utilities/interface/StreamID.h"
# include "Geometry/CaloTopology/interface/HGCalTopology.h"
# include "Geometry/Records/interface/IdealGeometryRecord.h"
# include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

# include <CLHEP/Vector/LorentzVector.h>

# include <algorithm>
# include <iostream>
# include <map>
# include <stdlib.h>
# include <string>
# include <type_traits>
# include <utility>
# include <vector>



class HGCalAlgoSuperClusterHoverE
{
    private:
    
    bool debug_;
    
    double coneDR_;
    double coneDR2_;
    double minClusE_;
    double minClusET_;
    
    
    public:
    
    HGCalAlgoSuperClusterHoverE(
        const edm::ParameterSet&
    );
    
    double getClusterBasedHoverE(
        const reco::SuperClusterRef&,
        const std::vector <reco::CaloCluster>&
    );
};


# endif
