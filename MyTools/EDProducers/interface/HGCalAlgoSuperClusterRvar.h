# ifndef HGCalAlgoSuperClusterRvar_H
# define HGCalAlgoSuperClusterRvar_H


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

# include "MyTools/EDProducers/interface/CommonUtilities.h"


class HGCalAlgoSuperClusterRvar
{
    private:
    
    bool debug_;
    
    int nLayer_;
    
    double cylinderR_;
    double cylinderR2_;
    
    double minHitE_;
    double minHitET_;
    double minHitET2_;
    
    hgcal::RecHitTools recHitTools_;
    
    
    public:
    
    HGCalAlgoSuperClusterRvar(
        const edm::ParameterSet&
    );
    
    double getRvar(
        const edm::EventSetup&,
        const reco::SuperClusterRef&,
        const std::vector <reco::PFRecHit>&
    );
};


# endif
