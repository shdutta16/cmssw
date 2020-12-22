# include "MyTools/EDProducers/interface/HGCalAlgoSuperClusterRvar.h"


HGCalAlgoSuperClusterRvar::HGCalAlgoSuperClusterRvar(
    const edm::ParameterSet &iConfig
)
{
    nLayer_ = iConfig.getParameter <int>("nLayer");
    assert(nLayer_ <= 28);
    
    cylinderR_ = iConfig.getParameter <double>("cylinderR");
    cylinderR2_ = cylinderR_*cylinderR_;
    
    minHitE_ = iConfig.getParameter <double>("minHitE");
    minHitET_ = iConfig.getParameter <double>("minHitET");
    minHitET2_ = minHitET_*minHitET_;
    
    debug_ = iConfig.getParameter <bool>("debug");
}


double HGCalAlgoSuperClusterRvar::getRvar(
    const edm::EventSetup &iSetup,
    const reco::SuperClusterRef &superClus,
    const std::vector <reco::PFRecHit> &recHits
)
{
    double superClus_E = superClus->energy();
    
    if(!superClus_E)
    {
        return std::numeric_limits<double>::max();
    }
    
    auto m_recHitPtr = CommonUtilities::getPFRecHitPtrMap(recHits);
    
    double totalE = 0;
    
    std::vector <std::pair <DetId, float> > v_superClus_HandF = superClus->hitsAndFractions();
    
    std::vector <std::pair <DetId, float> > v_superClus_HandF_selected;
    std::vector <double> v_superClus_hitE_selected;
    
    std::vector <double> v_layerEnergy(nLayer_, 0.0);
    std::vector <double> v_layerEnergyInR(nLayer_, 0.0);
    std::vector <double> v_layerRvar(nLayer_, 0.0);
    
    std::vector <ROOT::Math::XYZVector> v_layerCentroid(nLayer_);
    
    //for(int iLayer = 0; iLayer < nLayer_; iLayer++)
    //{
    //    ROOT::Math::XYZVector xyz_temp(0, 0, 0);
    //    
    //    v_layerCentroid.push_back(xyz_temp);
    //}
    
    //edm::ESHandle <CaloGeometry> geom;
    //iSetup.get<CaloGeometryRecord>().get(geom);
    //recHitTools_.setGeometry(*(geom.product()));
    
    CommonUtilities::initRecHitTools(recHitTools_, &iSetup);
    
    // Compute the centroid per layer
    for(auto HandF : v_superClus_HandF)
    {
        DetId hitId = HandF.first;
        double hitEfrac = HandF.second;
        
        int hitLayer = recHitTools_.getLayer(hitId) - 1;
        
        if(hitLayer > nLayer_-1)
        {
            continue;
        }
        
        if(hitId.det() != DetId::HGCalEE)
        {
            continue;
        }
        
        if(m_recHitPtr.find(hitId) == m_recHitPtr.end())
        {
            continue;
        }
        
        reco::PFRecHit recHit = *m_recHitPtr[hitId];
        
        if(recHit.energy() < minHitE_)
        {
            continue;
        }
        
        if(recHit.pt2() < minHitET2_)
        {
            continue;
        }
        
        
        // Fill the vector of selected hits
        v_superClus_HandF_selected.push_back(HandF);
        
        
        double hitE = recHit.energy() * hitEfrac;
        v_superClus_hitE_selected.push_back(hitE);
        
        totalE += hitE;
        
        auto hitPos = recHitTools_.getPosition(hitId);
        ROOT::Math::XYZVector hit_xyz(hitPos.x(), hitPos.y(), hitPos.z());
        
        v_layerEnergy.at(hitLayer) += hitE;
        v_layerCentroid.at(hitLayer) += hitE * hit_xyz;
    }
    
    for(int iLayer = 0; iLayer < nLayer_; iLayer++)
    {
        if(v_layerEnergy.at(iLayer))
        {
            v_layerCentroid.at(iLayer) /= v_layerEnergy.at(iLayer);
        }
        
        //printf("Layer %d: energy %0.2f \n", iLayer+1, v_layerEnergy.at(iLayer));
    }
    
    
    // Compute Rvar in a cylinder around the centroids
    int iHit = -1;
    double Rvar = 0;
    
    for(auto HandF : v_superClus_HandF_selected)
    {
        iHit++;
        
        DetId hitId = HandF.first;
        
        int hitLayer = recHitTools_.getLayer(hitId) - 1;
        
        auto hitPos = recHitTools_.getPosition(hitId);
        ROOT::Math::XYZVector hit_xyz(hitPos.x(), hitPos.y(), hitPos.z());
        
        auto dist_xyz = hit_xyz - v_layerCentroid.at(hitLayer);
        
        double R = std::sqrt(dist_xyz.x()*dist_xyz.x() + dist_xyz.y()*dist_xyz.y());
        double cellSize = CommonUtilities::getCellSize(hitId, &recHitTools_);
        
        // Including the cellSize seems to make the variable less sensitive to HD/LD the transition region
        if(R > cylinderR_+cellSize)
        {
            continue;
        }
        
        double hitE = v_superClus_hitE_selected.at(iHit);
        
        Rvar += hitE;
        //v_layerEnergyInR.at(hitLayer) += hitE;
    }
    
    
    Rvar /= superClus_E;
    
    
    return Rvar;
}
