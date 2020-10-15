# include "MyTools/EDProducers/interface/CommonUtilities.h"



namespace CommonUtilities
{
    void initRecHitTools(
        hgcal::RecHitTools &recHitTools,
        const edm::EventSetup *iSetup
    )
    {
        // In 11_1
        if(RecHitTools_has_getEventSetup<hgcal::RecHitTools>::value)
        {
            RecHitTools_has_getEventSetup<hgcal::RecHitTools>::eval(recHitTools, iSetup);
        }
        
        // In 11_2
        else if(RecHitTools_has_setGeometry<hgcal::RecHitTools>::value)
        {
            RecHitTools_has_setGeometry<hgcal::RecHitTools>::eval(recHitTools, iSetup);
        }
        
        else
        {
            printf("Error in CommonUtilities::initRecHitTools(...): cannot initialize. \n");
            
            fflush(stderr);
            fflush(stdout);
            
            exit(EXIT_FAILURE);
        }
    }
    
    
    std::map <DetId, int> getPFRecHitIndexMap(edm::Handle <std::vector <reco::PFRecHit> > v_recHit)
    {
        std::map <DetId, int> m_recHitIdx;
        
        int nHit = v_recHit->size();
        
        for(int iHit = 0; iHit < nHit; iHit++)
        {
            const reco::PFRecHit recHit = v_recHit->at(iHit);
            
            DetId hitId(recHit.detId());
            
            m_recHitIdx[hitId] = iHit;
        }
        
        return m_recHitIdx;
    }
    
    
    std::map <DetId, const reco::PFRecHit*> getPFRecHitPtrMap(
        const std::vector <reco::PFRecHit> &v_recHit
    )
    {
        std::map <DetId, const reco::PFRecHit*> m_recHitPtr;
        
        for(auto &recHit : v_recHit)
        {
            m_recHitPtr[recHit.detId()] = &recHit;
        }
        
        return m_recHitPtr;
    }
    
    
    std::map <DetId, const HGCRecHit*> getHGCRecHitPtrMap(
        edm::Handle <edm::SortedCollection <HGCRecHit,edm::StrictWeakOrdering <HGCRecHit> > > v_HGCEERecHit,
        edm::Handle <edm::SortedCollection <HGCRecHit,edm::StrictWeakOrdering <HGCRecHit> > > v_HGCHEFRecHit,
        edm::Handle <edm::SortedCollection <HGCRecHit,edm::StrictWeakOrdering <HGCRecHit> > > v_HGCHEBRecHit
    )
    {
        std::map <DetId, const HGCRecHit*> m_recHit;
        
        int nHGCEERecHit = v_HGCEERecHit->size();
        
        for(int iRecHit = 0; iRecHit < nHGCEERecHit; iRecHit++)
        {
            const HGCRecHit *recHit = &(*v_HGCEERecHit)[iRecHit];
            
            m_recHit[recHit->id()] = recHit;
        }
        
        
        //
        int nHGCHEFRecHit = v_HGCHEFRecHit->size();
        
        for(int iRecHit = 0; iRecHit < nHGCHEFRecHit; iRecHit++)
        {
            const HGCRecHit *recHit = &(*v_HGCHEFRecHit)[iRecHit];
            
            m_recHit[recHit->id()] = recHit;
        }
        
        
        //
        int nHGCHEBRecHit = v_HGCHEBRecHit->size();
        
        for(int iRecHit = 0; iRecHit < nHGCHEBRecHit; iRecHit++)
        {
            const HGCRecHit *recHit = &(*v_HGCHEBRecHit)[iRecHit];
            
            m_recHit[recHit->id()] = recHit;
        }
        
        
        return m_recHit;
    }
    
    
    double getCellSize(DetId detId, const hgcal::RecHitTools *recHitTools)
    {
        double SiThickness = recHitTools->getSiThickness(detId);
        
        // HD wafers
        if(SiThickness < 150)
        {
            return 0.465;
        }
        
        // LD wafers
        else
        {
            return 0.698;
        }
    }
}
