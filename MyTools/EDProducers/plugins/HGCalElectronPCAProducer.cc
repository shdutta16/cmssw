// -*- C++ -*-
//
// Package:    MyTools/HGCalElectronPCAProducer
// Class:      HGCalElectronPCAProducer
//
/**\class HGCalElectronPCAProducer HGCalElectronPCAProducer.cc MyTools/HGCalElectronPCAProducer/plugins/HGCalElectronPCAProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Soham Bhattacharya
//         Created:  Wed, 16 Sep 2020 14:18:20 GMT
//
//

// system include files
#include <memory>

// user include files
# include "DataFormats/CaloRecHit/interface/CaloCluster.h"
# include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
# include "DataFormats/FWLite/interface/ESHandle.h"
# include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
# include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
# include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
# include "DataFormats/HGCRecHit/interface/HGCRecHit.h"
# include "DataFormats/Math/interface/LorentzVector.h"
# include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
# include "DataFormats/ParticleFlowReco/interface/PFRecHitFraction.h"
# include "DataFormats/TrackReco/interface/Track.h"
# include "DataFormats/TrackReco/interface/TrackFwd.h"
# include "FWCore/Framework/interface/Event.h"
# include "FWCore/Framework/interface/Frameworkfwd.h"
# include "FWCore/Framework/interface/MakerMacros.h"
# include "FWCore/Framework/interface/stream/EDProducer.h"
# include "FWCore/ParameterSet/interface/ParameterSet.h"
# include "FWCore/Utilities/interface/StreamID.h"
# include "Geometry/CaloTopology/interface/HGCalTopology.h"
# include "Geometry/Records/interface/IdealGeometryRecord.h"
# include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
# include "RecoParticleFlow/PFClusterProducer/interface/InitialClusteringStepBase.h"

# include <CLHEP/Vector/LorentzVector.h>
# include <Math/VectorUtil.h>

# include <algorithm>
# include <iostream>
# include <map>
# include <stdlib.h>
# include <string>
# include <type_traits>
# include <utility>
# include <vector>

# include "MyTools/EDProducers/interface/CommonUtilities.h"


//
// class declaration
//

class HGCalElectronPCAProducer : public edm::stream::EDProducer<>
{
    public:
    
    explicit HGCalElectronPCAProducer(const edm::ParameterSet&);
    ~HGCalElectronPCAProducer();
    
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    
    
    private:
    
    void beginStream(edm::StreamID) override;
    void produce(edm::Event&, const edm::EventSetup&) override;
    void endStream() override;
    
    //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    
    // ----------member data ---------------------------
    
    std::string instanceName_;
    std::string instanceName_UU_;
    std::string instanceName_VV_;
    std::string instanceName_WW_;
    
    bool debug_;
    
    int nLayer_;
    double cylinderR_;
    double minHitE_;
    double minHitET_;
    
    edm::EDGetTokenT <std::vector <reco::GsfElectron> > tok_electron_;
    
    edm::EDGetTokenT <std::vector <reco::PFRecHit> > tok_PFRecHit_;
    
    edm::EDGetTokenT <edm::SortedCollection <HGCRecHit,edm::StrictWeakOrdering <HGCRecHit> > > tok_HGCEERecHit_;
    edm::EDGetTokenT <edm::SortedCollection <HGCRecHit,edm::StrictWeakOrdering <HGCRecHit> > > tok_HGCHEFRecHit_;
    edm::EDGetTokenT <edm::SortedCollection <HGCRecHit,edm::StrictWeakOrdering <HGCRecHit> > > tok_HGCHEBRecHit_;
    
    hgcal::RecHitTools recHitTools_;
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
HGCalElectronPCAProducer::HGCalElectronPCAProducer(const edm::ParameterSet& iConfig)
{
    //register your products
    /* Examples
    produces<ExampleData2>();
    
    //if do put with a label
    produces<ExampleData2>("label");
    
    //if you want to put into the Run
    produces<ExampleData2,InRun>();
    */
    //now do what ever other initialization is needed
    
    instanceName_ = iConfig.getParameter <std::string>("instanceName");
    
    tok_electron_ = consumes <std::vector <reco::GsfElectron> >(iConfig.getParameter <edm::InputTag>("electrons"));
    
    tok_PFRecHit_ = consumes <std::vector <reco::PFRecHit> >(iConfig.getParameter <edm::InputTag>("PFRecHits"));
    
    tok_HGCEERecHit_ = consumes <edm::SortedCollection <HGCRecHit,edm::StrictWeakOrdering <HGCRecHit> > >(iConfig.getUntrackedParameter <edm::InputTag>("HGCEERecHits"));
    tok_HGCHEFRecHit_ = consumes <edm::SortedCollection <HGCRecHit,edm::StrictWeakOrdering <HGCRecHit> > >(iConfig.getUntrackedParameter <edm::InputTag>("HGCHEFRecHits"));
    tok_HGCHEBRecHit_ = consumes <edm::SortedCollection <HGCRecHit,edm::StrictWeakOrdering <HGCRecHit> > >(iConfig.getUntrackedParameter <edm::InputTag>("HGCHEBRecHits"));
    
    nLayer_ = iConfig.getParameter <int>("nLayer");
    cylinderR_ = iConfig.getParameter <double>("cylinderR");
    
    minHitE_ = iConfig.getParameter <double>("minHitE");
    minHitET_ = iConfig.getParameter <double>("minHitET");
    
    debug_ = iConfig.getParameter <bool>("debug");
    
    
    instanceName_UU_ = instanceName_ + "Sigma2UU";
    instanceName_VV_ = instanceName_ + "Sigma2VV";
    instanceName_WW_ = instanceName_ + "Sigma2WW";
    
    produces <std::vector <double> > (instanceName_UU_);
    produces <std::vector <double> > (instanceName_VV_);
    produces <std::vector <double> > (instanceName_WW_);
}

HGCalElectronPCAProducer::~HGCalElectronPCAProducer() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void HGCalElectronPCAProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    /* This is an event example
    //Read 'ExampleData' from the Event
    ExampleData const& in = iEvent.get(inToken_);
    
    //Use the ExampleData to create an ExampleData2 which 
    // is put into the Event
    iEvent.put(std::make_unique<ExampleData2>(in));
    */
    
    /* this is an EventSetup example
    //Read SetupData from the SetupRecord in the EventSetup
    SetupData& setup = iSetup.getData(setupToken_);
    */
    
    CommonUtilities::initRecHitTools(recHitTools_, &iSetup);
    
    
    edm::Handle <std::vector <reco::PFRecHit> > v_PFRecHit;
    iEvent.getByToken(tok_PFRecHit_, v_PFRecHit);
    
    //edm::Handle <edm::SortedCollection <HGCRecHit,edm::StrictWeakOrdering <HGCRecHit> > > v_HGCEERecHit;
    //iEvent.getByToken(tok_HGCEERecHit_, v_HGCEERecHit);
    //
    //edm::Handle <edm::SortedCollection <HGCRecHit,edm::StrictWeakOrdering <HGCRecHit> > > v_HGCHEFRecHit;
    //iEvent.getByToken(tok_HGCHEFRecHit_, v_HGCHEFRecHit);
    //
    //edm::Handle <edm::SortedCollection <HGCRecHit,edm::StrictWeakOrdering <HGCRecHit> > > v_HGCHEBRecHit;
    //iEvent.getByToken(tok_HGCHEBRecHit_, v_HGCHEBRecHit);
    
    std::map <DetId, int> m_recHitIdx = CommonUtilities::getPFRecHitIndexMap(v_PFRecHit);
    //std::map <DetId, const HGCRecHit*> m_recHitPtr = CommonUtilities::getHGCRecHitPtrMap(v_HGCEERecHit, v_HGCHEFRecHit, v_HGCHEBRecHit);
    
    
    edm::Handle <std::vector <reco::GsfElectron> > v_electron;
    iEvent.getByToken(tok_electron_, v_electron);
    
    int nEle = v_electron->size();
    
    std::vector <double> v_sigma2UU;
    std::vector <double> v_sigma2VV;
    std::vector <double> v_sigma2WW;
    
    for(int iEle = 0; iEle < nEle; iEle++)
    {
        reco::GsfElectron ele = v_electron->at(iEle);
        
        std::vector <std::pair <DetId, float> > v_superClus_HandF = ele.superCluster()->hitsAndFractions();
        
        std::vector <double> v_layerEnergy(nLayer_, 0.0);
        std::vector <double> v_layerEnergyInR(nLayer_, 0.0);
        std::vector <double> v_layerRvar(nLayer_, 0.0);
        
        std::vector <double> v_superClus_isHitValid(v_superClus_HandF.size(), false);
        
        ROOT::Math::XYZVector centroid_xyz(0, 0, 0);
        std::vector <ROOT::Math::XYZVector> v_layerCentroid;
        
        for(int iLayer = 0; iLayer < nLayer_; iLayer++)
        {
            ROOT::Math::XYZVector xyz_temp(0, 0, 0);
            
            v_layerCentroid.push_back(xyz_temp);
        }
        
        
        double totalE = 0;
        
        // Compute the centroid per layer
        for(int iHit = 0; iHit < (int) v_superClus_HandF.size(); iHit++)
        {
            DetId hitId = v_superClus_HandF.at(iHit).first;
            double hitEfrac = v_superClus_HandF.at(iHit).second;
            
            int hitLayer = recHitTools_.getLayer(hitId) - 1;
            
            if(hitId.det() != DetId::HGCalEE)
            {
                //printf("Det %d, layer %d \n", (int) hitId.det(), hitLayer);
                continue;
            }
            
            if(hitLayer+1 > nLayer_)
            {
                continue;
            }
            
            //if(m_recHitPtr.find(hitId) == m_recHitPtr.end())
            //{
            //    continue;
            //}
            
            int hitIdx = m_recHitIdx.at(hitId);
            
            if(v_PFRecHit->at(hitIdx).energy() < minHitE_)
            {
                continue;
            }
            
            if(sqrt(v_PFRecHit->at(hitIdx).pt2()) < minHitET_)
            {
                continue;
            }
            
            v_superClus_isHitValid.at(iHit) = true;
            
            double hitE = v_PFRecHit->at(hitIdx).energy() * hitEfrac;
            
            //const HGCRecHit *recHit = m_recHitPtr.at(hitId);
            //double hitE = recHit->energy() * hitEfrac;
            
            totalE += hitE;
            
            auto hitPos = recHitTools_.getPosition(hitId);
            ROOT::Math::XYZVector hit_xyz(hitPos.x(), hitPos.y(), hitPos.z());
            
            v_layerEnergy.at(hitLayer) += hitE;
            
            v_layerCentroid.at(hitLayer) += hitE * hit_xyz;
            
            centroid_xyz += hitE * hit_xyz;
        }
        
        for(int iLayer = 0; iLayer < nLayer_; iLayer++)
        {
            if(v_layerEnergy.at(iLayer))
            {
                v_layerCentroid.at(iLayer) /= v_layerEnergy.at(iLayer);
            }
            
            //printf("Layer %d: energy %0.2f \n", iLayer+1, v_layerEnergy.at(iLayer));
        }
        
        if(totalE)
        {
            centroid_xyz /= totalE;
        }
        
        
        TMatrixD mat_cov(3, 3);
        //TMatrixD mat_dist(3, nValidHit);
        
        double dxdx = 0;
        double dydy = 0;
        double dzdz = 0;

        double dxdy = 0;
        double dydz = 0;
        double dzdx = 0;
        
        double weightTot = 0;
        
        for(int iHit = 0; iHit < (int) v_superClus_HandF.size(); iHit++)
        {
            DetId hitId = v_superClus_HandF.at(iHit).first;
            double hitEfrac = v_superClus_HandF.at(iHit).second;
            
            if(!v_superClus_isHitValid.at(iHit))
            {
                continue;
            }
            
            int hitIdx = m_recHitIdx.at(hitId);
            int hitLayer = recHitTools_.getLayer(hitId) - 1;
            
            double hitE = v_PFRecHit->at(hitIdx).energy() * hitEfrac;
            
            auto hitPos = recHitTools_.getPosition(hitId);
            ROOT::Math::XYZVector hit_xyz(hitPos.x(), hitPos.y(), hitPos.z());
            
            ROOT::Math::XYZVector distCyl_xyz = hit_xyz - v_layerCentroid.at(hitLayer);
            
            double rCyl = std::sqrt(distCyl_xyz.x()*distCyl_xyz.x() + distCyl_xyz.y()*distCyl_xyz.y());
            
            //double cellSize = CommonUtilities::getCellSize(hitId, &recHitTools_);
            
            // Compute in a cylinder
            //if(rCyl > cylinderR_+cellSize)
            if(rCyl > cylinderR_)
            {
                continue;
            }
            
            ROOT::Math::XYZVector dist_xyz = hit_xyz - centroid_xyz;
            
            double weight = hitE;
            weightTot += weight;
            
            dxdx += weight * dist_xyz.x() * dist_xyz.x();
            dydy += weight * dist_xyz.y() * dist_xyz.y();
            dzdz += weight * dist_xyz.z() * dist_xyz.z();
            
            dxdy += weight * dist_xyz.x() * dist_xyz.y();
            dydz += weight * dist_xyz.y() * dist_xyz.z();
            dzdx += weight * dist_xyz.z() * dist_xyz.x();
            
            
            //mat_dist(0, hitCount) = dist_xyz.x();
            //mat_dist(1, hitCount) = dist_xyz.y();
            //mat_dist(2, hitCount) = dist_xyz.z();
            //
            //
            //hitCount++;
        }
        
        dxdx /= weightTot;
        dydy /= weightTot;
        dzdz /= weightTot;
        
        dxdy /= weightTot;
        dydz /= weightTot;
        dzdx /= weightTot;
        
        mat_cov(0, 0) = dxdx;
        mat_cov(1, 1) = dydy;
        mat_cov(2, 2) = dzdz;
        
        mat_cov(0, 1) = mat_cov(1, 0) = dxdy;
        mat_cov(0, 2) = mat_cov(2, 0) = dzdx;
        mat_cov(1, 2) = mat_cov(2, 1) = dydz;
        
        // Get eigen values and vectors
        TVectorD v_eigVal(3);
        TMatrixD mat_eigVec(3, 3);
        
        if(weightTot > 0)// && ele.superCluster()->energy() > 0)
        {
            try
            {
                mat_eigVec = mat_cov.EigenVectors(v_eigVal);
            }
            
            catch(...)
            {
                printf("Warning in HGCalElectronPCAProducer::produce(...): Cannot get eigen values and vectors. \n");
                
                printf("Cov. matrix: \n");
                mat_cov.Print();
                
                fflush(stdout);
                fflush(stderr);
            }
        }
        
        if(debug_)
        {
            printf(
                "In HGCalElectronPCAProducer::produce(...) --> Ele %d/%d: eig vals (%0.4f, %0.4f, %0.4f) "
                "\n",
                
                iEle+1, nEle,
                v_eigVal(0),
                v_eigVal(1),
                v_eigVal(2)
            );
            
            printf("Cov. matrix: \n");
            mat_cov.Print();
        }
        
        //v_eigenVal.at(0) = v_eigVal(0);
        //v_eigenVal.at(1) = v_eigVal(1);
        //v_eigenVal.at(2) = v_eigVal(2);
        
        
        v_sigma2UU.push_back(v_eigVal(1));
        v_sigma2VV.push_back(v_eigVal(2));
        v_sigma2WW.push_back(v_eigVal(0));
    }
    
    
    iEvent.put(
        std::make_unique <std::vector <double> >(v_sigma2UU),
        instanceName_UU_
    );
    
    iEvent.put(
        std::make_unique <std::vector <double> >(v_sigma2VV),
        instanceName_VV_
    );
    
    iEvent.put(
        std::make_unique <std::vector <double> >(v_sigma2WW),
        instanceName_WW_
    );
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void HGCalElectronPCAProducer::beginStream(edm::StreamID) {
  // please remove this method if not needed
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void HGCalElectronPCAProducer::endStream() {
  // please remove this method if not needed
}

// ------------ method called when starting to processes a run  ------------
/*
void
HGCalElectronPCAProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
HGCalElectronPCAProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
HGCalElectronPCAProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
HGCalElectronPCAProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HGCalElectronPCAProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCalElectronPCAProducer);
