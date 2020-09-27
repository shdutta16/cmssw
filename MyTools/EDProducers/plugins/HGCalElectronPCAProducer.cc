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

# include "MyTools/EDProducers/plugins/CommonUtilities.h"


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
    
    std::string _instanceName;
    std::string _instanceName_UU;
    std::string _instanceName_VV;
    std::string _instanceName_WW;
    
    bool _debug;
    
    int _nLayer;
    double _cylinderR;
    double _minHitE;
    double _minHitET;
    
    edm::EDGetTokenT <std::vector <reco::GsfElectron> > _tok_electron;
    
    edm::EDGetTokenT <std::vector <reco::PFRecHit> > _tok_PFRecHit;
    
    edm::EDGetTokenT <edm::SortedCollection <HGCRecHit,edm::StrictWeakOrdering <HGCRecHit> > > _tok_HGCEERecHit;
    edm::EDGetTokenT <edm::SortedCollection <HGCRecHit,edm::StrictWeakOrdering <HGCRecHit> > > _tok_HGCHEFRecHit;
    edm::EDGetTokenT <edm::SortedCollection <HGCRecHit,edm::StrictWeakOrdering <HGCRecHit> > > _tok_HGCHEBRecHit;
    
    hgcal::RecHitTools recHitTools;
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
    
    _instanceName = iConfig.getParameter <std::string>("instanceName");
    
    _tok_electron = consumes <std::vector <reco::GsfElectron> >(iConfig.getParameter <edm::InputTag>("electrons"));
    
    _tok_PFRecHit = consumes <std::vector <reco::PFRecHit> >(iConfig.getParameter <edm::InputTag>("PFRecHits"));
    
    _tok_HGCEERecHit = consumes <edm::SortedCollection <HGCRecHit,edm::StrictWeakOrdering <HGCRecHit> > >(iConfig.getUntrackedParameter <edm::InputTag>("HGCEERecHits"));
    _tok_HGCHEFRecHit = consumes <edm::SortedCollection <HGCRecHit,edm::StrictWeakOrdering <HGCRecHit> > >(iConfig.getUntrackedParameter <edm::InputTag>("HGCHEFRecHits"));
    _tok_HGCHEBRecHit = consumes <edm::SortedCollection <HGCRecHit,edm::StrictWeakOrdering <HGCRecHit> > >(iConfig.getUntrackedParameter <edm::InputTag>("HGCHEBRecHits"));
    
    _nLayer = iConfig.getParameter <int>("nLayer");
    _cylinderR = iConfig.getParameter <double>("cylinderR");
    
    _minHitE = iConfig.getParameter <double>("minHitE");
    _minHitET = iConfig.getParameter <double>("minHitET");
    
    _debug = iConfig.getParameter <bool>("debug");
    
    
    _instanceName_UU = _instanceName + "Sigma2UU";
    _instanceName_VV = _instanceName + "Sigma2VV";
    _instanceName_WW = _instanceName + "Sigma2WW";
    
    produces <std::vector <double> > (_instanceName_UU);
    produces <std::vector <double> > (_instanceName_VV);
    produces <std::vector <double> > (_instanceName_WW);
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
    
    CommonUtilities::initRecHitTools(recHitTools, &iSetup);
    
    
    edm::Handle <std::vector <reco::PFRecHit> > v_PFRecHit;
    iEvent.getByToken(_tok_PFRecHit, v_PFRecHit);
    
    //edm::Handle <edm::SortedCollection <HGCRecHit,edm::StrictWeakOrdering <HGCRecHit> > > v_HGCEERecHit;
    //iEvent.getByToken(_tok_HGCEERecHit, v_HGCEERecHit);
    //
    //edm::Handle <edm::SortedCollection <HGCRecHit,edm::StrictWeakOrdering <HGCRecHit> > > v_HGCHEFRecHit;
    //iEvent.getByToken(_tok_HGCHEFRecHit, v_HGCHEFRecHit);
    //
    //edm::Handle <edm::SortedCollection <HGCRecHit,edm::StrictWeakOrdering <HGCRecHit> > > v_HGCHEBRecHit;
    //iEvent.getByToken(_tok_HGCHEBRecHit, v_HGCHEBRecHit);
    
    std::map <DetId, int> m_recHitIdx = CommonUtilities::getPFRecHitIndexMap(v_PFRecHit);
    //std::map <DetId, const HGCRecHit*> m_recHitPtr = CommonUtilities::getHGCRecHitPtrMap(v_HGCEERecHit, v_HGCHEFRecHit, v_HGCHEBRecHit);
    
    
    edm::Handle <std::vector <reco::GsfElectron> > v_electron;
    iEvent.getByToken(_tok_electron, v_electron);
    
    int nEle = v_electron->size();
    
    std::vector <double> v_sigma2UU;
    std::vector <double> v_sigma2VV;
    std::vector <double> v_sigma2WW;
    
    for(int iEle = 0; iEle < nEle; iEle++)
    {
        reco::GsfElectron ele = v_electron->at(iEle);
        
        std::vector <std::pair <DetId, float> > v_superClus_HandF = ele.superCluster()->hitsAndFractions();
        
        std::vector <double> v_layerEnergy(_nLayer, 0.0);
        std::vector <double> v_layerEnergyInR(_nLayer, 0.0);
        std::vector <double> v_layerRvar(_nLayer, 0.0);
        
        std::vector <double> v_superClus_isHitValid(v_superClus_HandF.size(), false);
        
        ROOT::Math::XYZVector centroid_xyz(0, 0, 0);
        std::vector <ROOT::Math::XYZVector> v_layerCentroid;
        
        for(int iLayer = 0; iLayer < _nLayer; iLayer++)
        {
            ROOT::Math::XYZVector xyz_temp(0, 0, 0);
            
            v_layerCentroid.push_back(xyz_temp);
        }
        
        
        double totalE = 0;
        
        // Compute the centroid per layer
        for(int iHit = 0; iHit < (int) v_superClus_HandF.size(); iHit++)
        {
            DetId hitId = v_superClus_HandF.at(iHit).first;
            DetId hitEfrac = v_superClus_HandF.at(iHit).second;
            
            int hitLayer = recHitTools.getLayer(hitId) - 1;
            
            if(hitId.det() != DetId::HGCalEE)
            {
                //printf("Det %d, layer %d \n", (int) hitId.det(), hitLayer);
                continue;
            }
            
            if(hitLayer+1 > _nLayer)
            {
                continue;
            }
            
            //if(m_recHitPtr.find(hitId) == m_recHitPtr.end())
            //{
            //    continue;
            //}
            
            int hitIdx = m_recHitIdx.at(hitId);
            
            if(v_PFRecHit->at(hitIdx).energy() < _minHitE)
            {
                continue;
            }
            
            if(sqrt(v_PFRecHit->at(hitIdx).pt2()) < _minHitET)
            {
                continue;
            }
            
            v_superClus_isHitValid.at(iHit) = true;
            
            double hitE = v_PFRecHit->at(hitIdx).energy() * hitEfrac;
            
            //const HGCRecHit *recHit = m_recHitPtr.at(hitId);
            //double hitE = recHit->energy() * hitEfrac;
            
            totalE += hitE;
            
            auto hitPos = recHitTools.getPosition(hitId);
            ROOT::Math::XYZVector hit_xyz(hitPos.x(), hitPos.y(), hitPos.z());
            
            v_layerEnergy.at(hitLayer) += hitE;
            
            v_layerCentroid.at(hitLayer) += hitE * hit_xyz;
            
            centroid_xyz += hitE * hit_xyz;
        }
        
        for(int iLayer = 0; iLayer < _nLayer; iLayer++)
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
            DetId hitEfrac = v_superClus_HandF.at(iHit).second;
            
            if(!v_superClus_isHitValid.at(iHit))
            {
                continue;
            }
            
            int hitIdx = m_recHitIdx.at(hitId);
            int hitLayer = recHitTools.getLayer(hitId) - 1;
            
            double hitE = v_PFRecHit->at(hitIdx).energy() * hitEfrac;
            
            auto hitPos = recHitTools.getPosition(hitId);
            ROOT::Math::XYZVector hit_xyz(hitPos.x(), hitPos.y(), hitPos.z());
            
            ROOT::Math::XYZVector distCyl_xyz = hit_xyz - v_layerCentroid.at(hitLayer);
            
            double rCyl = std::sqrt(distCyl_xyz.x()*distCyl_xyz.x() + distCyl_xyz.y()*distCyl_xyz.y());
            
            //double cellSize = CommonUtilities::getCellSize(hitId, &recHitTools);
            
            // Compute in a cylinder
            //if(rCyl > _cylinderR+cellSize)
            if(rCyl > _cylinderR)
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
        
        if(_debug)
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
        _instanceName_UU
    );
    
    iEvent.put(
        std::make_unique <std::vector <double> >(v_sigma2VV),
        _instanceName_VV
    );
    
    iEvent.put(
        std::make_unique <std::vector <double> >(v_sigma2WW),
        _instanceName_WW
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
