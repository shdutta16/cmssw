// -*- C++ -*-
//
// Package:    MyTools/HGCalPhotonPCAProducer
// Class:      HGCalPhotonPCAProducer
//
/**\class HGCalPhotonPCAProducer HGCalPhotonPCAProducer.cc MyTools/HGCalPhotonPCAProducer/plugins/HGCalPhotonPCAProducer.cc

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
# include "DataFormats/EgammaCandidates/interface/Photon.h"
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

//# include "MyTools/EDProducers/interface/CommonUtilities.h"
# include "RecoEgamma/EgammaTools/interface/HGCalShowerShapeHelper.h"


//
// class declaration
//

class HGCalPhotonPCAProducer : public edm::stream::EDProducer<>
{
    public:
    
    explicit HGCalPhotonPCAProducer(const edm::ParameterSet&);
    ~HGCalPhotonPCAProducer();
    
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
    
    HGCalShowerShapeHelper showerShapeHelper_;
    
    std::string instanceName_;
    std::string instanceName_UU_;
    std::string instanceName_VV_;
    std::string instanceName_WW_;
    
    bool debug_;
    
    double cylinderR_;
    double minHitE_;
    double minHitET_;
    
    double minPt_;
    
    edm::EDGetTokenT <std::vector <reco::Photon> > tok_photon_;
    
    edm::EDGetTokenT <std::vector <reco::PFRecHit> > tok_PFRecHit_;
    
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
HGCalPhotonPCAProducer::HGCalPhotonPCAProducer(const edm::ParameterSet& iConfig)
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
    
    tok_photon_ = consumes <std::vector <reco::Photon> >(iConfig.getParameter <edm::InputTag>("photons"));
    
    tok_PFRecHit_ = consumes <std::vector <reco::PFRecHit> >(iConfig.getParameter <edm::InputTag>("PFRecHits"));
    
    cylinderR_ = iConfig.getParameter <double>("cylinderR");
    minHitE_ = iConfig.getParameter <double>("minHitE");
    minHitET_ = iConfig.getParameter <double>("minHitET");
    
    minPt_ = iConfig.getParameter <double>("minPt");
    
    debug_ = iConfig.getParameter <bool>("debug");
    
    
    instanceName_UU_ = instanceName_ + "Sigma2UU";
    instanceName_VV_ = instanceName_ + "Sigma2VV";
    instanceName_WW_ = instanceName_ + "Sigma2WW";
    
    produces <std::vector <double> > (instanceName_UU_);
    produces <std::vector <double> > (instanceName_VV_);
    produces <std::vector <double> > (instanceName_WW_);
}

HGCalPhotonPCAProducer::~HGCalPhotonPCAProducer() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void HGCalPhotonPCAProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
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
    
    //CommonUtilities::initRecHitTools(recHitTools_, &iSetup);
    
    
    edm::Handle <std::vector <reco::PFRecHit> > h_PFRecHit;
    iEvent.getByToken(tok_PFRecHit_, h_PFRecHit);
    auto recHits = *h_PFRecHit;
    
    showerShapeHelper_.initPerEvent(
        iSetup,
        recHits
    );
    
    edm::Handle <std::vector <reco::Photon> > v_photon;
    iEvent.getByToken(tok_photon_, v_photon);
    
    int nPho = v_photon->size();
    
    std::vector <double> v_sigma2UU;
    std::vector <double> v_sigma2VV;
    std::vector <double> v_sigma2WW;
    
    for(int iPho = 0; iPho < nPho; iPho++)
    {
        reco::Photon pho = v_photon->at(iPho);
        
        if(pho.pt() < minPt_)
        {
            v_sigma2UU.push_back(-99);
            v_sigma2VV.push_back(-99);
            v_sigma2WW.push_back(-99);
            continue;
        }
        
        showerShapeHelper_.initPerObject(
            pho.superCluster()->hitsAndFractions(),
            minHitE_,
            minHitET_
        );
        
        auto pcaWidths = showerShapeHelper_.getPCAWidths(
            cylinderR_
        );
        
        if(debug_)
        {
            printf("In HGCalPhotonPCAProducer --> Pho %d/%d: s2uu %0.4f, s2vv %0.4f, s2ww %0.4f \n", iPho+1, nPho, pcaWidths.sigma2uu, pcaWidths.sigma2vv, pcaWidths.sigma2ww);
        }
        
        
        v_sigma2UU.push_back(pcaWidths.sigma2uu);
        v_sigma2VV.push_back(pcaWidths.sigma2vv);
        v_sigma2WW.push_back(pcaWidths.sigma2ww);
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
void HGCalPhotonPCAProducer::beginStream(edm::StreamID) {
  // please remove this method if not needed
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void HGCalPhotonPCAProducer::endStream() {
  // please remove this method if not needed
}

// ------------ method called when starting to processes a run  ------------
/*
void
HGCalPhotonPCAProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
HGCalPhotonPCAProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
HGCalPhotonPCAProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
HGCalPhotonPCAProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HGCalPhotonPCAProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCalPhotonPCAProducer);
