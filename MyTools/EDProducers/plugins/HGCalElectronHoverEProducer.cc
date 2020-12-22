// -*- C++ -*-
//
// Package:    MyTools/HGCalElectronHoverEProducer
// Class:      HGCalElectronHoverEProducer
//
/**\class HGCalElectronHoverEProducer HGCalElectronHoverEProducer.cc MyTools/HGCalElectronHoverEProducer/plugins/HGCalElectronHoverEProducer.cc

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

# include "RecoEgamma/EgammaTools/interface/HGCalClusterTools.h"

//
// class declaration
//

class HGCalElectronHoverEProducer : public edm::stream::EDProducer<>
{
    public:
    
    explicit HGCalElectronHoverEProducer(const edm::ParameterSet&);
    ~HGCalElectronHoverEProducer();
    
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
    
    HGCalClusterTools algo_HoverE_;
    
    std::string instanceName_;
    
    bool debug_;
    
    edm::EDGetTokenT <std::vector <reco::GsfElectron> > tok_electron_;
    edm::EDGetTokenT <std::vector <reco::CaloCluster> > tok_layerCluster_;
    
    double coneDR_;
    double minClusE_;
    double minClusET_;
    
    double minPt_;
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
HGCalElectronHoverEProducer::HGCalElectronHoverEProducer(const edm::ParameterSet& iConfig)
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
    tok_layerCluster_ = consumes <std::vector <reco::CaloCluster> >(iConfig.getParameter <edm::InputTag>("layerClusters"));
    
    minPt_ = iConfig.getParameter <double>("minPt");
    
    coneDR_ = iConfig.getParameter <double>("coneDR");
    minClusE_ = iConfig.getParameter <double>("minClusE");
    minClusET_ = iConfig.getParameter <double>("minClusET");
    
    debug_ = iConfig.getParameter <bool>("debug");
    
    
    produces <std::vector <double> > (instanceName_);
}

HGCalElectronHoverEProducer::~HGCalElectronHoverEProducer() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void HGCalElectronHoverEProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
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
    
    edm::Handle <std::vector <reco::GsfElectron> > h_electron;
    iEvent.getByToken(tok_electron_, h_electron);
    auto electrons = *h_electron;
    
    edm::Handle <std::vector <reco::CaloCluster> > h_layerCluster;
    iEvent.getByToken(tok_layerCluster_, h_layerCluster);
    auto layerClusters = *h_layerCluster;
    
    int iEle = -1;
    int nEle = h_electron->size();
    
    std::vector <double> v_HoverE;
    
    for(auto &ele : electrons)
    {
        iEle++;
        
        if(ele.pt() < minPt_)
        {
            v_HoverE.push_back(-99);
            continue;
        }
        
        double HoverE = algo_HoverE_.hadEnergyInCone(
            ele.superCluster()->eta(),
            ele.superCluster()->phi(),
            layerClusters,
            0.0,
            coneDR_,
            minClusET_,
            minClusE_,
            HGCalClusterTools::EType::ENERGY
        );
        
        HoverE /= ele.superCluster()->energy();
        
        if(debug_)
        {
            printf("In HGCalElectronHoverEProducer --> Ele %d/%d: H/E %0.4f \n", iEle+1, nEle, HoverE);
        }
        
        v_HoverE.push_back(HoverE);
    }
    
    
    iEvent.put(
        std::make_unique <std::vector <double> >(v_HoverE),
        instanceName_
    );
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void HGCalElectronHoverEProducer::beginStream(edm::StreamID) {
  // please remove this method if not needed
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void HGCalElectronHoverEProducer::endStream() {
  // please remove this method if not needed
}

// ------------ method called when starting to processes a run  ------------
/*
void
HGCalElectronHoverEProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
HGCalElectronHoverEProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
HGCalElectronHoverEProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
HGCalElectronHoverEProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HGCalElectronHoverEProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCalElectronHoverEProducer);
