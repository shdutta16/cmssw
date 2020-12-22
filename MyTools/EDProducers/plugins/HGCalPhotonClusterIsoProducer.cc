// -*- C++ -*-
//
// Package:    MyTools/HGCalPhotonClusterIsoProducer
// Class:      HGCalPhotonClusterIsoProducer
//
/**\class HGCalPhotonClusterIsoProducer HGCalPhotonClusterIsoProducer.cc MyTools/HGCalPhotonClusterIsoProducer/plugins/HGCalPhotonClusterIsoProducer.cc

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

class HGCalPhotonClusterIsoProducer : public edm::stream::EDProducer<>
{
    public:
    
    explicit HGCalPhotonClusterIsoProducer(const edm::ParameterSet&);
    ~HGCalPhotonClusterIsoProducer();
    
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
    
    HGCalClusterTools algo_iso_;
    
    std::string instanceName_;
    
    bool debug_;
    
    edm::EDGetTokenT <std::vector <reco::Photon> > tok_photon_;
    edm::EDGetTokenT <std::vector <reco::CaloCluster> > tok_layerCluster_;
    
    double coneDRmin_;
    double coneDRmax_;
    double minEmClusE_;
    double minEmClusET_;
    double minHadClusE_;
    double minHadClusET_;
    
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
HGCalPhotonClusterIsoProducer::HGCalPhotonClusterIsoProducer(const edm::ParameterSet& iConfig)
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
    tok_layerCluster_ = consumes <std::vector <reco::CaloCluster> >(iConfig.getParameter <edm::InputTag>("layerClusters"));
    
    minPt_ = iConfig.getParameter <double>("minPt");
    
    coneDRmin_ = iConfig.getParameter <double>("coneDRmin");
    coneDRmax_ = iConfig.getParameter <double>("coneDRmax");
    minEmClusE_ = iConfig.getParameter <double>("minEmClusE");
    minEmClusET_ = iConfig.getParameter <double>("minEmClusET");
    minHadClusE_ = iConfig.getParameter <double>("minHadClusE");
    minHadClusET_ = iConfig.getParameter <double>("minHadClusET");
    
    debug_ = iConfig.getParameter <bool>("debug");
    
    
    produces <std::vector <double> > (instanceName_);
}

HGCalPhotonClusterIsoProducer::~HGCalPhotonClusterIsoProducer() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void HGCalPhotonClusterIsoProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
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
    
    edm::Handle <std::vector <reco::Photon> > h_photon;
    iEvent.getByToken(tok_photon_, h_photon);
    auto photons = *h_photon;
    
    edm::Handle <std::vector <reco::CaloCluster> > h_layerCluster;
    iEvent.getByToken(tok_layerCluster_, h_layerCluster);
    auto layerClusters = *h_layerCluster;
    
    int iPho = -1;
    int nPho = h_photon->size();
    
    std::vector <double> v_clusIso;
    
    for(auto &pho : photons)
    {
        iPho++;
        
        if(pho.pt() < minPt_)
        {
            v_clusIso.push_back(-99);
            continue;
        }
        
        double emEnergy = algo_iso_.emEnergyInCone(
            pho.superCluster()->eta(),
            pho.superCluster()->phi(),
            layerClusters,
            coneDRmin_,
            coneDRmax_,
            minEmClusET_,
            minEmClusE_,
            HGCalClusterTools::EType::ENERGY
        );
        
        double hadEnergy = algo_iso_.hadEnergyInCone(
            pho.superCluster()->eta(),
            pho.superCluster()->phi(),
            layerClusters,
            0.0,
            coneDRmax_,
            minHadClusET_,
            minHadClusE_,
            HGCalClusterTools::EType::ENERGY
        );
        
        double energyInCone = emEnergy+hadEnergy;
        
        if(debug_)
        {
            printf("In HGCalPhotonClusterIsoProducer --> Pho %d/%d: iso %0.4f \n", iPho+1, nPho, energyInCone);
        }
        
        v_clusIso.push_back(energyInCone);
    }
    
    
    iEvent.put(
        std::make_unique <std::vector <double> >(v_clusIso),
        instanceName_
    );
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void HGCalPhotonClusterIsoProducer::beginStream(edm::StreamID) {
  // please remove this method if not needed
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void HGCalPhotonClusterIsoProducer::endStream() {
  // please remove this method if not needed
}

// ------------ method called when starting to processes a run  ------------
/*
void
HGCalPhotonClusterIsoProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
HGCalPhotonClusterIsoProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
HGCalPhotonClusterIsoProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
HGCalPhotonClusterIsoProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HGCalPhotonClusterIsoProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCalPhotonClusterIsoProducer);
