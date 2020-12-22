// -*- C++ -*-
//
// Package:    MyTools/HGCalElectronClusterIsoProducer
// Class:      HGCalElectronClusterIsoProducer
//
/**\class HGCalElectronClusterIsoProducer HGCalElectronClusterIsoProducer.cc MyTools/HGCalElectronClusterIsoProducer/plugins/HGCalElectronClusterIsoProducer.cc

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

class HGCalElectronClusterIsoProducer : public edm::stream::EDProducer<>
{
    public:
    
    explicit HGCalElectronClusterIsoProducer(const edm::ParameterSet&);
    ~HGCalElectronClusterIsoProducer();
    
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
    
    edm::EDGetTokenT <std::vector <reco::GsfElectron> > tok_electron_;
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
HGCalElectronClusterIsoProducer::HGCalElectronClusterIsoProducer(const edm::ParameterSet& iConfig)
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
    
    coneDRmin_ = iConfig.getParameter <double>("coneDRmin");
    coneDRmax_ = iConfig.getParameter <double>("coneDRmax");
    minEmClusE_ = iConfig.getParameter <double>("minEmClusE");
    minEmClusET_ = iConfig.getParameter <double>("minEmClusET");
    minHadClusE_ = iConfig.getParameter <double>("minHadClusE");
    minHadClusET_ = iConfig.getParameter <double>("minHadClusET");
    
    debug_ = iConfig.getParameter <bool>("debug");
    
    
    produces <std::vector <double> > (instanceName_);
}

HGCalElectronClusterIsoProducer::~HGCalElectronClusterIsoProducer() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void HGCalElectronClusterIsoProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
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
    
    std::vector <double> v_clusIso;
    
    for(auto &ele : electrons)
    {
        iEle++;
        
        if(ele.pt() < minPt_)
        {
            v_clusIso.push_back(-99);
            continue;
        }
        
        double emEnergy = algo_iso_.emEnergyInCone(
            ele.superCluster()->eta(),
            ele.superCluster()->phi(),
            layerClusters,
            coneDRmin_,
            coneDRmax_,
            minEmClusET_,
            minEmClusE_,
            HGCalClusterTools::EType::ENERGY
        );
        
        double hadEnergy = algo_iso_.hadEnergyInCone(
            ele.superCluster()->eta(),
            ele.superCluster()->phi(),
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
            printf("In HGCalElectronClusterIsoProducer --> Ele %d/%d: iso %0.4f \n", iEle+1, nEle, energyInCone);
        }
        
        v_clusIso.push_back(energyInCone);
    }
    
    
    iEvent.put(
        std::make_unique <std::vector <double> >(v_clusIso),
        instanceName_
    );
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void HGCalElectronClusterIsoProducer::beginStream(edm::StreamID) {
  // please remove this method if not needed
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void HGCalElectronClusterIsoProducer::endStream() {
  // please remove this method if not needed
}

// ------------ method called when starting to processes a run  ------------
/*
void
HGCalElectronClusterIsoProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
HGCalElectronClusterIsoProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
HGCalElectronClusterIsoProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
HGCalElectronClusterIsoProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HGCalElectronClusterIsoProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCalElectronClusterIsoProducer);
