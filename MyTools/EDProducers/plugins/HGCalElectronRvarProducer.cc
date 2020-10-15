// -*- C++ -*-
//
// Package:    MyTools/HGCalElectronRvarProducer
// Class:      HGCalElectronRvarProducer
//
/**\class HGCalElectronRvarProducer HGCalElectronRvarProducer.cc MyTools/HGCalElectronRvarProducer/plugins/HGCalElectronRvarProducer.cc

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

# include "MyTools/EDProducers/interface/HGCalAlgoSuperClusterRvar.h"


//
// class declaration
//

class HGCalElectronRvarProducer : public edm::stream::EDProducer<>
{
    public:
    
    explicit HGCalElectronRvarProducer(const edm::ParameterSet&);
    ~HGCalElectronRvarProducer();
    
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
    
    HGCalAlgoSuperClusterRvar algo_Rvar_;
    
    std::string instanceName_;
    
    bool debug_;
    
    edm::EDGetTokenT <std::vector <reco::GsfElectron> > tok_electron_;
    edm::EDGetTokenT <std::vector <reco::PFRecHit> > tok_PFRecHit_;
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
HGCalElectronRvarProducer::HGCalElectronRvarProducer(const edm::ParameterSet& iConfig) :
    algo_Rvar_(iConfig)
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
    
    debug_ = iConfig.getParameter <bool>("debug");
    
    
    produces <std::vector <double> > (instanceName_);
}

HGCalElectronRvarProducer::~HGCalElectronRvarProducer() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void HGCalElectronRvarProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
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
    
    edm::Handle <std::vector <reco::PFRecHit> > h_PFRecHit;
    iEvent.getByToken(tok_PFRecHit_, h_PFRecHit);
    auto recHits = *h_PFRecHit;
    
    edm::Handle <std::vector <reco::GsfElectron> > h_electron;
    iEvent.getByToken(tok_electron_, h_electron);
    auto electrons = *h_electron;
    
    int iEle = -1;
    int nEle = h_electron->size();
    
    std::vector <double> v_Rvar;
    
    for(auto &ele : electrons)
    {
        iEle++;
        
        double Rvar = algo_Rvar_.getRvar(
            iSetup,
            ele.superCluster(),
            recHits
        );
        
        if(debug_)
        {
            printf("In HGCalElectronRvarProducer --> Ele %d/%d: Rvar %0.4f \n", iEle+1, nEle, Rvar);
        }
        
        v_Rvar.push_back(Rvar);
    }
    
    
    iEvent.put(
        std::make_unique <std::vector <double> >(v_Rvar),
        instanceName_
    );
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void HGCalElectronRvarProducer::beginStream(edm::StreamID) {
  // please remove this method if not needed
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void HGCalElectronRvarProducer::endStream() {
  // please remove this method if not needed
}

// ------------ method called when starting to processes a run  ------------
/*
void
HGCalElectronRvarProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
HGCalElectronRvarProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
HGCalElectronRvarProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
HGCalElectronRvarProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HGCalElectronRvarProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCalElectronRvarProducer);
