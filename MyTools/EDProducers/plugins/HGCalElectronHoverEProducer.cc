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
    
    std::string _instanceName;
    
    bool _debug;
    
    double _coneDR;
    double _minClusE;
    double _minClusET;
    
    edm::EDGetTokenT <std::vector <reco::GsfElectron> > _tok_electron;
    edm::EDGetTokenT <std::vector <reco::CaloCluster> > _tok_layerCluster;
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
    
    _instanceName = iConfig.getParameter <std::string>("instanceName");
    
    _tok_electron = consumes <std::vector <reco::GsfElectron> >(iConfig.getParameter <edm::InputTag>("electrons"));
    _tok_layerCluster = consumes <std::vector <reco::CaloCluster> >(iConfig.getParameter <edm::InputTag>("layerClusters"));
    
    _coneDR = iConfig.getParameter <double>("coneDR");
    
    _minClusE = iConfig.getParameter <double>("minClusE");
    _minClusET = iConfig.getParameter <double>("minClusET");
    
    _debug = iConfig.getParameter <bool>("debug");
    
    
    produces <std::vector <double> > (_instanceName);
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
    
    edm::Handle <std::vector <reco::GsfElectron> > v_electron;
    iEvent.getByToken(_tok_electron, v_electron);
    
    edm::Handle <std::vector <reco::CaloCluster> > v_layerCluster;
    iEvent.getByToken(_tok_layerCluster, v_layerCluster);
    
    
    int nEle = v_electron->size();
    int nLayerClus = v_layerCluster->size();
    
    std::vector <double> v_HoverE;
    
    for(int iEle = 0; iEle < nEle; iEle++)
    {
        reco::GsfElectron ele = v_electron->at(iEle);
        
        CLHEP::HepLorentzVector ele_4mom;
        ele_4mom.setT(ele.energy());
        ele_4mom.setX(ele.px());
        ele_4mom.setY(ele.py());
        ele_4mom.setZ(ele.pz());
        
        double HoverE = 0;
        
        for(int iClus = 0; iClus < nLayerClus; iClus++)
        {
            reco::CaloCluster cluster = v_layerCluster->at(iClus);
            
            // E cut
            if(cluster.energy() < _minClusE)
            {
                continue;
            }
            
            // ET cut
            double clusET = cluster.energy() * std::sin(cluster.position().theta());
            
            if(clusET < _minClusET)
            {
                continue;
            }
            
            CLHEP::Hep3Vector cluster_3vec(
                cluster.x(),
                cluster.y(),
                cluster.z()
            );
            
            // dR cut
            double dR = cluster_3vec.deltaR(ele_4mom.v());
            
            if(dR > _coneDR)
            {
                continue;
            }
            
            
            if(cluster.seed().det() == DetId::HGCalHSi || cluster.seed().det() == DetId::HGCalHSc)
            {
                HoverE += cluster.energy();
            }
        }
        
        
        HoverE /= ele.energy();
        
        if(_debug)
        {
            printf("In HGCalElectronHoverEProducer --> Ele %d/%d: H/E %0.4f \n", iEle+1, nEle, HoverE);
        }
        
        v_HoverE.push_back(HoverE);
    }
    
    
    iEvent.put(
        std::make_unique <std::vector <double> >(v_HoverE),
        _instanceName
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
