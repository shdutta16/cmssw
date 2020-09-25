// -*- C++ -*-
//
// Package:    MyTools/MapProducer
// Class:      MapProducer
//
/**\class MapProducer MapProducer.cc MyTools/MapProducer/plugins/MapProducer.cc

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
# include "DataFormats/Common/interface/MapOfVectors.h"
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

//# include "MyTools/EDProducers/plugins/CommonUtilities.h"



//
// class declaration
//

class MapProducer : public edm::stream::EDProducer<>
{
    public:
    
    explicit MapProducer(const edm::ParameterSet&);
    ~MapProducer();
    
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
    bool _useProcessName;
    
    int _nLayer;
    double _cylinderR;
    double _minHitE;
    double _minHitET;
    
    std::vector <edm::InputTag> _v_inputTag;
    
    std::vector <edm::EDGetTokenT <std::vector <double> > > _v_tok_collection;
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
MapProducer::MapProducer(const edm::ParameterSet& iConfig)
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
    
    _debug = iConfig.getParameter <bool>("debug");
    _useProcessName = iConfig.getParameter <bool>("useProcessName");
    
    _v_inputTag = iConfig.getParameter <std::vector <edm::InputTag> >("collections");
    
    for(int iTag = 0; iTag < (int) _v_inputTag.size(); iTag++)
    {
        edm::InputTag tag = _v_inputTag.at(iTag);
        
        _v_tok_collection.push_back(
            consumes <std::vector <double> >(tag)
        );
    }
    
    
    produces <edm::MapOfVectors <std::string, double> > (_instanceName);
}

MapProducer::~MapProducer() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void MapProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
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
    
    std::map <std::string, std::vector <double> > m_collection;
    
    for(int iTag = 0; iTag < (int) _v_inputTag.size(); iTag++)
    {
        edm::Handle <std::vector <double> > v_value;
        iEvent.getByToken(_v_tok_collection.at(iTag), v_value);
        
        edm::InputTag tag = _v_inputTag.at(iTag);
        
        std::string key = tag.label();
        
        if(!tag.instance().empty())
        {
            key += "_" + tag.instance();
        }
        
        if(!tag.process().empty() && _useProcessName)
        {
            key += "_" + tag.process();
        }
        
        if(_debug)
        {
            printf("In MapProducer::produce(...): Adding collection \"%s\" to map. \n", key.c_str());
            
            printf("Values: ");
            for(int idx = 0; idx < (int) v_value->size(); idx++)
            {
                printf("%0.4e, ", v_value->at(idx));
            }
            printf("\n");
        }
        
        
        m_collection[key] = *v_value;
    }
    
    
    edm::MapOfVectors <std::string, double> outputMap(m_collection);
    
    
    iEvent.put(
        std::make_unique <edm::MapOfVectors <std::string, double> >(outputMap),
        _instanceName
    );
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void MapProducer::beginStream(edm::StreamID) {
  // please remove this method if not needed
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void MapProducer::endStream() {
  // please remove this method if not needed
}

// ------------ method called when starting to processes a run  ------------
/*
void
MapProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
MapProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
MapProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
MapProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void MapProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MapProducer);
