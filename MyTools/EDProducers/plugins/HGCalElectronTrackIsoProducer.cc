// -*- C++ -*-
//
// Package:    MyTools/HGCalElectronTrackIsoProducer
// Class:      HGCalElectronTrackIsoProducer
//
/**\class HGCalElectronTrackIsoProducer HGCalElectronTrackIsoProducer.cc MyTools/HGCalElectronTrackIsoProducer/plugins/HGCalElectronTrackIsoProducer.cc

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
# include "DataFormats/Math/interface/LorentzVector.h"
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


//
// class declaration
//

class HGCalElectronTrackIsoProducer : public edm::stream::EDProducer<>
{
    public:
    
    explicit HGCalElectronTrackIsoProducer(const edm::ParameterSet&);
    ~HGCalElectronTrackIsoProducer();
    
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
    
    double _isoConeDR;
    double _vetoConeDR;
    double _vetoPhiStripDeta;
    
    double _minTrackPt;
    double _maxTrackEleDz;
    
    edm::EDGetTokenT <std::vector <reco::GsfElectron> > _tok_electron;
    edm::EDGetTokenT <std::vector <reco::Track> > _tok_track;
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
HGCalElectronTrackIsoProducer::HGCalElectronTrackIsoProducer(const edm::ParameterSet& iConfig)
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
    _tok_track = consumes <std::vector <reco::Track> >(iConfig.getParameter <edm::InputTag>("tracks"));
    
    _isoConeDR = iConfig.getParameter <double>("isoConeDR");
    _vetoConeDR = iConfig.getParameter <double>("vetoConeDR");
    _vetoPhiStripDeta = iConfig.getParameter <double>("vetoPhiStripDeta");
    
    _minTrackPt = iConfig.getParameter <double>("minTrackPt");
    _maxTrackEleDz = iConfig.getParameter <double>("maxTrackEleDz");
    
    _debug = iConfig.getParameter <bool>("debug");
    
    
    produces <std::vector <double> > (_instanceName);
}

HGCalElectronTrackIsoProducer::~HGCalElectronTrackIsoProducer() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void HGCalElectronTrackIsoProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
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
    
    edm::Handle <std::vector <reco::Track> > v_track;
    iEvent.getByToken(_tok_track, v_track);
    
    
    int nEle = v_electron->size();
    int nTrack = v_track->size();
    
    std::vector <double> v_trackIso;
    
    for(int iEle = 0; iEle < nEle; iEle++)
    {
        reco::GsfElectron ele = v_electron->at(iEle);
        
        CLHEP::HepLorentzVector ele_4mom;
        ele_4mom.setT(ele.energy());
        ele_4mom.setX(ele.px());
        ele_4mom.setY(ele.py());
        ele_4mom.setZ(ele.pz());
        
        double trackIso = 0;
        
        const reco::Track *ele_track = &*ele.gsfTrack();
        
        for(int iTrack = 0; iTrack < nTrack; iTrack++)
        {
            reco::Track track = v_track->at(iTrack);
            
            // pT cut
            if(track.pt() < _minTrackPt)
            {
                continue;
            }
            
            // dz cut
            double trkEleDz = std::fabs(track.vz() - ele_track->vz());
            
            if(trkEleDz > _maxTrackEleDz)
            {
                continue;
            }
            
            // phi-strip veto
            double dEta = std::fabs(track.eta() - ele_track->eta());
            
            if(dEta < _vetoPhiStripDeta)
            {
                continue;
            }
            
            // signal and isolation cone veto
            double dR = ROOT::Math::VectorUtil::DeltaR(track.momentum(), ele_track->momentum());
            
            if(dR > _isoConeDR || dR < _vetoConeDR)
            {
                continue;
            }
            
            
            trackIso += track.pt();
        }
        
        
        if(_debug)
        {
            printf("In HGCalElectronTrackIsoProducer --> Ele %d/%d: isoTrack %0.4f \n", iEle+1, nEle, trackIso);
        }
        
        v_trackIso.push_back(trackIso);
    }
    
    
    iEvent.put(
        std::make_unique <std::vector <double> >(v_trackIso),
        _instanceName
    );
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void HGCalElectronTrackIsoProducer::beginStream(edm::StreamID) {
  // please remove this method if not needed
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void HGCalElectronTrackIsoProducer::endStream() {
  // please remove this method if not needed
}

// ------------ method called when starting to processes a run  ------------
/*
void
HGCalElectronTrackIsoProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
HGCalElectronTrackIsoProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
HGCalElectronTrackIsoProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
HGCalElectronTrackIsoProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HGCalElectronTrackIsoProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCalElectronTrackIsoProducer);
