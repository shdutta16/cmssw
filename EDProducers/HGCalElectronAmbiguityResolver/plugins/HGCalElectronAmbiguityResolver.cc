// -*- C++ -*-
//
// Package:    EDProducers/HGCalElectronAmbiguityResolver
// Class:      HGCalElectronAmbiguityResolver
// 
/**\class HGCalElectronAmbiguityResolver HGCalElectronAmbiguityResolver.cc EDProducers/HGCalElectronAmbiguityResolver/plugins/HGCalElectronAmbiguityResolver.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Soham Bhattacharya
//         Created:  Fri, 31 May 2019 19:55:49 GMT
//
//


// system include files
#include <memory>

// user include files
//# include "CommonTools/UtilAlgos/interface/TFileService.h"
//# include "DataFormats/CaloTowers/interface/CaloTowerDefs.h"
//# include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
//# include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
//# include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
//# include "DataFormats/HGCRecHit/interface/HGCRecHit.h"
//# include "DataFormats/HepMCCandidate/interface/GenParticle.h"
//# include "DataFormats/JetReco/interface/PFJet.h"
//# include "DataFormats/Math/interface/LorentzVector.h"
//# include "DataFormats/ParticleFlowReco/interface/HGCalMultiCluster.h"
//# include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
//# include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
//# include "DataFormats/TrackReco/interface/Track.h"
//# include "DataFormats/TrackReco/interface/TrackFwd.h"
//# include "DataFormats/VertexReco/interface/Vertex.h"
//# include "FWCore/Framework/interface/stream/EDProducer.h"
//# include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
//# include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
//# include "FWCore/ServiceRegistry/interface/Service.h"
//# include "FWCore/Utilities/interface/InputTag.h"
//# include "FWCore/Utilities/interface/StreamID.h"
//# include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
//# include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
//# include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
//# include "SimDataFormats/CaloHit/interface/PCaloHit.h"
//# include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
//# include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

# include "FWCore/Framework/interface/Frameworkfwd.h"
# include "FWCore/Framework/interface/stream/EDProducer.h"
# include "FWCore/Framework/interface/Event.h"
# include "FWCore/Framework/interface/MakerMacros.h"
# include "FWCore/ParameterSet/interface/ParameterSet.h"
# include "FWCore/Utilities/interface/StreamID.h"
# include "DataFormats/ParticleFlowReco/interface/HGCalMultiCluster.h"

# include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
# include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
# include "DataFormats/HGCRecHit/interface/HGCRecHit.h"
# include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
# include "DataFormats/Common/interface/SortedCollection.h"

//# include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
# include "RecoEgamma/EgammaElectronAlgos/src/EgAmbiguityTools.cc"

# include <CLHEP/Matrix/Matrix.h>
# include <CLHEP/Vector/ThreeVector.h>
# include <CLHEP/Vector/ThreeVector.h>

# include <Compression.h>
# include <TH1F.h>
# include <TH2F.h>
# include <TMatrixD.h>
# include <TTree.h> 
# include <TVectorD.h> 

//
// class declaration
//

class HGCalElectronAmbiguityResolver : public edm::stream::EDProducer<>
{
    public:
    
    explicit HGCalElectronAmbiguityResolver(const edm::ParameterSet&);
    ~HGCalElectronAmbiguityResolver();
    
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    
    private:
    
    virtual void beginStream(edm::StreamID) override;
    virtual void produce(edm::Event&, const edm::EventSetup&) override;
    virtual void endStream() override;
    
    //float HGCalElectronAmbiguityResolver::sharedEnergy(
    //    CaloCluster const& clu1,
    //    CaloCluster const& clu2,
    //    edm::SortedCollection <HGCRecHit,edm::StrictWeakOrdering <HGCRecHit> > const& endcapRecHits
    //);
    //
    //float HGCalElectronAmbiguityResolver::sharedEnergy(
    //    SuperClusterRef const& sc1,
    //    SuperClusterRef const& sc2,
    //    edm::SortedCollection <HGCRecHit,edm::StrictWeakOrdering <HGCRecHit> > const& endcapRecHits
    //);
    
    //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    
    // ----------member data ---------------------------
    
    hgcal::RecHitTools recHitTools;
    
    
    // From config file
    bool debug;
    
    std::string instanceName;
    
    edm::EDGetTokenT <std::vector <reco::GsfElectron> > tok_electron;
    
    edm::EDGetTokenT <edm::SortedCollection <HGCRecHit,edm::StrictWeakOrdering <HGCRecHit> > > tok_HGCEERecHit;
    edm::EDGetTokenT <edm::SortedCollection <HGCRecHit,edm::StrictWeakOrdering <HGCRecHit> > > tok_HGCHEFRecHit;
    edm::EDGetTokenT <edm::SortedCollection <HGCRecHit,edm::StrictWeakOrdering <HGCRecHit> > > tok_HGCHEBRecHit;
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
HGCalElectronAmbiguityResolver::HGCalElectronAmbiguityResolver(const edm::ParameterSet& iConfig)
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
    
    
    // From config file
    debug = iConfig.getParameter <bool>("debug");
    
    instanceName = iConfig.getParameter <std::string>("instanceName");
    
    tok_electron = consumes <std::vector <reco::GsfElectron> >(iConfig.getUntrackedParameter <edm::InputTag>("label_electron"));
    
    tok_HGCEERecHit = consumes <edm::SortedCollection <HGCRecHit,edm::StrictWeakOrdering <HGCRecHit> > >(iConfig.getUntrackedParameter <edm::InputTag>("label_HGCEERecHit"));
    tok_HGCHEFRecHit = consumes <edm::SortedCollection <HGCRecHit,edm::StrictWeakOrdering <HGCRecHit> > >(iConfig.getUntrackedParameter <edm::InputTag>("label_HGCHEFRecHit"));
    tok_HGCHEBRecHit = consumes <edm::SortedCollection <HGCRecHit,edm::StrictWeakOrdering <HGCRecHit> > >(iConfig.getUntrackedParameter <edm::InputTag>("label_HGCHEBRecHit"));
    
    // Produces
    instanceName = instanceName;
    
    produces <std::vector <reco::GsfElectron>> (instanceName);
}


HGCalElectronAmbiguityResolver::~HGCalElectronAmbiguityResolver()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//


//float HGCalElectronAmbiguityResolver::sharedEnergy(
//    CaloCluster const& clu1,
//    CaloCluster const& clu2,
//    edm::SortedCollection <HGCRecHit,edm::StrictWeakOrdering <HGCRecHit> > const& endcapRecHits
//)
//{
//    double fractionShared = 0;
//    
//    for(auto const& h1 : clu1.hitsAndFractions()) 
//    {
//        for(auto const& h2 : clu2.hitsAndFractions())
//        {
//            if (h1.first != h2.first)
//            {
//                continue;
//            }
//            
//        // here we have common Xtal id
//    
//        EcalRecHitCollection::const_iterator itt;
//        
//        if ((itt = endcapRecHits.find(h1.first)) != endcapRecHits.end())
//        {
//            fractionShared += itt->energy();
//        }
//    }
//    
//    //std::cout << "[sharedEnergy] shared energy /min(energy1,energy2) " << fractionShared << std::endl;
//    
//    return fractionShared;
//}
//
//
//float HGCalElectronAmbiguityResolver::sharedEnergy(
//    SuperClusterRef const& sc1,
//    SuperClusterRef const& sc2,
//    edm::SortedCollection <HGCRecHit,edm::StrictWeakOrdering <HGCRecHit> > const& endcapRecHits
//)
//{
//  double energyShared = 0;
//  
//    for(CaloCluster_iterator icl1 = sc1->clustersBegin(); icl1 != sc1->clustersEnd(); icl1++) 
//    {
//        for(CaloCluster_iterator icl2 = sc2->clustersBegin(); icl2 != sc2->clustersEnd(); icl2++)
//        {
//            energyShared += sharedEnergy(**icl1, **icl2, barrelRecHits, endcapRecHits);
//        }
//    }
//  
//  return energyShared;
//}


// ------------ method called to produce the data  ------------
void HGCalElectronAmbiguityResolver::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    
    /* This is an event example
    //Read 'ExampleData' from the Event
    Handle<ExampleData> pIn;
    iEvent.getByLabel("example",pIn);
    
    //Use the ExampleData to create an ExampleData2 which 
    // is put into the Event
    iEvent.put(std::make_unique<ExampleData2>(*pIn));
    */
    
    /* this is an EventSetup example
    //Read SetupData from the SetupRecord in the EventSetup
    ESHandle<SetupData> pSetup;
    iSetup.get<SetupRecord>().get(pSetup);
    */
    
    recHitTools.getEventSetup(iSetup);
    
    
    // RecHit dictionary
    edm::Handle <edm::SortedCollection <HGCRecHit,edm::StrictWeakOrdering <HGCRecHit> > > v_HGCEERecHit;
    iEvent.getByToken(tok_HGCEERecHit, v_HGCEERecHit);
    
    edm::Handle <edm::SortedCollection <HGCRecHit,edm::StrictWeakOrdering <HGCRecHit> > > v_HGCHEFRecHit;
    iEvent.getByToken(tok_HGCHEFRecHit, v_HGCHEFRecHit);
    
    edm::Handle <edm::SortedCollection <HGCRecHit,edm::StrictWeakOrdering <HGCRecHit> > > v_HGCHEBRecHit;
    iEvent.getByToken(tok_HGCHEBRecHit, v_HGCHEBRecHit);
    
    std::map <DetId, const HGCRecHit*> m_recHit;
    
    //
    int nHGCEERecHit = v_HGCEERecHit->size();
    
    for(int iRecHit = 0; iRecHit < nHGCEERecHit; iRecHit++)
    {
        const HGCRecHit *recHit = &(*v_HGCEERecHit)[iRecHit];
        
        m_recHit[recHit->id()] = recHit;
    }
    
    
    //
    int nHGCHEFRecHit = v_HGCHEFRecHit->size();
    
    for(int iRecHit = 0; iRecHit < nHGCHEFRecHit; iRecHit++)
    {
        const HGCRecHit *recHit = &(*v_HGCHEFRecHit)[iRecHit];
        
        m_recHit[recHit->id()] = recHit;
    }
    
    
    //
    int nHGCHEBRecHit = v_HGCHEBRecHit->size();
    
    for(int iRecHit = 0; iRecHit < nHGCHEBRecHit; iRecHit++)
    {
        const HGCRecHit *recHit = &(*v_HGCHEBRecHit)[iRecHit];
        
        m_recHit[recHit->id()] = recHit;
    }
    
    
    // Electrons
    edm::Handle <std::vector <reco::GsfElectron> > v_electron;
    iEvent.getByToken(tok_electron, v_electron);
    
    
    //std::sort(v_electron.begin(), v_electron.end(), egamma::isBetterElectron);
    //std::sort(v_electron->begin(), v_electron->end(), egamma::isInnermostElectron);
    
    // From https://cmssdt.cern.ch/lxr/source/RecoEgamma/EgammaElectronProducers/plugins/GsfElectronProducer.cc#0428
    
    // init
    //for(auto& electron : v_electron)
    //{
    //    electron.clearAmbiguousGsfTracks();
    //    electron.setAmbiguous(false);
    //}
    //
    //for(auto e1 = v_electron.begin(); e1 != v_electron.end(); ++e1) {
    //    if (e1->ambiguous())
    //        continue;
    //    if (ignoreNotPreselected && !isPreselected(*e1))
    //        continue;
    //        
    //    SuperClusterRef scRef1 = e1->superCluster();
    //    CaloClusterPtr eleClu1 = e1->electronCluster();
    //    LogDebug("GsfElectronAlgo") << "Blessing electron with E/P " << e1->eSuperClusterOverP() << ", cluster "
    //                                << scRef1.get() << " & track " << e1->gsfTrack().get();
    //
    //    for(auto e2 = e1 + 1; e2 != v_electron.end(); ++e2) {
    //        if (e2->ambiguous())
    //        continue;
    //        if (ignoreNotPreselected && !isPreselected(*e2))
    //        continue;
    //        
    //        SuperClusterRef scRef2 = e2->superCluster();
    //        CaloClusterPtr eleClu2 = e2->electronCluster();
    //
    //        // search if same cluster
    //        //bool sameCluster = false;
    //        //if (strategyCfg_.ambClustersOverlapStrategy == 0) {
    //        //sameCluster = (scRef1 == scRef2);
    //        //} else if (strategyCfg_.ambClustersOverlapStrategy == 1) {
    //        //float eMin = 1.;
    //        //float threshold = eMin * cosh(EleRelPoint(scRef1->position(), beamspot.position()).eta());
    //        //using egamma::sharedEnergy;
    //        ////sameCluster = ((sharedEnergy(*eleClu1, *eleClu2, barrelRecHits, endcapRecHits) >= threshold) ||
    //        ////                (sharedEnergy(*scRef1->seed(), *eleClu2, barrelRecHits, endcapRecHits) >= threshold) ||
    //        ////                (sharedEnergy(*eleClu1, *scRef2->seed(), barrelRecHits, endcapRecHits) >= threshold) ||
    //        ////                (sharedEnergy(*scRef1->seed(), *scRef2->seed(), barrelRecHits, endcapRecHits) >= threshold));
    //        //sameCluster = (sharedEnergy(*scRef1->seed(), *scRef2->seed(), barrelRecHits, endcapRecHits) >= threshold);
    //        //} else {
    //        //throw cms::Exception("GsfElectronAlgo|UnknownAmbiguityClustersOverlapStrategy")
    //        //    << "value of strategyCfg_.ambClustersOverlapStrategy is : " << strategyCfg_.ambClustersOverlapStrategy;
    //        //}
    //        
    //        bool sameCluster = (scRef1 == scRef2);
    //        
    //        // main instructions
    //        if (sameCluster) {
    //        LogDebug("GsfElectronAlgo") << "Discarding electron with E/P " << e2->eSuperClusterOverP() << ", cluster "
    //                                    << scRef2.get() << " and track " << e2->gsfTrack().get();
    //        e1->addAmbiguousGsfTrack(e2->gsfTrack());
    //        e2->setAmbiguous(true);
    //        } else if (e1->gsfTrack() == e2->gsfTrack()) {
    //        edm::LogWarning("GsfElectronAlgo") << "Forgetting electron with E/P " << e2->eSuperClusterOverP()
    //                                            << ", cluster " << scRef2.get() << " and track " << e2->gsfTrack().get();
    //        e2->setAmbiguous(true);
    //        }
    //    }
    //}
    //
    //
    //v_electron.erase(
    //    std::remove_if(
    //        v_electron.begin(),
    //        v_electron.end(),
    //        std::mem_fn(&reco::GsfElectron::ambiguous)
    //    ),
    //    v_electron.end()
    //);
    //
    //
    ////printf("multiClus_totE        %0.2f \n", multiClus_totE);
    ////printf("multiClus_recHit_totE %0.2f \n", multiClus_recHit_totE);
    ////printf("ratio %0.2f \n", multiClus_totE/multiClus_recHit_totE);
    //
    //// Put the collection in the event
    //iEvent.put(std::move(v_electron), instanceName);
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
HGCalElectronAmbiguityResolver::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
HGCalElectronAmbiguityResolver::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
HGCalElectronAmbiguityResolver::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
HGCalElectronAmbiguityResolver::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
HGCalElectronAmbiguityResolver::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
HGCalElectronAmbiguityResolver::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HGCalElectronAmbiguityResolver::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    //desc.setUnknown();
    desc.setAllowAnything();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCalElectronAmbiguityResolver);
