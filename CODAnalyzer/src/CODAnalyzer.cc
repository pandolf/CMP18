// -*- C++ -*-
//
// Package:    CODAnalyzer
// Class:      CODAnalyzer
// 
/**\class CODAnalyzer CODAnalyzer.cc CMP18/CODAnalyzer/src/CODAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Francesco Pandolfi,32 3-B02,+41227676027,
//         Created:  Wed Oct 18 10:51:59 CEST 2017
// $Id$
//
//



#include "../interface/CODAnalyzer.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

#include <cmath>


//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
CODAnalyzer::CODAnalyzer(const edm::ParameterSet& iConfig)

{

   nPhotonMax = 4;
   nPartonMax = 10;
   nPFCandMax = 500;
   nCaloTowerMax = 2000;
   

}


CODAnalyzer::~CODAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
CODAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   using namespace edm;
   using namespace reco;


   Handle<double> rhoHandle;
   iEvent.getByLabel("fixedGridRhoFastjetAll", rhoHandle);
   rho = (float) *rhoHandle;

   Handle< std::vector<Vertex> > vertexCollection;
   iEvent.getByLabel("offlinePrimaryVerticesWithBS", vertexCollection);

   Handle<std::vector<PileupSummaryInfo>> pupInfo;
   iEvent.getByLabel("addPileupInfo", pupInfo);


   event =   (int) iEvent.id().event();
   run   =   (int) iEvent.id().run();
   lumi  =   (int) iEvent.id().luminosityBlock();


   nVert = vertexCollection->size();
   nPU   = getPileUp(pupInfo);






  

}




void CODAnalyzer::computeQGvars( float sum_weight, float sum_pt, float sum_deta, float sum_dphi, float sum_deta2, float sum_dphi2, float sum_detadphi, float& a_axis1, float& a_axis2, float& a_ptD ) {

  float a = 0., b = 0., c = 0.;
  float ave_deta = 0., ave_dphi = 0., ave_deta2 = 0., ave_dphi2 = 0.;
  if(sum_weight > 0){
    ave_deta  = sum_deta/sum_weight;
    ave_dphi  = sum_dphi/sum_weight;
    ave_deta2 = sum_deta2/sum_weight;
    ave_dphi2 = sum_dphi2/sum_weight;
    a         = ave_deta2 - ave_deta*ave_deta;                          
    b         = ave_dphi2 - ave_dphi*ave_dphi;                          
    c         = -(sum_detadphi/sum_weight - ave_deta*ave_dphi);                
  }
  float delta = sqrt(fabs((a-b)*(a-b)+4*c*c));
  a_axis2 = (a+b-delta > 0 ?  sqrt(0.5*(a+b-delta)) : 0);
  a_axis1 = (a+b+delta > 0 ?  sqrt(0.5*(a+b+delta)) : 0);
  a_ptD   = (sum_weight > 0 ? sqrt(sum_weight)/sum_pt : 0);

  a_axis2 = -log(a_axis2);
  a_axis1 = -log(a_axis1);

}


int CODAnalyzer::getPileUp( edm::Handle<std::vector<PileupSummaryInfo>>& pupInfo ) {

  if(!pupInfo.isValid()) return -1;
  auto PVI = pupInfo->begin();
  while(PVI->getBunchCrossing() != 0 && PVI != pupInfo->end()) ++PVI;
  if(PVI != pupInfo->end()) return PVI->getPU_NumInteractions();
  else return -1;

} 


float CODAnalyzer::computeTau21( const std::vector< fastjet::PseudoJet >& newparts ) {
 
  fastjet::ClusterSequenceArea* thisClustering = new fastjet::ClusterSequenceArea(newparts, *jetDef, *fjAreaDefinition);
  std::vector<fastjet::PseudoJet> out_jets = sorted_by_pt(thisClustering->inclusive_jets(0.01));        
      
  return (*nSub2KT)(out_jets[0])/(*nSub1KT)(out_jets[0]);

}


//reco::GenParticleCollection::const_iterator CODAnalyzer::getMatchedGenParticle(const TLorentzVector& jet, edm::Handle<reco::GenParticleCollection>& genParticles ) {
//
//  float deltaRmin = 999.;
//  auto matchedGenParticle = genParticles->end();
//
//  for ( auto genParticle = genParticles->begin(); genParticle != genParticles->end(); ++genParticle ) {
//
//    if( !genParticle->isHardProcess() ) continue; // This status flag is exactly the pythia8 status-23 we need (i.e. the same as genParticles->status() == 23), probably also ok to use for other generators
//    if( abs(genParticle->pdgId()) > 5 && abs(genParticle->pdgId() != 21) ) continue; // only udscb quarks and gluons
//
//    float thisDeltaR = reco::deltaR(*genParticle, *jet);
//    if(thisDeltaR < deltaRmin && thisDeltaR < deltaRcut){
//      deltaRmin = thisDeltaR;
//      matchedGenParticle = genParticle;
//    }
//  }
//  return matchedGenParticle;
//}




// ------------ method called once each job just before starting event loop  ------------
void 
CODAnalyzer::beginJob()
{


  //file = TFile::Open("CODAnalyzer.root", "recreate" );
  //file->cd();

  //tree = new TTree( "CODAnalyzer", "" );

  tree = fs->make<TTree>("CODAnalyzer","CODAnalyzer");
  tree->Branch("event" , &event, "event/I");
  tree->Branch("run"   , &run  , "run/I");
  tree->Branch("lumi"  , &lumi , "lumi /I");
  tree->Branch("rho"   , &rho  , "rho/F");
  tree->Branch("nVert" , &nVert, "nVert/I");
  tree->Branch("nPU"   , &nPU  , "nPU/I");
  tree->Branch("pt"    , &pt   , "pt/F");
  tree->Branch("eta"   , &eta  , "eta/F");
  tree->Branch("phi"   , &phi  , "phi/F");
  tree->Branch("mass"  , &mass , "mass/F");
  tree->Branch("axis1" , &axis1, "axis1/F");
  tree->Branch("axis2" , &axis2, "axis2/F");
  tree->Branch("ptD"   , &ptD  , "ptD/F");
  tree->Branch("tau21" , &tau21, "tau21/F");
  tree->Branch("ptGen"    , &ptGen   , "ptGen/F");
  tree->Branch("etaGen"   , &etaGen  , "etaGen/F");
  tree->Branch("phiGen"   , &phiGen  , "phiGen/F");
  tree->Branch("massGen"  , &massGen , "massGen/F");
  tree->Branch("axis1Gen" , &axis1Gen, "axis1Gen/F");
  tree->Branch("axis2Gen" , &axis2Gen, "axis2Gen/F");
  tree->Branch("ptDGen"   , &ptDGen  , "ptDGen/F");
  tree->Branch("tau21Gen" , &tau21Gen, "tau21Gen/F");
  tree->Branch("btag"  , &btag , "btag/F");
  tree->Branch("partonId"  , &partonId , "partonId/I");
  tree->Branch("jetIdLevel"  , &jetIdLevel , "jetIdLevel/I");
  tree->Branch("pixelSize", &pixelSize, "pixelSize/F");
  tree->Branch("drMax", &drMax, "drMax/F");
  tree->Branch("nPix", &nPix, "nPix/I");
  tree->Branch("jetImageReco",jetImageReco, "jetImageReco[nPix]/F");
  tree->Branch("jetImageGen" ,jetImageGen , "jetImageGen[nPix]/F");
  tree->Branch("frac_pt" , &frac_pt , "frac_pt/F");



}

// ------------ method called once each job just after ending the event loop  ------------
void 
CODAnalyzer::endJob() 
{

  //file->cd();
  //tree->Write();
  //file->Close();

}

// ------------ method called when starting to processes a run  ------------
void 
CODAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
CODAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
CODAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
CODAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CODAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(CODAnalyzer);