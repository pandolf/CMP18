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

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

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

   Handle<PhotonCollection> photons_h;
   iEvent.getByLabel( "photons", photons_h);
   const PhotonCollection* photons = photons_h.product();

   //edm::EDGetTokenT<EcalRecHitCollection> ebReducedRecHitCollection_;
   //ebReducedRecHitCollection_        = mayConsume<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>

   //edm::EDGetTokenT<edm::View<reco::Candidate> > pfCandidatesToken_;
   //pfCandidatesToken_        = mayConsume< edm::View<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("pfCandidates")); 
   //edm::Handle< edm::View<reco::Candidate> > pfCandidatesHandle;
   //iEvent.getByToken(pfCandidatesToken_, pfCandidatesHandle);





   event =   (int) iEvent.id().event();
   run   =   (int) iEvent.id().run();
   lumi  =   (int) iEvent.id().luminosityBlock();


   nVert = vertexCollection->size();
   nPU   = getPileUp(pupInfo);



   nPhoton = photons->size();

   for( int iPhot=0; iPhot<nPhoton; ++iPhot ) {

      pt_phot  [iPhot] = photons->at(iPhot).pt();
      eta_phot [iPhot] = photons->at(iPhot).eta();
      phi_phot [iPhot] = photons->at(iPhot).phi();
      eRaw_phot[iPhot] = photons->at(iPhot).superCluster()->rawEnergy();
      eSC_phot [iPhot] = photons->at(iPhot).superCluster()->energy();

      e1x5_phot[iPhot] = photons->at(iPhot).e1x5();
      e2x5_phot[iPhot] = photons->at(iPhot).e2x5();
      e3x3_phot[iPhot] = photons->at(iPhot).e3x3();
      e5x5_phot[iPhot] = photons->at(iPhot).e5x5();
      maxEnergyXtal_phot[iPhot] = photons->at(iPhot).maxEnergyXtal();
      sigmaEtaEta_phot[iPhot] = photons->at(iPhot).sigmaEtaEta();
      sigmaIetaIeta_phot[iPhot] = photons->at(iPhot).sigmaIetaIeta();
      //sigmaIphiIphi_phot[iPhot] = photons->at(iPhot).sigmaIphiIphi();
      r1x5_phot[iPhot] = photons->at(iPhot).r1x5();
      r2x5_phot[iPhot] = photons->at(iPhot).r2x5();
      r9_phot[iPhot] = photons->at(iPhot).r9();

      //full5x5_e1x5_phot[iPhot] = photons->at(iPhot).full5x5_e1x5();
      //full5x5_e2x5_phot[iPhot] = photons->at(iPhot).full5x5_e2x5();
      //full5x5_e3x3_phot[iPhot] = photons->at(iPhot).full5x5_e3x3();
      //full5x5_e5x5_phot[iPhot] = photons->at(iPhot).full5x5_e5x5();
      //full5x5_maxEnergyXtal_phot[iPhot] = photons->at(iPhot).full5x5_maxEnergyXtal();
      //full5x5_sigmaEtaEta_phot[iPhot] = photons->at(iPhot).full5x5_sigmaEtaEta();
      //full5x5_sigmaIetaIeta_phot[iPhot] = photons->at(iPhot).full5x5_sigmaIetaIeta();
      //full5x5_sigmaIphiIphi_phot[iPhot] = photons->at(iPhot).full5x5_sigmaIphiIphi();
      //full5x5_r1x5_phot[iPhot] = photons->at(iPhot).full5x5_r1x5();
      //full5x5_r2x5_phot[iPhot] = photons->at(iPhot).full5x5_r2x5();
      //full5x5_r9_phot[iPhot] = photons->at(iPhot).full5x5_r9();

      hOverE_phot[iPhot]  = photons->at(iPhot).hadronicOverEm();
      hTowOverE_phot[iPhot]  = photons->at(iPhot).hadTowOverEm();
      chIso_phot[iPhot] = photons->at(iPhot).chargedHadronIso();
      nhIso_phot[iPhot] = photons->at(iPhot).neutralHadronIso();
      phIso_phot[iPhot] = photons->at(iPhot).photonIso();
      //convSafeEleVeto_phot = photons[i]->
      //pixelSeedVeto_phot = photons[i]->
     
   } // for photons
  

  tree->Fill();

}




//void CODAnalyzer::computeQGvars( float sum_weight, float sum_pt, float sum_deta, float sum_dphi, float sum_deta2, float sum_dphi2, float sum_detadphi, float& a_axis1, float& a_axis2, float& a_ptD ) {
//
//  float a = 0., b = 0., c = 0.;
//  float ave_deta = 0., ave_dphi = 0., ave_deta2 = 0., ave_dphi2 = 0.;
//  if(sum_weight > 0){
//    ave_deta  = sum_deta/sum_weight;
//    ave_dphi  = sum_dphi/sum_weight;
//    ave_deta2 = sum_deta2/sum_weight;
//    ave_dphi2 = sum_dphi2/sum_weight;
//    a         = ave_deta2 - ave_deta*ave_deta;                          
//    b         = ave_dphi2 - ave_dphi*ave_dphi;                          
//    c         = -(sum_detadphi/sum_weight - ave_deta*ave_dphi);                
//  }
//  float delta = sqrt(fabs((a-b)*(a-b)+4*c*c));
//  a_axis2 = (a+b-delta > 0 ?  sqrt(0.5*(a+b-delta)) : 0);
//  a_axis1 = (a+b+delta > 0 ?  sqrt(0.5*(a+b+delta)) : 0);
//  a_ptD   = (sum_weight > 0 ? sqrt(sum_weight)/sum_pt : 0);
//
//  a_axis2 = -log(a_axis2);
//  a_axis1 = -log(a_axis1);
//
//}


int CODAnalyzer::getPileUp( edm::Handle<std::vector<PileupSummaryInfo>>& pupInfo ) {

  if(!pupInfo.isValid()) return -1;
  auto PVI = pupInfo->begin();
  while(PVI->getBunchCrossing() != 0 && PVI != pupInfo->end()) ++PVI;
  if(PVI != pupInfo->end()) return PVI->getPU_NumInteractions();
  else return -1;

} 


//float CODAnalyzer::computeTau21( const std::vector< fastjet::PseudoJet >& newparts ) {
// 
//  fastjet::ClusterSequenceArea* thisClustering = new fastjet::ClusterSequenceArea(newparts, *jetDef, *fjAreaDefinition);
//  std::vector<fastjet::PseudoJet> out_jets = sorted_by_pt(thisClustering->inclusive_jets(0.01));        
//      
//  return (*nSub2KT)(out_jets[0])/(*nSub1KT)(out_jets[0]);
//
//}


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

  tree->Branch("nPhoton"  , &nPhoton , "nPhoton/I" );
  tree->Branch("pt_phot"  , pt_phot  , "pt_phot[nPhoton]/F");
  tree->Branch("eta_phot" , eta_phot , "eta_phot[nPhoton]/F");
  tree->Branch("phi_phot" , phi_phot , "phi_phot[nPhoton]/F");
  tree->Branch("eRaw_phot", eRaw_phot, "eRaw_phot[nPhoton]/F");
  tree->Branch("eSC_phot" , eSC_phot , "eSC_phot[nPhoton]/F");

  tree->Branch("e1x5_phot" , e1x5_phot , "e1x5_phot[nPhoton]/F");
  tree->Branch("e2x5_phot" , e2x5_phot , "e2x5_phot[nPhoton]/F");
  tree->Branch("e3x3_phot" , e3x3_phot , "e3x3_phot[nPhoton]/F");
  tree->Branch("e5x5_phot" , e5x5_phot , "e5x5_phot[nPhoton]/F");
  tree->Branch("maxEnergyXtal_phot" , maxEnergyXtal_phot , "maxEnergyXtal_phot[nPhoton]/F");
  tree->Branch("sigmaEtaEta_phot" , sigmaEtaEta_phot , "sigmaEtaEta_phot[nPhoton]/F");
  tree->Branch("sigmaIetaIeta_phot" , sigmaIetaIeta_phot , "sigmaIetaIeta_phot[nPhoton]/F");
  tree->Branch("r1x5_phot" , r1x5_phot , "r1x5_phot[nPhoton]/F");
  tree->Branch("r2x5_phot" , r2x5_phot , "r2x5_phot[nPhoton]/F");
  tree->Branch("r9_phot" , r9_phot , "r9_phot[nPhoton]/F");

  //tree->Branch(full5x5_e1x5_phot[nPhoton];
  //tree->Branch(full5x5_e2x5_phot[nPhoton];
  //tree->Branch(full5x5_e3x3_phot[nPhoton];
  //tree->Branch(full5x5_e5x5_phot[nPhoton];
  //tree->Branch(full5x5_maxEnergyXtal_phot[nPhoton];
  //tree->Branch(full5x5_sigmaEtaEta_phot[nPhoton];
  //tree->Branch(full5x5_sigmaIetaIeta_phot[nPhoton];
  //tree->Branch(full5x5_sigmaIphiIphi_phot[nPhoton];
  //tree->Branch(full5x5_r1x5_phot[nPhoton];
  //tree->Branch(full5x5_r2x5_phot[nPhoton];
  //tree->Branch(full5x5_r9_phot[nPhoton];

  tree->Branch("hOverE_phot" , hOverE_phot , "hOverE_phot[nPhoton]/F");
  tree->Branch("hTowOverE_phot" , hTowOverE_phot , "hTowOverE_phot[nPhoton]/F");
  tree->Branch("chIso_phot" , chIso_phot , "chIso_phot[nPhoton]/F");
  tree->Branch("nhIso_phot" , nhIso_phot , "nhIso_phot[nPhoton]/F");
  tree->Branch("phIso_phot" , phIso_phot , "phIso_phot[nPhoton]/F");

  //tree->Branch("eSC_phot" , eSC_phot , "eSC_phot[nPhoton]/F");
  //tree->Branch("eSC_phot" , eSC_phot , "eSC_phot[nPhoton]/F");
  //tree->Branch(hOverE_phot  [nPhoton];
  //tree->Branch(hTowOverE_phot  [nPhoton];
  //tree->Branch(chIso_phot  [nPhoton];
  //tree->Branch(nhIso_phot  [nPhoton];
  //tree->Branch(phIso_phot  [nPhoton];
  //tree->Branch(convSafeEleVeto_phot  [nPhoton];
  //tree->Branch(pixelSeedVeto_phot  [nPhoton];



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
