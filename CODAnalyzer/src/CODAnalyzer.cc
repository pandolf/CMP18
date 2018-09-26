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

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"

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

   Handle< EcalRecHitCollection > ebReducedRecHitCollection_h;
   iEvent.getByLabel( "reducedEcalRecHitsEB", ebReducedRecHitCollection_h);
   const EcalRecHitCollection* ebRecHits = ebReducedRecHitCollection_h.product();

   Handle< PFCandidateCollection > pfCandidateCollection_h;
   iEvent.getByLabel( "particleFlow", pfCandidateCollection_h);
   const PFCandidateCollection* pfCands = pfCandidateCollection_h.product();

   Handle< CaloTowerCollection > caloTowerCollection_h;
   iEvent.getByLabel( "towerMaker", caloTowerCollection_h);
   const CaloTowerCollection* caloTowers = caloTowerCollection_h.product();
 
   Handle< GenParticleCollection > GenParticleCollection_h;
   iEvent.getByLabel( "genParticles", GenParticleCollection_h);
   const GenParticleCollection* genParticles = GenParticleCollection_h.product();
 


   event =   (int) iEvent.id().event();
   run   =   (int) iEvent.id().run();
   lumi  =   (int) iEvent.id().luminosityBlock();


   nVert = vertexCollection->size();
   nPU   = getPileUp(pupInfo);

   
   
   
   
   
   // BEGIN GENPARTICLES -------------------------------------

   int iGenPart_ind = 0;

   for( GenParticleCollection::const_iterator iGenPart = genParticles->begin(); iGenPart!=genParticles->end(); ++iGenPart ) {

     if( iGenPart->status()==3 ) {

       pt_genPart   [iGenPart_ind] = (float)(iGenPart->pt());
       eta_genPart  [iGenPart_ind] = (float)(iGenPart->eta());
       phi_genPart  [iGenPart_ind] = (float)(iGenPart->phi());
       mass_genPart [iGenPart_ind] = (float)(iGenPart->mass());
       pdgId_genPart[iGenPart_ind] = (float)(iGenPart->pdgId());

       iGenPart_ind++;

     } // if status 3

   } // for genparticles

   nGenPart = iGenPart_ind;





   // BEGIN PHOTONS -------------------------------------

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

      ecalRecHitSumEtConeDR04_phot[iPhot] = photons->at(iPhot).ecalRecHitSumEtConeDR04();
      hcalTowerSumEtConeDR04_phot[iPhot] = photons->at(iPhot).hcalTowerSumEtConeDR04();
      hcalTowerSumEtBcConeDR04_phot[iPhot] = photons->at(iPhot).hcalTowerSumEtBcConeDR04();
      trkSumPtSolidConeDR04_phot[iPhot] = photons->at(iPhot).trkSumPtSolidConeDR04();
      trkSumPtHollowConeDR04_phot[iPhot] = photons->at(iPhot).trkSumPtHollowConeDR04();
      nTrkSolidConeDR04_phot[iPhot] = photons->at(iPhot).nTrkSolidConeDR04();
      nTrkHollowConeDR04_phot[iPhot] = photons->at(iPhot).nTrkHollowConeDR04();

      ecalRecHitSumEtConeDR03_phot[iPhot] = photons->at(iPhot).ecalRecHitSumEtConeDR03();
      hcalTowerSumEtConeDR03_phot[iPhot] = photons->at(iPhot).hcalTowerSumEtConeDR03();
      hcalTowerSumEtBcConeDR03_phot[iPhot] = photons->at(iPhot).hcalTowerSumEtBcConeDR03();
      trkSumPtSolidConeDR03_phot[iPhot] = photons->at(iPhot).trkSumPtSolidConeDR03();
      trkSumPtHollowConeDR03_phot[iPhot] = photons->at(iPhot).trkSumPtHollowConeDR03();
      nTrkSolidConeDR03_phot[iPhot] = photons->at(iPhot).nTrkSolidConeDR03();
      nTrkHollowConeDR03_phot[iPhot] = photons->at(iPhot).nTrkHollowConeDR03();

      chIsoPF_phot[iPhot] = photons->at(iPhot).chargedHadronIso();
      nhIsoPF_phot[iPhot] = photons->at(iPhot).neutralHadronIso();
      phIsoPF_phot[iPhot] = photons->at(iPhot).photonIso();
      
      //convSafeEleVeto_phot = photons[i]->
      hasPixelSeed_phot[iPhot] = photons->at(iPhot).hasPixelSeed();
     
   } // for photons
  



   // BEGIN EB RECHITS ----------------------------------------

   // initialize:
   for( int i=0; i<NECALRECHITMAX; ++i ) {
     energy_ecalRecHit[i] = 0.;
   }

   nEcalRecHit = ebRecHits->size();
 
   int irechit_ind=0;
 
   for( EcalRecHitCollection::const_iterator iRecHit = ebRecHits->begin(); iRecHit != ebRecHits->end(); ++iRecHit ) {
 
     EBDetId edDetId(iRecHit->id());
     
     iEta_ecalRecHit[irechit_ind] = edDetId.ieta();
     iPhi_ecalRecHit[irechit_ind] = edDetId.iphi();
     energy_ecalRecHit[irechit_ind] = iRecHit->energy();
 
     irechit_ind++;
 
   }  // for recHits
 
 
 
   // BEGIN PF CANDIDATES -------------------------------------

   for( int i=0; i<NPFCANDMAX; ++i ) {
     pt_pfCand[i] = 0.;
     eta_pfCand[i] = 0.;
     phi_pfCand[i] = 0.;
     mass_pfCand[i] = 0.;
     pdgId_pfCand[i] = 0;
   }



 
   int iCand_ind = 0;
 
   for( PFCandidateCollection::const_iterator iCand=pfCands->begin(); iCand!=pfCands->end(); ++iCand ) {
 
     if( fabs(iCand->eta())<2.5 ) { //&& iCand->pdgId()!=0 ) {
 
       float mass = -1.;

       if( abs(iCand->pdgId())==22 )
         mass = 0.;
       else if( abs(iCand->pdgId())==130 ) 
         mass = 0.497648;
       else if( abs(iCand->pdgId())==13 || abs(iCand->pdgId())==11 || abs(iCand->pdgId())==211 )
         mass = (float)(iCand->mass());
       else {
         std::cout << "Unknown PDG ID: " << iCand->pdgId() << std::endl;
         continue;
       }


       if( mass >= 0. ) {
      
         pt_pfCand   [iCand_ind] = (float)(iCand->pt());
         eta_pfCand  [iCand_ind] = (float)(iCand->eta());
         phi_pfCand  [iCand_ind] = (float)(iCand->phi());
         pdgId_pfCand[iCand_ind] = (int)  (iCand->pdgId());
         mass_pfCand [iCand_ind] = mass;
 
         iCand_ind++;
 
       } // if mass ok
 
     } // if eta 2.5
     
   } // for pf cands
 
   nPFCand = iCand_ind;



   // BEGIN CALOTOWERS -------------------------------

   int iCaloTower_ind = 0;

   for( CaloTowerCollection::const_iterator iCaloTower=caloTowers->begin(); iCaloTower!=caloTowers->end(); iCaloTower++ ) {

     pt_caloTower   [iCaloTower_ind] = (float)(iCaloTower->pt());
     eta_caloTower  [iCaloTower_ind] = (float)(iCaloTower->eta());
     phi_caloTower  [iCaloTower_ind] = (float)(iCaloTower->phi());
     mass_caloTower [iCaloTower_ind] = (float)(iCaloTower->mass());

     energy_caloTower      [iCaloTower_ind] = (float)(iCaloTower->energy());
     emEnergy_caloTower    [iCaloTower_ind] = (float)(iCaloTower->emEnergy());
     hadEnergy_caloTower   [iCaloTower_ind] = (float)(iCaloTower->hadEnergy());
     outerEnergy_caloTower [iCaloTower_ind] = (float)(iCaloTower->outerEnergy());

     iCaloTower_ind++;

   } // for calotowers

   nCaloTower = iCaloTower_ind;



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


  tree = fs->make<TTree>("codtree","");
  tree->Branch("event" , &event, "event/I");
  tree->Branch("run"   , &run  , "run/I");
  tree->Branch("lumi"  , &lumi , "lumi /I");
  tree->Branch("rho"   , &rho  , "rho/F");
  tree->Branch("nVert" , &nVert, "nVert/I");
  tree->Branch("nPU"   , &nPU  , "nPU/I");

  tree->Branch("nGenPart" , &nGenPart , "nGenPart/I");
  tree->Branch("pt_genPart" , pt_genPart , "pt_genPart[nGenPart]/F");
  tree->Branch("eta_genPart" , eta_genPart , "eta_genPart[nGenPart]/F");
  tree->Branch("phi_genPart" , phi_genPart , "phi_genPart[nGenPart]/F");
  tree->Branch("mass_genPart" , mass_genPart , "mass_genPart[nGenPart]/F");
  tree->Branch("pdgId_genPart" , pdgId_genPart , "pdgId_genPart[nGenPart]/I");

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

  tree->Branch("ecalRecHitSumEtConeDR04_phot" , ecalRecHitSumEtConeDR04_phot , "ecalRecHitSumEtConeDR04_phot[nPhoton]/F");
  tree->Branch("hcalTowerSumEtConeDR04_phot" , hcalTowerSumEtConeDR04_phot , "hcalTowerSumEtConeDR04_phot[nPhoton]/F");
  tree->Branch("hcalTowerSumEtBcConeDR04_phot" , hcalTowerSumEtBcConeDR04_phot , "hcalTowerSumEtBcConeDR04_phot[nPhoton]/F");
  tree->Branch("trkSumPtSolidConeDR04_phot" , trkSumPtSolidConeDR04_phot , "trkSumPtSolidConeDR04_phot[nPhoton]/F");
  tree->Branch("trkSumPtHollowConeDR04_phot" , trkSumPtHollowConeDR04_phot , "trkSumPtHollowConeDR04_phot[nPhoton]/F");
  tree->Branch("nTrkSolidConeDR04_phot" , nTrkSolidConeDR04_phot , "nTrkSolidConeDR04_phot[nPhoton]/I");
  tree->Branch("nTrkHollowConeDR04_phot" , nTrkHollowConeDR04_phot , "nTrkHollowConeDR04_phot[nPhoton]/I");

  tree->Branch("ecalRecHitSumEtConeDR03_phot" , ecalRecHitSumEtConeDR03_phot , "ecalRecHitSumEtConeDR03_phot[nPhoton]/F");
  tree->Branch("hcalTowerSumEtConeDR03_phot" , hcalTowerSumEtConeDR03_phot , "hcalTowerSumEtConeDR03_phot[nPhoton]/F");
  tree->Branch("hcalTowerSumEtBcConeDR03_phot" , hcalTowerSumEtBcConeDR03_phot , "hcalTowerSumEtBcConeDR03_phot[nPhoton]/F");
  tree->Branch("trkSumPtSolidConeDR03_phot" , trkSumPtSolidConeDR03_phot , "trkSumPtSolidConeDR03_phot[nPhoton]/F");
  tree->Branch("trkSumPtHollowConeDR03_phot" , trkSumPtHollowConeDR03_phot , "trkSumPtHollowConeDR03_phot[nPhoton]/F");
  tree->Branch("nTrkSolidConeDR03_phot" , nTrkSolidConeDR03_phot , "nTrkSolidConeDR03_phot[nPhoton]/I");
  tree->Branch("nTrkHollowConeDR03_phot" , nTrkHollowConeDR03_phot , "nTrkHollowConeDR03_phot[nPhoton]/I");

  tree->Branch("chIsoPF_phot" , chIsoPF_phot , "chIsoPF_phot[nPhoton]/F");
  tree->Branch("nhIsoPF_phot" , nhIsoPF_phot , "nhIsoPF_phot[nPhoton]/F");
  tree->Branch("phIsoPF_phot" , phIsoPF_phot , "phIsoPF_phot[nPhoton]/F");

  tree->Branch("hasPixelSeed_phot" , hasPixelSeed_phot , "hasPixelSeed_phot[nPhoton]/O");

  tree->Branch("nEcalRecHit" , &nEcalRecHit , "nEcalRecHit/I");
  tree->Branch("iPhi_ecalRecHit" , iPhi_ecalRecHit , "iPhi_ecalRecHit[nEcalRecHit]/I");
  tree->Branch("iEta_ecalRecHit" , iEta_ecalRecHit , "iEta_ecalRecHit[nEcalRecHit]/I");
  tree->Branch("energy_ecalRecHit" , energy_ecalRecHit , "energy_ecalRecHit[nEcalRecHit]/F");

  tree->Branch("nPFCand" , &nPFCand , "nPFCand/I");
  tree->Branch("pt_pfCand" , pt_pfCand , "pt_pfCand[nPFCand]/F");
  tree->Branch("eta_pfCand" , eta_pfCand , "eta_pfCand[nPFCand]/F");
  tree->Branch("phi_pfCand" , phi_pfCand , "phi_pfCand[nPFCand]/F");
  tree->Branch("mass_pfCand" , mass_pfCand , "mass_pfCand[nPFCand]/F");
  tree->Branch("pdgId_pfCand" , pdgId_pfCand , "pdgId_pfCand[nPFCand]/I");

  tree->Branch("nCaloTower" , &nCaloTower , "nCaloTower/I");
  tree->Branch("pt_caloTower" , pt_caloTower , "pt_caloTower[nCaloTower]/F");
  tree->Branch("eta_caloTower" , eta_caloTower , "eta_caloTower[nCaloTower]/F");
  tree->Branch("phi_caloTower" , phi_caloTower , "phi_caloTower[nCaloTower]/F");
  tree->Branch("mass_caloTower" , mass_caloTower , "mass_caloTower[nCaloTower]/F");
  tree->Branch("energy_caloTower" , energy_caloTower , "energy_caloTower[nCaloTower]/F");
  tree->Branch("emEnergy_caloTower" , emEnergy_caloTower , "emEnergy_caloTower[nCaloTower]/F");
  tree->Branch("hadEnergy_caloTower" , hadEnergy_caloTower , "hadEnergy_caloTower[nCaloTower]/F");
  tree->Branch("outerEnergy_caloTower" , outerEnergy_caloTower , "outerEnergy_caloTower[nCaloTower]/F");

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
