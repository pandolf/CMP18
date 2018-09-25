#ifndef CODAnalyzer_h
#define CODAnalyzer_h



#include <memory>


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"

#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/NjettinessPlugin.hh"



#define nIMGMAX 99999

class TTree;
class PileupSummaryInfo;


class CODAnalyzer : public edm::EDAnalyzer {

   public:

      explicit CODAnalyzer(const edm::ParameterSet&);
      ~CODAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:

      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      void fillImage( float ptRatio, float dEta, float dPhi, int nPix_1D, float pixelSize, float* image );
      void computeQGvars( float sum_weight, float sum_pt, float sum_deta, float sum_dphi, float sum_deta2, float sum_dphi2, float sum_detadphi, float& a_axis1, float& a_axis2, float& a_ptD );
      float computeTau21( const std::vector< fastjet::PseudoJet >& newparts );
      int getPileUp( edm::Handle<std::vector<PileupSummaryInfo>>& pupInfo );


      // ----------member data ---------------------------

      //const JetCorrector *JEC;
      edm::Service<TFileService> fs;
      TTree* tree;

      float rho;
      int event, run, lumi, nVert, nPU;

      float ecalRecHits[61000];

      int nPhoton;
      float pt_phot   [nPhotonMax];
      float eta_phot  [nPhotonMax];
      float phi_phot  [nPhotonMax];
      float eRaw_phot [nPhotonMax];
      float eSC_phot  [nPhotonMax];

      float e1x5_phot[nPhotonMax];
      float e2x5_phot[nPhotonMax];
      float e3x3_phot[nPhotonMax];
      float e5x5_phot[nPhotonMax];
      float maxEnergyXtal_phot[nPhotonMax];
      float sigmaEtaEta_phot[nPhotonMax];
      float sigmaIetaIeta_phot[nPhotonMax];
      float sigmaIphiIphi_phot[nPhotonMax];
      float r1x5_phot[nPhotonMax];
      float r2x5_phot[nPhotonMax];
      float r9_phot[nPhotonMax];

      float full5x5_e1x5_phot[nPhotonMax];
      float full5x5_e2x5_phot[nPhotonMax];
      float full5x5_e3x3_phot[nPhotonMax];
      float full5x5_e5x5_phot[nPhotonMax];
      float full5x5_maxEnergyXtal_phot[nPhotonMax];
      float full5x5_sigmaEtaEta_phot[nPhotonMax];
      float full5x5_sigmaIetaIeta_phot[nPhotonMax];
      float full5x5_sigmaIphiIphi_phot[nPhotonMax];
      float full5x5_r1x5_phot[nPhotonMax];
      float full5x5_r2x5_phot[nPhotonMax];
      float full5x5_r9_phot[nPhotonMax];

      float hOverE_phot  [nPhotonMax];
      float chIso_phot  [nPhotonMax];
      float nhIso_phot  [nPhotonMax];
      float phIso_phot  [nPhotonMax];
      float convSafeEleVeto_phot  [nPhotonMax];
      float pixelSeedVeto_phot  [nPhotonMax];


      int nParton;
      float pt_parton   [nPartonMax];
      float eta_parton  [nPartonMax];
      float phi_parton  [nPartonMax];
      float mass_parton [nPartonMax];
      int   pdgId_parton[nPartonMax];
   

      int nPFCand;
      float pt_pfCand   [nPFCandMax];
      float eta_pfCand  [nPFCandMax];
      float phi_pfCand  [nPFCandMax];
      float mass_pfCand [nPFCandMax];
      int   pdgId_pfCand[nPFCandMax];

      int nCaloTower;
      float pt_caloTower   [nCaloTowerMax];
      float eta_caloTower  [nCaloTowerMax];
      float phi_caloTower  [nCaloTowerMax];
      float mass_caloTower [nCaloTowerMax];
      float emf_caloTower  [nCaloTowerMax];



};


#endif
