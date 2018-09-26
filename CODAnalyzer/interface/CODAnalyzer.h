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



#define NPHOTONMAX 6
#define NPARTONMAX 10
#define NPFCANDMAX 2000
#define NCALOTOWERMAX 4000
#define NECALRECHITMAX 2000


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

      //void fillImage( float ptRatio, float dEta, float dPhi, int nPix_1D, float pixelSize, float* image );
      //void computeQGvars( float sum_weight, float sum_pt, float sum_deta, float sum_dphi, float sum_deta2, float sum_dphi2, float sum_detadphi, float& a_axis1, float& a_axis2, float& a_ptD );
      //float computeTau21( const std::vector< fastjet::PseudoJet >& newparts );
      int getPileUp( edm::Handle<std::vector<PileupSummaryInfo>>& pupInfo );


      // ----------member data ---------------------------

      edm::Service<TFileService> fs;
      TTree* tree;

      float rho;
      int event, run, lumi, nVert, nPU;

      int nEcalRecHit;
      int iEta_ecalRecHit[NECALRECHITMAX];
      int iPhi_ecalRecHit[NECALRECHITMAX];
      float energy_ecalRecHit[NECALRECHITMAX];


      int nPhoton;
      float pt_phot  [NPHOTONMAX];
      float eta_phot [NPHOTONMAX];
      float phi_phot [NPHOTONMAX];
      float eRaw_phot[NPHOTONMAX];
      float eSC_phot [NPHOTONMAX];

      float e1x5_phot[NPHOTONMAX];
      float e2x5_phot[NPHOTONMAX];
      float e3x3_phot[NPHOTONMAX];
      float e5x5_phot[NPHOTONMAX];
      float maxEnergyXtal_phot[NPHOTONMAX];
      float sigmaEtaEta_phot[NPHOTONMAX];
      float sigmaIetaIeta_phot[NPHOTONMAX];
      float r1x5_phot[NPHOTONMAX];
      float r2x5_phot[NPHOTONMAX];
      float r9_phot[NPHOTONMAX];

      float ecalRecHitSumEtConeDR04_phot[NPHOTONMAX];
      float hcalTowerSumEtConeDR04_phot[NPHOTONMAX];
      float hcalTowerSumEtBcConeDR04_phot[NPHOTONMAX];
      float trkSumPtSolidConeDR04_phot[NPHOTONMAX];
      float trkSumPtHollowConeDR04_phot[NPHOTONMAX];
      int   nTrkSolidConeDR04_phot[NPHOTONMAX];
      int   nTrkHollowConeDR04_phot[NPHOTONMAX];

      float ecalRecHitSumEtConeDR03_phot[NPHOTONMAX];
      float hcalTowerSumEtConeDR03_phot[NPHOTONMAX];
      float hcalTowerSumEtBcConeDR03_phot[NPHOTONMAX];
      float trkSumPtSolidConeDR03_phot[NPHOTONMAX];
      float trkSumPtHollowConeDR03_phot[NPHOTONMAX];
      int   nTrkSolidConeDR03_phot[NPHOTONMAX];
      int   nTrkHollowConeDR03_phot[NPHOTONMAX];

      //float full5x5_e1x5_phot[NPHOTONMAX];
      //float full5x5_e2x5_phot[NPHOTONMAX];
      //float full5x5_e3x3_phot[NPHOTONMAX];
      //float full5x5_e5x5_phot[NPHOTONMAX];
      //float full5x5_maxEnergyXtal_phot[NPHOTONMAX];
      //float full5x5_sigmaEtaEta_phot[NPHOTONMAX];
      //float full5x5_sigmaIetaIeta_phot[NPHOTONMAX];
      //float full5x5_sigmaIphiIphi_phot[NPHOTONMAX];
      //float full5x5_r1x5_phot[NPHOTONMAX];
      //float full5x5_r2x5_phot[NPHOTONMAX];
      //float full5x5_r9_phot[NPHOTONMAX];

      float hOverE_phot  [NPHOTONMAX];
      float hTowOverE_phot  [NPHOTONMAX];
      float chIsoPF_phot  [NPHOTONMAX];
      float nhIsoPF_phot  [NPHOTONMAX];
      float phIsoPF_phot  [NPHOTONMAX];

      //float convSafeEleVeto_phot  [NPHOTONMAX];
      bool hasPixelSeed_phot  [NPHOTONMAX];


      int nParton;
      float pt_parton   [NPARTONMAX];
      float eta_parton  [NPARTONMAX];
      float phi_parton  [NPARTONMAX];
      float mass_parton [NPARTONMAX];
      int   pdgId_parton[NPARTONMAX];
   

      int nPFCand;
      float pt_pfCand   [NPFCANDMAX];
      float eta_pfCand  [NPFCANDMAX];
      float phi_pfCand  [NPFCANDMAX];
      float mass_pfCand [NPFCANDMAX];
      int   pdgId_pfCand[NPFCANDMAX];

      int nCaloTower;
      float pt_caloTower   [NCALOTOWERMAX];
      float eta_caloTower  [NCALOTOWERMAX];
      float phi_caloTower  [NCALOTOWERMAX];
      float mass_caloTower [NCALOTOWERMAX];

      float energy_caloTower      [NCALOTOWERMAX];
      float emEnergy_caloTower    [NCALOTOWERMAX];
      float hadEnergy_caloTower   [NCALOTOWERMAX];
      float outerEnergy_caloTower [NCALOTOWERMAX];



};


#endif
