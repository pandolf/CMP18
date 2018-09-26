#include <iostream>
#include <stdlib.h>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TLorentzVector.h"





int main( int argc, char* argv[] ) {


  TFile* file = TFile::Open( "codTree.root" );
  TTree* tree = (TTree*)file->Get("codana/codtree");


  int entry = 0;
  if( argc > 1 ) {
    entry = atoi(argv[1]);
    std::cout << "Will look at entry n." << entry << std::endl;
  }


  int nPFCand;
  tree->SetBranchAddress( "nPFCand", &nPFCand );
  float pt_pfCand[4000];
  tree->SetBranchAddress( "pt_pfCand", pt_pfCand );
  float eta_pfCand[4000];
  tree->SetBranchAddress( "eta_pfCand", eta_pfCand );
  float phi_pfCand[4000];
  tree->SetBranchAddress( "phi_pfCand", phi_pfCand );
  float energy_pfCand[4000];
  tree->SetBranchAddress( "energy_pfCand", energy_pfCand );

  int nCaloTower;
  tree->SetBranchAddress( "nCaloTower", &nCaloTower );
  float pt_caloTower[4000];
  tree->SetBranchAddress( "pt_caloTower", pt_caloTower );
  float eta_caloTower[4000];
  tree->SetBranchAddress( "eta_caloTower", eta_caloTower );
  float phi_caloTower[4000];
  tree->SetBranchAddress( "phi_caloTower", phi_caloTower );
  float energy_caloTower[4000];
  tree->SetBranchAddress( "energy_caloTower", energy_caloTower );


  tree->GetEntry(entry);


  TH2D* h2_pfEtaPhi = new TH2D( "pfEtaPhi", "",  250, -2.5, 2.5, 360, -3.1416, 3.1416 );
  h2_pfEtaPhi->SetXTitle( "Eta" );
  h2_pfEtaPhi->SetYTitle( "Phi" );
  

  for( int i=0; i<nPFCand; ++i ) {

    int ieta = h2_pfEtaPhi->GetXaxis()->FindBin( eta_pfCand[i] );
    int iphi = h2_pfEtaPhi->GetYaxis()->FindBin( phi_pfCand[i] );

    h2_pfEtaPhi->SetBinContent( ieta, iphi, energy_pfCand[i] );

  }




  TH2D* h2_ctEtaPhi = new TH2D( "ctEtaPhi", "",  60, -3., 3., 63, -3.1416, 3.1416 );
  h2_ctEtaPhi->SetXTitle( "Eta" );
  h2_ctEtaPhi->SetYTitle( "Phi" );
  

  for( int i=0; i<nCaloTower; ++i ) {

    int ieta = h2_ctEtaPhi->GetXaxis()->FindBin( eta_caloTower[i] );
    int iphi = h2_ctEtaPhi->GetYaxis()->FindBin( phi_caloTower[i] );

    h2_ctEtaPhi->SetBinContent( ieta, iphi, energy_caloTower[i] );


  }


  TFile* outfile = TFile::Open( "etaPhiJet.root", "recreate" );
  outfile->cd();

  h2_pfEtaPhi->Write();
  h2_ctEtaPhi->Write();

  outfile->Close();

  return 0;

}
