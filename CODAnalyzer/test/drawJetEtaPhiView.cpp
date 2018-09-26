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
  float mass_pfCand[4000];
  tree->SetBranchAddress( "mass_pfCand", mass_pfCand );

  tree->GetEntry(entry);


  TH2D* h2_pfEtaPhi = new TH2D( "pfEtaPhi", "",  250, -2.5, 2.5, 360, 0., 360.00001 );
  

  for( int i=0; i<nPFCand; ++i ) {

    TLorentzVector p;
    p.SetPtEtaPhiM( pt_pfCand[i], eta_pfCand[i], phi_pfCand[i], mass_pfCand[i] );

    int ieta = h2_pfEtaPhi->GetXaxis()->FindBin( p.Eta() );
    int iphi = h2_pfEtaPhi->GetYaxis()->FindBin( p.Phi() );

    h2_pfEtaPhi->SetBinContent( ieta, iphi, p.Energy() );

  }


  TFile* outfile = TFile::Open( "etaPhiJet.root", "recreate" );
  outfile->cd();

  h2_pfEtaPhi->Write();

  outfile->Close();

  return 0;

}
