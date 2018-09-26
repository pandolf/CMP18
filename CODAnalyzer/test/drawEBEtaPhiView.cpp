#include <iostream>
#include <stdlib.h>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2D.h"





int main( int argc, char* argv[] ) {


  TFile* file = TFile::Open( "codTree.root" );
  TTree* tree = (TTree*)file->Get("codana/codtree");

  int entry = 0;
  if( argc > 1 ) {
    entry = atoi(argv[1]);
    std::cout << "Will look at entry n." << entry << std::endl;
  }


  int nEcalRecHit;
  tree->SetBranchAddress( "nEcalRecHit", &nEcalRecHit );
  int iEta_ecalRecHit[4000];
  tree->SetBranchAddress( "iEta_ecalRecHit", iEta_ecalRecHit );
  int iPhi_ecalRecHit[4000];
  tree->SetBranchAddress( "iPhi_ecalRecHit", iPhi_ecalRecHit );
  float energy_ecalRecHit[4000];
  tree->SetBranchAddress( "energy_ecalRecHit", energy_ecalRecHit );

  tree->GetEntry(entry);

  TH2D* h2_ebEtaPhi = new TH2D( "ebEtaPhi", "",  170, -85., 85.000001, 360, 0., 360.00001 );
  

  for( int i=0; i<nEcalRecHit; ++i ) 
    h2_ebEtaPhi->SetBinContent( iEta_ecalRecHit[i], iPhi_ecalRecHit[i], energy_ecalRecHit[i] );

  TFile* outfile = TFile::Open( "prova.root", "recreate" );
  outfile->cd();

  h2_ebEtaPhi->Write();

  outfile->Close();

  return 0;

}
