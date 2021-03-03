#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <array>
#include <stdlib.h>
#include <functional>
#include <time.h>
#include <math.h>
#include <numeric>

#include "TROOT.h"

#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TLorentzVector.h"

#include "ExRootAnalysis/ExRootProgressBar.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"

using namespace std;


//---------------------------------------------------------------------------


static int NMAX = 9000;

struct PFCand
{
  float pt = 0;
  float eta = 0;
  float phi = 0;
  float e = 0;
  float puppi = 1;
  float pdgid = 0;
  float hardfrac = 1;  
  float cluster_idx = -1;
  float cluster_hardch_pt = 0;
  float cluster_puch_pt = 0;
  float cluster_r = 0;
  float vtxid = -1;
  float npv = 0;
};




//---------------------------------------------------------------------------

int main(int argc, char *argv[])
{

  srand(time(NULL));

  if(argc < 3) {
    cout << " Usage: " << "PapuDelphes" << " input_file"
         << " output_file" << endl;
    cout << " input_file - input file in ROOT format," << endl;
    cout << " output_file - output file in ROOT format" << endl;
    return 1;
  }

  // figure out how to read the file here 
  //

  TFile* ifile = TFile::Open(argv[1], "READ");
  TTree* itree = (TTree*)ifile->Get("Delphes;1");

  auto* fout = TFile::Open(argv[2], "RECREATE");
  auto* tout = new TTree("events", "events");

  unsigned int nevt = itree->GetEntries();
  TBranch* pfbranch = (TBranch*)itree->GetBranch("ParticleFlowCandidate");
  TBranch* genjetbranch = (TBranch*)itree->GetBranch("GenJet");
  std::cout << "NEVT: " << nevt << std::endl;
  vector<PFCand> input_particles;

  vector<PFCand> output_particles;
  output_particles.reserve(NMAX);

  vector<float> vpt, veta, vphi, ve, vpuppi, vpdgid, vhardfrac, vcluster_idx, vvtxid, vcluster_r, vcluster_hardch_pt, vcluster_puch_pt, vnpv;
  vpt.reserve(NMAX); veta.reserve(NMAX); vphi.reserve(NMAX); 
  ve.reserve(NMAX); vpuppi.reserve(NMAX); vpdgid.reserve(NMAX); 
  vhardfrac.reserve(NMAX); vcluster_idx.reserve(NMAX); vvtxid.reserve(NMAX);
  vcluster_r.reserve(NMAX); vcluster_hardch_pt.reserve(NMAX); vcluster_puch_pt.reserve(NMAX);
  vnpv.reserve(NMAX);

  float genmet=-99., genmetphi=-99.;
  float genjet1pt=-99., genjet1eta=-99., genjet1phi=-99., genjet1e=-99.;
  float genjet2pt=-99., genjet2eta=-99., genjet2phi=-99., genjet2e=-99.;
  
  tout->Branch("pt", &vpt);
  tout->Branch("eta", &veta);
  tout->Branch("phi", &vphi);
  tout->Branch("e", &ve);
  tout->Branch("puppi", &vpuppi);
  tout->Branch("pdgid", &vpdgid);
  tout->Branch("hardfrac", &vhardfrac);
  tout->Branch("cluster_idx", &vcluster_idx);
  tout->Branch("cluster_r", &vcluster_r);
  tout->Branch("cluster_hardch_pt", &vcluster_hardch_pt);
  tout->Branch("cluster_puch_pt", &vcluster_puch_pt);
  tout->Branch("vtxid", &vvtxid);
  tout->Branch("npv", &vnpv);
  TBranch* b_genmet = tout->Branch("genmet",&genmet, "genmet/F");
  TBranch* b_genmetphi = tout->Branch("genmetphi",&genmetphi, "genmetphi/F");
  TBranch* b_genjet1pt = tout->Branch("genjet1pt",&genjet1pt, "genjet1pt/F");
  TBranch* b_genjet1eta = tout->Branch("genjet1eta",&genjet1eta, "genjet1eta/F");
  TBranch* b_genjet1phi = tout->Branch("genjet1phi",&genjet1phi, "genjet1phi/F");
  TBranch* b_genjet1e = tout->Branch("genjet1e",&genjet1e, "genjet1e/F");
  TBranch* b_genjet2pt = tout->Branch("genjet2pt",&genjet2pt, "genjet2pt/F");
  TBranch* b_genjet2eta = tout->Branch("genjet2eta",&genjet2eta, "genjet2eta/F");
  TBranch* b_genjet2phi = tout->Branch("genjet2phi",&genjet2phi, "genjet2phi/F");
  TBranch* b_genjet2e = tout->Branch("genjet2e",&genjet2e, "genjet2e/F");


  ExRootProgressBar progressBar(nevt);
  
  auto comp_pt = [](auto &a, auto &b) { return a.sum_pt() > b.sum_pt(); };
  auto comp_p4 = [](auto &a, auto &b) { return a.pt > b.pt; };

  for (unsigned int k=0; k<nevt; k++){
    itree->GetEntry(k);
    
    float npv = itree->GetLeaf("Vertex_size")->GetValue(0);
    genmet = itree->GetLeaf("GenMissingET.MET")->GetValue(0);
    genmetphi = itree->GetLeaf("GenMissingET.Phi")->GetValue(0);

    unsigned int ngenjets = genjetbranch->GetEntries();
    ngenjets = itree->GetLeaf("GenJet_size")->GetValue(0);
    for (unsigned int j=0; j<ngenjets; j++){
      if (j>1)
	break;
      TLorentzVector tmpjet;
      tmpjet.SetPtEtaPhiM(itree->GetLeaf("GenJet.PT")->GetValue(j),itree->GetLeaf("GenJet.Eta")->GetValue(j),itree->GetLeaf("GenJet.Phi")->GetValue(j),itree->GetLeaf("GenJet.Mass")->GetValue(j));
      if (j==0){
	genjet1pt = tmpjet.Pt();
	genjet1eta = tmpjet.Eta();
	genjet1phi = tmpjet.Phi();
	genjet1e = tmpjet.E();
      }
      if (j==1){
	genjet2pt = tmpjet.Pt();
	genjet2eta = tmpjet.Eta();
	genjet2phi = tmpjet.Phi();
	genjet2e = tmpjet.E();
      }      
    }

    
    tout->Fill();

    progressBar.Update(k, k);

  }

  fout->Write();
  fout->Close();

}
