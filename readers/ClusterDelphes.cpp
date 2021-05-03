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

#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TLorentzVector.h"

#include "ExRootAnalysis/ExRootProgressBar.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"

#include "external/fastjet/contribs/EnergyCorrelator/EnergyCorrelator.hh"
//#include "fastjet/contrib/EnergyCorrelator.hh"
//#include "modules/EnergyCorrelator.hh"

#include "fastjet/plugins/CDFCones/fastjet/CDFJetCluPlugin.hh"
#include "fastjet/plugins/CDFCones/fastjet/CDFMidPointPlugin.hh"
#include "fastjet/plugins/SISCone/fastjet/SISConePlugin.hh"

#include "fastjet/contribs/Nsubjettiness/ExtraRecombiners.hh"
#include "fastjet/contribs/Nsubjettiness/Njettiness.hh"
#include "fastjet/contribs/Nsubjettiness/NjettinessPlugin.hh"
#include "fastjet/contribs/Nsubjettiness/XConePlugin.hh"
#include "fastjet/contribs/Nsubjettiness/Nsubjettiness.hh"

#include "fastjet/contribs/ValenciaPlugin/ValenciaPlugin.hh"

#include "fastjet/contribs/RecursiveTools/SoftDrop.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"

#include "modules/EnergyCorrelations.h"

using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;

//---------------------------------------------------------------------------

struct PFCand
{
  float pt = 0;
  float eta = 0;
  float phi = 0;
  float e = 0;
  float puppi = 1;
  float vtxid = -1;
};


struct puppiak8
{
  float pt = 0;
  float eta = 0;
  float phi = 0;
  float msd = 0;
  float n2 = 1;
};

struct chsak4
{
  float pt = 0;
  float eta = 0;
  float phi = 0;
  float mass = 0;
};

struct xcone4
{
  float pt = 0;
  float eta = 0;
  float phi = 0;
  float mass = 0;
};


template <typename T>
void 
fill(vector<float> &vattr, vector<puppiak8> &fatjets, T fn_attr)
{
  vattr.clear();
  for (auto& p : fatjets)
    vattr.push_back(fn_attr(p));
}

template <typename T>
void 
fillak4(vector<float> &vattr, vector<chsak4> &jets, T fn_attr)
{
  vattr.clear();
  for (auto& p : jets)
    vattr.push_back(fn_attr(p));
}

template <typename T>
void 
fillxcone4(vector<float> &vattr, vector<xcone4> &jets, T fn_attr)
{
  vattr.clear();
  for (auto& p : jets)
    vattr.push_back(fn_attr(p));
}


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
  std::cout << "NEVT: " << nevt << std::endl;
  vector<PFCand> input_particles;
  vector<puppiak8> output_fatjets;
  vector<chsak4> output_jets;
  vector<xcone4> output_xjets;
  output_fatjets.reserve(10);
  output_jets.reserve(30);
  output_xjets.reserve(3);

  std::vector<float> puppijetpt, puppijeteta, puppijetphi, puppijetm, puppijetn2;
  puppijetpt.reserve(10); puppijeteta.reserve(10); puppijetphi.reserve(10); puppijetm.reserve(10); puppijetn2.reserve(10);
  std::vector<float> chsjetpt, chsjeteta, chsjetphi, chsjetm;
  chsjetpt.reserve(30); chsjeteta.reserve(30); chsjetphi.reserve(30); chsjetm.reserve(30); 
  std::vector<float> xjetpt, xjeteta, xjetphi, xjetm;
  xjetpt.reserve(30); xjeteta.reserve(30); xjetphi.reserve(30); xjetm.reserve(30); 

  float ht=0;

  TBranch* b_ht = tout->Branch("HT",&ht, "HT/F");

  tout->Branch("FatJet_pt", &puppijetpt);
  tout->Branch("FatJet_eta", &puppijeteta);
  tout->Branch("FatJet_phi", &puppijetphi);
  tout->Branch("FatJet_msoftdrop", &puppijetm);
  tout->Branch("FatJet_n2b1", &puppijetn2);
  tout->Branch("Jet_pt", &chsjetpt);
  tout->Branch("Jet_eta", &chsjeteta);
  tout->Branch("Jet_phi", &chsjetphi);
  tout->Branch("Jet_mass", &chsjetm);
  tout->Branch("XJet_pt", &xjetpt);
  tout->Branch("XJet_eta", &xjeteta);
  tout->Branch("XJet_phi", &xjetphi);
  tout->Branch("XJet_mass", &xjetm);

  ExRootProgressBar progressBar(nevt);
  
  auto comp_p4 = [](auto &a, auto &b) { return a.pt > b.pt; };

  fastjet::JetDefinition *jetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm,0.8);
  fastjet::JetDefinition *jetDef_ak4 = new fastjet::JetDefinition(fastjet::antikt_algorithm,0.4);

  // Soft drop -- only needed for large radius jets
  double radius = 0.8;
  double sdZcut = 0.1;
  double sdBeta = 0.;
  
  fastjet::contrib::SoftDrop softDrop = fastjet::contrib::SoftDrop(sdBeta,sdZcut,radius);

  TH1D* hDTotalMCWeight = new TH1D("hDTotalMCWeight","hDTotalMCWeight",1,0,2);

  for (unsigned int k=0; k<nevt; k++){
    itree->GetEntry(k);

    hDTotalMCWeight->SetBinContent(1,hDTotalMCWeight->GetBinContent(1)+1);

    if (k%100==0)
      std::cout << k << " / " << nevt << std::endl;

    input_particles.clear();
    unsigned int npfs = pfbranch->GetEntries();
    npfs = itree->GetLeaf("ParticleFlowCandidate_size")->GetValue(0);

    for (unsigned int j=0; j<npfs; j++){
      PFCand tmppf;
      tmppf.pt = itree->GetLeaf("ParticleFlowCandidate.PT")->GetValue(j);
      tmppf.eta = itree->GetLeaf("ParticleFlowCandidate.Eta")->GetValue(j);
      tmppf.phi = itree->GetLeaf("ParticleFlowCandidate.Phi")->GetValue(j);
      tmppf.e = itree->GetLeaf("ParticleFlowCandidate.E")->GetValue(j);
      tmppf.puppi = itree->GetLeaf("ParticleFlowCandidate.PuppiW")->GetValue(j);
      if (itree->GetLeaf("ParticleFlowCandidate.Charge")->GetValue(j)!=0){
        if (itree->GetLeaf("ParticleFlowCandidate.hardfrac")->GetValue(j)==1)
          tmppf.vtxid = 0;
        else
	  continue; // Charged Hadron Subtraction
      }
      else
        tmppf.vtxid = -1;
      input_particles.push_back(tmppf);
    }

    // sorting input particles by pT
    sort(input_particles.begin(), input_particles.end(), comp_p4);


    vector<fastjet::PseudoJet> finalStates_ak8;
    for(auto &p : input_particles){
      if (p.vtxid==1)
	continue; // Charged Hadron Subtraction
      TLorentzVector tmp;
      tmp.SetPtEtaPhiE(p.pt,p.eta,p.phi,p.e);

      // Puppi
      if (p.puppi>0.){
	fastjet::PseudoJet curjet_puppi(p.puppi*tmp.Px(), p.puppi*tmp.Py(), p.puppi*tmp.Pz(), p.puppi*tmp.E());
	finalStates_ak8.emplace_back(curjet_puppi);
      }
    }

    fastjet::ClusterSequence seq_puppi(sorted_by_pt(finalStates_ak8), *jetDef);

    float minpt = 150;

    vector<fastjet::PseudoJet> allJets_puppi(sorted_by_pt(seq_puppi.inclusive_jets(minpt)));

    output_fatjets.clear();

    for (auto& puppiJet : allJets_puppi) {
      if (puppiJet.perp() < minpt)
	break;

      fastjet::PseudoJet sdJet = (softDrop)(puppiJet);
      float beta = 1.0;
      EnergyCorrelatorN2 N2(beta);

      puppiak8 tmpak8;
      tmpak8.pt = puppiJet.perp();
      tmpak8.eta = puppiJet.eta();
      tmpak8.phi = puppiJet.phi();
      tmpak8.msd = sdJet.m();
      tmpak8.n2 = N2(sdJet);

      output_fatjets.push_back(tmpak8);

    }

    vector<fastjet::PseudoJet> finalStates_ak4;
    for(auto &p : input_particles){
      if (p.vtxid==1)
	continue; // Charged Hadron Subtraction
      TLorentzVector tmp;
      tmp.SetPtEtaPhiE(p.pt,p.eta,p.phi,p.e);

      fastjet::PseudoJet curjet_chs(tmp.Px(),tmp.Py(),tmp.Pz(),tmp.E());
      finalStates_ak4.emplace_back(curjet_chs);
    }


    fastjet::ClusterSequence seq_chs(sorted_by_pt(finalStates_ak4), *jetDef_ak4);

    minpt = 30;

    vector<fastjet::PseudoJet> allJets_chs(sorted_by_pt(seq_chs.inclusive_jets(minpt)));

    output_jets.clear();
    for (auto& chsJet : allJets_chs) {
      
      ht += chsJet.perp();
      
      if (chsJet.perp() < minpt)
	break;

      chsak4 tmpak4;
      tmpak4.pt = chsJet.perp();
      tmpak4.eta = chsJet.eta();
      tmpak4.phi = chsJet.phi();
      tmpak4.mass = chsJet.m();

      output_jets.push_back(tmpak4);
    }


    XConePlugin xcone_plugin3(3, 0.4, 2.0);
    JetDefinition xcone_jetDef3(&xcone_plugin3);
    ClusterSequence xcone_seq3(finalStates_ak4, xcone_jetDef3);
    vector<PseudoJet> xcone_jets3 = xcone_seq3.inclusive_jets();

    output_xjets.clear();

    for(auto & j : sorted_by_pt(xcone_jets3)){
      xcone4 tmpx;
      tmpx.pt = j.pt();
      tmpx.eta = j.eta();
      tmpx.phi = j.phi();
      tmpx.mass = j.m();
      //cout << j.pt() << endl;
      output_xjets.push_back(tmpx);
    }

    //output_fatjets.resize(2);

    fill(puppijetpt, output_fatjets, [](puppiak8& p) { return p.pt; }); 
    fill(puppijeteta, output_fatjets, [](puppiak8& p) { return p.eta; }); 
    fill(puppijetphi, output_fatjets, [](puppiak8& p) { return p.phi; }); 
    fill(puppijetm, output_fatjets, [](puppiak8& p) { return p.msd; }); 
    fill(puppijetn2, output_fatjets, [](puppiak8& p) { return p.n2; }); 
    fillak4(chsjetpt, output_jets, [](chsak4& p) { return p.pt; }); 
    fillak4(chsjeteta, output_jets, [](chsak4& p) { return p.eta; }); 
    fillak4(chsjetphi, output_jets, [](chsak4& p) { return p.phi; }); 
    fillak4(chsjetm, output_jets, [](chsak4& p) { return p.mass; }); 
    fillxcone4(xjetpt, output_xjets, [](xcone4& p) { return p.pt; }); 
    fillxcone4(xjeteta, output_xjets, [](xcone4& p) { return p.eta; }); 
    fillxcone4(xjetphi, output_xjets, [](xcone4& p) { return p.phi; }); 
    fillxcone4(xjetm, output_xjets, [](xcone4& p) { return p.mass; }); 

    tout->Fill();

    progressBar.Update(k, k);

  }

  hDTotalMCWeight->SetDirectory(0);
  fout->WriteTObject(hDTotalMCWeight); delete hDTotalMCWeight; hDTotalMCWeight = nullptr;
  fout->Write();
  fout->Close();

}
