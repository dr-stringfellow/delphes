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

#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"

#include "fastjet/plugins/CDFCones/fastjet/CDFJetCluPlugin.hh"
#include "fastjet/plugins/CDFCones/fastjet/CDFMidPointPlugin.hh"
#include "fastjet/plugins/SISCone/fastjet/SISConePlugin.hh"

#include "fastjet/contribs/Nsubjettiness/ExtraRecombiners.hh"
#include "fastjet/contribs/Nsubjettiness/Njettiness.hh"
#include "fastjet/contribs/Nsubjettiness/NjettinessPlugin.hh"
#include "fastjet/contribs/Nsubjettiness/Nsubjettiness.hh"

#include "fastjet/contribs/ValenciaPlugin/ValenciaPlugin.hh"

#include "fastjet/contribs/RecursiveTools/SoftDrop.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"

#include "modules/EnergyCorrelations.h"

using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;

static int NMAX = 1000;
const std::vector<double> betas = {0.5, 1.0, 2.0, 4.0};
const std::vector<int> ibetas = {0,1,2,3};
const std::vector<int> Ns = {1,2,3,4};
const std::vector<int> orders = {1,2,3};

//---------------------------------------------------------------------------


TString makeECFString(svj::ECFParams p) {
  return TString::Format("ECFN_%i_%i_%.2i",p.order,p.N,int(10*betas.at(p.ibeta)));
}

struct PFCand
{
  float pt = 0;
  float eta = 0;
  float phi = 0;
  float e = 0;
  float puppi = 1;
  float charge = 1;
  float hardfrac = 1;  
  float vtxid = -1;
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
  TBranch* partbranch = (TBranch*)itree->GetBranch("Particle");
  TBranch* pfbranch = (TBranch*)itree->GetBranch("ParticleFlowCandidate");
  TBranch* genjetbranch = (TBranch*)itree->GetBranch("GenJet");
  std::cout << "NEVT: " << nevt << std::endl;
  vector<PFCand> input_particles;
  vector<PFCand> output_particles;
  output_particles.reserve(NMAX);

  vector<float> vpt, veta, vphi, ve, vpuppi, vpdgid, vhardfrac;
  vpt.reserve(NMAX); veta.reserve(NMAX); vphi.reserve(NMAX); 
  ve.reserve(NMAX); vpuppi.reserve(NMAX); vpdgid.reserve(NMAX); 
  vhardfrac.reserve(NMAX);

  float genjetpt=-99., genjeteta=-99., genjetphi=-99., genjetm=-99.;
  float puppijetpt=-99., puppijeteta=-99., puppijetphi=-99., puppijetm=-99.;

  // jet branches
  TBranch* b_genjetpt = tout->Branch("genjetpt",&genjetpt, "genjetpt/F");
  TBranch* b_genjeteta = tout->Branch("genjeteta",&genjeteta, "genjeteta/F");
  TBranch* b_genjetphi = tout->Branch("genjetphi",&genjetphi, "genjetphi/F");
  TBranch* b_genjetm = tout->Branch("genjetm",&genjetm, "genjetm/F");
  TBranch* b_puppijetpt = tout->Branch("puppijetpt",&puppijetpt, "puppijetpt/F");
  TBranch* b_puppijeteta = tout->Branch("puppijeteta",&puppijeteta, "puppijeteta/F");
  TBranch* b_puppijetphi = tout->Branch("puppijetphi",&puppijetphi, "puppijetphi/F");
  TBranch* b_puppijetm = tout->Branch("puppijetm",&puppijetm, "puppijetm/F");

  //for (auto p : ecfParams) {
  //TString ecfn(makeECFString(p));
  //tout->Branch("fj"+ecfn,&(fjECFNs[p][0]),"fj"+ecfn+"/F");
  //}

  ExRootProgressBar progressBar(nevt);
  
  auto comp_p4 = [](auto &a, auto &b) { return a.pt > b.pt; };

  fastjet::JetDefinition *jetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm,0.8);

  // Soft drop -- only needed for large radius jets
  double radius = 0.8;
  double sdZcut = 0.1;
  double sdBeta = 0.;
  
  fastjet::contrib::SoftDrop softDrop = fastjet::contrib::SoftDrop(sdBeta,sdZcut,radius);

  for (unsigned int k=0; k<nevt; k++){
    itree->GetEntry(k);

    if (k%100==0)
      std::cout << k << " / " << nevt << std::endl;

    unsigned int ngenjets = genjetbranch->GetEntries();
    ngenjets = itree->GetLeaf("GenJet_size")->GetValue(0);
    TLorentzVector genjet;
    for (unsigned int j=0; j<ngenjets; j++){
      genjet.SetPtEtaPhiM(itree->GetLeaf("GenJet.PT")->GetValue(j),itree->GetLeaf("GenJet.Eta")->GetValue(j),itree->GetLeaf("GenJet.Phi")->GetValue(j),itree->GetLeaf("GenJet.Mass")->GetValue(j));
      genjetpt = genjet.Pt();
      genjeteta = genjet.Eta();
      genjetphi = genjet.Phi();
      genjetm = genjet.M();
      break;// only get leading gen jet
    }

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
      tmppf.hardfrac = itree->GetLeaf("ParticleFlowCandidate.hardfrac")->GetValue(j);
      if (itree->GetLeaf("ParticleFlowCandidate.Charge")->GetValue(j)!=0){
        if (itree->GetLeaf("ParticleFlowCandidate.hardfrac")->GetValue(j)==1)
          tmppf.vtxid = 0;
        else
          tmppf.vtxid = 1;
      }
      else
        tmppf.vtxid = -1;
      input_particles.push_back(tmppf);
    }

    // sorting input particles by pT
    sort(input_particles.begin(), input_particles.end(), comp_p4);

    vector<fastjet::PseudoJet> finalStates_puppi;
    vector<fastjet::PseudoJet> finalStates_truth;
    vector<fastjet::PseudoJet> finalStates_pf;
    int pfid = 0;
    for(auto &p : input_particles){
      if (p.vtxid==1)
	continue; // Charged Hadron Subtraction
      TLorentzVector tmp;
      tmp.SetPtEtaPhiE(p.pt,p.eta,p.phi,p.e);

      // Puppi
      if (p.puppi>0.){
	fastjet::PseudoJet curjet_puppi(p.puppi*tmp.Px(), p.puppi*tmp.Py(), p.puppi*tmp.Pz(), p.puppi*tmp.E());
	curjet_puppi.set_user_index(pfid);
	finalStates_puppi.emplace_back(curjet_puppi);
      }
      pfid += 1;
    }

    fastjet::ClusterSequence seq_puppi(sorted_by_pt(finalStates_puppi), *jetDef);

    float minpt = 400;

    vector<fastjet::PseudoJet> allJets_puppi(sorted_by_pt(seq_puppi.inclusive_jets(minpt)));


    for (auto& puppiJet : allJets_puppi) {
      if (puppiJet.perp() < minpt)
	break;

      std::unique_ptr<svj::ECFCalculator>ecfcalc{nullptr};
      ecfcalc.reset(new svj::ECFCalculator());

      fastjet::PseudoJet sdJet = (softDrop)(puppiJet);
      vector<fastjet::PseudoJet> sdConstituents = sorted_by_pt(sdJet.constituents());

      TLorentzVector tmp;
      tmp.SetPtEtaPhiM(puppiJet.perp(),puppiJet.eta(),puppiJet.phi(),puppiJet.m());
      if (tmp.DeltaR(genjet)<0.4){
	puppijetpt = tmp.Pt();
	puppijeteta = tmp.Eta();
	puppijetphi = tmp.Phi();
	puppijetm = sdJet.m();
	break;
      }

      // get ecfs                                                                                                                                 
      unsigned nFilter = min(30, (int)sdConstituents.size());
      vector<fastjet::PseudoJet> sdConstsFiltered(sdConstituents.begin(), sdConstituents.begin() + nFilter);

      svj::ECFParams ep;
      ecfcalc->calculate(sdConstsFiltered);
      for (auto iter = ecfcalc->begin(); iter != ecfcalc->end(); ++iter) {
	int N = iter.get<svj::ECFCalculator::nP>();
	int o = iter.get<svj::ECFCalculator::oP>();
	int beta = iter.get<svj::ECFCalculator::bP>();
	float ecf(iter.get<svj::ECFCalculator::ecfP>());
	//genJetInfo.ecfs[o][N][beta] = ecf;
	std::cout << "ECF: " << ecf << " beta: " << beta << " o: " << o << " N: " << N << std::endl;
	ep.order = o + 1; ep.N = N + 1, ep.ibeta = beta;
	//gt.fjECFNs[ep][0] = ecf;
      }

    }


    tout->Fill();

    progressBar.Update(k, k);

  }

  fout->Write();
  fout->Close();

}
