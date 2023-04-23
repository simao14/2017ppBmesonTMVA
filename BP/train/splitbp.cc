#include <string>
#include <vector>
#include <iostream>
#include "TROOT.h"
#include "TString.h"
#include "TTree.h"
#include "TFile.h"
#include "TEntryList.h"
#include "TH1D.h"

// estimate the number of events after precuts
// save the skimmed output tree for mva training

int splitbp()
{
  gROOT->SetBatch(kTRUE);
  TString inputSname ="~/dat/presel/BPMC.root";
  TString inputBname ="~/dat/presel/new/BPData4.root";

  // signal
  TFile* inputS = TFile::Open(inputSname);
  TTree* signal = (TTree*)inputS->Get("Bfinder/ntKp");
  signal->AddFriend("hltanalysis/HltTree");
  signal->AddFriend("hiEvtAnalyzer/HiTree");
  signal->AddFriend("skimanalysis/HltTree");
  TTree* sghlt = (TTree*)inputS->Get("hltanalysis/HltTree");
  TTree* sghie = (TTree*)inputS->Get("hiEvtAnalyzer/HiTree");
  TTree* sgska = (TTree*)inputS->Get("skimanalysis/HltTree");
  TTree* sggen = (TTree*)inputS->Get("Bfinder/ntGen");
  TTree* sgroot = (TTree*)inputS->Get("Bfinder/root");

  signal->AddFriend(sggen);
  signal->AddFriend(sgroot);

  // background
  TFile* inputB = TFile::Open(inputBname);
  TTree* background = (TTree*)inputB->Get("Bfinder/ntKp");
  background->AddFriend("hltanalysis/HltTree");
  background->AddFriend("hiEvtAnalyzer/HiTree");
  background->AddFriend("skimanalysis/HltTree");
  TTree* bghlt = (TTree*)inputB->Get("hltanalysis/HltTree");
  TTree* bghie = (TTree*)inputB->Get("hiEvtAnalyzer/HiTree");
  TTree* bgska = (TTree*)inputB->Get("skimanalysis/HltTree");
  TTree* bgroot = (TTree*)inputB->Get("Bfinder/root");
  background->AddFriend(bgroot);


  // std::vector<double> ptbins = {5, 7, 10, 15, 20, 50, 60};
  // std::vector<double> ptbins = {7, 10, 15, 20, 50};
  std::vector<double> ptbins = {5, 60};

  ulong max = TVirtualTreePlayer::kMaxEntries;



  for (unsigned i = 0; i < ptbins.size() - 1; ++i) {
    std::cout << "Selecting pT " << ptbins[i] << "\n";

    TString sel = TString::Format("Bpt > %f && Bpt < %f", ptbins[i], ptbins[i + 1]);
    TString gsel = TString::Format("Gpt > %f && Gpt < %f", ptbins[i], ptbins[i + 1]);
    TString sgsel = TString::Format("(%s) || (%s)", sel.Data(), gsel.Data());
    auto outname = TString::Format("../sample/BPMC_%.0f_%.0f.root", ptbins[i], ptbins[i+1]);

    signal->Draw(">>slist", sel, "entrylist", max);
    TEntryList *slist = (TEntryList*)gDirectory->Get("slist");

    TFile* outputS = new TFile(outname, "recreate");
    outputS->mkdir("Bfinder");
    outputS->mkdir("hltanalysis");
    outputS->mkdir("hiEvtAnalyzer");
    outputS->mkdir("skimanalysis");
    outputS->cd("Bfinder");
    TTree* outStree = signal->CloneTree(0);
    TTree* outGen = sggen->CloneTree(0);
    TTree* outSroot = sgroot->CloneTree(0);
    outputS->cd("hltanalysis");
    TTree* outShlt = sghlt->CloneTree(0);
    outputS->cd("hiEvtAnalyzer");
    TTree* outShie = sghie->CloneTree(0);
    outputS->cd("skimanalysis");
    TTree* outSska = sgska->CloneTree(0);

    for(auto ient = 0; ient < slist->GetN(); ++ient){
      auto entry = slist->GetEntry(ient);
      signal->GetEntry(entry);
      outStree->Fill();
      outSroot->Fill();
      outShlt->Fill();
      outShie->Fill();
      outSska->Fill();

      sggen->GetEntry(entry);
      outGen->Fill();
    }
    outputS->Write();
    outputS->Close();

    background->Draw(">>blist", sel, "entrylist", max);
    TEntryList *blist = (TEntryList*)gDirectory->Get("blist");

    outname = TString::Format("../sample/BPData_%.0f_%.0f.root", ptbins[i], ptbins[i+1]);
    TFile* outputB = new TFile(outname, "recreate");
    outputB->mkdir("Bfinder");
    outputB->mkdir("hltanalysis");
    outputB->mkdir("hiEvtAnalyzer");
    outputB->mkdir("skimanalysis");
    outputB->cd("Bfinder");
    TTree* outBtree = background->CloneTree(0);
    TTree* outBroot = bgroot->CloneTree(0);
    outputB->cd("hltanalysis");
    TTree* outBhlt = bghlt->CloneTree(0);
    outputB->cd("hiEvtAnalyzer");
    TTree* outBhie = bghie->CloneTree(0);
    outputB->cd("skimanalysis");
    TTree* outBska = bgska->CloneTree(0);

    for(auto ient = 0; ient < blist->GetN(); ++ient){
      auto entry = blist->GetEntry(ient);
      // std::cout << entry << "\n";

      background->GetEntry(entry);
      outBtree->Fill();
      outBroot->Fill();
      outBhlt->Fill();
      outBhie->Fill();
      outBska->Fill();
      }
    outputB->Write();
    outputB->Close();

  }

  return 0;

}
