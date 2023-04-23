#include <string>
#include <vector>
#include <iostream>
#include "TROOT.h"
#include "TString.h"
#include "TTree.h"
#include "TFile.h"
#include "TEntryListArray.h"
#include "TH1D.h"

// estimate the number of events after precuts
// save the skimmed output tree for mva training

int precut(std::string inputSname, std::string inputBname, std::string mycuts, std::string mycutb,
           std::string outputname)
{
  bool massPlot = false;

  gROOT->SetBatch(kTRUE);

  // signal
  TFile* inputS = TFile::Open(inputSname.c_str());
  TTree* signal = (TTree*)inputS->Get("Bfinder/ntphi");
  signal->AddFriend("hltanalysis/HltTree");
  signal->AddFriend("hiEvtAnalyzer/HiTree");
  signal->AddFriend("skimanalysis/HltTree");
  TTree* sghlt = (TTree*)inputS->Get("hltanalysis/HltTree");
  TTree* sghie = (TTree*)inputS->Get("hiEvtAnalyzer/HiTree");
  TTree* sgska = (TTree*)inputS->Get("skimanalysis/HltTree");

  // background
  TFile* inputB = TFile::Open(inputBname.c_str());
  TTree* background = (TTree*)inputB->Get("Bfinder/ntphi");
  background->AddFriend("hltanalysis/HltTree");
  background->AddFriend("hiEvtAnalyzer/HiTree");
  background->AddFriend("skimanalysis/HltTree");
  TTree* bghlt = (TTree*)inputB->Get("hltanalysis/HltTree");
  TTree* bghie = (TTree*)inputB->Get("hiEvtAnalyzer/HiTree");
  TTree* bgska = (TTree*)inputB->Get("skimanalysis/HltTree");

  // // TFile* outputB = new TFile(outputname.c_str(), "recreate");
  // TFile* outputB = new TFile(outputname.c_str(), "recreate");
  // outputB->mkdir("Bfinder");
  // outputB->mkdir("hltanalysis");
  // outputB->mkdir("hiEvtAnalyzer");
  // outputB->mkdir("skimanalysis");
  // outputB->cd("Bfinder");
  // TTree* outBtree = background->CloneTree(0);
  // outputB->cd("hltanalysis");
  // TTree* outBhlt = bghlt->CloneTree(0);
  // outputB->cd("hiEvtAnalyzer");
  // TTree* outBhie = bghie->CloneTree(0);
  // outputB->cd("skimanalysis");
  // TTree* outBska = bgska->CloneTree(0);

  TString cut = "(abs(Btktkmass-1.019455)<0.015)&&(((((abs(Btktkmass-1.019455)<0.015)&& TMath::Abs(Bmumumass-3.096916)<0.15 && Bpt > 0 && Bpt < 5 && (abs(Btrk1Eta)<2.4 && abs(Btrk2Eta)<2.4 && Btrk1Pt>0.0 && Btrk2Pt>0.0) && Btrk1Pt > 0.2 && Btrk2Pt > 0.2  && Bchi2cl > 0.05 && BsvpvDistance/BsvpvDisErr > 2.0)  && ( (Bpt < 2 && Bpt > 0) || (Bpt < 3 && Bpt > 2) || (Bpt < 5 && Bpt > 3)  )))  ||  ( Bpt > 3 && ((HBHENoiseFilterResult == 1 && pPAprimaryVertexFilter == 1 && pBeamScrapingFilter == 1 && HLT_HIL1DoubleMu0_v1 == 1 && (abs(PVz)<15))  &&  (Bmu1isTriggered == 1 && Bmu2isTriggered == 1 ) &&  (Bchi2cl > 0.05 && BsvpvDistance/BsvpvDisErr > 2.0)    && (TMath::Abs(By)<2.4&&TMath::Abs(Bmumumass-3.096916)<0.15&&((abs(Bmu1eta)<1.2&&Bmu1pt>3.5)||(abs(Bmu1eta)>1.2&&abs(Bmu1eta)<2.1&&Bmu1pt>(5.47-1.89*abs(Bmu1eta)))||(abs(Bmu1eta)>2.1&&abs(Bmu1eta)<2.4&&Bmu1pt>1.5))&&((abs(Bmu2eta)<1.2&&Bmu2pt>3.5)||(abs(Bmu2eta)>1.2&&abs(Bmu2eta)<2.1&&Bmu2pt>(5.47-1.89*abs(Bmu2eta)))||(abs(Bmu2eta)>2.1&&abs(Bmu2eta)<2.4&&Bmu2pt>1.5))&&Bmu1InPixelLayer>0&&(Bmu1InPixelLayer+Bmu1InStripLayer)>5&&Bmu2InPixelLayer>0&&(Bmu2InPixelLayer+Bmu2InStripLayer)>5&&Bmu1dxyPV<0.3&&Bmu2dxyPV<0.3&&Bmu1dzPV<20&&Bmu2dzPV<20&&Bmu1isTrackerMuon&&Bmu2isTrackerMuon&&Bmu1isGlobalMuon&&Bmu2isGlobalMuon)  && ( Btrk1Pt > 0.2 && Btrk2Pt > 0.2 && abs(Btrk1Eta-0.0) < 2.4 && abs(Btrk2Eta-0.0) < 2.4  && Btrk1highPurity  && Btrk2highPurity  && Btrk1PixelHit + Btrk1StripHit > 10  && Btrk2PixelHit + Btrk2StripHit > 10) &&  (Btrk1PtErr/Btrk1Pt < 0.1)  &&  (Btrk2PtErr/Btrk2Pt < 0.1)    && Btrk1Chi2ndf/(Btrk1nStripLayer+Btrk1nPixelLayer) < 0.18   && Btrk2Chi2ndf/(Btrk2nStripLayer+Btrk2nPixelLayer) < 0.18 )))";

  std::vector<double> ptbins = {7, 10, 15, 20, 30, 50};

  ulong max = TVirtualTreePlayer::kMaxEntries;

  outputS->cd();
  // TTree* slisttree = new TTree("slisttree", "signal entry list");
  signal->Draw(">>slist", cut, "entrylistarray", max);
  TEntryListArray *slist = (TEntryListArray*)gDirectory->Get("slist");
  // slisttree->Branch("slist", slist);
  // slisttree->Fill();


  signal->SetEntryList(slist);
  for (unsigned i = 0; i < ptbins.size() - 1; ++i) {
    auto outname = TString::Format("BsMC_%.0f_%.0f.root", ptbins[i], ptbins[i+1]);
    // massBefore[i] = new TH1D(iname, "B mass", 50, 5.36682 - 0.3, 5.36682 + 0.3);
    TString sel = TString::Format("Bpt > %f && Bpt << %f", ptbins[i], ptbins[i + 1]);
    // std::cout << sel << ", " <<  tevt << "\n";

    // massBefore[i]->Write();
    TFile* outputS = new TFile(outname, ), "recreate");
    outputS->mkdir("Bfinder");
    outputS->mkdir("hltanalysis");
    outputS->mkdir("hiEvtAnalyzer");
    outputS->mkdir("skimanalysis");
    outputS->cd("Bfinder");
    TTree* outStree = signal->CloneTree(0);
    outputS->cd("hltanalysis");
    TTree* outShlt = sghlt->CloneTree(0);
    outputS->cd("hiEvtAnalyzer");
    TTree* outShie = sghie->CloneTree(0);
    outputS->cd("skimanalysis");
    TTree* outSska = sgska->CloneTree(0);

    for(auto i = 0; i < slist->GetN(); ++i){
      auto entry = slist->GetEntry(i);
      signal->GetEntry(entry);
      outStree->Fill();
      outShlt->Fill();
      outShie->Fill();
      outSska->Fill();
      }
    slist->Delete();
  }

  return 0;

  // outputB->cd();
  // TFile* precutB = new TFile("precutB.root", "recreate");

  // background->Draw(">>elist", cut, "entrylistarray", max);
  // TEntryListArray *list = (TEntryListArray*)gDirectory->Get("elist");
  // // try {
  // //   list->Print();
  // //   list->Write();
  // //   // list->Print("all");
  // //   // list->GetSubListForEntry(1098)->Print("all");
  // //   precutB->Close();
  // // } catch (const std::exception& e) {
  // //   std::cout << " a standard exception was caught, with message '"
  // //             << e.what() << "'\n";
  // // }


  // // evaluate eff in each pT bin
  // if (massPlot) {
  //   int nevt = background->GetEntries();
  //   int tevt;
  //   std::vector<TH1D*> massBefore(ptbins.size() - 1);
  //   std::vector<TH1D*> massAfter(ptbins.size() - 1);

  //   TFile* outplot = new TFile("massplot.root", "recreate");
  //   for (unsigned i = 0; i < ptbins.size() - 1; ++i) {
  //     auto iname = TString::Format("massb_%.0f_%.0f", ptbins[i], ptbins[i+1]);
  //     massBefore[i] = new TH1D(iname, "B mass", 50, 5.36682 - 0.3, 5.36682 + 0.3);
  //     TString sel = TString::Format("Bpt > %f && Bpt << %f", ptbins[i], ptbins[i + 1]);
  //     tevt = background->Draw("Bmass >> " + iname, sel, "", max);
  //     std::cout << sel << ", " <<  tevt << "\n";
  //     massBefore[i]->Write();
  //   }
  //   background->SetEntryList(list);
  //   for (unsigned i = 0; i < ptbins.size() - 1; ++i) {
  //     auto fname = TString::Format("massa_%.0f_%.0f", ptbins[i], ptbins[i+1]);
  //     massAfter[i] = new TH1D(fname, "B mass", 50, 5.36682 - 0.3, 5.36682 + 0.3);
  //     TString sel = TString::Format("Bpt > %f && Bpt << %f", ptbins[i], ptbins[i + 1]);
  //     tevt = background->Draw("Bmass >> " + fname, sel, "", max);
  //     std::cout << sel << ", " <<  tevt << "\n";
  //     massAfter[i]->Write();
  //   }
  //   outplot->Close();
  // }

  // for(auto i = 0; i < list->GetN(); ++i){
  //   auto entry = list->GetEntry(i);
  //   background->GetEntry(entry);
  //   outBtree->Fill();
  //   outBhlt->Fill();
  //   outBhie->Fill();
  //   outBska->Fill();
  //   }
  // list->Delete();

  // outputB->Write();
  // // outputB->ls();
  // outputB->Close();

  // return 0;
}
