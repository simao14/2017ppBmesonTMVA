#include "TROOT.h"
#include "TH1.h"
#include "TTree.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TSystem.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom.h"
#include <TH1D.h>
#include <RtypesCore.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>

#include <algorithm>
#include <utility>
#include <vector>


constexpr bool noBDT = false;
// constexpr bool noBDT = true;
constexpr bool filter = false;
constexpr bool jpsi = false;
//constexpr bool jpsi = true;


// constexpr std::array<double, 5> nomBDT = {0.23, 0.20, 0.12, 0.11, 0.08};
// constexpr std::array<double, 5> preBDT = {0.51, 0.45, 0.38, 0.30, 0.25};
constexpr std::array<double, 6> nomBDT = {0.09, 0.08, 0.07, 0.0, 0.02, 0.04};

//#include "his.h"
using namespace std;

using std::cout;
using std::endl;

enum Tracking{
  loose = 0,
  standard = 1,
  tight = 2,
};

std::map<Tracking, double> ptErr{
  {Tracking::standard, 0.1},
  {Tracking::tight, 0.05},
  {Tracking::loose, 0.15}
};

std::map<Tracking, double> chi2Nlayer{
  {Tracking::standard, 0.18},
  {Tracking::tight, 0.15},
  {Tracking::loose, 0.18}
};

void BPPreSel(int doMC, bool nominal){
/*
  TString inList="~/dat/Unskimmed/BPData.root";
  TFile* flist = new TFile(inList.Data());
  TTree* tlist = flist->Get("Bfinder/root");
  int EvtNo;
  int flumi;
  tlist->SetBranchStatus("*", 0);
  tlist->SetBranchStatus("EvtInfo.EvtNo", 1);
  tlist->SetBranchStatus("EvtInfo.LumiNo", 1);
  tlist->SetBranchAddress("EvtInfo.EvtNo", &EvtNo);
  tlist->SetBranchAddress("EvtInfo.LumiNo", &flumi);
  unsigned umax = -1;
  std::cout << "filter list contains " << tlist->GetEntries() << " entries" << "\n";
  using is_pair = std::pair<Int_t, UShort_t>;
  // std::vector<char> evtList(umax, 0);
  std::vector<is_pair> filterList(tlist->GetEntries());
  // std::vector<is_pair> filterList(100000);

  // std::vector< std::vector<UShort_t> > lumiList(1, std::vector<UShort_t>(umax, 0) );
  cout << "initialized" << "\n";

  if (filter) {
    for (auto i = 0; i < tlist->GetEntries(); ++i) {
    // for (auto i = 0; i < 100000; ++i) {
      // if (i % 100000 == 0) {
      //   std::cout << "\r entry: " << i;
      // }
      tlist->GetEntry(i);
      unsigned uevtno = EvtNo;
      UShort_t clumi = (UShort_t) flumi;
      filterList[i] = std::make_pair(EvtNo, clumi);
      // evtList[uevtno] = evtList[uevtno] + 1;
      // lumiList[uevtno].push_back(clumi);
      // if (evtList[uevtno] > 3) {
      //   cout << "warning: more than 3 list needed" << "\n";
      //   cout << "evtno: " << EvtNo << ", no:" << (int) evtList[uevtno] << "\n";
      // }
    }
    // std::cout << "\n built list" << std::endl;
    // auto max = std::max_element(evtList.begin(), evtList.end());
    // cout << "max repeated evt:" << (int) *max << "\n";

    // sorting
    std::sort(filterList.begin(), filterList.end());
  }
  flist->Close();
  // for (auto p : filterList) {
  //   cout << p.first << ", " << p.second << "\n";
  // }
  // return;
*/

	int EvtNo;

	TString FileName;


	// if(doMC == 0) FileName = "~/dat/presel/BPData.root";
	// if(doMC == 1) FileName = "~/dat/presel/BPMC.root";

	// if(doMC == 0) FileName = "~/dat/presel/output/BP_Data_BDTs_pre.root";
	// if(doMC == 1) FileName = "~/dat/presel/output/BP_MC_BDTs_pre.root";
  // if (jpsi) {
  //   FileName = "~/dat/presel/output/jpsinp_pre.root";
  // }

  // Use nominal BDT score regardless of preselection
  if(doMC == 0) FileName = "../output/BP_Data_BDTs_nom.root";
  // if(doMC == 1) FileName = "../output/BP_MC_BDTs_nom.root";
  if(doMC == 1) FileName = "../output/BP_MC_BDTs_pre_tnp.root";
  if (jpsi) {
    // FileName = "../output/jpsinp_nom.root";
    FileName = "../output/jpsinp_pre.root";
  }

  if (nominal) {
    if(doMC == 0) FileName = "../output/BP_Data_BDTs_nom.root";
    // if(doMC == 1) FileName = "../output/BP_MC_BDTs_nom.root";
    if(doMC == 1) FileName = "../output/BP_MC_BDTs_nom.root";
    if (jpsi) {
      // FileName = "../output/jpsinp_nom.root";
      FileName = "../output/jpsinp_nom.root";
    }
  }


	//	TString FileName2 = "../../FirstSamples/WithoutBDT/BPMC.root";
	//	FileName="/data/szhaozho/Bmeson2017pp/BDTOutput/AllMerge/BPMCAllBDT.root";
	//FileName="/data/szhaozho/ppNewTMVA/CMSSW_10_3_2/src/Bs/SamplesNew2/BsMCAllNew.root"; 
	//FileName="/data/szhaozho/ppNewTMVA/CMSSW_10_3_2/src/Bs/SamplesNew/BsMCAllBDT.root"; 
	//FileName ="BsMC.root";

	TFile * fin = new TFile(FileName.Data());

	TTree * root = (TTree * ) fin->Get("Bfinder/root");

	TTree * ntKp = (TTree * ) fin->Get("Bfinder/ntKp");
	//	TTree * ntHi = (TTree * ) fin->Get("hiEvtAnalyzer/HiTree");
	TTree * ntSkim = (TTree * ) fin->Get("skimanalysis/HltTree");
	TTree * ntHlt = (TTree *) fin->Get("hltanalysis/HltTree");

	//	TFile * fin2 = new TFile(FileName2.Data());
	TTree * ntHi = (TTree * ) fin->Get("hiEvtAnalyzer/HiTree");
	
	TTree * BDT1 = (TTree *) fin->Get("BDT_pt_3_5");
    TTree * BDT2 = (TTree *) fin->Get("BDT_pt_5_7");
    TTree * BDT3 = (TTree *) fin->Get("BDT_pt_7_10");
	TTree * BDT4 = (TTree *) fin->Get("BDT_pt_10_15");
	TTree * BDT5 = (TTree *) fin->Get("BDT_pt_15_20");
	TTree * BDT6 = (TTree *) fin->Get("BDT_pt_20_50");
	// TTree * BDT7 = (TTree *) fin->Get("BDT_pt_50_60");
	// TTree * BDT8 = (TTree *) fin->Get("BDT_pt_1_2");
	// TTree * BDT9 = (TTree *) fin->Get("BDT_pt_50_100");

	/*


	   TF1 * func;

	   if(CentMin == 0 && CentMax == 90) func = new TF1("func","1/0.144708+TMath::Exp(-1.035696*(x-15.321432))+TMath::Exp(-0.204131*(x-30.289313))",5,50);

	   if(CentMin == 0 && CentMax == 30) func = new TF1("func","1/0.131733+TMath::Exp(-1.051514*(x-15.468517))+TMath::Exp(-0.204272*(x-31.002324))",5,50);

	//	if(CentMin == 30 && CentMax == 40) func = new TF1("func","1/0.131733+TMath::Exp(-1.051514*(x-15.468517))+TMath::Exp(-0.204272*(x-31.002324))",5,50);
	//	if(CentMin == 40 && CentMax == 90) func = new TF1("func","1/0.131733+TMath::Exp(-1.051514*(x-15.468517))+TMath::Exp(-0.204272*(x-31.002324))",5,50);


	if(CentMin == 30 && CentMax == 90) func = new TF1("func","1/0.204307+TMath::Exp(-1.041731*(x-14.608514))+TMath::Exp(-0.206513*(x-27.599694))",5,50);
	*/




	const int NCand = 8000;

	int run;
	int lumi;
	ULong64_t evt;
	int hiBin;
	Float_t PVz;
	Int_t pBeamScrapingFilter;
	Int_t pPAprimaryVertexFilter;
	Int_t HBHENoiseFilterResult;









	//finSystWeight->Close();




	Int_t   Bsize;
	Float_t Btrk1Pt[NCand];
	Float_t Btrk2Pt[NCand];

	Float_t Btrk1PtErr[NCand];
	Float_t Btrk2PtErr[NCand];


	Float_t Bchi2cl[NCand];
	Float_t BsvpvDistance[NCand];
	Float_t BsvpvDisErr[NCand];
	Float_t Bpt[NCand];
	Float_t Btrk1Eta[NCand];
	Float_t Btrk2Eta[NCand];
	Float_t By[NCand];

	Float_t BsvpvDistance_2D[NCand];
	Float_t BsvpvDisErr_2D[NCand];

	Float_t Bmass[NCand];


	Float_t Bmumumass[NCand];
	Float_t Bmu1eta[NCand];
	Float_t Bmu1pt[NCand];
	Float_t Bmu2eta[NCand];
	Float_t Bmu2pt[NCand];

	Bool_t Bmu1TMOneStationTight[NCand];
	Int_t Bmu1InPixelLayer[NCand];
	Int_t Bmu1InStripLayer[NCand];

	Bool_t Bmu2TMOneStationTight[NCand];	
	Int_t Bmu2InPixelLayer[NCand];
	Int_t Bmu2InStripLayer[NCand];


	Bool_t Bmu1isGlobalMuon[NCand];
	Bool_t Bmu2isGlobalMuon[NCand];


	Bool_t Bmu1isTrackerMuon[NCand];
	Bool_t Bmu2isTrackerMuon[NCand];

	Float_t Bmu1dxyPV[NCand];
	Float_t Bmu2dxyPV[NCand];

	Float_t Bmu1dzPV[NCand];
	Float_t Bmu2dzPV[NCand];

	Bool_t Btrk1highPurity[NCand];
	Bool_t Btrk2highPurity[NCand];

	Float_t Btktkmass[NCand];

	Float_t Btrk1PixelHit[NCand];
	Float_t Btrk2PixelHit[NCand];

	Float_t Btrk1StripHit[NCand];
	Float_t Btrk2StripHit[NCand];

	Float_t Btrk1Chi2ndf[NCand];
	Float_t Btrk2Chi2ndf[NCand];


	Float_t Btrk1nStripLayer[NCand];
	Float_t Btrk2nStripLayer[NCand];

	Float_t Btrk1nPixelLayer[NCand];
	Float_t Btrk2nPixelLayer[NCand];


	Double_t BDT_pt_3_5[NCand];
	Double_t BDT_pt_5_7[NCand];
	Double_t BDT_pt_7_10[NCand];
	Double_t BDT_pt_10_15[NCand];
	Double_t BDT_pt_15_20[NCand];
	Double_t BDT_pt_20_50[NCand];
	Double_t BDT_pt_50_60[NCand];
	// Float_t BDT_pt_2_3[NCand];
	// Float_t BDT_pt_1_2[NCand];
	// Double_t BDT_pt_50_100[NCand];

	Float_t Bgen[NCand];
	Float_t Bgenpt[NCand];

	//Additional Variables//

	Float_t Bd0[NCand];
	Float_t Bd0Err[NCand];
	Float_t Bdtheta[NCand];
	Float_t Balpha[NCand];
	Float_t Btrk1Dz1[NCand];
	Float_t Btrk2Dz1[NCand];

	Float_t Btrk1DzError1[NCand];
	Float_t Btrk2DzError1[NCand];

	Float_t Btrk1Dxy1[NCand];
	Float_t Btrk2Dxy1[NCand];

	Float_t Btrk1DxyError1[NCand];
	Float_t Btrk2DxyError1[NCand];


	Float_t Bmumueta[NCand];
	Float_t Bmumuphi[NCand];
	Float_t Bmumupt[NCand];

	//DONE Additional Variables//



	Float_t Bmu1TrgMatchFilterE[NCand];
	Float_t Bmu2TrgMatchFilterE[NCand];




	Float_t DataWeight_5_10[NCand];
	Float_t DataWeight_10_15[NCand];
	Float_t DataWeight_15_20[NCand];
	Float_t DataWeight_20_50[NCand];



	Int_t HLT_HIL1DoubleMu0_v1;
	ntHlt->SetBranchAddress("HLT_HIL1DoubleMu0_v1",&HLT_HIL1DoubleMu0_v1);

	Float_t pthat;
	Float_t weight;


	ntSkim->SetBranchAddress("HBHENoiseFilterResult",&HBHENoiseFilterResult);
	ntSkim->SetBranchAddress("pPAprimaryVertexFilter",&pPAprimaryVertexFilter);
	ntSkim->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter);

	Int_t nMult;


   if (!jpsi) {
    root->SetBranchAddress("EvtInfo.nMult",&nMult);
    root->SetBranchAddress("EvtInfo.LumiNo", &lumi);
   }

	ntKp->SetBranchAddress("Bsize",&Bsize);
	ntKp->SetBranchAddress("PVz",&PVz);

	ntKp->SetBranchAddress("Btrk1Pt",Btrk1Pt);
	ntKp->SetBranchAddress("Btrk2Pt",Btrk2Pt);
	ntKp->SetBranchAddress("Btrk1PtErr",Btrk1PtErr);
	ntKp->SetBranchAddress("Btrk2PtErr",Btrk2PtErr);
	ntKp->SetBranchAddress("Bchi2cl",Bchi2cl);
	ntKp->SetBranchAddress("BsvpvDistance",BsvpvDistance);
	ntKp->SetBranchAddress("BsvpvDisErr",BsvpvDisErr);
	ntKp->SetBranchAddress("Bpt",Bpt);
	ntKp->SetBranchAddress("By",By);
	ntKp->SetBranchAddress("Btrk1Eta",Btrk1Eta);
	ntKp->SetBranchAddress("Btrk2Eta",Btrk2Eta);
	ntKp->SetBranchAddress("Bmass",Bmass);
	ntKp->SetBranchAddress("Bgen",Bgen);
	ntKp->SetBranchAddress("Bgenpt",Bgenpt);

	ntKp->SetBranchAddress("BsvpvDistance_2D",BsvpvDistance_2D);
	ntKp->SetBranchAddress("BsvpvDisErr_2D",BsvpvDisErr_2D);



	ntKp->SetBranchAddress("Bmu1TrgMatchFilterE",Bmu1TrgMatchFilterE);
	ntKp->SetBranchAddress("Bmu2TrgMatchFilterE",Bmu2TrgMatchFilterE);

	ntKp->SetBranchAddress("Bmumumass",Bmumumass);
	ntKp->SetBranchAddress("Bmu1eta",Bmu1eta);
	ntKp->SetBranchAddress("Bmu2eta",Bmu2eta);
	ntKp->SetBranchAddress("Bmu1pt",Bmu1pt);
	ntKp->SetBranchAddress("Bmu2pt",Bmu2pt);

	ntKp->SetBranchAddress("Bmu1TMOneStationTight",Bmu1TMOneStationTight);
	ntKp->SetBranchAddress("Bmu1InPixelLayer",Bmu1InPixelLayer);
	ntKp->SetBranchAddress("Bmu1InStripLayer",Bmu1InStripLayer);

	ntKp->SetBranchAddress("Bmu2TMOneStationTight",Bmu2TMOneStationTight);
	ntKp->SetBranchAddress("Bmu2InPixelLayer",Bmu2InPixelLayer);
	ntKp->SetBranchAddress("Bmu2InStripLayer",Bmu2InStripLayer);


	ntKp->SetBranchAddress("Bmu1isGlobalMuon",Bmu1isGlobalMuon);
	ntKp->SetBranchAddress("Bmu2isGlobalMuon",Bmu2isGlobalMuon);

	ntKp->SetBranchAddress("Bmu1isTrackerMuon",Bmu1isTrackerMuon);
	ntKp->SetBranchAddress("Bmu2isTrackerMuon",Bmu2isTrackerMuon);


	ntKp->SetBranchAddress("Bmu1dxyPV",Bmu1dxyPV);
	ntKp->SetBranchAddress("Bmu2dxyPV",Bmu2dxyPV);
	ntKp->SetBranchAddress("Bmu1dzPV",Bmu1dzPV);
	ntKp->SetBranchAddress("Bmu2dzPV",Bmu2dzPV);


	ntKp->SetBranchAddress("Btrk1highPurity",Btrk1highPurity);
	ntKp->SetBranchAddress("Btrk2highPurity",Btrk2highPurity);

	ntKp->SetBranchAddress("Btktkmass",Btktkmass);


	ntKp->SetBranchAddress("Btrk1PixelHit",Btrk1PixelHit);
	ntKp->SetBranchAddress("Btrk2PixelHit",Btrk2PixelHit);
	ntKp->SetBranchAddress("Btrk1StripHit",Btrk1StripHit);
	ntKp->SetBranchAddress("Btrk2StripHit",Btrk2StripHit);



	ntKp->SetBranchAddress("Btrk1Chi2ndf",Btrk1Chi2ndf);
	ntKp->SetBranchAddress("Btrk2Chi2ndf",Btrk2Chi2ndf);



	ntKp->SetBranchAddress("evt",&evt);	



	ntKp->SetBranchAddress("Btrk1nStripLayer",Btrk1nStripLayer);	
	ntKp->SetBranchAddress("Btrk2nStripLayer",Btrk2nStripLayer);
	ntKp->SetBranchAddress("Btrk1nPixelLayer",Btrk1nPixelLayer);
	ntKp->SetBranchAddress("Btrk2nPixelLayer",Btrk2nPixelLayer);

	//Additional//




	//Additional Variables//

	ntKp->SetBranchAddress("Bd0",Bd0);
	ntKp->SetBranchAddress("Bd0Err",Bd0Err);
	ntKp->SetBranchAddress("Bdtheta",Bdtheta);
	ntKp->SetBranchAddress("Balpha",Balpha);
	ntKp->SetBranchAddress("Btrk1Dz1",Btrk1Dz1);

	ntKp->SetBranchAddress("Btrk2Dz1",Btrk2Dz1);
	ntKp->SetBranchAddress("Btrk1DzError1",Btrk1DzError1);
	ntKp->SetBranchAddress("Btrk2DzError1",Btrk2DzError1);
	ntKp->SetBranchAddress("Btrk1Dxy1",Btrk1Dxy1);
	ntKp->SetBranchAddress("Btrk2Dxy1",Btrk2Dxy1);
	ntKp->SetBranchAddress("Btrk1DxyError1",Btrk1DxyError1);
	ntKp->SetBranchAddress("Btrk2DxyError1",Btrk2DxyError1);
	ntKp->SetBranchAddress("Bmumueta",Bmumueta);
	ntKp->SetBranchAddress("Bmumuphi",Bmumuphi);
	ntKp->SetBranchAddress("Bmumupt",Bmumupt);


	ntHi->SetBranchAddress("pthat",&pthat);
	ntHi->SetBranchAddress("weight",&weight);


	
	BDT1->SetBranchAddress("BDT_pt_3_5",BDT_pt_3_5);
	BDT2->SetBranchAddress("BDT_pt_5_7",BDT_pt_5_7);
	BDT3->SetBranchAddress("BDT_pt_7_10",BDT_pt_7_10);
	BDT4->SetBranchAddress("BDT_pt_10_15",BDT_pt_10_15);
	BDT5->SetBranchAddress("BDT_pt_15_20",BDT_pt_15_20);
	BDT6->SetBranchAddress("BDT_pt_20_50",BDT_pt_20_50);
	// BDT7->SetBranchAddress("BDT_pt_50_60",BDT_pt_50_60);
	// BDT8->SetBranchAddress("BDT_pt_1_2",BDT_pt_1_2);
	// BDT9->SetBranchAddress("BDT_pt_50_100",BDT_pt_50_100);

	//BDT->SetBranchAddress("BDT_pt_2_3",BDT_pt_2_3);

	//More Debug Variables//
	Float_t PVx; 
	Float_t PVy; 
	//	Float_t PVz; 

	Float_t PVxE; 
	Float_t PVyE; 
	Float_t PVzE; 

	Float_t BvtxX[NCand]; 
	Float_t BvtxY[NCand]; 
	Float_t BvtxZ[NCand]; 
	Float_t Bdxyz[NCand]; 


	ntKp->SetBranchAddress("PVx",&PVx);
	ntKp->SetBranchAddress("PVy",&PVy);
	//	ntKp->SetBranchAddress("PVz",&PVz);

	ntKp->SetBranchAddress("PVxE",&PVxE);
	ntKp->SetBranchAddress("PVyE",&PVyE);
	ntKp->SetBranchAddress("PVzE",&PVzE);


	ntKp->SetBranchAddress("BvtxX",BvtxX);
	ntKp->SetBranchAddress("BvtxY",BvtxY);
	ntKp->SetBranchAddress("Bdxyz",Bdxyz);


	Float_t  BSx;
	Float_t  BSy;
	Float_t  BSz;

	Float_t  BSxErr;
	Float_t  BSyErr;
	Float_t  BSzErr;

	Float_t  BSdxdz;
	Float_t  BSdydz;


	Float_t  BSdxdzErr;
	Float_t  BSdydzErr;	



	Float_t BSWidthX;
	Float_t BSWidthXErr;


	Float_t BSWidthY;
	Float_t BSWidthYErr;


	ntKp->SetBranchAddress("BSx",&BSx);
	ntKp->SetBranchAddress("BSy",&BSy);
	ntKp->SetBranchAddress("BSz",&BSz);


	ntKp->SetBranchAddress("BSxErr",&BSxErr);
	ntKp->SetBranchAddress("BSyErr",&BSyErr);
	ntKp->SetBranchAddress("BSzErr",&BSzErr);


	ntKp->SetBranchAddress("BSdxdz",&BSdxdz);
	ntKp->SetBranchAddress("BSdydz",&BSdydz);


	ntKp->SetBranchAddress("BSdxdzErr",&BSdxdzErr);
	ntKp->SetBranchAddress("BSdydzErr",&BSdydzErr);

	ntKp->SetBranchAddress("BSWidthX",&BSWidthX);
	ntKp->SetBranchAddress("BSWidthXErr",&BSWidthXErr);

	ntKp->SetBranchAddress("BSWidthY",&BSWidthY);
	ntKp->SetBranchAddress("BSWidthYErr",&BSWidthYErr);



	//DONE Additional Variables//


	TFile * fout;
  TString foutstr;

  if(doMC == 0) 	foutstr = "/data3/smcosta/data/BPData_presel.root";
  if(doMC == 1) 	foutstr = "/data3/smcosta/data/BPMC_presel.root";
  if (filter) {
    if(doMC == 0) 	foutstr = "/data3/smcosta/data/BPData_presel_sort.root";
    if(doMC == 1) 	foutstr = "/data3/smcosta/data/BPMC_presel_sort.root";
  }
  if (jpsi) {
    foutstr = "/data3/smcosta/data/jpsinp_presel.root";
  }
  if (nominal) {
    if(doMC == 0) 	foutstr = "/data3/smcosta/data/BPData_nom.root";
    if(doMC == 1) 	foutstr = "/data3/smcosta/data/BPMC_nom.root";
    if (jpsi) {
      foutstr = "/data3/smcosta/data/jpsinp_nom.root";
    }
  }
  if (noBDT) {
    foutstr = foutstr(0, foutstr.First(".")) + "_noBDT.root";
  }
	fout = new TFile(foutstr,"RECREATE");

	fout->cd();
	TTree* EffInfoTree = new TTree("EffInfoTree","EffInfoTree");	

	Int_t BsizeNew;
	Int_t runNew;
	Int_t lumiNew;
	Int_t evtNew;
	Float_t BmassNew[NCand];
	Float_t BptNew[NCand];
	Float_t ByNew[NCand];
	Float_t BEff[NCand];
	Float_t BEffInv[NCand];
	Float_t BEffInvErr[NCand];

	Float_t BEffInvUp[NCand];
	Float_t BEffInvErrUp[NCand];
	Float_t BEffInvDown[NCand];
	Float_t BEffInvErrDown[NCand];



	Float_t BEffInvBDTWeighted[NCand];
	Float_t BEffInvErrBDTWeighted[NCand];


	Float_t BEff1D[NCand];
	Float_t BEffInv1D[NCand];
	Float_t BEffInvErr1D[NCand];

	Float_t BEffInvFit[NCand];
	Float_t BEffInvErrFit[NCand];

	Float_t BEffInvCent[NCand];
	Float_t BEffInvErrCent[NCand];



	Float_t TnPScale[NCand];
	Float_t TnPErrUp[NCand];
	Float_t TnPErrDown[NCand];

    Float_t BDT_pt_3_5_New[NCand];
	Float_t BDT_pt_5_7_New[NCand];
	Float_t BDT_pt_7_10_New[NCand];
	Float_t BDT_pt_10_15_New[NCand];
	Float_t BDT_pt_15_20_New[NCand];
	Float_t BDT_pt_20_50_New[NCand];
	Float_t BDT_pt_50_60_New[NCand];
	// Float_t BDT_pt_50_100_New[NCand];


	//Additional Variables//

	Float_t Btrk1PtNew[NCand];
	Float_t Btrk2PtNew[NCand];

	Float_t Btrk1PtErrNew[NCand];
	Float_t Btrk2PtErrNew[NCand];


	Float_t Bchi2clNew[NCand];
	Float_t BsvpvDistanceNew[NCand];
	Float_t BsvpvDisErrNew[NCand];
	Float_t BsvpvDistRatioNew[NCand];
	Float_t Btrk1EtaNew[NCand];
	Float_t Btrk2EtaNew[NCand];


	Float_t BsvpvDistance_2DNew[NCand];
	Float_t BsvpvDisErr_2DNew[NCand];


	Float_t BmumumassNew[NCand];
	Float_t Bmu1etaNew[NCand];
	Float_t Bmu1ptNew[NCand];
	Float_t Bmu2etaNew[NCand];
	Float_t Bmu2ptNew[NCand];



	Float_t Bmu1dxyPVNew[NCand];
	Float_t Bmu2dxyPVNew[NCand];

	Float_t Bmu1dzPVNew[NCand];
	Float_t Bmu2dzPVNew[NCand];

	Float_t BtktkmassNew[NCand];


	Float_t Bd0New[NCand];
	Float_t Bd0ErrNew[NCand];
	Float_t BdthetaNew[NCand];
	Float_t BalphaNew[NCand];
	Float_t Btrk1Dz1New[NCand];
	Float_t Btrk2Dz1New[NCand];

	Float_t Btrk1DzError1New[NCand];
	Float_t Btrk2DzError1New[NCand];

	Float_t Btrk1Dxy1New[NCand];
	Float_t Btrk2Dxy1New[NCand];

	Float_t Btrk1DxyError1New[NCand];
	Float_t Btrk2DxyError1New[NCand];


	Float_t BmumuetaNew[NCand];
	Float_t BmumuphiNew[NCand];
	Float_t BmumuptNew[NCand];

	Float_t pthatNew;
	Float_t weightNew;

	Int_t HLT_HIL1DoubleMu0_v1New;

	//More Debug Variables//

	Float_t PVxNew; 
	Float_t PVyNew; 
	Float_t PVzNew; 

	Float_t PVxENew; 
	Float_t PVyENew; 
	Float_t PVzENew; 

	Float_t BvtxXNew[NCand]; 
	Float_t BvtxYNew[NCand]; 
	Float_t BvtxZNew[NCand]; 
	Float_t BdxyzNew[NCand]; 



	EffInfoTree->Branch("PVx",&PVxNew,"PVx/F");
	EffInfoTree->Branch("PVy",&PVyNew,"PVy/F");
	EffInfoTree->Branch("PVz",&PVzNew,"PVz/F");

	EffInfoTree->Branch("PVxE",&PVxENew,"PVxE/F");
	EffInfoTree->Branch("PVyE",&PVyENew,"PVyE/F");
	EffInfoTree->Branch("PVzE",&PVzENew,"PVzE/F");


	EffInfoTree->Branch("BvtxX",BvtxXNew,"BvtxX/F");
	EffInfoTree->Branch("BvtxY",BvtxYNew,"BvtxY/F");
	EffInfoTree->Branch("BvtxZtoPVZ",BvtxZNew,"BvtxZtoPVZ/F");

	Float_t  BSxNew;
	Float_t  BSyNew;
	Float_t  BSzNew;

	Float_t  BSxErrNew;
	Float_t  BSyErrNew;
	Float_t  BSzErrNew;

	Float_t  BSdxdzNew;
	Float_t  BSdydzNew;


	Float_t  BSdxdzErrNew;
	Float_t  BSdydzErrNew;


	Float_t BSWidthXNew;
	Float_t BSWidthXErrNew;


	Float_t BSWidthYNew;
	Float_t BSWidthYErrNew;


	EffInfoTree->Branch("BSx",&BSxNew,"BSx/F");
	EffInfoTree->Branch("BSy",&BSyNew,"BSy/F");
	EffInfoTree->Branch("BSz",&BSzNew,"BSz/F");


	EffInfoTree->Branch("BSxErr",&BSxErrNew,"BSxErr/F");
	EffInfoTree->Branch("BSyErr",&BSyErrNew,"BSyErr/F");
	EffInfoTree->Branch("BSzErr",&BSzErrNew,"BSzErr/F");

	EffInfoTree->Branch("BSdxdz",&BSdxdzNew,"BSdxdz/F");
	EffInfoTree->Branch("BSdydz",&BSdydzNew,"BSdydz/F");

	EffInfoTree->Branch("BSdxdzErr",&BSdxdzErrNew,"BSdxdzErr/F");
	EffInfoTree->Branch("BSdydzErr",&BSdydzErrNew,"BSdydzErr/F");

	EffInfoTree->Branch("BSWidthX",&BSWidthXNew,"BSWidthX/F");
	EffInfoTree->Branch("BSWidthXErr",&BSWidthXErrNew,"BSWidthXErr/F");

	EffInfoTree->Branch("BSWidthY",&BSWidthYNew,"BSWidthY/F");
	EffInfoTree->Branch("BSWidthYErr",&BSWidthYErrNew,"BSWidthYErr/F");


	//DONE Variable//




	Float_t BgenNew[NCand];
	Float_t BgenptNew[NCand];


	int NumCand = 0;
	int BSizeCount;
	int iPass;

	int XLoc;
	int YLoc;

	int XLocTnP;
	int YLocTnP;

	int XLoc1D;
	int XLocCent;
	
	Int_t nMultNew;
	
	EffInfoTree->Branch("nMult",&nMultNew,"nMult/I");


	EffInfoTree->Branch("Bsize",&BsizeNew,"Bsize/I");
	EffInfoTree->Branch("HLT_HIL1DoubleMu0_v1",&HLT_HIL1DoubleMu0_v1New,"HLT_HIL1DoubleMu0_v1/I");



	EffInfoTree->Branch("Bmass",BmassNew,"Bmass/F");
	EffInfoTree->Branch("By",ByNew,"By/F");
	EffInfoTree->Branch("Bpt",BptNew,"Bpt/F");
	EffInfoTree->Branch("Bgen",BgenNew,"Bgen/F");
	EffInfoTree->Branch("Bgenpt",BgenptNew,"Bgenpt/F");


	//Additional Variables//
	EffInfoTree->Branch("Btrk1Pt",Btrk1PtNew,"Btrk1Pt/F");
	EffInfoTree->Branch("Btrk1Eta",Btrk1EtaNew,"Btrk1Eta/F");

	EffInfoTree->Branch("Btktkmass",BtktkmassNew,"Btktkmass/F");


	EffInfoTree->Branch("Btrk1PtErr",Btrk1PtErrNew,"Btrk1PtErr/F");

	EffInfoTree->Branch("Bchi2cl",Bchi2clNew,"Bchi2cl/F");
	EffInfoTree->Branch("BsvpvDistance",BsvpvDistanceNew,"BsvpvDistance/F");
	EffInfoTree->Branch("BsvpvDisErr",BsvpvDisErrNew,"BsvpvDisErr/F");
	EffInfoTree->Branch("BsvpvDistRatio", BsvpvDistRatioNew, "BsvpvDistRatio/F");

	EffInfoTree->Branch("BsvpvDistance_2D",BsvpvDistance_2DNew,"BsvpvDistance_2D/F");
	EffInfoTree->Branch("BsvpvDisErr_2D",BsvpvDisErr_2DNew,"BsvpvDisErr_2D/F");


	EffInfoTree->Branch("Bmumumass",BmumumassNew,"Bmumumass/F");
	EffInfoTree->Branch("Bmu1eta",Bmu1etaNew,"Bmu1eta/F");
	EffInfoTree->Branch("Bmu2eta",Bmu2etaNew,"Bmu2eta/F");
	EffInfoTree->Branch("Bmu1pt",Bmu1ptNew,"Bmu1pt/F");
	EffInfoTree->Branch("Bmu2pt",Bmu2ptNew,"Bmu2pt/F");

	EffInfoTree->Branch("Bmu1dxyPV",Bmu1dxyPVNew,"Bmu1dxyPV/F");
	EffInfoTree->Branch("Bmu2dxyPV",Bmu2dxyPVNew,"Bmu2dxyPV/F");
	EffInfoTree->Branch("Bmu1dzPV",Bmu1dzPVNew,"Bmu1dzPV/F");
	EffInfoTree->Branch("Bmu2dzPV",Bmu2dzPVNew,"Bmu2dzPV/F");




	EffInfoTree->Branch("Bd0",Bd0New,"Bd0/F");
	EffInfoTree->Branch("Bd0Err",Bd0ErrNew,"Bd0Err/F");
	EffInfoTree->Branch("Bdtheta",BdthetaNew,"Bdtheta/F");
	EffInfoTree->Branch("Balpha",BalphaNew,"Balpha/F");
	EffInfoTree->Branch("Btrk1Dz1",Btrk1Dz1New,"Btrk1Dz1/F");
	EffInfoTree->Branch("Btrk1DzError1",Btrk1DzError1New,"Btrk1DzError1/F");
	EffInfoTree->Branch("Btrk1Dxy1",Btrk1Dxy1New,"Btrk1Dxy1/F");
	EffInfoTree->Branch("Btrk1DxyError1",Btrk1DxyError1New,"Btrk1DxyError1/F");
	EffInfoTree->Branch("Bmumueta",BmumuetaNew,"Bmumueta/F");
	EffInfoTree->Branch("Bmumuphi",BmumuphiNew,"Bmumuphi/F");
	EffInfoTree->Branch("Bmumupt",BmumuptNew,"Bmumupt/F");

	EffInfoTree->Branch("pthat",&pthatNew,"pthat/F");
	EffInfoTree->Branch("weight",&weightNew,"weight/F");



	Bool_t Bmu1isTriggered[NCand];
	Bool_t Bmu2isTriggered[NCand];
	ntKp->SetBranchAddress("Bmu1isTriggered",Bmu1isTriggered);
	ntKp->SetBranchAddress("Bmu2isTriggered",Bmu2isTriggered);

	
	   EffInfoTree->Branch("BDT_pt_3_5",BDT_pt_3_5_New,"BDT_pt_3_5/F");
	   EffInfoTree->Branch("BDT_pt_5_7",BDT_pt_5_7_New,"BDT_pt_5_7/F");
	   EffInfoTree->Branch("BDT_pt_7_10",BDT_pt_7_10_New,"BDT_pt_7_10/F");
	   EffInfoTree->Branch("BDT_pt_10_15",BDT_pt_10_15_New,"BDT_pt_10_15/F");
	   EffInfoTree->Branch("BDT_pt_15_20",BDT_pt_15_20_New,"BDT_pt_15_20/F");
	   EffInfoTree->Branch("BDT_pt_20_50",BDT_pt_20_50_New,"BDT_pt_20_50/F");
	   EffInfoTree->Branch("BDT_pt_50_60",BDT_pt_50_60_New,"BDT_pt_50_60/F");
	
	   // EffInfoTree->Branch("BDT_pt_50_100",BDT_pt_50_100_New,"BDT_pt_50_100/F");

	   // Float_t BDT_pt_1_2_New[NCand];


	   // EffInfoTree->Branch("BDT_pt_1_2",BDT_pt_1_2_New,"BDT_pt_1_2/F");


	//Float_t BDT_pt_2_3_New[NCand];

	//EffInfoTree->Branch("BDT_pt_2_3",BDT_pt_2_3_New,"BDT_pt_2_3/F");


	



	//Done Additional Variables//



	//More//


	//Additional Variables for TMVA

	Int_t Bmu1InPixelLayerNew[NCand]; 
	Int_t Bmu2InPixelLayerNew[NCand]; 
	Int_t Bmu1InStripLayerNew[NCand]; 
	Int_t Bmu2InStripLayerNew[NCand]; 
	Int_t Bmu1isTrackerMuonNew[NCand]; 
	Int_t Bmu2isTrackerMuonNew[NCand]; 
	Int_t Bmu1isGlobalMuonNew[NCand]; 
	Int_t Bmu2isGlobalMuonNew[NCand]; 
	Int_t Btrk1highPurityNew[NCand]; 
	Int_t Btrk2highPurityNew[NCand]; 
	Int_t Btrk1PixelHitNew[NCand];
	Int_t Btrk2PixelHitNew[NCand];
	Int_t Btrk1StripHitNew[NCand];
	Int_t Btrk2StripHitNew[NCand];
	Float_t Btrk1Chi2ndfNew[NCand];
	Float_t Btrk2Chi2ndfNew[NCand];
	Int_t Btrk1nStripLayerNew[NCand];
	Int_t Btrk2nStripLayerNew[NCand];
	Int_t Btrk1nPixelLayerNew[NCand];
	Int_t Btrk2nPixelLayerNew[NCand];

	Int_t Bmu1isTriggeredNew[NCand];
	Int_t Bmu2isTriggeredNew[NCand];


	EffInfoTree->Branch("Bmu1isTriggered",Bmu1isTriggeredNew,"Bmu1isTriggered/I");
	EffInfoTree->Branch("Bmu2isTriggered",Bmu2isTriggeredNew,"Bmu2isTriggered/I");


	Int_t HBHENoiseFilterResultNew;
	Int_t pPAprimaryVertexFilterNew;
	Int_t pBeamScrapingFilterNew;
//	Int_t HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1New;

	EffInfoTree->Branch("Bmu1InPixelLayer",Bmu1InPixelLayerNew,"Bmu1InPixelLayer/I");
	EffInfoTree->Branch("Bmu2InPixelLayer",Bmu2InPixelLayerNew,"Bmu2InPixelLayer/I");
	EffInfoTree->Branch("Bmu1InStripLayer",Bmu1InStripLayerNew,"Bmu1InStripLayer/I");
	EffInfoTree->Branch("Bmu2InStripLayer",Bmu2InStripLayerNew,"Bmu2InStripLayer/I");

	EffInfoTree->Branch("Bmu1isTrackerMuon",Bmu1isTrackerMuonNew,"Bmu1isTrackerMuon/I");
	EffInfoTree->Branch("Bmu2isTrackerMuon",Bmu2isTrackerMuonNew,"Bmu2isTrackerMuon/I");

	EffInfoTree->Branch("Bmu1isGlobalMuon",Bmu1isGlobalMuonNew,"Bmu1isGlobalMuon/I");
	EffInfoTree->Branch("Bmu2isGlobalMuon",Bmu2isGlobalMuonNew,"Bmu2isGlobalMuon/I");

	EffInfoTree->Branch("Btrk1highPurity",Btrk1highPurityNew,"Btrk1highPurity/I");
	EffInfoTree->Branch("Btrk2highPurity",Btrk2highPurityNew,"Btrk2highPurity/I");

	EffInfoTree->Branch("HBHENoiseFilterResult",&HBHENoiseFilterResultNew,"HBHENoiseFilterResult/I");
	EffInfoTree->Branch("pPAprimaryVertexFilter",&pPAprimaryVertexFilterNew,"pPAprimaryVertexFilter/I");
	EffInfoTree->Branch("pBeamScrapingFilter",&pBeamScrapingFilterNew,"pBeamScrapingFilter/I");
//	EffInfoTree->Branch("HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1",&HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1New,"HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1/I");


	EffInfoTree->Branch("Btrk1PixelHit",Btrk1PixelHitNew,"Btrk1PixelHit/I");
	EffInfoTree->Branch("Btrk2PixelHit",Btrk2PixelHitNew,"Btrk2PixelHit/I");

	EffInfoTree->Branch("Btrk1StripHit",Btrk1StripHitNew,"Btrk1StripHit/I");
	EffInfoTree->Branch("Btrk2StripHit",Btrk2StripHitNew,"Btrk2StripHit/I");

	EffInfoTree->Branch("Btrk1Chi2ndf",Btrk1Chi2ndfNew,"Btrk1Chi2ndf/F");
	EffInfoTree->Branch("Btrk2Chi2ndf",Btrk2Chi2ndfNew,"Btrk2Chi2ndf/F");


	EffInfoTree->Branch("Btrk1nStripLayer",Btrk1nStripLayerNew,"Btrk1nStripLayer/I");
	EffInfoTree->Branch("Btrk2nStripLayer",Btrk2nStripLayerNew,"Btrk2nStripLayer/I");

	EffInfoTree->Branch("Btrk1nPixelLayer",Btrk1nPixelLayerNew,"Btrk1nPixelLayer/I");
	EffInfoTree->Branch("Btrk2nPixelLayer",Btrk2nPixelLayerNew,"Btrk2nPixelLayer/I");



	Int_t MuonAcc[NCand];
	EffInfoTree->Branch("MuonAcc",MuonAcc,"MuonAcc/I");



	int NEvents = ntKp->GetEntries();


	int LocWeightX;


	int LocWeightY;

	double BDTSystemWeight;

  int trackSelection = -1;
	EffInfoTree->Branch("track", &trackSelection,"track/I");

	//Copy a new tree//

	cout << "Copying Fit Tree for Framework" << endl;  

	TTree* EffInfoTree_NewFit;
	EffInfoTree_NewFit = EffInfoTree->CloneTree(0);
	EffInfoTree_NewFit->SetObject("ntKp","ntKp");

	//NEvents = 1000000;

	//Done Copying//

  if (!jpsi) {
    root->SetBranchAddress("EvtInfo.EvtNo", &EvtNo);
  }

	for(int i = 0; i < NEvents; i++){

		if(i%100000==0) std::cout<<std::setiosflags(std::ios::left)<<"  [ \033[1;36m"<<std::setw(10)<<i<<"\033[0m"<<" / "<<std::setw(10)<<ntKp->GetEntries()<<" ] "<<"\033[1;36m"<<Form("%.0f",100.*i/ntKp->GetEntries())<<"%\033[0m"<<"\r"<<std::flush;


		ntKp->GetEntry(i);
		ntSkim->GetEntry(i);
		ntHi->GetEntry(i);
		ntHlt->GetEntry(i);

		//BDT//
		
		   // BDT1->GetEntry(i);
		   BDT2->GetEntry(i);
		   BDT3->GetEntry(i);
		   BDT4->GetEntry(i);
		   BDT5->GetEntry(i);
		   BDT6->GetEntry(i);
		   // BDT7->GetEntry(i);
		   // BDT8->GetEntry(i);
		   // BDT9->GetEntry(i);
		   
	//	BDT->GetEntry(i);


       if (!jpsi) {
         root->GetEntry(i);
       }

    // if (filter && !evtList[uevt]) {
    //   continue;
    // } else if (clumi != lumiList[uevt]) {
    //   cout << "discarding evt " << EvtNo << "\n";
    //   cout << "nom lumi:" << (unsigned) lumiList[uevt] <<
    //     ", extra lumi:" << (unsigned) clumi <<
    //     ", int lumi:" << lumi <<
    //     "\n";
    //   continue;
    // }
	/*
    if (filter) {
      // unsigned uevt = EvtNo;
      UShort_t clumi = (UShort_t) lumi;
      auto thispair = std::make_pair(EvtNo, clumi);
      if (!std::binary_search(filterList.begin(), filterList.end(), thispair)) {
          cout << "discarding extra evt " << EvtNo << "\n";
          cout << "lumi: " << lumi << "\n";
        continue;
      }
    }
	*/

		BSizeCount = 0;
		iPass = 0;

		runNew = run;
		lumiNew = lumi;
		evtNew = evt;

		double lowstBDTCut = 0.32;

    // std::array<double, 5> bdtCut = (nominal)? nomBDT : preBDT;
    std::array<double, 6> bdtCut = nomBDT;

		for(int j =0; j < Bsize; j++){

      bool passBDT = ( (Bpt[j] > 3 && Bpt[j] < 5 && BDT_pt_3_5[j] > bdtCut[0]) ||
					   (Bpt[j] > 5 && Bpt[j] < 7 && BDT_pt_5_7[j] > bdtCut[1]) ||
                       (Bpt[j] > 7 && Bpt[j] < 10 && BDT_pt_7_10[j] > bdtCut[2]) ||
                       (Bpt[j] > 10 && Bpt[j] < 15 && BDT_pt_10_15[j] > bdtCut[3]) ||
                       (Bpt[j] > 15 && Bpt[j] < 20 && BDT_pt_15_20[j] > bdtCut[4])  ||
                       (Bpt[j] > 20 && Bpt[j] < 50 && BDT_pt_20_50[j] > bdtCut[5]) ||
                       (Bpt[j] > 50 && Bpt[j] < 60));

      bool passPre = BsvpvDistance[j]/BsvpvDisErr[j] > 2.0 & Bchi2cl[j] > 0.05 &&
        TMath::Abs(Bmumumass[j]-3.096916)<0.15;
      if ((pPAprimaryVertexFilter == 1 && pBeamScrapingFilter == 1 && HLT_HIL1DoubleMu0_v1 == 1)  &&  (Bmu1isTriggered[j] == 1 && Bmu2isTriggered[j] == 1 ) && (Btrk1Pt[j] > 0.5 && Bpt[j] > 2 && abs(Btrk1Eta[j]-0.0) < 2.4  && (TMath::Abs(By[j])<2.4&&((abs(Bmu1eta[j])<1.2&&Bmu1pt[j]>3.5)||(abs(Bmu1eta[j])>1.2&&abs(Bmu1eta[j])<2.1&&Bmu1pt[j]>(5.47-1.89*abs(Bmu1eta[j])))||(abs(Bmu1eta[j])>2.1&&abs(Bmu1eta[j])<2.4&&Bmu1pt[j]>1.5))&&((abs(Bmu2eta[j])<1.2&&Bmu2pt[j]>3.5)||(abs(Bmu2eta[j])>1.2&&abs(Bmu2eta[j])<2.1&&Bmu2pt[j]>(5.47-1.89*abs(Bmu2eta[j])))||(abs(Bmu2eta[j])>2.1&&abs(Bmu2eta[j])<2.4&&Bmu2pt[j]>1.5))&&Bmu1TMOneStationTight[j]&&Bmu2TMOneStationTight[j]&&Bmu1InPixelLayer[j]>0&&(Bmu1InPixelLayer[j]+Bmu1InStripLayer[j])>5&&Bmu2InPixelLayer[j]>0&&(Bmu2InPixelLayer[j]+Bmu2InStripLayer[j])>5&&Bmu1dxyPV[j]<0.3&&Bmu2dxyPV[j]<0.3&&Bmu1dzPV[j]<20&&Bmu2dzPV[j]<20&&Bmu1isTrackerMuon[j]&&Bmu2isTrackerMuon[j]&&Bmu1isGlobalMuon[j]&&Bmu2isGlobalMuon[j]&&Btrk1highPurity[j]&&abs(Btrk1Eta[j])<2.4&&Btrk1Pt[j]>0.5)  && (Btrk1PixelHit[j] + Btrk1StripHit[j] > 10) &&  (Btrk1PtErr[j]/Btrk1Pt[j] < ptErr[Tracking::loose]) && Btrk1Chi2ndf[j]/(Btrk1nStripLayer[j]+Btrk1nPixelLayer[j]) < chi2Nlayer[Tracking::loose]   && (abs(PVz)<15))){

        if (nominal && !passPre) {
          continue;
        }
        if (!nominal && Bchi2cl[j] <= 0.003) {
          continue;
        }
        // if (!noBDT && !passBDT && !passPre) {
        if (!noBDT && !passBDT) {
          continue;
        }

        if (Bpt[j] < 5 && Bpt[j] > 60) {
          continue;
        }

				if(doMC==1 && !jpsi){

					if(Bgen[j] !=23333) continue;
				}

        // tracking selection flag
        auto passTracking = [=] (Tracking cut) {
          return (Btrk1PtErr[j]/Btrk1Pt[j] < ptErr[cut] &&
                  Btrk1Chi2ndf[j]/(Btrk1nStripLayer[j]+Btrk1nPixelLayer[j])
                  < chi2Nlayer[cut]);
        };

        if (passTracking(Tracking::tight)) {
          trackSelection = Tracking::tight;
        } else if (passTracking(Tracking::standard)) {
          trackSelection = Tracking::standard;
        } else {
          trackSelection = Tracking::loose;
        }


				nMultNew = nMult;

				MuonAcc[iPass] = 0;

				if(((abs(Bmu1eta[j])<1.2 && Bmu1pt[j]>3.5) || (abs(Bmu1eta[j])>1.2 && abs(Bmu1eta[j])<2.1 && Bmu1pt[j]>(5.77-1.8*abs(Bmu1eta[j]))) || (abs(Bmu1eta[j])>2.1 && abs(Bmu1eta[j])<2.4 && Bmu1pt[j]>1.8)) && ((abs(Bmu2eta[j])<1.2 && Bmu2pt[j]>3.5) || (abs(Bmu2eta[j])>1.2 && abs(Bmu2eta[j])<2.1 && Bmu2pt[j]>(5.77-1.8*abs(Bmu2eta[j]))) || (abs(Bmu2eta[j])>2.1 && abs(Bmu2eta[j])<2.4 && Bmu2pt[j]>1.8)) && Bmu1TMOneStationTight[j] && Bmu2TMOneStationTight[j] && Bmu1InPixelLayer[j] > 0 && (Bmu1InPixelLayer[j]+Bmu1InStripLayer[j]) > 5 && Bmu2InPixelLayer[j] > 0 && (Bmu2InPixelLayer[j]+Bmu2InStripLayer[j]) > 5 && Bmu1dxyPV[j]< 0.3 && Bmu2dxyPV[j]< 0.3 && Bmu1dzPV[j]<20 && Bmu2dzPV[j]<20 && Bmu1isGlobalMuon[j] && Bmu2isGlobalMuon[j]) MuonAcc[iPass] = 1;

				BmassNew[iPass] =  Bmass[j];
				BptNew[iPass] =  Bpt[j];
				ByNew[iPass] =  By[j];
				BgenNew[iPass] = Bgen[j];
				BgenptNew[iPass] = Bgenpt[j];



				HLT_HIL1DoubleMu0_v1New = HLT_HIL1DoubleMu0_v1;
				HBHENoiseFilterResultNew = HBHENoiseFilterResult;
				pPAprimaryVertexFilterNew =  pPAprimaryVertexFilter;
				pBeamScrapingFilterNew = pBeamScrapingFilter;
				
			    BDT_pt_3_5_New[iPass] = BDT_pt_3_5[j];
				BDT_pt_5_7_New[iPass] = BDT_pt_5_7[j];
				BDT_pt_7_10_New[iPass] = BDT_pt_7_10[j];
				BDT_pt_10_15_New[iPass] = BDT_pt_10_15[j];
				BDT_pt_15_20_New[iPass] = BDT_pt_15_20[j];
				BDT_pt_20_50_New[iPass] = BDT_pt_20_50[j];
				BDT_pt_50_60_New[iPass] = BDT_pt_50_60[j];
				// BDT_pt_2_3_New[iPass] = BDT_pt_2_3[j];
				// BDT_pt_50_100_New[iPass] = BDT_pt_50_100[j];



				//Additional Variables//

				Btrk1PtNew[iPass] = Btrk1Pt[j];
				Btrk1EtaNew[iPass] = Btrk1Eta[j];

				Btrk1PtErrNew[iPass] = Btrk1PtErr[j];

				Bchi2clNew[iPass] = Bchi2cl[j];
				BsvpvDistanceNew[iPass] = BsvpvDistance[j];
				BsvpvDisErrNew[iPass] = BsvpvDisErr[j];
				BsvpvDistRatioNew[iPass] = BsvpvDistance[j] / BsvpvDisErr[j];
				BmumumassNew[iPass] = Bmumumass[j];
				Bmu1etaNew[iPass] = Bmu1eta[j];
				Bmu2etaNew[iPass] = Bmu2eta[j];
				Bmu1ptNew[iPass] = Bmu1pt[j];
				Bmu2ptNew[iPass] = Bmu2pt[j];

				Bmu1dxyPVNew[iPass] = Bmu1dxyPV[j];
				Bmu2dxyPVNew[iPass] = Bmu2dxyPV[j];
				Bmu1dzPVNew[iPass] = Bmu1dzPV[j];
				Bmu2dzPVNew[iPass] = Bmu2dzPV[j];


				BsvpvDistance_2DNew[iPass] = BsvpvDistance_2D[j];
				BsvpvDisErr_2DNew[iPass] = BsvpvDisErr_2D[j];


				Bd0New[iPass] = Bd0[j];
				Bd0ErrNew[iPass] = Bd0Err[j];
				BdthetaNew[iPass] = Bdtheta[j];
				BalphaNew[iPass] =  Balpha[j];
				Btrk1Dz1New[iPass] = Btrk1Dz1[j];
				Btrk1DzError1[iPass] = Btrk1DzError1[j];


				Btrk1Dxy1New[iPass] = Btrk1Dxy1[j];
				Btrk1DxyError1New[iPass] = Btrk1DxyError1[j];


				BmumuetaNew[iPass] = Bmumueta[j];
				BmumuptNew[iPass] = Bmumupt[j];
				BmumuphiNew[iPass] = Bmumuphi[j];


				pthatNew = pthat;
				weightNew = weight;
				//DONE Additional Variables

				PVxNew = PVx;
				PVyNew = PVy;
				PVzNew = PVz;


				PVxENew = PVxE;
				PVyENew = PVyE;
				PVzENew = PVzE;

				BvtxXNew[iPass] = BvtxX[j];
				BvtxYNew[iPass] = BvtxY[j];
				BvtxZNew[iPass] = sqrt(Bdxyz[j] * Bdxyz[j] - (BvtxX[j] - PVx) * (BvtxX[j] - PVx) - (BvtxY[j] - PVy) * (BvtxY[j] - PVy));


				BSxNew = BSx;
				BSyNew = BSy;
				BSzNew = BSz;

				BSxErrNew = BSxErr;
				BSyErrNew = BSyErr;
				BSzErrNew = BSzErr;


				BSdxdzNew = BSdxdz;
				BSdydzNew = BSdydz;



				BSdxdzErrNew = BSdxdzErr;
				BSdydzErrNew = BSdydzErr;


				BSWidthXNew = BSWidthX;
				BSWidthXErrNew = BSWidthXErr;


				BSWidthYNew = BSWidthY;
				BSWidthYErrNew = BSWidthYErr;

				//More Muon Info



				Bmu1InPixelLayerNew[iPass] = Bmu1InPixelLayer[j]; 
				Bmu2InPixelLayerNew[iPass] = Bmu2InPixelLayer[j]; 
				Bmu1InStripLayerNew[iPass] = Bmu1InStripLayer[j]; 
				Bmu2InStripLayerNew[iPass] = Bmu2InStripLayer[j]; 
				Bmu1isTrackerMuonNew[iPass] = Bmu1isTrackerMuon[j]; 
				Bmu2isTrackerMuonNew[iPass] = Bmu2isTrackerMuon[j]; 
				Bmu1isGlobalMuonNew[iPass] = Bmu1isGlobalMuon[j]; 
				Bmu2isGlobalMuonNew[iPass] = Bmu2isGlobalMuon[j]; 
				Btrk1highPurityNew[iPass] = Btrk1highPurity[j]; 
				Btrk2highPurityNew[iPass] = Btrk2highPurity[j]; 
				Btrk1PixelHitNew[iPass] = Btrk1PixelHit[j];
				Btrk2PixelHitNew[iPass] = Btrk2PixelHit[j];
				Btrk1StripHitNew[iPass] = Btrk1StripHit[j];
				Btrk2StripHitNew[iPass] = Btrk2StripHit[j];

				Btrk1Chi2ndfNew[iPass] = Btrk1Chi2ndf[j];
				Btrk2Chi2ndfNew[iPass] = Btrk2Chi2ndf[j];
				Btrk1nStripLayerNew[iPass] = Btrk1nStripLayer[j];
				Btrk2nStripLayerNew[iPass] = Btrk2nStripLayer[j];
				Btrk1nPixelLayerNew[iPass] = Btrk1nPixelLayer[j];
				Btrk2nPixelLayerNew[iPass] = Btrk2nPixelLayer[j];

	
				BtktkmassNew[iPass] = Btktkmass[j];


				//More Debug Variables DONE//



				NumCand = NumCand + 1;
				BSizeCount = BSizeCount + 1;
			//	iPass = iPass + 1;
				BsizeNew = BSizeCount;
				if(BsizeNew > 0) EffInfoTree->Fill();
				if(BsizeNew > 0) EffInfoTree_NewFit->Fill();
			}



			}



		}

		fout->Write();
		fout->Close();


		}



