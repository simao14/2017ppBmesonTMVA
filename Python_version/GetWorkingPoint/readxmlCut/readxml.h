#ifndef _READXML_H_
#define _READXML_H_

#include "../includes/uti.h"
#include "../includes/TMVAClassification.h"
#include "../includes/cfout.h"

const int NmaxFonll = 1198; //fonll data points number
float fcentral[NmaxFonll],fpt[NmaxFonll];
const int NEff = 100;
const int nMVA = 200;
Float_t minMVA = -1;
Float_t maxMVA = 1;
Float_t widMVA = (maxMVA-minMVA)/nMVA;
Float_t gmvaBins[nMVA+1];

Double_t norm = 0.401; // FF of c->D0, i.e., D fraction
Double_t BR = 6.02061e-5;
Double_t deltapt = 0.25;

Float_t massD = 5.27934;
Float_t dmassDsignal = 0.08;
Float_t dmassDsidbandL = 0.20;
Float_t dmassDsidbandH = 0.25;
Float_t massDsignalL = massD - dmassDsignal;
Float_t massDsignalH = massD + dmassDsignal;

Float_t massDsidebandRH = massD + dmassDsidbandH;
Float_t massDsidebandRL = massD + dmassDsidbandL;
Float_t massDsidebandLH = massD - dmassDsidbandL;
Float_t massDsidebandLL = massD - dmassDsidbandH;

TString outputresult = "results/fresult";
TString outputfonll = "plots/cfonll";
TString outputmvadis = "plots/cmvadistribution";
TString outputefficiency = "plots/cefficiency";
TString outputsignificance = "plots/csignificance";




int calRatio(TTree* signal, TTree* background, TTree* generated, 
		TString cuts, TString cutb, TString cutg, TString weightfunctiongen, TString weightfunctionreco,
		Float_t ptmin, Float_t ptmax,
		TString fonllfile, Float_t lumi, Float_t raa, TString outputname,
		Float_t* results)
{
	// background sideband -> signal region
	signal->AddFriend("hiEvtAnalyzer/HiTree");
	generated->AddFriend("hiEvtAnalyzer/HiTree");
	
	cout << "pass  6 "<< endl;


	TH1D* hMassD = new TH1D("hMassD",";B mass (GeV/c^{2});Background Entries",50,5.0,6.0);
	background->Project("hMassD","Bmass",cutb);
//	Float_t nB = hMassD->Integral()*dmassDsignal/(dmassDsidbandH-dmassDsidbandL);

        int massDsidebandRHBin = hMassD->GetXaxis()->FindBin(massDsidebandRH+0.00001);
        int massDsidebandRLBin = hMassD->GetXaxis()->FindBin(massDsidebandRL -0.00001 );


        int massDsidebandLHBin = hMassD->GetXaxis()->FindBin(massDsidebandLH+ 0.00001);
        int massDsidebandLLBin = hMassD->GetXaxis()->FindBin(massDsidebandLL -0.00001);

        Float_t nB = (hMassD->Integral(massDsidebandRLBin,massDsidebandRHBin) + hMassD->Integral(massDsidebandLLBin,massDsidebandLHBin))/((massDsidebandRH - massDsidebandRL) + (massDsidebandLH - massDsidebandLL))*(massDsignalH- massDsignalL);

	cout << "weightfunctiongen = " << weightfunctiongen.Data() << endl;
	// fonll predictions
	std::ifstream getdata(fonllfile);
	cout << "FONLL = " << fonllfile.Data() << endl;
	if(!getdata.is_open()) 
	{
		std::cout<<"Opening fonll file fails"<<std::endl;
		return 1;
	}
	cout << "ptmin = " << ptmin << "  ptmax = " << ptmax  << endl;
	Float_t tem;
	Int_t nbin=0;
	while (!getdata.eof())
	{
	//	getdata>>fpt[nbin]>>fcentral[nbin]>>tem>>tem>>tem>>tem>>tem>>tem>>tem>>tem;
		getdata>>fpt[nbin]>>tem>>tem>>fcentral[nbin]>>tem>>tem>>tem>>tem>>tem>>tem; //use max FONLL as central
		if(fpt[nbin]>=ptmin&&fpt[nbin]<=ptmax) nbin++;
		//  cout << "nbin = " << nbin << " pT = " << fpt[nbin-1] << "  x-section = " << fcentral[nbin-1] << endl;
	}

	cout << "pass  7 "<< endl;



	TH1D* hfonll = new TH1D("hfonll",";B_{s} p_{T} (GeV/c);FONLL differential xsection",nbin-1,fpt);
	for(int i=0;i<nbin;i++)
	{
		hfonll->SetBinContent(i,fcentral[i]);
	}
	TCanvas* cfonll = new TCanvas("cfonll","",600,600);
	hfonll->Draw();
	cfonll->SaveAs(Form("%s_%s.pdf",outputfonll.Data(),outputname.Data()));
	cout << "nbin = " << nbin << endl;
	hMassD->Draw();
	cfonll->SaveAs("plots/bkgdmass.png");

	TH1D* hrec = new TH1D("hrec",";B_{s} p_{T} (GeV/c);Signal reco entries",nbin-1,fpt);
	TH1D* hgen = new TH1D("hgen",";B_{s} p_{T} (GeV/c);Generated entries",nbin-1,fpt);
	TH1D* heff = new TH1D("heff",";B_{s} p_{T} (GeV/c);Prefilter efficiency",nbin-1,fpt);
	signal->Project("hrec","Bpt",Form("%s*(%s)",weightfunctionreco.Data(),cuts.Data()));
	generated->Project("hgen","Gpt",Form("%s*(%s)",weightfunctiongen.Data(),cutg.Data()));
	cout << "Gen Cut = " << Form("%s*(%s)",weightfunctiongen.Data(),cutg.Data()) << endl;
	cout << "Total Signal = " << hrec->Integral() << endl;
	cout << "Total Generated = " << hgen->Integral() << endl;

	heff->Divide(hrec,hgen,1.,1.,"B");
	cout << "norm = " <<  norm << endl; 
	heff->Draw();
	cfonll->SaveAs(Form("plots/Efficiency_%d_%d.png",ptmin,ptmax));
	TH1D* htheoryreco = new TH1D("htheoryreco","",nbin-1,fpt);
	htheoryreco->Multiply(heff,hfonll,1,1,"B");
	Float_t nS = htheoryreco->Integral()*norm*BR*deltapt*lumi*raa*2; // x2
//	Float_t nS = hfonll->Integral()*hrec->Integral()/hgen->Integral()*norm*BR*deltapt*lumi*raa*2;
	results[0] = nB;
	results[1] = nS;
	
cout << "htheoryreco->Integral() = " << htheoryreco->Integral() << endl;

cout << "nS = " <<nS << endl; 
	 cout << "nB = " <<nB << endl;
	return 0;
}

#endif

