#include "uti.h"

//#include "TMVA_BDT_PbPb_pt_5p0_10p0.class.C"
#include "readxml/weights/TMVAClassification_BDTs.class.C"

#define MAX_XB       20000
using namespace std;
constexpr bool debug = false;

void BDT(TString infname, TString ofname, int ptmin, int ptmax)
{
	TFile* inf = new TFile(infname);
	TTree* t = (TTree*)inf->Get("Bfinder/ntKp");
	//TTree* t = (TTree*)inf->Get("ntphi");
	//	t->AddFriend("mvaTree");

	Int_t           Bsize;
	Int_t           RunNo;
	Int_t           EvtNo;
	Int_t           LumiNo;
	Float_t         Bchi2cl[MAX_XB];   //[Bsize]
	Float_t         BsvpvDistance[MAX_XB];   //[Bsize]
	Float_t         BsvpvDisErr[MAX_XB];   //[Bsize]
	Float_t         Balpha[MAX_XB];   //[Bsize]   
	Float_t         Btrk1Pt[MAX_XB];   //[Bsize]   
	Float_t         Btrk2Pt[MAX_XB];   //[Bsize]   
	Float_t         Btrk1Eta[MAX_XB];   //[Bsize]   
	Float_t         Btrk2Eta[MAX_XB];   //[Bsize]   
	Float_t         Btrk1Dxy1[MAX_XB];   //[Bsize]   
	Float_t         Btrk2Dxy1[MAX_XB];   //[Bsize]   
	Float_t         Btrk1DxyError1[MAX_XB];   //[Bsize]   
	Float_t         Btrk2DxyError1[MAX_XB];   //[Bsize]   
	Float_t         Btrk1Dz1[MAX_XB];   //[Bsize]   
	Float_t         Btrk2Dz1[MAX_XB];   //[Bsize]   
	Float_t         Btrk1DzError1[MAX_XB];   //[Bsize]   
	Float_t         Btrk2DzError1[MAX_XB];   //[Bsize]   

	Float_t         Btrk1D0Err[MAX_XB];   //[Bsize]   
	Float_t         Btrk2D0Err[MAX_XB];   //[Bsize]   
	Float_t         Btktkmass[MAX_XB]; 
	Float_t         Bd0[MAX_XB]; 
	Float_t         Bd0Err[MAX_XB]; 
	Float_t         Bdtheta[MAX_XB]; 
	Float_t         Bpt[MAX_XB]; 

	t->SetBranchAddress("RunNo", &RunNo);
	t->SetBranchAddress("EvtNo", &EvtNo);
	t->SetBranchAddress("LumiNo", &LumiNo);


	t->SetBranchAddress("Bsize", &Bsize);
	t->SetBranchAddress("Balpha", Balpha);
	t->SetBranchAddress("Bchi2cl", Bchi2cl);
	t->SetBranchAddress("BsvpvDistance", BsvpvDistance);
	t->SetBranchAddress("BsvpvDisErr", BsvpvDisErr);
	//  t->SetBranchAddress("MVA", MVA);

	t->SetBranchAddress("Bd0Err", Bd0Err);
	t->SetBranchAddress("Btrk1Pt", Btrk1Pt);
	t->SetBranchAddress("Bd0", Bd0);
	//t->SetBranchAddress("Bd0Err", Bd0Err);
	t->SetBranchAddress("Bdtheta", Bdtheta);
	t->SetBranchAddress("Btktkmass", Btktkmass);
	t->SetBranchAddress("Btrk1Dxy1", Btrk1Dxy1);
	t->SetBranchAddress("Btrk1Eta", Btrk1Eta);

	t->SetBranchAddress("Btrk1DxyError1", Btrk1DxyError1);
	t->SetBranchAddress("Bpt", Bpt);
	t->SetBranchAddress("Btrk1Dz1", Btrk1Dz1);
	t->SetBranchAddress("Btrk1DzError1", Btrk1DzError1);



	std::vector<std::string> theInputVars;

	std::string a1 = "Btrk1Pt";
	std::string a2 = "abs(Btrk1Dz1/Btrk1DzError1)";
	std::string a3 = "abs(Btrk1Dxy1/Btrk1DxyError1)";
	std::string a4 = "BsvpvDistance/BsvpvDisErr";
	std::string a5 = "Balpha";
	std::string a6 = "Bchi2cl";

	theInputVars.push_back(a1);
	theInputVars.push_back(a2);
	theInputVars.push_back(a3);
	theInputVars.push_back(a4);
	theInputVars.push_back(a5);
	theInputVars.push_back(a6);



	std::vector<double> inputValues;
	ReadBDTs mva(theInputVars);
	int run;
	int evt;
	int lumi;



	TFile* outf = new TFile(ofname,"recreate");
	outf->cd();
	TTree* mvaTree = new TTree(Form("BDT_pt_%d_%d",ptmin,ptmax),"BDT");	

	double BDT[MAX_XB];
	mvaTree->Branch("Bsize",&Bsize,"Bsize/I");
	mvaTree->Branch(Form("BDT_pt_%d_%d",ptmin,ptmax),BDT,Form("BDT_pt_%d_%d[Bsize]/D",ptmin,ptmax));
	mvaTree->Branch("run",&run,"run/I");
	mvaTree->Branch("evt",&evt,"evt/I");
	mvaTree->Branch("lumi",&lumi,"lumi/I");
	std::cout<<std::endl;
	std::cout<<"  Input file: "<<infname<<std::endl;
	std::cout<<"  Calculating MVA values ..."<<std::endl;
	for(int i=0;i< t->GetEntries();i++)
	{
	
		t->GetEntry(i);
		run = RunNo;
		evt = EvtNo;
		lumi = LumiNo;

		if(i%1000==0) std::cout<<std::setiosflags(std::ios::left)<<"  [ \033[1;36m"<<std::setw(10)<<i<<"\033[0m"<<" / "<<std::setw(10)<<t->GetEntries()<<" ] "<<"\033[1;36m"<<Form("%.0f",100.*i/t->GetEntries())<<"%\033[0m"<<"\r"<<std::flush;

		for(int j=0;j<Bsize;j++)
		{

			inputValues.clear();

			inputValues.push_back(Btrk1Pt[j]);
			inputValues.push_back(abs(Btrk1Dz1[j]/Btrk1DzError1[j]));
			inputValues.push_back(abs(Btrk1Dxy1[j]/Btrk1DxyError1[j]));
			inputValues.push_back(BsvpvDistance[j]/BsvpvDisErr[j]);
			inputValues.push_back(Balpha[j]);
			inputValues.push_back(Bchi2cl[j]);

			BDT[j] = mva.GetMvaValue(inputValues);  

      // debugging
      if (debug) {
        std::vector<double> pt1list =
          {10.498739,
           10.499157,
           10.498937,
           10.497864,
           10.504761,
           10.499465,
           10.496091,
           10.501258,
           10.501850};

        std::vector<double> bdtlist =
          { 0.6392934,
            0.6215225,
            0.6559369,
            0.5754145,
            0.6784788,
            0.4103814,
            0.6560600,
            0.6033133,
            0.2050227 };

        for (auto ipt1 = 0; ipt1 < pt1list.size(); ++ipt1) {
          auto pt1 = pt1list[ipt1];

          if (abs(Btrk1Pt[j] - pt1) < 1e-6) {
            cout << "evt: " << evt << "\n";
            for (auto i : inputValues) {
              cout << i << " : ";
            }
            cout << "BDT: " << BDT[j] << ", " << bdtlist[ipt1] << "\n";
          }
        }
      }
    }

    mvaTree->Fill();
	}
	std::cout<<std::endl<<"  Processed "<<"\033[1;31m"<<t->GetEntries()<<"\033[0m event(s)."<<std::endl;
	outf->Write();
	outf->Close();
}



int main(int argc, char* argv[])
{
	if(argc==5)
	{
		BDT(argv[1],argv[2],atof(argv[3]), atof(argv[4]));
		return 0;
	}
	else
	{
		std::cout<<"  Error: invalid argument number - BDT()"<<std::endl;
		return 1;
	}
}

