#ifndef _VARCOMP_NTUPLE_H_
#define _VARCOMP_NTUPLE_H_

#include <TTree.h>
#include <TMath.h>

#ifndef MAX_XB
#define MAX_XB     20000
#endif

#ifndef MASS_JPSI
#define MASS_JPSI  3.096900
#endif

namespace mytmva
{
  class ntuple
  {
  public:
    ntuple(TTree* nt) : fnt(nt) { fweight = nt->FindBranch("pthatweight"); setbranchaddress(); }
    ~ntuple() { fnt = 0; }
    bool passedpre(int j);
    bool passedsig(int j) { return Bgen[j]==23333; }
    bool passedbkg(int j) { return ((Bmass[j]>4.98&&Bmass[j]<5.08) || (Bmass[j]>5.48&&Bmass[j]<5.58)); }
    bool isweight() { return fweight; }

    TTree* getnt() { return fnt; }

    float  pthatweight;
    int    Bsize;
    float  Bgen[MAX_XB];
    float  Bpt[MAX_XB];
    float  Balpha[MAX_XB];
    float  Btrk1Pt[MAX_XB];
    float  Btrk2Pt[MAX_XB];
    float  Bmass[MAX_XB];
    float  Btktkmass[MAX_XB];
    float  Bujpt[MAX_XB];
    float  Bchi2cl[MAX_XB];
    float  Bdtheta[MAX_XB];
    float  Bujphi[MAX_XB];
    float  Btrk1Phi[MAX_XB];
    float  Btrk2Phi[MAX_XB];
    float  BsvpvDistance[MAX_XB];
    float  BsvpvDisErr[MAX_XB];
    float  BsvpvDistance_2D[MAX_XB];
    float  BsvpvDisErr_2D[MAX_XB];
    float  Bujeta[MAX_XB];
    float  Btrk1Eta[MAX_XB];
    float  Btrk2Eta[MAX_XB];
    float  Btrk1Dxy[MAX_XB];
    float  Btrk2Dxy[MAX_XB];
    float  Btrk1D0Err[MAX_XB];
    float  Btrk2D0Err[MAX_XB];
    float  By[MAX_XB];

  private:
    TTree* fnt;
    bool   fweight;
    void setbranchaddress();
    //
    bool   Bmu1TMOneStationTight[MAX_XB];
    int    Bmu1InPixelLayer[MAX_XB];
    int    Bmu1InStripLayer[MAX_XB];
    float  Bmu1dxyPV[MAX_XB];
    float  Bmu1dzPV[MAX_XB];
    bool   Bmu1isGlobalMuon[MAX_XB];
    float  Bmu1eta[MAX_XB];
    float  Bmu1pt[MAX_XB];
    bool   Bmu2TMOneStationTight[MAX_XB];
    int    Bmu2InPixelLayer[MAX_XB];
    int    Bmu2InStripLayer[MAX_XB];
    float  Bmu2dxyPV[MAX_XB];
    float  Bmu2dzPV[MAX_XB];
    bool   Bmu2isGlobalMuon[MAX_XB];
    float  Bmu2eta[MAX_XB];
    float  Bmu2pt[MAX_XB];
    float  Bmumumass[MAX_XB];
    bool   Btrk1highPurity[MAX_XB];
    float  Btrk1PixelHit[MAX_XB];
    float  Btrk1StripHit[MAX_XB];
    float  Btrk1PtErr[MAX_XB];
    bool   Btrk2highPurity[MAX_XB];
    float  Btrk2PixelHit[MAX_XB];
    float  Btrk2StripHit[MAX_XB];
    float  Btrk2PtErr[MAX_XB];
  };
}

bool mytmva::ntuple::passedpre(int j)
{
  if( Bmu1TMOneStationTight[j] && Bmu1InPixelLayer[j]>0 && (Bmu1InPixelLayer[j]+Bmu1InStripLayer[j])>5 && Bmu1dxyPV[j]<0.3 && Bmu1dzPV[j]<20 && Bmu1isGlobalMuon[j] && ((TMath::Abs(Bmu1eta[j])<1.2 && Bmu1pt[j]>3.5) || (TMath::Abs(Bmu1eta[j])>1.2 && TMath::Abs(Bmu1eta[j])<2.1 && Bmu1pt[j]>(5.77-1.8*TMath::Abs(Bmu1eta[j]))) || (TMath::Abs(Bmu1eta[j])>2.1 && TMath::Abs(Bmu1eta[j])<2.4 && Bmu1pt[j]>1.8)) && //mu1
Bmu2TMOneStationTight[j] && Bmu2InPixelLayer[j]>0 && (Bmu2InPixelLayer[j]+Bmu2InStripLayer[j])>5 && Bmu2dxyPV[j]<0.3 && Bmu2dzPV[j]<20 && Bmu2isGlobalMuon[j] && ((TMath::Abs(Bmu2eta[j])<1.2 && Bmu2pt[j]>3.5) || (TMath::Abs(Bmu2eta[j])>1.2 && TMath::Abs(Bmu2eta[j])<2.1 && Bmu2pt[j]>(5.77-1.8*TMath::Abs(Bmu2eta[j]))) || (TMath::Abs(Bmu2eta[j])>2.1 && TMath::Abs(Bmu2eta[j])<2.4 && Bmu2pt[j]>1.8)) && // mu2
TMath::Abs(Bmumumass[j]-3.096900)<0.15 && // jpsi
Btrk1highPurity[j] && TMath::Abs(Btrk1Eta[j])<2.4 && Btrk1Pt[j]>1.0 && TMath::Abs(Btrk1PtErr[j]/Btrk1Pt[j])<0.3 && // trk1
TMath::Abs(By[j])<2.4 && Bchi2cl[j]>0.05 // B
      ) return true;
  return false;
  /*
  */
}

void mytmva::ntuple::setbranchaddress()
{
  // public:
  if(fweight) { fnt->SetBranchAddress("pthatweight", &pthatweight); }

  fnt->SetBranchAddress("Bsize", &Bsize);
  fnt->SetBranchAddress("Bgen", Bgen);
  fnt->SetBranchAddress("Bpt", Bpt);
  fnt->SetBranchAddress("Balpha", Balpha);
  fnt->SetBranchAddress("Btrk1Pt", Btrk1Pt);
  fnt->SetBranchAddress("Btrk2Pt", Btrk2Pt);
  fnt->SetBranchAddress("Bmass", Bmass);
  fnt->SetBranchAddress("Btktkmass", Btktkmass);
  fnt->SetBranchAddress("Bujpt", Bujpt);
  fnt->SetBranchAddress("Bchi2cl", Bchi2cl);
  fnt->SetBranchAddress("Bdtheta", Bdtheta);
  fnt->SetBranchAddress("Bujphi", Bujphi);
  fnt->SetBranchAddress("Btrk1Phi", Btrk1Phi);
  fnt->SetBranchAddress("Btrk2Phi", Btrk2Phi);
  fnt->SetBranchAddress("BsvpvDistance", BsvpvDistance);
  fnt->SetBranchAddress("BsvpvDisErr", BsvpvDisErr);
  fnt->SetBranchAddress("BsvpvDistance_2D", BsvpvDistance_2D);
  fnt->SetBranchAddress("BsvpvDisErr_2D", BsvpvDisErr_2D);
  fnt->SetBranchAddress("Bujeta", Bujeta);
  fnt->SetBranchAddress("Btrk1Eta", Btrk1Eta);
  fnt->SetBranchAddress("Btrk2Eta", Btrk2Eta);
  fnt->SetBranchAddress("Btrk1Dxy", Btrk1Dxy);
  fnt->SetBranchAddress("Btrk2Dxy", Btrk2Dxy);
  fnt->SetBranchAddress("Btrk1D0Err", Btrk1D0Err);
  fnt->SetBranchAddress("Btrk2D0Err", Btrk2D0Err);
  // private
  fnt->SetBranchAddress("Bmu1TMOneStationTight", Bmu1TMOneStationTight);
  fnt->SetBranchAddress("Bmu1InPixelLayer", Bmu1InPixelLayer);
  fnt->SetBranchAddress("Bmu1InStripLayer", Bmu1InStripLayer);
  fnt->SetBranchAddress("Bmu1dxyPV", Bmu1dxyPV);
  fnt->SetBranchAddress("Bmu1dzPV", Bmu1dzPV);
  fnt->SetBranchAddress("Bmu1isGlobalMuon", Bmu1isGlobalMuon);
  fnt->SetBranchAddress("Bmu1eta", Bmu1eta);
  fnt->SetBranchAddress("Bmu1pt", Bmu1pt);
  fnt->SetBranchAddress("Bmu2TMOneStationTight", Bmu2TMOneStationTight);
  fnt->SetBranchAddress("Bmu2InPixelLayer", Bmu2InPixelLayer);
  fnt->SetBranchAddress("Bmu2InStripLayer", Bmu2InStripLayer);
  fnt->SetBranchAddress("Bmu2dxyPV", Bmu2dxyPV);
  fnt->SetBranchAddress("Bmu2dzPV", Bmu2dzPV);
  fnt->SetBranchAddress("Bmu2isGlobalMuon", Bmu2isGlobalMuon);
  fnt->SetBranchAddress("Bmu2eta", Bmu2eta);
  fnt->SetBranchAddress("Bmu2pt", Bmu2pt);
  fnt->SetBranchAddress("Bmumumass", Bmumumass);
  fnt->SetBranchAddress("Btrk1highPurity", Btrk1highPurity);
  fnt->SetBranchAddress("Btrk1PixelHit", Btrk1PixelHit);
  fnt->SetBranchAddress("Btrk1StripHit", Btrk1StripHit);
  fnt->SetBranchAddress("Btrk1PtErr", Btrk1PtErr);
  fnt->SetBranchAddress("Btrk2highPurity", Btrk2highPurity);
  fnt->SetBranchAddress("Btrk2PixelHit", Btrk2PixelHit);
  fnt->SetBranchAddress("Btrk2StripHit", Btrk2StripHit);
  fnt->SetBranchAddress("Btrk2PtErr", Btrk2PtErr);
  fnt->SetBranchAddress("By", By);
  fnt->SetBranchAddress("Btktkmass", Btktkmass);
}

#endif
