// https://root.cern.ch/doc/v610/TMVAGui_8cxx_source.html
// https://root.cern.ch/doc/v608/efficiencies_8cxx_source.html
#include <iostream>
#include <string>

#include <TKey.h>
#include <TList.h>
#include "TMVA/efficiencies.h"
#include "TMVA/mvas.h"
#include "TMVA/correlations.h"

#include "mvaeffs.h"
#include "xjjcuti.h"

namespace mytmva
{
  void guiefficiencies();
  void efficiencies(std::string outfname);
}

void mytmva::guiefficiencies()
{
  std::string outfname = "dataset/results/rootfiles/TMVA_pytorch_10.0_15.0_0_2_4_7_8_11.root";
  mytmva::efficiencies(outfname);
}

void mytmva::efficiencies(std::string outfname)
{
  TString dataset("");
  std::string outputstr = xjjc::str_replaceallspecial(outfname);

  // set up dataset
  TFile* file = TFile::Open( outfname.c_str() );
  if(!file)
    { std::cout << "==> Abort " << __FUNCTION__ << ", please verify filename." << std::endl; return; }
  if(file->GetListOfKeys()->GetEntries()<=0)
    { std::cout << "==> Abort " << __FUNCTION__ << ", please verify if dataset exist." << std::endl; return; }
  // --->>
  if( (dataset==""||dataset.IsWhitespace()) && (file->GetListOfKeys()->GetEntries()==1))
    { dataset = ((TKey*)file->GetListOfKeys()->At(0))->GetName(); }
  // <<---
  else if((dataset==""||dataset.IsWhitespace()) && (file->GetListOfKeys()->GetEntries()>=1))
    {
      std::cout << "==> Warning " << __FUNCTION__ << ", more than 1 dataset." << std::endl;
      file->ls(); return;
    }
  else { return; }

  // TMVA::efficiencies(dataset.Data(), outfname.c_str(), 1);
  TMVA::efficiencies(dataset.Data(), outfname.c_str(), 2);
  // TMVA::efficiencies(dataset.Data(), outfname.c_str(), 3);
  TMVA::mvas(dataset.Data(), outfname.c_str(), TMVA::kCompareType);
  TMVA::correlations(dataset.Data(), outfname.c_str());
  //mytmva::mvaeffs(dataset.Data(), outfname.c_str());
  //mytmva::mvaeffs(dataset.Data(), outfname.c_str(), 1.e+3, 1.e+5);
  mytmva::mvaeffs(dataset.Data(), outfname.c_str(), 101.953617, 13730.981893); //S^prime and B^prime values, respectively
  //S^prime 101.953617   194.055714   498.886649  208.675714  221.456861 102.726613 9.017866
  //B^prime 13730.981893 12019.550583 6573.084962 1073.808659 325.775993 41.322664  10.236871

  //gSystem->Exec(Form("rm %s/plots/*.png", dataset.Data()));
  gSystem->Exec(Form("mkdir -p %s/plots/%s", dataset.Data(), outputstr.c_str()));
  gSystem->Exec(Form("mv %s/plots/*.png %s/plots/%s/", dataset.Data(), dataset.Data(), outputstr.c_str()));
}

int main()
{

  mytmva::guiefficiencies(); return 0; 
  
}
