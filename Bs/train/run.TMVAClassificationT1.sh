#!/bin/bash
##
ptmin=$3; ptmax=$4 ;

#ptmin=5; ptmax=10 ;


#inputs="/export/d00/scratch/zzshi/CMSSW_7_5_8_patch3/Merge/2018Ana/BsRAA2015RunII/PthatService/output/PthatMergedreweighted.root" 
#inputs="/export/d00/scratch/zzshi/CMSSW_7_5_8_patch3/Merge/2018Ana/BsTMVA/TMVAHigh3BinsRedo/MC_Bs_PbPb_TMVA_BDT_pp_CentWighted.root"
#inputs="/export/d00/scratch/zzshi/CMSSW_7_5_8_patch3/Merge/2018Ana/BsRAA2015RunII/PthatService/output/PthatMergedCentreweight.root"
inputs="/mnt/submit-hi2/d00/scratch/zzshi/CMSSW_7_5_8_patch3/Merge/2018Ana/BsRAA2015RunII/PthatService/output/OffialBdPbPbMCAll.root"
#inputs=/mnt/hadoop/cms/store/user/gwangjun/Bp2018/MC/crab_Bfinder_20190221_Pythia8_BuToJpsiK_Bpt0p0_1032_NoJSON_pthatweight_hardcut_v2.root ;
## inputs=/mnt/hadoop/cms/store/user/gwangjun/Bp2018/MC/Bntuple20160816_Bpt7svpv5p5Bpt10svpv3p5_BfinderMC_PbPb_Pythia8_BuToJpsiK_TuneCUEP8M1_20160816_bPt5jpsiPt0tkPt0p8_Bp_pthatweight_JingBDT.root ;
output=rootfiles/TMVA_X ;
inputb="/mnt/submit-hi2/d00/scratch/zzshi/CMSSW_7_5_8_patch3/Merge/ToCopy/Data/New/crab_Bfinder_20190513_HIDoubleMuon__PsiPeri__HIRun2018A_04Apr2019_v1_1033p1_GoldenJSON_skimhltBsize_ntphi.root";

inputm=$inputs ;
# outputmva=/mnt/hadoop/cms/store/user/gwangjun/tmva ;
outputmva=/mnt/submit-hi2/d00/scratch/zzshi/CMSSW_10_3_2/src/JingTMVA/TMVA/BmesonRun2018/tmva/train;

# prefilter
#OLD CUTS
#cut="Btrk1Pt > 1.0 && Btrk2Pt > 1.0 && Bchi2cl > 0.05 && BsvpvDistance/BsvpvDisErr > 2.2 && Bpt > 5 && abs(Btrk1Eta-0.0) < 2.4  && abs(Btrk2Eta-0.0) < 2.4 && (TMath::Abs(By)<2.4&&TMath::Abs(Bmumumass-3.096916)<0.15&&((abs(Bmu1eta)<1.2&&Bmu1pt>3.5)||(abs(Bmu1eta)>1.2&&abs(Bmu1eta)<2.1&&Bmu1pt>(5.77-1.8*abs(Bmu1eta)))||(abs(Bmu1eta)>2.1&&abs(Bmu1eta)<2.4&&Bmu1pt>1.8))&&((abs(Bmu2eta)<1.2&&Bmu2pt>3.5)||(abs(Bmu2eta)>1.2&&abs(Bmu2eta)<2.1&&Bmu2pt>(5.77-1.8*abs(Bmu2eta)))||(abs(Bmu2eta)>2.1&&abs(Bmu2eta)<2.4&&Bmu2pt>1.8))&&Bmu1TMOneStationTight&&Bmu2TMOneStationTight&&Bmu1InPixelLayer>0&&(Bmu1InPixelLayer+Bmu1InStripLayer)>5&&Bmu2InPixelLayer>0&&(Bmu2InPixelLayer+Bmu2InStripLayer)>5&&Bmu1dxyPV<0.3&&Bmu2dxyPV<0.3&&Bmu1dzPV<20&&Bmu2dzPV<20&&Bmu1isGlobalMuon&&Bmu2isGlobalMuon&&Btrk1highPurity&&Btrk2highPurity&&abs(Btrk1Eta)<2.4&&abs(Btrk2Eta)<2.4&&Btrk1Pt>1.&&Btrk2Pt>1.&&abs(Btktkmass-1.019455)<0.015) && (abs(PVz)<15&&pclusterCompatibilityFilter&&pprimaryVertexFilter) && (Btrk1PixelHit + Btrk1StripHit > 10) && (Btrk2PixelHit + Btrk2StripHit > 10) && (Btrk1PtErr/Btrk1Pt < 0.1)&& (Btrk2PtErr/Btrk2Pt < 0.1) && Btrk1Chi2ndf/(Btrk1nStripLayer+Btrk1nPixelLayer) < 0.18 && Btrk2Chi2ndf/(Btrk2nStripLayer+Btrk2nPixelLayer) < 0.18"

cut="(hiBin < 181) && Btrk1Pt > 0.9 && Btrk2Pt > 0.9 && Bchi2cl > 0.05 && BsvpvDistance/BsvpvDisErr > 2.0 && Bpt > 5 && abs(Btrk1Eta-0.0) < 2.4  && abs(Btrk2Eta-0.0) < 2.4 && (TMath::Abs(By)<2.4&&TMath::Abs(Bmumumass-3.096916)<0.15&&((abs(Bmu1eta)<1.2&&Bmu1pt>3.5)||(abs(Bmu1eta)>1.2&&abs(Bmu1eta)<2.1&&Bmu1pt>(5.77-1.89*abs(Bmu1eta)))||(abs(Bmu1eta)>2.1&&abs(Bmu1eta)<2.4&&Bmu1pt>1.5))&&((abs(Bmu2eta)<1.2&&Bmu2pt>3.5)||(abs(Bmu2eta)>1.2&&abs(Bmu2eta)<2.1&&Bmu2pt>(5.77-1.89*abs(Bmu2eta)))||(abs(Bmu2eta)>2.1&&abs(Bmu2eta)<2.4&&Bmu2pt>1.5))&&Bmu1TMOneStationTight&&Bmu2TMOneStationTight&&Bmu1InPixelLayer>0&&(Bmu1InPixelLayer+Bmu1InStripLayer)>5&&Bmu2InPixelLayer>0&&(Bmu2InPixelLayer+Bmu2InStripLayer)>5&&Bmu1dxyPV<0.3&&Bmu2dxyPV<0.3&&Bmu1dzPV<20&&Bmu2dzPV<20&&Bmu1isGlobalMuon&&Bmu2isGlobalMuon&&Btrk1highPurity&&Btrk2highPurity&&abs(Btrk1Eta)<2.4&&abs(Btrk2Eta)<2.4&&Btrk1Pt>1.&&Btrk2Pt>1.&&abs(Btktkmass-1.019455)<0.015) && (abs(PVz)<15&&pclusterCompatibilityFilter&&pprimaryVertexFilter) && (Btrk1PixelHit + Btrk1StripHit > 10) && (Btrk2PixelHit + Btrk2StripHit > 10) && (Btrk1PtErr/Btrk1Pt < 0.1)&& (Btrk2PtErr/Btrk2Pt < 0.1) && Btrk1Chi2ndf/(Btrk1nStripLayer+Btrk1nPixelLayer) < 0.18 && Btrk2Chi2ndf/(Btrk2nStripLayer+Btrk2nPixelLayer) < 0.18 && (HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1) && (abs(PVz)<15&&pclusterCompatibilityFilter&&pprimaryVertexFilter && phfCoincFilter2Th4)"


#algo="CutsGA,CutsSA"
#algo="BDT,BDTG,CutsGA,CutsSA"

#algo="BDT,BDTG,CutsGA,CutsSA"

algo="BDT,MLP,CutsSA,DNN"


stages="0,1,2,3,4,5,6,7,8,9,10,11,12,13" # see definition below

sequence=1 # if it's on, remove stages 1 by 1. 

## ===== do not change the lines below =====

varlist=(
    '#  0  :  ("Btrk1Pt"  , "Btrk1Pt"                                                                                        , "FMax")  #' 
    '#  1  :  ("Btrk2Pt"  , "Btrk2Pt"                                                                                        , "FMin")  #' 
	'#  2  :  ("Trk1DCAz" , "Trk1DCAz := abs(Btrk1Dz1/Btrk1DzError1)"                                                                                       , "FMin")  #' 
	'#  3  :  ("Trk2DCAz" , "Trk2DCAz := abs(Btrk2Dz1/Btrk2DzError1)"                                                                                       , "FMin")  #' 
	'#  4  :  ("Trk1DCAxy"  , "Trk1DCAxy := abs(Btrk1Dxy1/Btrk1DxyError1)"                                                                 , "FMin")  #' 
	'#  5  :  ("Trk2DCAxy"  , "Trk2DCAxy := abs(Btrk2Dxy1/Btrk2DxyError1)"                                                                 , "FMin")  #' 
    '#  6  :  ("MassDis"  , "MassDis := abs(Btktkmass-1.019455)"                                                             , "FMax")  #' 
    '#  7  :  ("dls"      , "dls := BsvpvDistance/BsvpvDisErr"                                                                 , "FMax")  #' 
    '#  8  :  ("Balpha"   , "Balpha"                                                                                          , "FMax")  #' 
    '#  9  :  ("dls2D"    , "dls2D := Bd0"                                                                             , "FMin")  #' 
    '#  10 :  ("cos(Bdtheta)", "cos(Bdtheta)"                                                                                 , "FMax")  #' 
    '#  11 :  ("Bchi2cl",    "Bchi2cl"                                                                                             , ""    )  #'
	'#  12 :  ("Btrk1Eta", "Btrk1Eta"                                                                                 , "FMin")  #' 
	'#  13 :  ("Btrk2Eta", "Btrk2Eta"                                                                                 , "FMin")  #' 
)


cuts="${cut} && Bgen==23333"
cutb="${cut} &&  (abs(Bmass - 5.36682 ) > 0.20 &&  abs(Bmass - 5.36682) < 0.30)" # Bmass_pdg=5.28GeV, sideband 0.2 ~ 0.3 for each side.
# rootfiles=rootfiles

## ===== do not do not do not change the lines below =====
IFS=','; allstages=($stages); unset IFS;
echo -e '

####################################################################################################################################
#                                                Variables \e[1;32m(To be used)\e[0m                                                            #
####################################################################################################################################
#                                                                                                                                  #'
vv=0
while(( $vv < ${#varlist[@]})) ; do
    for ss in ${allstages[@]} ; do [[ $ss == $vv ]] && { echo -en "\e[1;32m" ; break ; } ; done ;
    echo -e "${varlist[vv]}\e[0m"
    vv=$((vv+1))
done
echo '#                                                                                                                                  #
####################################################################################################################################

'
##
mkdir -p $output ; rm -r $output ;
tmp=$(date +%y%m%d%H%M%S)

##
[[ $# -eq 0 ]] && echo "usage: ./run_TMVAClassification.sh [train] [draw curves] [create BDT tree]"

g++ TMVAClassification.C $(root-config --libs --cflags) -lTMVA -lTMVAGui -g -o TMVAClassification_${tmp}.exe || exit 1
g++ guivariables.C $(root-config --libs --cflags) -lTMVA -lTMVAGui -g -o guivariables_${tmp}.exe || { rm *_${tmp}.exe ; exit 1 ; }
g++ guiefficiencies.C $(root-config --libs --cflags) -lTMVA -lTMVAGui -g -o guiefficiencies_${tmp}.exe || { rm *_${tmp}.exe ; exit 1 ; }
g++ guieffvar.C $(root-config --libs --cflags) -lTMVA -lTMVAGui -g -o guieffvar_${tmp}.exe || { rm *_${tmp}.exe ; exit 1 ; }
g++ mvaprod.C $(root-config --libs --cflags) -lTMVA -lTMVAGui -g -o mvaprod_${tmp}.exe || { rm *_${tmp}.exe ; exit 1 ; }

[[ ${1:-0} -eq 1 ]] && {
    conf=
    echo -e "\e[2m==> Do you really want to run\e[0m \e[1mTMVAClassification.C\e[0m \e[2m(it might be very slow)?\e[0m [y/n]"
    read conf
    while [[ $conf != 'y' && $conf != 'n' ]] ; do { echo "warning: input [y/n]" ; read conf ; } ; done ;
    [[ $conf == 'n' ]] && { rm *_${tmp}.exe ; exit ; }
}

stage=$stages
while [[ $stage == *,* ]]
do
# train
    [[ ${1:-0} -eq 1 ]] && { ./TMVAClassification_${tmp}.exe $inputs $inputb "$cuts" "$cutb" $output $ptmin $ptmax "$algo" "$stage"; } &
    [[ $sequence -eq 0 ]] && break;
    while [[ $stage != *, ]] ; do stage=${stage%%[0-9]} ; done ;
    stage=${stage%%,}
done
wait

# draw curves
[[ ${2:-0} -eq 1 ]] && { 
    ./guivariables_${tmp}.exe $output $ptmin $ptmax "$algo" "$stages"
    ./guiefficiencies_${tmp}.exe $output $ptmin $ptmax "$algo" "$stages"
}

# draw curve vs. var
[[ ${2:-0} -eq 1 && $sequence -eq 1 ]] && ./guieffvar_${tmp}.exe $output $ptmin $ptmax "$algo" "$stages"

# produce mva values
[[ ${3:-0} -eq 1 ]] && ./mvaprod_${tmp}.exe $inputm "Bfinder/ntphi" $output $outputmva $ptmin $ptmax "$algo" "${stages}"

##
rm *_${tmp}.exe
