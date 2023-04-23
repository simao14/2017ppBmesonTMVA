i=$1
N=1
ptmin=( 7 10 15 20  )
ptmax=( 10 15 20 50 )

doMVA=1
doMerge=0
doReadxml=0

mvatype="BDTs"

# InputB="../sample/BsData_${ptmin[ i ]}_${ptmax[ i ]}.root"
# InputB="/home/tasheng/braa/Unskimmed/BsData.root"
InputB="../sample/BsData_7_50.root"
OutputB="../output/Bs_Data_${mvatype}_nom_pt_${ptmin[ i ]}_${ptmax[ i ]}.root"

# InputS="~/braa/Unskimmed/NewOfficialMC/BsMC.root"
# OutputS="~/dat/Unskimmed_gen/Bs_MC_${mvatype}_trk5_pt_${ptmin[ i ]}_${ptmax[ i ]}.root"

# InputS="../sample/BsMC_7_50.root"
# OutputS="../output/Bs_MC_${mvatype}_trk5_pt_${ptmin[ i ]}_${ptmax[ i ]}.root"

InputS="~/dat/presel/BsMC.root"
OutputS="../output/Bs_MC_${mvatype}_nom_pt_${ptmin[ i ]}_${ptmax[ i ]}.root"

# MergedB="../sample/BsData_a_${ptmin[ i ]}_${ptmax[ i ]}.root"
# MergedS="../sample/BsMC_a_${ptmin[ i ]}_${ptmax[ i ]}.root"
MergedB="../output/Bs_Data_BDTs_nom.root"
MergedS="../output/Bs_MC_BDTs_nom.root"

#OUTPUTNAME="BP_BDT_Training_${ptmin}_${ptmax}"
# echo $InputS; exit


collisionsyst="pp"

lumi=303.25
#lumi=76.0

weightfunctiongen="weight"
weightfunctionreco="weight"






cut="(abs(Btktkmass-1.019455)<0.015)&&(((((abs(Btktkmass-1.019455)<0.015)&& TMath::Abs(Bmumumass-3.096916)<0.15 && Bpt > 0 && Bpt < 5 && (abs(Btrk1Eta)<2.4 && abs(Btrk2Eta)<2.4 && Btrk1Pt>0.0 && Btrk2Pt>0.0) && Btrk1Pt > 0.5 && Btrk2Pt > 0.5  && Bchi2cl > 0.05 && BsvpvDistance/BsvpvDisErr > 2.0)  && ( (Bpt < 2 && Bpt > 0) || (Bpt < 3 && Bpt > 2) || (Bpt < 5 && Bpt > 3)  )))  ||  ( Bpt > 3 && ((pPAprimaryVertexFilter == 1 && pBeamScrapingFilter == 1 && HLT_HIL1DoubleMu0_v1 == 1 && (abs(PVz)<15))  &&  (Bmu1isTriggered == 1 && Bmu2isTriggered == 1 ) &&  (Bchi2cl > 0.05 && BsvpvDistance/BsvpvDisErr > 2.0)    && (TMath::Abs(By)<2.4&&TMath::Abs(Bmumumass-3.096916)<0.15&&((abs(Bmu1eta)<1.2&&Bmu1pt>3.5)||(abs(Bmu1eta)>1.2&&abs(Bmu1eta)<2.1&&Bmu1pt>(5.47-1.89*abs(Bmu1eta)))||(abs(Bmu1eta)>2.1&&abs(Bmu1eta)<2.4&&Bmu1pt>1.5))&&((abs(Bmu2eta)<1.2&&Bmu2pt>3.5)||(abs(Bmu2eta)>1.2&&abs(Bmu2eta)<2.1&&Bmu2pt>(5.47-1.89*abs(Bmu2eta)))||(abs(Bmu2eta)>2.1&&abs(Bmu2eta)<2.4&&Bmu2pt>1.5))&&Bmu1InPixelLayer>0&&(Bmu1InPixelLayer+Bmu1InStripLayer)>5&&Bmu2InPixelLayer>0&&(Bmu2InPixelLayer+Bmu2InStripLayer)>5&&Bmu1dxyPV<0.3&&Bmu2dxyPV<0.3&&Bmu1dzPV<20&&Bmu2dzPV<20&&Bmu1isTrackerMuon&&Bmu2isTrackerMuon&&Bmu1isGlobalMuon&&Bmu2isGlobalMuon)  && ( Btrk1Pt > 0.5 && Btrk2Pt > 0.5 && abs(Btrk1Eta-0.0) < 2.4 && abs(Btrk2Eta-0.0) < 2.4  && Btrk1highPurity  && Btrk2highPurity  && Btrk1PixelHit + Btrk1StripHit > 10  && Btrk2PixelHit + Btrk2StripHit > 10) &&  (Btrk1PtErr/Btrk1Pt < 0.1)  &&  (Btrk2PtErr/Btrk2Pt < 0.1)    && Btrk1Chi2ndf/(Btrk1nStripLayer+Btrk1nPixelLayer) < 0.18   && Btrk2Chi2ndf/(Btrk2nStripLayer+Btrk2nPixelLayer) < 0.18 )))"

# cut="(abs(Btktkmass-1.019455)<0.015)&&(((((abs(Btktkmass-1.019455)<0.015) && Bpt > 0 && Bpt < 5 && (abs(Btrk1Eta)<2.4 && abs(Btrk2Eta)<2.4 && Btrk1Pt>0.0 && Btrk2Pt>0.0) && Btrk1Pt > 0.5 && Btrk2Pt > 0.5  )  && ( (Bpt < 2 && Bpt > 0) || (Bpt < 3 && Bpt > 2) || (Bpt < 5 && Bpt > 3)  )))  ||  ( Bpt > 3 && ((pPAprimaryVertexFilter == 1 && pBeamScrapingFilter == 1 && HLT_HIL1DoubleMu0_v1 == 1 && (abs(PVz)<15))  &&  (Bmu1isTriggered == 1 && Bmu2isTriggered == 1 ) &&  (TMath::Abs(By)<2.4&&((abs(Bmu1eta)<1.2&&Bmu1pt>3.5)||(abs(Bmu1eta)>1.2&&abs(Bmu1eta)<2.1&&Bmu1pt>(5.47-1.89*abs(Bmu1eta)))||(abs(Bmu1eta)>2.1&&abs(Bmu1eta)<2.4&&Bmu1pt>1.5))&&((abs(Bmu2eta)<1.2&&Bmu2pt>3.5)||(abs(Bmu2eta)>1.2&&abs(Bmu2eta)<2.1&&Bmu2pt>(5.47-1.89*abs(Bmu2eta)))||(abs(Bmu2eta)>2.1&&abs(Bmu2eta)<2.4&&Bmu2pt>1.5))&&Bmu1InPixelLayer>0&&(Bmu1InPixelLayer+Bmu1InStripLayer)>5&&Bmu2InPixelLayer>0&&(Bmu2InPixelLayer+Bmu2InStripLayer)>5&&Bmu1dxyPV<0.3&&Bmu2dxyPV<0.3&&Bmu1dzPV<20&&Bmu2dzPV<20&&Bmu1isTrackerMuon&&Bmu2isTrackerMuon&&Bmu1isGlobalMuon&&Bmu2isGlobalMuon)  && ( Btrk1Pt > 0.5 && Btrk2Pt > 0.5 && abs(Btrk1Eta-0.0) < 2.4 && abs(Btrk2Eta-0.0) < 2.4  && Btrk1highPurity  && Btrk2highPurity  && Btrk1PixelHit + Btrk1StripHit > 10  && Btrk2PixelHit + Btrk2StripHit > 10) &&  (Btrk1PtErr/Btrk1Pt < 0.1)  &&  (Btrk2PtErr/Btrk2Pt < 0.1)    && Btrk1Chi2ndf/(Btrk1nStripLayer+Btrk1nPixelLayer) < 0.18   && Btrk2Chi2ndf/(Btrk2nStripLayer+Btrk2nPixelLayer) < 0.18 )))"

# CUT="${cut} &&(BsvpvDistance_2D > 0) && (BsvpvDisErr_2D>0) && (BsvpvDistance_2D/BsvpvDisErr_2D > 2)"
CUT="${cut}"


#CUT="((abs(PVz)<15))"

MYCUTS=("(${CUT}) && Bgen==23333")

#MYCUTB=("${CUT}&& ((Bmass - 5.27929 ) > 0.20 &&  (Bmass - 5.27929) < 0.25)")

MYCUTB=("${CUT[0]}&& abs(Bmass - 5.36682 ) > 0.20 &&  abs(Bmass - 5.36682) < 0.30")

#MYCUTG=("(Gy>-2.40&&Gy<2.40) && (TMath::Abs(GpdgId)==521 && GisSignal==1 && Collisions==0)");
#MYCUTG=("abs(Gy)<2.4&&abs(GpdgId)==511&&GisSignal>0");
MYCUTG="TMath::Abs(Gy)<2.4&&TMath::Abs(GpdgId)==531&&GisSignal>0"

# copy weight files
xmlfile="weights/TMVAClassification_${mvatype}.weights.xml"


if [ $doMVA -eq 1 ]; then   

    # weight_c="../train/dataset/weights/rootfiles_TMVA_B_s_BDT_BDTs_BDTs4_DNN_${ptmin[ i ]}p0_${ptmax[ i ]}p0_0-1-2-3-4-5-6-7-8-11_root/TMVAClassification_${mvatype}.class.C"
    # weight_xml="../train/dataset/weights/rootfiles_TMVA_B_s_BDT_BDTs_BDTs4_DNN_${ptmin[ i ]}p0_${ptmax[ i ]}p0_0-1-2-3-4-5-6-7-8-11_root/TMVAClassification_${mvatype}.weights.xml"
    weight_c="/home/tasheng/bmva_trk5/TMVA/Bs/train/dataset/weights/rootfiles_TMVA_B_s_BDTs_BDTs4_${ptmin[ i ]}p0_${ptmax[ i ]}p0_0-1-2-3-4-5-6-7-8-11_root/TMVAClassification_${mvatype}.class.C"
    weight_xml="/home/tasheng/bmva_trk5/TMVA/Bs/train/dataset/weights/rootfiles_TMVA_B_s_BDTs_BDTs4_${ptmin[ i ]}p0_${ptmax[ i ]}p0_0-1-2-3-4-5-6-7-8-11_root/TMVAClassification_${mvatype}.weights.xml"

    echo "copying weight files"
    cp $weight_c readxml/weights/
    cp $weight_xml readxml/weights/

    echo "calculating BDT cuts with pT ${ptmin[ i ]} - ${ptmax[ i ]}"
    echo "calculating BDT cuts with ${InputS} ${OutputS}"
	g++ ${mvatype}.C $(root-config --cflags --libs) -g -o  ${mvatype}.exe
	./${mvatype}.exe "$InputS" "$OutputS"  ${ptmin[ i ]} ${ptmax[ i ]} &
	# gdb -ex=r --args ${mvatype}.exe "$InputS" "$OutputS"  $ptmin $ptmax
	./${mvatype}.exe "$InputB" "$OutputB"  ${ptmin[ i ]} ${ptmax[ i ]} &
	# gdb -ex=r --args ${mvatype}.exe "$InputB" "$OutputB"  $ptmin $ptmax
  wait
	rm ${mvatype}.exe

fi


if [ $doMerge -eq 1 ]; then
	  hadd -f $MergedS $InputS ../output/Bs_MC_${mvatype}_nom_pt_*.root
	  hadd -f $MergedB $InputB ../output/Bs_Data_${mvatype}_nom_pt_*.root
fi



cd readxml/

if [ $doReadxml -eq 1 ]; then

    weight_c="../../train/dataset/weights/rootfiles_TMVA_B_s_BDTs_BDTs4_${ptmin[ i ]}p0_${ptmax[ i ]}p0_0-1-2-3-4-5-6-7-8-11_root/TMVAClassification_BDTs4.class.C"
    weight_xml="../../train/dataset/weights/rootfiles_TMVA_B_s_BDTs_BDTs4_${ptmin[ i ]}p0_${ptmax[ i ]}p0_0-1-2-3-4-5-6-7-8-11_root/TMVAClassification_BDTs4.weights.xml"
    echo "copying weight files"
    cp $weight_c weights/
    cp $weight_xml weights/

    echo "ptmin Now = " ${ptmin[i]}  "     ptmax Now =  "  ${ptmax[i]} 

    OUTPUTNAME="Bs_${mvatype}_Training_${ptmin[i]}_${ptmax[i]}"
    xmlfile="weights/TMVAClassification_${mvatype}.weights.xml"
    

    g++ readxml.cc $(root-config --cflags --libs) -l TMVA -l XMLIO -g -o readxml.exe

    ./readxml.exe "../$MergedS"  "../$MergedB"  "$OUTPUTNAME"  "$xmlfile"  "$collisionsyst"  "$MYCUTS"  "$MYCUTB"  "$MYCUTG"  "$weightfunctiongen"  "$weightfunctionreco"  "$mvatype"  ${ptmin[i]}  ${ptmax[i]}  $lumi  $R

fi	
