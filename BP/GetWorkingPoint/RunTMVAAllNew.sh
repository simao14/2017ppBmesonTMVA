echo "Note: using bmva_trk5 weights"
i=$1
N=1
# i=0
# N=6
ptmin=(5 7 10 15 20 50)
ptmax=(7 10 15 20 50 60)

# jpsi=true
jpsi=false

doMVA=1
doMerge=0
doReadxml=0

mvatype="BDTs"

# InputB="/data/szhaozho/Bmeson2017pp/BPData.root"
# OutputB="/data/szhaozho/Bmeson2017pp/${mvatype}Output/BP_Data_${mvatype}_pt_${ptmin}_${ptmax}.root"

InputB="../sample/BPData_5_60.root"
OutputB="../output/BP_Data_${mvatype}_nom_pt_${ptmin[ i ]}_${ptmax[ i ]}.root"


# InputS="~/braa/Unskimmed/NewOfficialMC/BPMC.root"
# OutputS="~/dat/Unskimmed_gen/BP_MC_${mvatype}_trk5_pt_${ptmin[ i ]}_${ptmax[ i ]}.root"

# InputS="../sample/BPMC_5_60.root"
# OutputS="../output/BP_MC_${mvatype}_trk5_pt_${ptmin[ i ]}_${ptmax[ i ]}.root"
InputS="~/dat/presel/BPMC.root"
OutputS="../output/BP_MC_${mvatype}_nom_pt_${ptmin[ i ]}_${ptmax[ i ]}.root"

MergedB="../output/BP_Data_BDTs_nom.root"
MergedS="../output/BP_MC_BDTs_nom.root"


if [ "$jpsi" = true ]; then
    InputS="../sample/jpsinp_5_60.root"
    OutputS="../output/jpsinp_nom_${ptmin[ i ]}_${ptmax[ i ]}.root"
    MergedS="../output/jpsinp_nom.root"
fi



collisionsyst="pp"

lumi=303.25
#lumi=76.0

weightfunctiongen="weight"
weightfunctionreco="weight"




RAA=1.0


# to match the skim cut
CUT="(pPAprimaryVertexFilter == 1 && pBeamScrapingFilter == 1 && HLT_HIL1DoubleMu0_v1 == 1)  &&  (Bmu1isTriggered == 1 && Bmu2isTriggered == 1 ) && (Btrk1Pt > 0.5 && Bchi2cl > 0.05 && BsvpvDistance/BsvpvDisErr > 2.0 && Bpt > 2 && abs(Btrk1Eta-0.0) < 2.4  && (TMath::Abs(By)<2.4&&TMath::Abs(Bmumumass-3.096916)<0.15&&((abs(Bmu1eta)<1.2&&Bmu1pt>3.5)||(abs(Bmu1eta)>1.2&&abs(Bmu1eta)<2.1&&Bmu1pt>(5.47-1.89*abs(Bmu1eta)))||(abs(Bmu1eta)>2.1&&abs(Bmu1eta)<2.4&&Bmu1pt>1.5))&&((abs(Bmu2eta)<1.2&&Bmu2pt>3.5)||(abs(Bmu2eta)>1.2&&abs(Bmu2eta)<2.1&&Bmu2pt>(5.47-1.89*abs(Bmu2eta)))||(abs(Bmu2eta)>2.1&&abs(Bmu2eta)<2.4&&Bmu2pt>1.5))&&Bmu1TMOneStationTight&&Bmu2TMOneStationTight&&Bmu1InPixelLayer>0&&(Bmu1InPixelLayer+Bmu1InStripLayer)>5&&Bmu2InPixelLayer>0&&(Bmu2InPixelLayer+Bmu2InStripLayer)>5&&Bmu1dxyPV<0.3&&Bmu2dxyPV<0.3&&Bmu1dzPV<20&&Bmu2dzPV<20&&Bmu1isTrackerMuon&&Bmu2isTrackerMuon&&Bmu1isGlobalMuon&&Bmu2isGlobalMuon&&Btrk1highPurity&&abs(Btrk1Eta)<2.4&&Btrk1Pt>0.5)  && (Btrk1PixelHit + Btrk1StripHit > 10) &&  (Btrk1PtErr/Btrk1Pt < 0.1)&& Btrk1Chi2ndf/(Btrk1nStripLayer+Btrk1nPixelLayer) < 0.18   && (abs(PVz)<15))"
# cut="(pPAprimaryVertexFilter == 1 && pBeamScrapingFilter == 1 && HLT_HIL1DoubleMu0_v1 == 1)  &&  (Bmu1isTriggered == 1 && Bmu2isTriggered == 1 ) && (Btrk1Pt > 0.5 && Bpt > 2 && abs(Btrk1Eta-0.0) < 2.4  && (TMath::Abs(By)<2.4 && ((abs(Bmu1eta)<1.2&&Bmu1pt>3.5)||(abs(Bmu1eta)>1.2&&abs(Bmu1eta)<2.1&&Bmu1pt>(5.47-1.89*abs(Bmu1eta)))||(abs(Bmu1eta)>2.1&&abs(Bmu1eta)<2.4&&Bmu1pt>1.5))&&((abs(Bmu2eta)<1.2&&Bmu2pt>3.5)||(abs(Bmu2eta)>1.2&&abs(Bmu2eta)<2.1&&Bmu2pt>(5.47-1.89*abs(Bmu2eta)))||(abs(Bmu2eta)>2.1&&abs(Bmu2eta)<2.4&&Bmu2pt>1.5))&&Bmu1TMOneStationTight&&Bmu2TMOneStationTight&&Bmu1InPixelLayer>0&&(Bmu1InPixelLayer+Bmu1InStripLayer)>5&&Bmu2InPixelLayer>0&&(Bmu2InPixelLayer+Bmu2InStripLayer)>5&&Bmu1dxyPV<0.3&&Bmu2dxyPV<0.3&&Bmu1dzPV<20&&Bmu2dzPV<20&&Bmu1isTrackerMuon&&Bmu2isTrackerMuon&&Bmu1isGlobalMuon&&Bmu2isGlobalMuon&&Btrk1highPurity&&abs(Btrk1Eta)<2.4&&Btrk1Pt>0.5)  && (Btrk1PixelHit + Btrk1StripHit > 10) &&  (Btrk1PtErr/Btrk1Pt < 0.1)&& Btrk1Chi2ndf/(Btrk1nStripLayer+Btrk1nPixelLayer) < 0.18   && (abs(PVz)<15))"



#CUT="((abs(PVz)<15))"

MYCUTS=("(${CUT}&&Bgen==23333)")

MYCUTB=("${CUT}&& ((Bmass - 5.27929 ) > 0.25 &&  (Bmass - 5.27929) < 0.30)")

#MYCUTG=("(Gy>-2.40&&Gy<2.40) && (TMath::Abs(GpdgId)==521 && GisSignal==1 && GcollisionId==0)");
MYCUTG=("abs(Gy)<2.4&&abs(GpdgId)==521&&GisSignal==1");



if [ $doMVA -eq 1 ]; then
weight_c="/home/tasheng/bmva_trk5/TMVA/BP/train/dataset/weights/rootfiles_TMVA_B_s_BDT_BDTs_${ptmin[ i ]}p0_${ptmax[ i ]}p0_0-2-4-7-8-11_root/TMVAClassification_BDTs.class.C"
weight_xml="/home/tasheng/bmva_trk5/TMVA/BP/train/dataset/weights/rootfiles_TMVA_B_s_BDT_BDTs_${ptmin[ i ]}p0_${ptmax[ i ]}p0_0-2-4-7-8-11_root/TMVAClassification_BDTs.weights.xml"
# weight_c="../train/dataset/weights/rootfiles_TMVA_B_s_BDT_BDTs_${ptmin[ i ]}p0_${ptmax[ i ]}p0_0-2-4-7-8-11_root/TMVAClassification_BDTs.class.C"
# weight_xml="../train/dataset/weights/rootfiles_TMVA_B_s_BDT_BDTs_${ptmin[ i ]}p0_${ptmax[ i ]}p0_0-2-4-7-8-11_root/TMVAClassification_BDTs.weights.xml"
    echo "copying weight files"
    cp $weight_c readxml/weights/
    cp $weight_xml readxml/weights/

    echo "calculating BDT cuts with files ${InputS} ${OutputS}"
    # echo "Only doing MC"
	g++ ${mvatype}.C $(root-config --cflags --libs) -g -o  ${mvatype}.exe
	./${mvatype}.exe "$InputS" "$OutputS"  ${ptmin[ i ]} ${ptmax[ i ]} &
  if [ "$jpsi" = false ]; then
      ./${mvatype}.exe "$InputB" "$OutputB"  ${ptmin[ i ]} ${ptmax[ i ]} &
  fi
  wait
	rm ${mvatype}.exe

fi




if [ $doMerge -eq 1 ]; then
    if [ "$jpsi" = true ]; then
        # hadd -f $MergedS $InputS ../output/jpsinp_nom_pt_*.root &
        hadd -f $MergedS $InputS ../output/jpsinp_nom_*.root &
    else
	      hadd -f $MergedS $InputS ../output/BP_MC_${mvatype}_nom_pt_*.root &
	      hadd -f $MergedB $InputB ../output/BP_Data_${mvatype}_nom_pt_*.root &
    fi
    wait
fi


cd readxml/


# while [ $i -lt $N ]
# do

	if [ $doReadxml -eq 1 ]; then  

      weight_c="../../train/dataset/weights/rootfiles_TMVA_B_s_BDT_BDTs_${ptmin[ i ]}p0_${ptmax[ i ]}p0_0-2-4-7-8-11_root/TMVAClassification_BDTs.class.C"
      weight_xml="../../train/dataset/weights/rootfiles_TMVA_B_s_BDT_BDTs_${ptmin[ i ]}p0_${ptmax[ i ]}p0_0-2-4-7-8-11_root/TMVAClassification_BDTs.weights.xml"
      echo "copying weight files"
      cp $weight_c weights/
      cp $weight_xml weights/
		echo "ptmin Now = " ${ptmin[i]}  "     ptmax Now =  "  ${ptmax[i]} 
	
		OUTPUTNAME="BP_${mvatype}_Training_${ptmin[i]}_${ptmax[i]}"
		xmlfile="weights/TMVAClassification_${mvatype}.weights.xml"
		

		g++ readxml.cc $(root-config --cflags --libs) -l TMVA -l XMLIO -g -o readxml.exe

		./readxml.exe "../$MergedS"  "../$MergedB"  "$OUTPUTNAME"  "$xmlfile"  "$collisionsyst"  "$MYCUTS"  "$MYCUTB"  "$MYCUTG"  "$weightfunctiongen"  "$weightfunctionreco"  "$mvatype"  ${ptmin[i]}  ${ptmax[i]}  $lumi  $RAA

		rm readxml.exe

	fi

	i=$(($i+1))
	

# done

