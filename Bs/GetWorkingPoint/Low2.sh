i=0
N=2
ptmin=( 1 2 )
ptmax=( 2 3 )

doMVA=0
doMerge=0
doReadxml=1

mvatype="BDT"

InputB="/data/szhaozho/ppNewTMVA/CMSSW_10_3_2/src/BZ/Samples/BZData.root"
OutputB="/data/szhaozho/Bmeson2017pp/${mvatype}Output/BP_Data_${mvatype}_pt_${ptmin}_${ptmax}.root"

InputS="/data/szhaozho/ppNewTMVA/CMSSW_10_3_2/src/BZ/Samples/BZMC.root"
OutputS="/data/szhaozho/Bmeson2017pp/${mvatype}Output/BP_MC_${mvatype}_pt_${ptmin}_${ptmax}.root"

#MergedB="/data/szhaozho/Bmeson2017ppNew/BDTOutput/AllMerge/BPDataAllBDT.root"
#MergedS="/data/szhaozho/Bmeson2017ppNew/BDTOutput/AllMerge/BPMCAllBDT.root"

#MergedB="/data/szhaozho/ppNewTMVA/CMSSW_10_3_2/src/BZ/Samples/Final/BZDataAllBDT.root"
#MergedS="/data/szhaozho/ppNewTMVA/CMSSW_10_3_2/src/BZ/Samples/Final/BZMCAllBDT.root"

MergedB="/data/szhaozho/ppNewTMVA/CMSSW_10_3_2/src/Bs/Samples/BsNewLow.root"
MergedS="/data/szhaozho/ppNewTMVA/CMSSW_10_3_2/src/Bs/Samples/BsNewLowMC.root"

#MergedB="/data/szhaozho/Bmeson2017pp/BDTOutput/Backup/BP_Data_${mvatype}_Merged_pt_${ptmin}_${ptmax}.root"
#MergedS="/data/szhaozho/Bmeson2017pp/BDTOutput/Backup/BP_MC_${mvatype}_Merged_pt_${ptmin}_${ptmax}.root"

#OUTPUTNAME="BP_BDT_Training_${ptmin}_${ptmax}"



collisionsyst="pp"

lumi=303.25
#lumi=76.0

weightfunctiongen="weight"
weightfunctionreco="weight"




RAA=0.3

#CUT="(Btrk1Pt > 0.2 && Btrk2Pt > 0.2 && Bchi2cl > 0.05 && BsvpvDistance/BsvpvDisErr > 2.0 && Bpt > 2 && abs(Btrk1Eta-0.0) < 2.4  && abs(Btrk2Eta-0.0) < 2.4 && (TMath::Abs(By)<2.4&&TMath::Abs(Bmumumass-3.096916)<0.15&&((abs(Bmu1eta)<1.2&&Bmu1pt>3.5)||(abs(Bmu1eta)>1.2&&abs(Bmu1eta)<2.1&&Bmu1pt>(5.47-1.89*abs(Bmu1eta)))||(abs(Bmu1eta)>2.1&&abs(Bmu1eta)<2.4&&Bmu1pt>1.5))&&((abs(Bmu2eta)<1.2&&Bmu2pt>3.5)||(abs(Bmu2eta)>1.2&&abs(Bmu2eta)<2.1&&Bmu2pt>(5.47-1.89*abs(Bmu2eta)))||(abs(Bmu2eta)>2.1&&abs(Bmu2eta)<2.4&&Bmu2pt>1.5))&&Bmu1TMOneStationTight&&Bmu2TMOneStationTight&&Bmu1InPixelLayer>0&&(Bmu1InPixelLayer+Bmu1InStripLayer)>5&&Bmu2InPixelLayer>0&&(Bmu2InPixelLayer+Bmu2InStripLayer)>5&&Bmu1dxyPV<0.3&&Bmu2dxyPV<0.3&&Bmu1dzPV<20&&Bmu2dzPV<20&&Bmu1isTrackerMuon&&Bmu2isTrackerMuon&&Bmu1isGlobalMuon&&Bmu2isGlobalMuon&&Btrk1highPurity&&Btrk2highPurity&&abs(Btrk1Eta)<2.4&&abs(Btrk2Eta)<2.4&&Btrk1Pt>1.&&Btrk2Pt>1.&&abs(Btktkmass-1.019455)<0.015)  && (Btrk1PixelHit + Btrk1StripHit > 10) && (Btrk2PixelHit + Btrk2StripHit > 10) && (Btrk1PtErr/Btrk1Pt < 0.1)&& (Btrk2PtErr/Btrk2Pt < 0.1) && Btrk1Chi2ndf/(Btrk1nStripLayer+Btrk1nPixelLayer) < 0.18 && Btrk2Chi2ndf/(Btrk2nStripLayer+Btrk2nPixelLayer) < 0.18  && (abs(PVz)<15))"

#CUT="(Btrk1Pt > 0.2 && Btrk2Pt > 0.2 && Bchi2cl > 0.05 && BsvpvDistance/BsvpvDisErr > 2.0 && Bpt > 2 && abs(Btrk1Eta-0.0) < 2.4  && abs(Btrk2Eta-0.0) < 2.4 && (TMath::Abs(By)<2.4&&TMath::Abs(Bmumumass-3.096916)<0.15&&((abs(Bmu1eta)<1.2&&Bmu1pt>3.5)||(abs(Bmu1eta)>1.2&&abs(Bmu1eta)<2.1&&Bmu1pt>(5.47-1.89*abs(Bmu1eta)))||(abs(Bmu1eta)>2.1&&abs(Bmu1eta)<2.4&&Bmu1pt>1.5))&&((abs(Bmu2eta)<1.2&&Bmu2pt>3.5)||(abs(Bmu2eta)>1.2&&abs(Bmu2eta)<2.1&&Bmu2pt>(5.47-1.89*abs(Bmu2eta)))||(abs(Bmu2eta)>2.1&&abs(Bmu2eta)<2.4&&Bmu2pt>1.5))&&Bmu1TMOneStationTight&&Bmu2TMOneStationTight&&Bmu1InPixelLayer>0&&(Bmu1InPixelLayer+Bmu1InStripLayer)>5&&Bmu2InPixelLayer>0&&(Bmu2InPixelLayer+Bmu2InStripLayer)>5&&Bmu1dxyPV<0.3&&Bmu2dxyPV<0.3&&Bmu1dzPV<20&&Bmu2dzPV<20&&Bmu1isTrackerMuon&&Bmu2isTrackerMuon&&Bmu1isGlobalMuon&&Bmu2isGlobalMuon&&Btrk1highPurity&&Btrk2highPurity&&abs(Btrk1Eta)<2.4&&abs(Btrk2Eta)<2.4&&Btrk1Pt>0.&&Btrk2Pt>0.&&abs(Btktkmass-1.019455)<0.015)  && (Btrk1PixelHit + Btrk1StripHit > 10) && (Btrk2PixelHit + Btrk2StripHit > 10) && (Btrk1PtErr/Btrk1Pt < 0.1)&& (Btrk2PtErr/Btrk2Pt < 0.1) && Btrk1Chi2ndf/(Btrk1nStripLayer+Btrk1nPixelLayer) < 0.18 && Btrk2Chi2ndf/(Btrk2nStripLayer+Btrk2nPixelLayer) < 0.18  && (abs(PVz)<15))"


#CUT="(Btrk1Pt > 0.2  && BsvpvDistance/BsvpvDisErr > 2.0 && Bpt > 2 && abs(Btrk1Eta-0.0) < 2.4 && (TMath::Abs(By)<2.4&&TMath::Abs(Bmumumass-3.096916)<0.15&&((abs(Bmu1eta)<1.2&&Bmu1pt>3.5)||(abs(Bmu1eta)>1.2&&abs(Bmu1eta)<2.1&&Bmu1pt>(5.47-1.89*abs(Bmu1eta)))||(abs(Bmu1eta)>2.1&&abs(Bmu1eta)<2.4&&Bmu1pt>1.5))&&((abs(Bmu2eta)<1.2&&Bmu2pt>3.5)||(abs(Bmu2eta)>1.2&&abs(Bmu2eta)<2.1&&Bmu2pt>(5.47-1.89*abs(Bmu2eta)))||(abs(Bmu2eta)>2.1&&abs(Bmu2eta)<2.4&&Bmu2pt>1.5))&&Bmu1TMOneStationTight&&Bmu2TMOneStationTight&&Bmu1InPixelLayer>0&&(Bmu1InPixelLayer+Bmu1InStripLayer)>5&&Bmu2InPixelLayer>0&&(Bmu2InPixelLayer+Bmu2InStripLayer)>5&&Bmu1dxyPV<0.3&&Bmu2dxyPV<0.3&&Bmu1dzPV<20&&Bmu2dzPV<20&&Bmu1isTrackerMuon&&Bmu2isTrackerMuon&&Bmu1isGlobalMuon&&Bmu2isGlobalMuon&&Btrk1highPurity &&Btrk1Pt>0.  && (Btrk1PixelHit + Btrk1StripHit > 10)  && (Btrk1PtErr/Btrk1Pt < 0.1)&& Btrk1Chi2ndf/(Btrk1nStripLayer+Btrk1nPixelLayer) < 0.18  && (abs(PVz)<15)))"
#CUT="(Bmu1isTriggered == 1 && Bmu2isTriggered ==1)&&(HBHENoiseFilterResult == 1 && pPAprimaryVertexFilter == 1 && pBeamScrapingFilter == 1)&&Btrk1Pt > 0.2 && Bchi2cl > 0.05 && BsvpvDistance/BsvpvDisErr > 2.0 && Bpt > 2 && abs(Btrk1Eta-0.0) < 2.4  && (TMath::Abs(By)<2.4&&TMath::Abs(Bmumumass-3.096916)<0.15&&((abs(Bmu1eta)<1.2&&Bmu1pt>3.5)||(abs(Bmu1eta)>1.2&&abs(Bmu1eta)<2.1&&Bmu1pt>(5.47-1.89*abs(Bmu1eta)))||(abs(Bmu1eta)>2.1&&abs(Bmu1eta)<2.4&&Bmu1pt>1.5))&&((abs(Bmu2eta)<1.2&&Bmu2pt>3.5)||(abs(Bmu2eta)>1.2&&abs(Bmu2eta)<2.1&&Bmu2pt>(5.47-1.89*abs(Bmu2eta)))||(abs(Bmu2eta)>2.1&&abs(Bmu2eta)<2.4&&Bmu2pt>1.5))&&Bmu1TMOneStationTight&&Bmu2TMOneStationTight&&Bmu1InPixelLayer>0&&(Bmu1InPixelLayer+Bmu1InStripLayer)>5&&Bmu2InPixelLayer>0&&(Bmu2InPixelLayer+Bmu2InStripLayer)>5&&Bmu1dxyPV<0.3&&Bmu2dxyPV<0.3&&Bmu1dzPV<20&&Bmu2dzPV<20&&Bmu1isTrackerMuon&&Bmu2isTrackerMuon&&Bmu1isGlobalMuon&&Bmu2isGlobalMuon&&Btrk1highPurity&&abs(Btrk1Eta)<2.4&&Btrk1Pt>0.2)  && (Btrk1PixelHit + Btrk1StripHit > 10) &&  (Btrk1PtErr/Btrk1Pt < 0.1)&& Btrk1Chi2ndf/(Btrk1nStripLayer+Btrk1nPixelLayer) < 0.18   && (abs(PVz)<15)"
#CUT="Btrk1Pt > 0.2 && Bchi2cl > 0.05 && BsvpvDistance/BsvpvDisErr > 2.0 && Bpt > 0 && abs(Btrk1Eta-0.0) < 2.4  && (TMath::Abs(By)<2.4&&TMath::Abs(Bmumumass-3.096916)<0.15&&((abs(Bmu1eta)<1.2&&Bmu1pt>3.5)||(abs(Bmu1eta)>1.2&&abs(Bmu1eta)<2.1&&Bmu1pt>(5.47-1.89*abs(Bmu1eta)))||(abs(Bmu1eta)>2.1&&abs(Bmu1eta)<2.4&&Bmu1pt>1.5))&&((abs(Bmu2eta)<1.2&&Bmu2pt>3.5)||(abs(Bmu2eta)>1.2&&abs(Bmu2eta)<2.1&&Bmu2pt>(5.47-1.89*abs(Bmu2eta)))||(abs(Bmu2eta)>2.1&&abs(Bmu2eta)<2.4&&Bmu2pt>1.5))&&Bmu1TMOneStationTight&&Bmu2TMOneStationTight&&Bmu1InPixelLayer>0&&(Bmu1InPixelLayer+Bmu1InStripLayer)>5&&Bmu2InPixelLayer>0&&(Bmu2InPixelLayer+Bmu2InStripLayer)>5&&Bmu1dxyPV<0.3&&Bmu2dxyPV<0.3&&Bmu1dzPV<20&&Bmu2dzPV<20&&Bmu1isTrackerMuon&&Bmu2isTrackerMuon&&Bmu1isGlobalMuon&&Bmu2isGlobalMuon&&Btrk1highPurity &&Btrk2highPurity&& abs(Btrk1Eta)<2.4 && abs(Btrk1Eta)<2.4 && Btrk1Pt>0.2 && Btrk2Pt>0.2)  && (Btrk1PixelHit + Btrk1StripHit > 10)  && (Btrk2PixelHit + Btrk2StripHit > 10) &&  (Btrk1PtErr/Btrk1Pt < 0.1) && (Btrk2PtErr/Btrk2Pt < 0.1)   && Btrk1Chi2ndf/(Btrk1nStripLayer+Btrk1nPixelLayer) < 0.18 && Btrk2Chi2ndf/(Btrk2nStripLayer+Btrk2nPixelLayer) < 0.18   && (abs(PVz)<15)"

#CUT="(abs(Btktkmass-1.019455)<0.028)"

CUT="(TMath::Abs(Bmumumass-3.096916)<0.14 && abs(Btrk1Eta)<2.4 && abs(Btrk2Eta)<2.4 && Btrk1Pt>0.0 && Btrk2Pt>0.0) && Btrk1Pt > 0.2 && Btrk2Pt> 0.2  && Bchi2cl > 0.12 && BsvpvDistance/BsvpvDisErr > 4.0 && (abs(Btktkmass-1.019455)<0.030)"

#CUT="((abs(PVz)<15))"

MYCUTS=("(${CUT}&&Bgen==23333)")

#MYCUTB=("${CUT}&& ((Bmass - 5.27929 ) > 0.20 &&  (Bmass - 5.27929) < 0.25)")

MYCUTB=("${CUT[0]}&& abs(Bmass - 5.36682 ) > 0.20 &&  abs(Bmass - 5.36682) < 0.30")

#MYCUTG=("(Gy>-2.40&&Gy<2.40) && (TMath::Abs(GpdgId)==521 && GisSignal==1 && GcollisionId==0)");
#MYCUTG=("abs(Gy)<2.4&&abs(GpdgId)==511&&GisSignal>0");
MYCUTG="TMath::Abs(Gy)<2.4&&TMath::Abs(GpdgId)==531&&GisSignal>0"


if [ $doMVA -eq 1 ]; then   

	g++ ${mvatype}.C $(root-config --cflags --libs) -g -o  ${mvatype}.exe 
	./${mvatype}.exe "$InputS" "$OutputS"  $ptmin $ptmax
	./${mvatype}.exe "$InputB" "$OutputB"  $ptmin $ptmax
	rm ${mvatype}.exe

fi




if [ $doMerge -eq 1 ]; then   

	hadd $MergedS $InputS $OutputS
	hadd $MergedB $InputB $OutputB

fi


cd readxml/


while [ $i -lt $N ]
do

	if [ $doReadxml -eq 1 ]; then  

		echo "ptmin Now = " ${ptmin[i]}  "     ptmax Now =  "  ${ptmax[i]} 
	
		OUTPUTNAME="BP_${mvatype}_Training_${ptmin[i]}_${ptmax[i]}"
		xmlfile="weights/TMVAClassification_${mvatype}.weights.xml"
		

		g++ readxml.cc $(root-config --cflags --libs) -l TMVA -l XMLIO -g -o readxml.exe

		./readxml.exe "$MergedS"  "$MergedB"  "$OUTPUTNAME"  "$xmlfile"  "$collisionsyst"  "$MYCUTS"  "$MYCUTB"  "$MYCUTG"  "$weightfunctiongen"  "$weightfunctionreco"  "$mvatype"  ${ptmin[i]}  ${ptmax[i]}  $lumi  $R

	fi	

   i=$(($i+1))
done
