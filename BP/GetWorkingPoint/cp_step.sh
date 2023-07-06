ptmin=5
ptmax=10
stages="0-2-4-7-8-11"
method="BDTs"
method2="BDTh"

cp ../train/dataset/weights/rootfiles_TMVA_B_s_${method}_${method2}_${ptmin}p0_${ptmax}p0_${stages}_root/TMVAClassification_${method}.class.C readxmlCut/weights

cp ../train/dataset/weights/rootfiles_TMVA_B_s_${method}_${method2}_${ptmin}p0_${ptmax}p0_${stages}_root/TMVAClassification_${method}.weights.xml readxmlCut/weights

cp ../train/dataset/weights/rootfiles_TMVA_B_s_${method}_${method2}_${ptmin}p0_${ptmax}p0_${stages}_root/TMVAClassification_${method2}.class.C readxmlCut/weights

cp ../train/dataset/weights/rootfiles_TMVA_B_s_${method}_${method2}_${ptmin}p0_${ptmax}p0_${stages}_root/TMVAClassification_${method2}.weights.xml readxmlCut/weights