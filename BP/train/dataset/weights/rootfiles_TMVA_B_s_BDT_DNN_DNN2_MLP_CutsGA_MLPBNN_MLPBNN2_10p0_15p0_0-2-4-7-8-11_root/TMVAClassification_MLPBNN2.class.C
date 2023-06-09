// Class: ReadMLPBNN2
// Automatically generated by MethodBase::MakeClass
//

/* configuration options =====================================================

#GEN -*-*-*-*-*-*-*-*-*-*-*- general info -*-*-*-*-*-*-*-*-*-*-*-

Method         : MLP::MLPBNN2
TMVA Release   : 4.2.1         [262657]
ROOT Release   : 6.10/09       [395785]
Creator        : tasheng
Date           : Mon Jan 17 19:36:33 2022
Host           : Linux cmsbuild64.cern.ch 2.6.32-696.30.1.el6.x86_64 #1 SMP Tue May 22 06:09:36 CEST 2018 x86_64 x86_64 x86_64 GNU/Linux
Dir            : /home/tasheng/bmva/TMVA/BP/train
Training events: 57215
Analysis type  : [Classification]


#OPT -*-*-*-*-*-*-*-*-*-*-*-*- options -*-*-*-*-*-*-*-*-*-*-*-*-

# Set by User:
NCycles: "600" [Number of training cycles]
HiddenLayers: "N+5,N" [Specification of hidden layer architecture]
NeuronType: "tanh" [Neuron activation function type]
V: "False" [Verbose output (short form of "VerbosityLevel" below - overrides the latter one)]
VarTransform: "N" [List of variable transformations performed before training, e.g., "D_Background,P_Signal,G,N_AllClasses" for: "Decorrelation, PCA-transformation, Gaussianisation, Normalisation, each for the given class of events ('AllClasses' denotes all events of all classes, if no class indication is given, 'All' is assumed)"]
H: "True" [Print method-specific help message]
TrainingMethod: "BFGS" [Train with Back-Propagation (BP), BFGS Algorithm (BFGS), or Genetic Algorithm (GA - slower and worse)]
TestRate: "5" [Test for overtraining performed at each #th epochs]
UseRegulator: "True" [Use regulator to avoid over-training]
# Default:
RandomSeed: "1" [Random seed for initial synapse weights (0 means unique seed for each run; default value '1')]
EstimatorType: "CE" [MSE (Mean Square Estimator) for Gaussian Likelihood or CE(Cross-Entropy) for Bernoulli Likelihood]
NeuronInputType: "sum" [Neuron input function type]
VerbosityLevel: "Default" [Verbosity level]
CreateMVAPdfs: "False" [Create PDFs for classifier outputs (signal and background)]
IgnoreNegWeightsInTraining: "False" [Events with negative weights are ignored in the training (but are included for testing and performance evaluation)]
LearningRate: "2.000000e-02" [ANN learning rate parameter]
DecayRate: "1.000000e-02" [Decay rate for learning parameter]
EpochMonitoring: "False" [Provide epoch-wise monitoring plots according to TestRate (caution: causes big ROOT output file!)]
Sampling: "1.000000e+00" [Only 'Sampling' (randomly selected) events are trained each epoch]
SamplingEpoch: "1.000000e+00" [Sampling is used for the first 'SamplingEpoch' epochs, afterwards, all events are taken for training]
SamplingImportance: "1.000000e+00" [ The sampling weights of events in epochs which successful (worse estimator than before) are multiplied with SamplingImportance, else they are divided.]
SamplingTraining: "True" [The training sample is sampled]
SamplingTesting: "False" [The testing sample is sampled]
ResetStep: "50" [How often BFGS should reset history]
Tau: "3.000000e+00" [LineSearch "size step"]
BPMode: "sequential" [Back-propagation learning mode: sequential or batch]
BatchSize: "-1" [Batch size: number of events/batch, only set if in Batch Mode, -1 for BatchSize=number_of_events]
ConvergenceImprove: "1.000000e-30" [Minimum improvement which counts as improvement (<0 means automatic convergence check is turned off)]
ConvergenceTests: "-1" [Number of steps (without improvement) required for convergence (<0 means automatic convergence check is turned off)]
UpdateLimit: "10000" [Maximum times of regulator update]
CalculateErrors: "False" [Calculates inverse Hessian matrix at the end of the training to be able to calculate the uncertainties of an MVA value]
WeightRange: "1.000000e+00" [Take the events for the estimator calculations from small deviations from the desired value to large deviations only over the weight range]
##


#VAR -*-*-*-*-*-*-*-*-*-*-*-* variables *-*-*-*-*-*-*-*-*-*-*-*-

NVar 6
Btrk1Pt                       Btrk1Pt                       Btrk1Pt                       Btrk1Pt                                                         'F'    [0.200013399124,10.0547304153]
abs(Btrk1Dz1/Btrk1DzError1)   Trk1DCAz                      Trk1DCAz                      Trk1DCAz                                                        'F'    [4.01805991714e-05,8233.625]
abs(Btrk1Dxy1/Btrk1DxyError1) Trk1DCAxy                     Trk1DCAxy                     Trk1DCAxy                                                       'F'    [9.54853530857e-05,227.588607788]
BsvpvDistance/BsvpvDisErr     dls                           dls                           dls                                                             'F'    [2.00039958954,9310.62890625]
Balpha                        Balpha                        Balpha                        Balpha                                                          'F'    [4.76321656606e-05,3.14021730423]
Bchi2cl                       Bchi2cl                       Bchi2cl                       Bchi2cl                                                         'F'    [0.0500016920269,0.999972939491]
NSpec 0


============================================================================ */

#include <vector>
#include <cmath>
#include <string>
#include <iostream>

#ifndef IClassifierReader__def
#define IClassifierReader__def

class IClassifierReader {

 public:

   // constructor
   IClassifierReader() : fStatusIsClean( true ) {}
   virtual ~IClassifierReader() {}

   // return classifier response
   virtual double GetMvaValue( const std::vector<double>& inputValues ) const = 0;

   // returns classifier status
   bool IsStatusClean() const { return fStatusIsClean; }

 protected:

   bool fStatusIsClean;
};

#endif

class ReadMLPBNN2 : public IClassifierReader {

 public:

   // constructor
   ReadMLPBNN2( std::vector<std::string>& theInputVars ) 
      : IClassifierReader(),
        fClassName( "ReadMLPBNN2" ),
        fNvars( 6 ),
        fIsNormalised( false )
   {      
      // the training input variables
      const char* inputVars[] = { "Btrk1Pt", "abs(Btrk1Dz1/Btrk1DzError1)", "abs(Btrk1Dxy1/Btrk1DxyError1)", "BsvpvDistance/BsvpvDisErr", "Balpha", "Bchi2cl" };

      // sanity checks
      if (theInputVars.size() <= 0) {
         std::cout << "Problem in class \"" << fClassName << "\": empty input vector" << std::endl;
         fStatusIsClean = false;
      }

      if (theInputVars.size() != fNvars) {
         std::cout << "Problem in class \"" << fClassName << "\": mismatch in number of input values: "
                   << theInputVars.size() << " != " << fNvars << std::endl;
         fStatusIsClean = false;
      }

      // validate input variables
      for (size_t ivar = 0; ivar < theInputVars.size(); ivar++) {
         if (theInputVars[ivar] != inputVars[ivar]) {
            std::cout << "Problem in class \"" << fClassName << "\": mismatch in input variable names" << std::endl
                      << " for variable [" << ivar << "]: " << theInputVars[ivar].c_str() << " != " << inputVars[ivar] << std::endl;
            fStatusIsClean = false;
         }
      }

      // initialize min and max vectors (for normalisation)
      fVmin[0] = -1;
      fVmax[0] = 1;
      fVmin[1] = -1;
      fVmax[1] = 1;
      fVmin[2] = -1;
      fVmax[2] = 1;
      fVmin[3] = -1;
      fVmax[3] = 1;
      fVmin[4] = -1;
      fVmax[4] = 1;
      fVmin[5] = -1;
      fVmax[5] = 1;

      // initialize input variable types
      fType[0] = 'F';
      fType[1] = 'F';
      fType[2] = 'F';
      fType[3] = 'F';
      fType[4] = 'F';
      fType[5] = 'F';

      // initialize constants
      Initialize();

      // initialize transformation
      InitTransform();
   }

   // destructor
   virtual ~ReadMLPBNN2() {
      Clear(); // method-specific
   }

   // the classifier response
   // "inputValues" is a vector of input values in the same order as the 
   // variables given to the constructor
   double GetMvaValue( const std::vector<double>& inputValues ) const;

 private:

   // method-specific destructor
   void Clear();

   // input variable transformation

   double fMin_1[3][6];
   double fMax_1[3][6];
   void InitTransform_1();
   void Transform_1( std::vector<double> & iv, int sigOrBgd ) const;
   void InitTransform();
   void Transform( std::vector<double> & iv, int sigOrBgd ) const;

   // common member variables
   const char* fClassName;

   const size_t fNvars;
   size_t GetNvar()           const { return fNvars; }
   char   GetType( int ivar ) const { return fType[ivar]; }

   // normalisation of input variables
   const bool fIsNormalised;
   bool IsNormalised() const { return fIsNormalised; }
   double fVmin[6];
   double fVmax[6];
   double NormVariable( double x, double xmin, double xmax ) const {
      // normalise to output range: [-1, 1]
      return 2*(x - xmin)/(xmax - xmin) - 1.0;
   }

   // type of input variable: 'F' or 'I'
   char   fType[6];

   // initialize internal variables
   void Initialize();
   double GetMvaValue__( const std::vector<double>& inputValues ) const;

   // private members (method specific)

   double ActivationFnc(double x) const;
   double OutputActivationFnc(double x) const;

   int fLayers;
   int fLayerSize[4];
   double fWeightMatrix0to1[12][7];   // weight matrix from layer 0 to 1
   double fWeightMatrix1to2[7][12];   // weight matrix from layer 1 to 2
   double fWeightMatrix2to3[1][7];   // weight matrix from layer 2 to 3

   double * fWeights[4];
};

inline void ReadMLPBNN2::Initialize()
{
   // build network structure
   fLayers = 4;
   fLayerSize[0] = 7; fWeights[0] = new double[7]; 
   fLayerSize[1] = 12; fWeights[1] = new double[12]; 
   fLayerSize[2] = 7; fWeights[2] = new double[7]; 
   fLayerSize[3] = 1; fWeights[3] = new double[1]; 
   // weight matrix from layer 0 to 1
   fWeightMatrix0to1[0][0] = -0.430055355739917;
   fWeightMatrix0to1[1][0] = 0.984089766562825;
   fWeightMatrix0to1[2][0] = 0.347869070111069;
   fWeightMatrix0to1[3][0] = 0.440193063302422;
   fWeightMatrix0to1[4][0] = -1.63106732613374;
   fWeightMatrix0to1[5][0] = 0.0517483414777738;
   fWeightMatrix0to1[6][0] = -0.544952578754601;
   fWeightMatrix0to1[7][0] = 0.476073140023487;
   fWeightMatrix0to1[8][0] = -0.192782715082236;
   fWeightMatrix0to1[9][0] = -0.49202112973683;
   fWeightMatrix0to1[10][0] = -2.15135955722438;
   fWeightMatrix0to1[0][1] = -0.629756144845196;
   fWeightMatrix0to1[1][1] = -0.355489418506103;
   fWeightMatrix0to1[2][1] = -0.681452974469694;
   fWeightMatrix0to1[3][1] = -1.22519259970004;
   fWeightMatrix0to1[4][1] = 0.368047825119326;
   fWeightMatrix0to1[5][1] = -0.568830932257474;
   fWeightMatrix0to1[6][1] = 1.55639082786623;
   fWeightMatrix0to1[7][1] = -0.424997972013558;
   fWeightMatrix0to1[8][1] = 1.9959966216849;
   fWeightMatrix0to1[9][1] = -0.456146623769775;
   fWeightMatrix0to1[10][1] = 0.434406595351388;
   fWeightMatrix0to1[0][2] = -0.248175668704181;
   fWeightMatrix0to1[1][2] = 3.99815102793814;
   fWeightMatrix0to1[2][2] = 0.270382904626342;
   fWeightMatrix0to1[3][2] = 0.282149387057768;
   fWeightMatrix0to1[4][2] = 1.31529308696511;
   fWeightMatrix0to1[5][2] = -4.61270758577886;
   fWeightMatrix0to1[6][2] = -0.985608562948398;
   fWeightMatrix0to1[7][2] = 2.31277650711206;
   fWeightMatrix0to1[8][2] = -1.56841138461125;
   fWeightMatrix0to1[9][2] = 1.18980627833826;
   fWeightMatrix0to1[10][2] = -0.0973110219351366;
   fWeightMatrix0to1[0][3] = -0.968244586188514;
   fWeightMatrix0to1[1][3] = 3.49455544011701;
   fWeightMatrix0to1[2][3] = 0.222902153847121;
   fWeightMatrix0to1[3][3] = -0.480675270980661;
   fWeightMatrix0to1[4][3] = 1.6307063746422;
   fWeightMatrix0to1[5][3] = -0.319089619498518;
   fWeightMatrix0to1[6][3] = 0.970452065416955;
   fWeightMatrix0to1[7][3] = 0.0618505806610365;
   fWeightMatrix0to1[8][3] = 2.15455406095873;
   fWeightMatrix0to1[9][3] = 0.847900567972977;
   fWeightMatrix0to1[10][3] = 0.0461513543044774;
   fWeightMatrix0to1[0][4] = -0.579895448633752;
   fWeightMatrix0to1[1][4] = -6.95611320604195;
   fWeightMatrix0to1[2][4] = 0.711313170167701;
   fWeightMatrix0to1[3][4] = -2.51878987687861;
   fWeightMatrix0to1[4][4] = 1.06507228721706;
   fWeightMatrix0to1[5][4] = -0.0586870247411154;
   fWeightMatrix0to1[6][4] = -0.00318328336651011;
   fWeightMatrix0to1[7][4] = -0.466042798973928;
   fWeightMatrix0to1[8][4] = -2.00150265189683;
   fWeightMatrix0to1[9][4] = 1.78643831887096;
   fWeightMatrix0to1[10][4] = -0.466027619332558;
   fWeightMatrix0to1[0][5] = 0.169406683737321;
   fWeightMatrix0to1[1][5] = 0.0544663965533545;
   fWeightMatrix0to1[2][5] = -0.651277657775991;
   fWeightMatrix0to1[3][5] = 0.25796606198843;
   fWeightMatrix0to1[4][5] = -0.0972893898276932;
   fWeightMatrix0to1[5][5] = -0.253807525104342;
   fWeightMatrix0to1[6][5] = 0.443561596632714;
   fWeightMatrix0to1[7][5] = -0.209885498538278;
   fWeightMatrix0to1[8][5] = -0.100883998308796;
   fWeightMatrix0to1[9][5] = 0.324278062104257;
   fWeightMatrix0to1[10][5] = -0.11203604317343;
   fWeightMatrix0to1[0][6] = 0.624543174531467;
   fWeightMatrix0to1[1][6] = -0.144078855161311;
   fWeightMatrix0to1[2][6] = 1.17651871708238;
   fWeightMatrix0to1[3][6] = 1.40935035813282;
   fWeightMatrix0to1[4][6] = -0.94380543179032;
   fWeightMatrix0to1[5][6] = -4.72199632752815;
   fWeightMatrix0to1[6][6] = 0.103997099405775;
   fWeightMatrix0to1[7][6] = 3.0291077621234;
   fWeightMatrix0to1[8][6] = 4.3292684841673;
   fWeightMatrix0to1[9][6] = -0.283351492676293;
   fWeightMatrix0to1[10][6] = -2.62320981877348;
   // weight matrix from layer 1 to 2
   fWeightMatrix1to2[0][0] = -1.53932800450246;
   fWeightMatrix1to2[1][0] = 1.35048054833054;
   fWeightMatrix1to2[2][0] = 0.913053569095608;
   fWeightMatrix1to2[3][0] = 2.42971336344971;
   fWeightMatrix1to2[4][0] = 0.440017745805155;
   fWeightMatrix1to2[5][0] = 1.68481210840193;
   fWeightMatrix1to2[0][1] = -0.63666635763826;
   fWeightMatrix1to2[1][1] = -2.12074407994695;
   fWeightMatrix1to2[2][1] = -0.446791621841251;
   fWeightMatrix1to2[3][1] = 3.10881734526515;
   fWeightMatrix1to2[4][1] = -1.02896275971907;
   fWeightMatrix1to2[5][1] = -2.27278214764847;
   fWeightMatrix1to2[0][2] = -1.40929780659678;
   fWeightMatrix1to2[1][2] = -0.0761837577940548;
   fWeightMatrix1to2[2][2] = 0.848142965694829;
   fWeightMatrix1to2[3][2] = 1.17891575687362;
   fWeightMatrix1to2[4][2] = -1.43006318106675;
   fWeightMatrix1to2[5][2] = 0.0872942902594484;
   fWeightMatrix1to2[0][3] = 1.57656661795336;
   fWeightMatrix1to2[1][3] = 0.135695644849228;
   fWeightMatrix1to2[2][3] = -0.812750852085143;
   fWeightMatrix1to2[3][3] = -0.706581955444624;
   fWeightMatrix1to2[4][3] = 0.121280820265326;
   fWeightMatrix1to2[5][3] = -1.63595094597265;
   fWeightMatrix1to2[0][4] = -1.70726802711202;
   fWeightMatrix1to2[1][4] = 0.269039768946155;
   fWeightMatrix1to2[2][4] = 1.27396511528195;
   fWeightMatrix1to2[3][4] = -2.07781561514283;
   fWeightMatrix1to2[4][4] = 1.7243129199504;
   fWeightMatrix1to2[5][4] = -1.52418741026943;
   fWeightMatrix1to2[0][5] = -0.776503711390803;
   fWeightMatrix1to2[1][5] = 3.08189883363424;
   fWeightMatrix1to2[2][5] = 1.54876909382356;
   fWeightMatrix1to2[3][5] = -4.46464160475799;
   fWeightMatrix1to2[4][5] = 0.492902274665917;
   fWeightMatrix1to2[5][5] = 1.2356099124259;
   fWeightMatrix1to2[0][6] = 0.215159177235193;
   fWeightMatrix1to2[1][6] = 0.930049065238399;
   fWeightMatrix1to2[2][6] = 0.913976996756227;
   fWeightMatrix1to2[3][6] = -0.685499893032718;
   fWeightMatrix1to2[4][6] = -0.708080805755583;
   fWeightMatrix1to2[5][6] = -0.0400648063286739;
   fWeightMatrix1to2[0][7] = -1.09624191169973;
   fWeightMatrix1to2[1][7] = -2.62866954196081;
   fWeightMatrix1to2[2][7] = -0.33730961567201;
   fWeightMatrix1to2[3][7] = 1.22736284506906;
   fWeightMatrix1to2[4][7] = 1.74419743127646;
   fWeightMatrix1to2[5][7] = 0.195371861966852;
   fWeightMatrix1to2[0][8] = 0.253476824417695;
   fWeightMatrix1to2[1][8] = -1.49622465176103;
   fWeightMatrix1to2[2][8] = -0.440729860959914;
   fWeightMatrix1to2[3][8] = -2.86923743945223;
   fWeightMatrix1to2[4][8] = 0.891023651409129;
   fWeightMatrix1to2[5][8] = 0.745701010121103;
   fWeightMatrix1to2[0][9] = -1.84735044615258;
   fWeightMatrix1to2[1][9] = -0.797669266128251;
   fWeightMatrix1to2[2][9] = -1.74310196858477;
   fWeightMatrix1to2[3][9] = -0.665517728681165;
   fWeightMatrix1to2[4][9] = 1.01016968627564;
   fWeightMatrix1to2[5][9] = 0.914535043370396;
   fWeightMatrix1to2[0][10] = -1.57556040216526;
   fWeightMatrix1to2[1][10] = 1.71308045009465;
   fWeightMatrix1to2[2][10] = 1.07755217241104;
   fWeightMatrix1to2[3][10] = 0.64324889049016;
   fWeightMatrix1to2[4][10] = 0.410157610836762;
   fWeightMatrix1to2[5][10] = -1.33341588998146;
   fWeightMatrix1to2[0][11] = 0.0353577453831073;
   fWeightMatrix1to2[1][11] = -1.09194030749053;
   fWeightMatrix1to2[2][11] = 1.63806683315204;
   fWeightMatrix1to2[3][11] = 2.13985008704683;
   fWeightMatrix1to2[4][11] = -0.663308554368684;
   fWeightMatrix1to2[5][11] = 0.296111423301027;
   // weight matrix from layer 2 to 3
   fWeightMatrix2to3[0][0] = 0.256081299761059;
   fWeightMatrix2to3[0][1] = -7.11243548947181;
   fWeightMatrix2to3[0][2] = -1.08099637294068;
   fWeightMatrix2to3[0][3] = 3.38578260372807;
   fWeightMatrix2to3[0][4] = 1.18638723823502;
   fWeightMatrix2to3[0][5] = -2.1889058151459;
   fWeightMatrix2to3[0][6] = 0.092235569427176;
}

inline double ReadMLPBNN2::GetMvaValue__( const std::vector<double>& inputValues ) const
{
   if (inputValues.size() != (unsigned int)fLayerSize[0]-1) {
      std::cout << "Input vector needs to be of size " << fLayerSize[0]-1 << std::endl;
      return 0;
   }

   for (int l=0; l<fLayers; l++)
      for (int i=0; i<fLayerSize[l]; i++) fWeights[l][i]=0;

   for (int l=0; l<fLayers-1; l++)
      fWeights[l][fLayerSize[l]-1]=1;

   for (int i=0; i<fLayerSize[0]-1; i++)
      fWeights[0][i]=inputValues[i];

   // layer 0 to 1
   for (int o=0; o<fLayerSize[1]-1; o++) {
      for (int i=0; i<fLayerSize[0]; i++) {
         double inputVal = fWeightMatrix0to1[o][i] * fWeights[0][i];
         fWeights[1][o] += inputVal;
      }
      fWeights[1][o] = ActivationFnc(fWeights[1][o]);
   }
   // layer 1 to 2
   for (int o=0; o<fLayerSize[2]-1; o++) {
      for (int i=0; i<fLayerSize[1]; i++) {
         double inputVal = fWeightMatrix1to2[o][i] * fWeights[1][i];
         fWeights[2][o] += inputVal;
      }
      fWeights[2][o] = ActivationFnc(fWeights[2][o]);
   }
   // layer 2 to 3
   for (int o=0; o<fLayerSize[3]; o++) {
      for (int i=0; i<fLayerSize[2]; i++) {
         double inputVal = fWeightMatrix2to3[o][i] * fWeights[2][i];
         fWeights[3][o] += inputVal;
      }
      fWeights[3][o] = OutputActivationFnc(fWeights[3][o]);
   }

   return fWeights[3][0];
}

double ReadMLPBNN2::ActivationFnc(double x) const {
   // hyperbolic tan
   return tanh(x);
}
double ReadMLPBNN2::OutputActivationFnc(double x) const {
   // sigmoid
   return 1.0/(1.0+exp(-x));
}
   
// Clean up
inline void ReadMLPBNN2::Clear() 
{
   // clean up the arrays
   for (int lIdx = 0; lIdx < 4; lIdx++) {
      delete[] fWeights[lIdx];
   }
}
   inline double ReadMLPBNN2::GetMvaValue( const std::vector<double>& inputValues ) const
   {
      // classifier response value
      double retval = 0;

      // classifier response, sanity check first
      if (!IsStatusClean()) {
         std::cout << "Problem in class \"" << fClassName << "\": cannot return classifier response"
                   << " because status is dirty" << std::endl;
         retval = 0;
      }
      else {
         if (IsNormalised()) {
            // normalise variables
            std::vector<double> iV;
            iV.reserve(inputValues.size());
            int ivar = 0;
            for (std::vector<double>::const_iterator varIt = inputValues.begin();
                 varIt != inputValues.end(); varIt++, ivar++) {
               iV.push_back(NormVariable( *varIt, fVmin[ivar], fVmax[ivar] ));
            }
            Transform( iV, -1 );
            retval = GetMvaValue__( iV );
         }
         else {
            std::vector<double> iV;
            int ivar = 0;
            for (std::vector<double>::const_iterator varIt = inputValues.begin();
                 varIt != inputValues.end(); varIt++, ivar++) {
               iV.push_back(*varIt);
            }
            Transform( iV, -1 );
            retval = GetMvaValue__( iV );
         }
      }

      return retval;
   }

//_______________________________________________________________________
inline void ReadMLPBNN2::InitTransform_1()
{
   // Normalization transformation, initialisation
   fMin_1[0][0] = 0.200558498502;
   fMax_1[0][0] = 10.0547304153;
   fMin_1[1][0] = 0.200013399124;
   fMax_1[1][0] = 9.56284618378;
   fMin_1[2][0] = 0.200013399124;
   fMax_1[2][0] = 10.0547304153;
   fMin_1[0][1] = 4.01805991714e-05;
   fMax_1[0][1] = 8233.625;
   fMin_1[1][1] = 0.000844108290039;
   fMax_1[1][1] = 3256.72607422;
   fMin_1[2][1] = 4.01805991714e-05;
   fMax_1[2][1] = 8233.625;
   fMin_1[0][2] = 0.000327922782162;
   fMax_1[0][2] = 227.588607788;
   fMin_1[1][2] = 9.54853530857e-05;
   fMax_1[1][2] = 31.6222686768;
   fMin_1[2][2] = 9.54853530857e-05;
   fMax_1[2][2] = 227.588607788;
   fMin_1[0][3] = 2.00099754333;
   fMax_1[0][3] = 9310.62890625;
   fMin_1[1][3] = 2.00039958954;
   fMax_1[1][3] = 7566.11328125;
   fMin_1[2][3] = 2.00039958954;
   fMax_1[2][3] = 9310.62890625;
   fMin_1[0][4] = 4.76321656606e-05;
   fMax_1[0][4] = 3.09399151802;
   fMin_1[1][4] = 0.000779491965659;
   fMax_1[1][4] = 3.14021730423;
   fMin_1[2][4] = 4.76321656606e-05;
   fMax_1[2][4] = 3.14021730423;
   fMin_1[0][5] = 0.0500016920269;
   fMax_1[0][5] = 0.999972939491;
   fMin_1[1][5] = 0.0500183738768;
   fMax_1[1][5] = 0.999963581562;
   fMin_1[2][5] = 0.0500016920269;
   fMax_1[2][5] = 0.999972939491;
}

//_______________________________________________________________________
inline void ReadMLPBNN2::Transform_1( std::vector<double>& iv, int cls) const
{
   // Normalization transformation
   if (cls < 0 || cls > 2) {
   if (2 > 1 ) cls = 2;
      else cls = 2;
   }
   const int nVar = 6;

   // get indices of used variables

   // define the indices of the variables which are transformed by this transformation
   static std::vector<int> indicesGet;
   static std::vector<int> indicesPut;

   if ( indicesGet.empty() ) { 
      indicesGet.reserve(fNvars);
      indicesGet.push_back( 0);
      indicesGet.push_back( 1);
      indicesGet.push_back( 2);
      indicesGet.push_back( 3);
      indicesGet.push_back( 4);
      indicesGet.push_back( 5);
   } 
   if ( indicesPut.empty() ) { 
      indicesPut.reserve(fNvars);
      indicesPut.push_back( 0);
      indicesPut.push_back( 1);
      indicesPut.push_back( 2);
      indicesPut.push_back( 3);
      indicesPut.push_back( 4);
      indicesPut.push_back( 5);
   } 

   static std::vector<double> dv;
   dv.resize(nVar);
   for (int ivar=0; ivar<nVar; ivar++) dv[ivar] = iv[indicesGet.at(ivar)];
   for (int ivar=0;ivar<6;ivar++) {
      double offset = fMin_1[cls][ivar];
      double scale  = 1.0/(fMax_1[cls][ivar]-fMin_1[cls][ivar]);
      iv[indicesPut.at(ivar)] = (dv[ivar]-offset)*scale * 2 - 1;
   }
}

//_______________________________________________________________________
inline void ReadMLPBNN2::InitTransform()
{
   InitTransform_1();
}

//_______________________________________________________________________
inline void ReadMLPBNN2::Transform( std::vector<double>& iv, int sigOrBgd ) const
{
   Transform_1( iv, sigOrBgd );
}
