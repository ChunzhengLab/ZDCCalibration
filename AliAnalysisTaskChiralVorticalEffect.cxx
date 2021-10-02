#include <iostream>
#include <cstdlib>
#include <sys/time.h>
#include <vector>
#include <algorithm>
// ROOT classes
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TList.h"
#include "TH1I.h"
#include "TH3I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TExMap.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TComplex.h"
#include "THnBase.h"
#include "THn.h"
// Alice analysis base class
#include "AliAnalysisTaskSE.h"
// Alice analysis additional classes
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
// Alice AOD classes
#include "AliAODInputHandler.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODVZERO.h"
#include "AliAODZDC.h"

#include "AliAODv0.h"
#include "AliAODTrack.h"
#include "AliAODHeader.h"

// Alice classes
#include "AliCentrality.h"
#include "AliEventplane.h"
//#include "AliEventCuts.h"
#include "AliAnalysisUtils.h"
// Alice MC classes
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliAODMCParticle.h"
// Alice "V" classes
#include "AliVParticle.h"
#include "AliVEvent.h"
#include "AliVVertex.h"
#include "AliVVZERO.h"
// Alice PID classes
#include "AliAODPid.h"
#include "AliAODpidUtil.h"
#include "AliPID.h"
#include "AliPIDCombined.h"
#include "AliPIDResponse.h"
//#include "AliMultSelection.h"

#include "AliAnalysisTaskChiralVorticalEffect.h"

class AliAnalysisTaskMyTask; 

using std::cout;
using std::endl;
using std::vector;

using namespace std;

ClassImp(AliAnalysisTaskChiralVorticalEffect);

//---------------------------------------------------
AliAnalysisTaskChiralVorticalEffect::AliAnalysisTaskChiralVorticalEffect() : AliAnalysisTaskSE(),
fAOD(0), fOutputList(0),
mHarmonic(2.),
fFilterBit(1),
fPtMin(0.2),
fPtMax(5.),
fEtaMax(0.8),
fNhitsMin(80),
fChi2Max(4.),
fDeDxMin(10),
fDcaXyMax(3.),
fDcaZMax(3.),

fV0CPAMin(0.995),
fV0DCAToPrimVtxMax(1.5),
fV0DecayLengthMax(100.),
fV0DecayLengthMin(3.),
fV0DcaBetweenDaughtersMax(1.),
fV0PtMin(0.5),
fV0RapidityMax(0.5),

fDaughtersPtMax(20.),
fDaughtersEtaMax(0.8),
fDaughtersTPCNclsMin(70),
fDaughtersDCAToPrimVtxMin(0.02),
fDaughtersNsigma(3.),

fMassMean(1.115683),
fLambdaMassCut(0.005)
{
  runNum       = -999;
  runNumBin    = -999;
  for (int i=0; i<3; ++i) vtx[i] = -999;
  vzBin        = -999;
  cent         = -999;
  centBin      = -999;

  fOutputList = NULL;
  // Event-wise
  hEvtCount  = NULL;
  hRunNumBin = NULL;
  hCent      = NULL;
  for (int i=0; i<2; ++i) hCentCorr[i] = NULL;
  for (int i=0; i<2; ++i) hVxy[i]      = NULL;
  for (int i=0; i<2; ++i) hVz[i]       = NULL;

  //Random Event Plane
  hPsi2RDM = NULL;
  hPsi3RDM = NULL;

  //TPC Event Plane
  hPsi2TPC = NULL;
  hPsi3TPC = NULL;

  // Track-wise
  for (int i=0; i<2; ++i) hPt[i]    = NULL;
  for (int i=0; i<2; ++i) hEta[i]   = NULL;
  for (int i=0; i<2; ++i) hPhi[i]   = NULL;
  for (int i=0; i<2; ++i) hNhits[i] = NULL;
  for (int i=0; i<2; ++i) hDcaXy[i] = NULL;
  for (int i=0; i<2; ++i) hDcaZ[i]  = NULL;
  hPDedx = NULL;

  //V0-wise
  hV0Pt = NULL;
  hV0Eta = NULL;
  hV0DcaToPrimVertex = NULL;
  hV0CPA = NULL;
  hV0DecayLength = NULL;

  for (int i=0; i<2; ++i) hLambdaPt[i] = NULL;
  for (int i=0; i<2; ++i) hLambdaEta[i] = NULL;
  for (int i=0; i<2; ++i) hLambdaDcaToPrimVertex[i] = NULL;
  for (int i=0; i<2; ++i) hLambdaCPA[i] = NULL;
  for (int i=0; i<2; ++i) hLambdaDecayLength[i] = NULL;
  for (int i=0; i<2; ++i) hLambdaMass[i] = NULL;
  for (int i=0; i<10; ++i) hLambdaMassCent[i] = NULL;

  for (int i=0; i<2; ++i) hAntiLambdaPt[i] = NULL;
  for (int i=0; i<2; ++i) hAntiLambdaEta[i] = NULL;
  for (int i=0; i<2; ++i) hAntiLambdaDcaToPrimVertex[i] = NULL;
  for (int i=0; i<2; ++i) hAntiLambdaCPA[i] = NULL;
  for (int i=0; i<2; ++i) hAntiLambdaDecayLength[i] = NULL;
  for (int i=0; i<2; ++i) hAntiLambdaMass[i] = NULL;
  for (int i=0; i<10; ++i) hAntiLambdaMassCent[i] = NULL;

  //Inclusive
  //Delta
  pDelta_hPos_hPos = NULL;
  pDelta_hNeg_hNeg = NULL;
  pDelta_hPos_hNeg = NULL;
  //Gamma
  //TPC Plane
  pGammaTPC_hPos_hPos = NULL;
  pGammaTPC_hNeg_hNeg = NULL;
  pGammaTPC_hPos_hNeg = NULL;
  pCosCosTPC_hPos_hPos = NULL;
  pCosCosTPC_hNeg_hNeg = NULL;
  pCosCosTPC_hPos_hNeg = NULL;
  pSinSinTPC_hPos_hPos = NULL;
  pSinSinTPC_hNeg_hNeg = NULL;
  pSinSinTPC_hPos_hNeg = NULL;
  //Random Plane 
  pGammaRDM_hPos_hPos = NULL;
  pGammaRDM_hNeg_hNeg = NULL;
  pGammaRDM_hPos_hNeg = NULL;
  pCosCosRDM_hPos_hPos = NULL;
  pCosCosRDM_hNeg_hNeg = NULL;
  pCosCosRDM_hPos_hNeg = NULL;
  pSinSinRDM_hPos_hPos = NULL;
  pSinSinRDM_hNeg_hNeg = NULL;
  pSinSinRDM_hPos_hNeg = NULL;
  //Diff
  for(int i=0;i<10;i++) pGammaTPCVsDeltaPt_hPos_hPos[i] = NULL;
  for(int i=0;i<10;i++) pGammaTPCVsDeltaPt_hPos_hNeg[i] = NULL;
  for(int i=0;i<10;i++) pGammaTPCVsDeltaPt_hNeg_hNeg[i] = NULL;

  for(int i=0;i<10;i++) pGammaTPVVsMeanPt_hPos_hPos[i] = NULL;
  for(int i=0;i<10;i++) pGammaTPVVsMeanPt_hPos_hNeg[i] = NULL;
  for(int i=0;i<10;i++) pGammaTPVVsMeanPt_hNeg_hNeg[i] = NULL;

  for(int i=0;i<10;i++) pGammaTPCVsDeltaEta_hPos_hPos[i] = NULL;
  for(int i=0;i<10;i++) pGammaTPCVsDeltaEta_hPos_hNeg[i] = NULL;
  for(int i=0;i<10;i++) pGammaTPCVsDeltaEta_hNeg_hNeg[i] = NULL;


  //pion - Inclusive
  // pi+ - h+
  pDelta_pion_hPos = NULL;
  pGammaTPC_pion_hPos = NULL;
  pGammaRDM_pion_hPos = NULL;
  for(int i=0;i<10;i++) pGammaTPCVsPt_pion_hPos[i] = NULL;
  // pi- - h+
  pDelta_antiPion_hPos = NULL;
  pGammaTPC_antiPion_hPos = NULL;
  pGammaRDM_antiPion_hPos = NULL;
  for(int i=0;i<10;i++) pGammaTPCVsPt_antiPion_hPos[i] = NULL;
  // pi+ - h-
  pDelta_pion_hNeg = NULL;
  pGammaTPC_pion_hNeg = NULL;
  pGammaRDM_pion_hNeg = NULL;
  for(int i=0;i<10;i++) pGammaTPCVsPt_pion_hNeg[i] = NULL;
  // pi- - h-
  pDelta_antiPion_hNeg = NULL;
  pGammaTPC_antiPion_hNeg = NULL;
  pGammaRDM_antiPion_hNeg = NULL;
  for(int i=0;i<10;i++) pGammaTPCVsPt_antiPion_hNeg[i] = NULL;
  //kaon - Inclusive
  // ka+ - h+
  pDelta_kaon_hPos = NULL;
  pGammaTPC_kaon_hPos = NULL;
  pGammaRDM_kaon_hPos = NULL;
  for(int i=0;i<10;i++) pGammaTPCVsPt_kaon_hPos[i] = NULL;
  // ka- - h+
  pDelta_antiKaon_hPos = NULL;
  pGammaTPC_antiKaon_hPos = NULL;
  pGammaRDM_antiKaon_hPos = NULL;
  for(int i=0;i<10;i++) pGammaTPCVsPt_antiKaon_hPos[i] = NULL;
  // ka+ - h-
  pDelta_kaon_hNeg = NULL;
  pGammaTPC_kaon_hNeg = NULL;
  pGammaRDM_kaon_hNeg = NULL;
  for(int i=0;i<10;i++) pGammaTPCVsPt_kaon_hNeg[i] = NULL;
  // ka- - h-
  pDelta_antiKaon_hNeg = NULL;
  pGammaTPC_antiKaon_hNeg = NULL;
  pGammaRDM_antiKaon_hNeg = NULL;
  for(int i=0;i<10;i++) pGammaTPCVsPt_antiKaon_hNeg[i] = NULL;
  //proton - Inclusive
  // p+ - h+
  pDelta_proton_hPos = NULL;
  pGammaTPC_proton_hPos = NULL;
  pGammaRDM_proton_hPos = NULL;
  for(int i=0;i<10;i++) pGammaTPCVsPt_proton_hPos[i] = NULL;
  // p- - h+
  pDelta_antiProton_hPos = NULL;
  pGammaTPC_antiProton_hPos = NULL;
  pGammaRDM_antiProton_hPos = NULL;
  for(int i=0;i<10;i++) pGammaTPCVsPt_antiProton_hPos[i] = NULL;
  // p+ - h-
  pDelta_proton_hNeg = NULL;
  pGammaTPC_proton_hNeg = NULL;
  pGammaRDM_proton_hNeg = NULL;
  for(int i=0;i<10;i++) pGammaTPCVsPt_proton_hNeg[i] = NULL;
  // p- - h-
  pDelta_antiProton_hNeg = NULL;
  pGammaTPC_antiProton_hNeg = NULL;
  pGammaRDM_antiProton_hNeg = NULL;
  for(int i=0;i<10;i++) pGammaTPCVsPt_antiProton_hNeg[i] = NULL;


  //pion - pion
  //pi+ - pi+
  pDelta_pion_pion = NULL;
  pGammaTPC_pion_pion = NULL;
  pGammaRDM_pion_pion = NULL;
  //pi- - pi-
  pDelta_antiPion_antiPion = NULL;
  pGammaTPC_antiPion_antiPion = NULL;
  pGammaRDM_antiPion_antiPion = NULL;
  //pi+ - pi-
  pDelta_pion_antiPion = NULL;
  pGammaTPC_pion_antiPion = NULL;
  pGammaRDM_pion_antiPion = NULL;
  //kaon - kaon
  //ka+ - ka+
  pDelta_kaon_kaon = NULL;
  pGammaTPC_kaon_kaon = NULL;
  pGammaRDM_kaon_kaon = NULL;
  //ka- - ka-
  pDelta_antiKaon_antiKaon = NULL;
  pGammaTPC_antiKaon_antiKaon = NULL;
  pGammaRDM_antiKaon_antiKaon = NULL;
  //ka+ - ka-
  pDelta_kaon_antiKaon = NULL;
  pGammaTPC_kaon_antiKaon = NULL;
  pGammaRDM_kaon_antiKaon = NULL;
  //proton - proton
  //p+ - p+
  pDelta_proton_proton = NULL;
  pGammaTPC_proton_proton = NULL;
  pGammaRDM_proton_proton = NULL;
  //p- - p-
  pDelta_antiProton_antiProton = NULL;
  pGammaTPC_antiProton_antiProton = NULL;
  pGammaRDM_antiProton_antiProton = NULL;
  //p+ - p-
  pDelta_proton_antiProton = NULL;
  pGammaTPC_proton_antiProton = NULL;
  pGammaRDM_proton_antiProton = NULL;


  //pion - kaon
  // pi+ - k+
  pDelta_pion_kaon = NULL;
  pGammaTPC_pion_kaon = NULL;
  pGammaRDM_pion_kaon = NULL;
  // pi- - k+
  pDelta_antiPion_kaon = NULL;
  pGammaTPC_antiPion_kaon = NULL;
  pGammaRDM_antiPion_kaon = NULL;
  // pi+ - k-
  pDelta_pion_antiKaon = NULL;
  pGammaTPC_pion_antiKaon = NULL;
  pGammaRDM_pion_antiKaon = NULL;
  // pi- - k-
  pDelta_antiPion_antiKaon = NULL;
  pGammaTPC_antiPion_antiKaon = NULL;
  pGammaRDM_antiPion_antiKaon = NULL;
  //kaon - proton
  // k+ - p+
  pDelta_kaon_proton = NULL;
  pGammaTPC_kaon_proton = NULL;
  pGammaRDM_kaon_proton = NULL;
  // k- - p+
  pDelta_antiKaon_proton = NULL;
  pGammaTPC_antiKaon_proton = NULL;
  pGammaRDM_antiKaon_proton = NULL;
  // k+ - p-
  pDelta_kaon_antiProton = NULL;
  pGammaTPC_kaon_antiProton = NULL;
  pGammaRDM_kaon_antiProton = NULL;
  // k- - p-
  pDelta_antiKaon_antiProton = NULL;
  pGammaTPC_antiKaon_antiProton = NULL;
  pGammaRDM_antiKaon_antiProton = NULL;
  //proton - pi
  // p+ - pi+
  pDelta_proton_pion = NULL;
  pGammaTPC_proton_pion = NULL;
  pGammaRDM_proton_pion = NULL;
  // p- - pi+
  pDelta_antiProton_pion = NULL;
  pGammaTPC_antiProton_pion = NULL;
  pGammaRDM_antiProton_pion = NULL;
  // p+ - pi-
  pDelta_proton_antiPion = NULL;
  pGammaTPC_proton_antiPion = NULL;
  pGammaRDM_proton_antiPion = NULL;
  // p- - pi-
  pDelta_antiProton_antiPion = NULL;
  pGammaTPC_antiProton_antiPion = NULL;
  pGammaRDM_antiProton_antiPion = NULL;

  //lambda - h
  //lambda - h+
  pDelta_lambda_hPos = NULL;
  pGammaTPC_lambda_hPos = NULL;
  pGammaRDM_lambda_hPos = NULL;
  //antiLambda - h+
  pDelta_antiLambda_hPos = NULL;
  pGammaTPC_antiLambda_hPos = NULL;
  pGammaRDM_antiLambda_hPos = NULL;
  //lambda - h-
  pDelta_lambda_hNeg = NULL;
  pGammaTPC_lambda_hNeg = NULL;
  pGammaRDM_lambda_hNeg = NULL;
  //antiLambda - h-
  pDelta_antiLambda_hNeg = NULL;
  pGammaTPC_antiLambda_hNeg = NULL;
  pGammaRDM_antiLambda_hNeg = NULL;

  //lambda - p
  //lambda - p+
  pDelta_lambda_proton = NULL;
  pGammaTPC_lambda_proton = NULL;
  pGammaRDM_lambda_proton = NULL;
  //antiLambda - p+
  pDelta_antiLambda_proton = NULL;
  pGammaTPC_antiLambda_proton = NULL;
  pGammaRDM_antiLambda_proton = NULL;
  //lambda - p-
  pDelta_lambda_antiProton = NULL;
  pGammaTPC_lambda_antiProton = NULL;
  pGammaRDM_lambda_antiProton = NULL;
  //antiLambda - p-
  pDelta_antiLambda_antiProton = NULL;
  pGammaTPC_antiLambda_antiProton = NULL;
  pGammaRDM_antiLambda_antiProton = NULL;
}

//---------------------------------------------------
AliAnalysisTaskChiralVorticalEffect::AliAnalysisTaskChiralVorticalEffect(const char *name) : AliAnalysisTaskSE(name),
fAOD(0), fOutputList(0),
mHarmonic(2.),
fFilterBit(1),
fPtMin(0.2),
fPtMax(5.),
fEtaMax(0.8),
fNhitsMin(80),
fChi2Max(4.),
fDeDxMin(10),
fDcaXyMax(3.),
fDcaZMax(3.),

fV0CPAMin(0.995),
fV0DCAToPrimVtxMax(1.5),
fV0DecayLengthMax(100.),
fV0DecayLengthMin(3.),
fV0DcaBetweenDaughtersMax(1.),
fV0PtMin(0.5),
fV0RapidityMax(0.5),

fDaughtersPtMax(20.),
fDaughtersEtaMax(0.8),
fDaughtersTPCNclsMin(70),
fDaughtersDCAToPrimVtxMin(0.02),
fDaughtersNsigma(3.),

fMassMean(1.115683),
fLambdaMassCut(0.005)
{
  runNum       = -999;
  runNumBin    = -999;
  for (int i=0; i<3; ++i) vtx[i] = -999;
  vzBin        = -999;
  cent         = -999;
  centBin      = -999;

  fOutputList = NULL;
  // Event-wise
  hEvtCount  = NULL;
  hRunNumBin = NULL;
  hCent      = NULL;
  for (int i=0; i<2; ++i) hCentCorr[i] = NULL;
  for (int i=0; i<2; ++i) hVxy[i]      = NULL;
  for (int i=0; i<2; ++i) hVz[i]       = NULL;

  //Random Event Plane
  hPsi2RDM = NULL;
  hPsi3RDM = NULL;

  //TPC Event Plane
  hPsi2TPC = NULL;
  hPsi3TPC = NULL;

  // Track-wise
  for (int i=0; i<2; ++i) hPt[i]    = NULL;
  for (int i=0; i<2; ++i) hEta[i]   = NULL;
  for (int i=0; i<2; ++i) hPhi[i]   = NULL;
  for (int i=0; i<2; ++i) hNhits[i] = NULL;
  for (int i=0; i<2; ++i) hDcaXy[i] = NULL;
  for (int i=0; i<2; ++i) hDcaZ[i]  = NULL;
  hPDedx = NULL;

  //V0-wise
  hV0Pt = NULL;
  hV0Eta = NULL;
  hV0DcaToPrimVertex = NULL;
  hV0CPA = NULL;
  hV0DecayLength = NULL;

  for (int i=0; i<2; ++i) hLambdaPt[i] = NULL;
  for (int i=0; i<2; ++i) hLambdaEta[i] = NULL;
  for (int i=0; i<2; ++i) hLambdaDcaToPrimVertex[i] = NULL;
  for (int i=0; i<2; ++i) hLambdaCPA[i] = NULL;
  for (int i=0; i<2; ++i) hLambdaDecayLength[i] = NULL;
  for (int i=0; i<2; ++i) hLambdaMass[i] = NULL;
  for (int i=0; i<10; ++i) hLambdaMassCent[i] = NULL;

  for (int i=0; i<2; ++i) hAntiLambdaPt[i] = NULL;
  for (int i=0; i<2; ++i) hAntiLambdaEta[i] = NULL;
  for (int i=0; i<2; ++i) hAntiLambdaDcaToPrimVertex[i] = NULL;
  for (int i=0; i<2; ++i) hAntiLambdaCPA[i] = NULL;
  for (int i=0; i<2; ++i) hAntiLambdaDecayLength[i] = NULL;
  for (int i=0; i<2; ++i) hAntiLambdaMass[i] = NULL;
  for (int i=0; i<10; ++i) hAntiLambdaMassCent[i] = NULL;

  //Inclusive
  //Delta
  pDelta_hPos_hPos = NULL;
  pDelta_hNeg_hNeg = NULL;
  pDelta_hPos_hNeg = NULL;
  //Gamma
  //TPC Plane
  pGammaTPC_hPos_hPos = NULL;
  pGammaTPC_hNeg_hNeg = NULL;
  pGammaTPC_hPos_hNeg = NULL;
  pCosCosTPC_hPos_hPos = NULL;
  pCosCosTPC_hNeg_hNeg = NULL;
  pCosCosTPC_hPos_hNeg = NULL;
  pSinSinTPC_hPos_hPos = NULL;
  pSinSinTPC_hNeg_hNeg = NULL;
  pSinSinTPC_hPos_hNeg = NULL;
  //Random Plane 
  pGammaRDM_hPos_hPos = NULL;
  pGammaRDM_hNeg_hNeg = NULL;
  pGammaRDM_hPos_hNeg = NULL;
  pCosCosRDM_hPos_hPos = NULL;
  pCosCosRDM_hNeg_hNeg = NULL;
  pCosCosRDM_hPos_hNeg = NULL;
  pSinSinRDM_hPos_hPos = NULL;
  pSinSinRDM_hNeg_hNeg = NULL;
  pSinSinRDM_hPos_hNeg = NULL;
  //Diff
  for(int i=0;i<10;i++) pGammaTPCVsDeltaPt_hPos_hPos[i] = NULL;
  for(int i=0;i<10;i++) pGammaTPCVsDeltaPt_hPos_hNeg[i] = NULL;
  for(int i=0;i<10;i++) pGammaTPCVsDeltaPt_hNeg_hNeg[i] = NULL;

  for(int i=0;i<10;i++) pGammaTPVVsMeanPt_hPos_hPos[i] = NULL;
  for(int i=0;i<10;i++) pGammaTPVVsMeanPt_hPos_hNeg[i] = NULL;
  for(int i=0;i<10;i++) pGammaTPVVsMeanPt_hNeg_hNeg[i] = NULL;

  for(int i=0;i<10;i++) pGammaTPCVsDeltaEta_hPos_hPos[i] = NULL;
  for(int i=0;i<10;i++) pGammaTPCVsDeltaEta_hPos_hNeg[i] = NULL;
  for(int i=0;i<10;i++) pGammaTPCVsDeltaEta_hNeg_hNeg[i] = NULL;


  //pion - Inclusive
  // pi+ - h+
  pDelta_pion_hPos = NULL;
  pGammaTPC_pion_hPos = NULL;
  pGammaRDM_pion_hPos = NULL;
  for(int i=0;i<10;i++) pGammaTPCVsPt_pion_hPos[i] = NULL;
  // pi- - h+
  pDelta_antiPion_hPos = NULL;
  pGammaTPC_antiPion_hPos = NULL;
  pGammaRDM_antiPion_hPos = NULL;
  for(int i=0;i<10;i++) pGammaTPCVsPt_antiPion_hPos[i] = NULL;
  // pi+ - h-
  pDelta_pion_hNeg = NULL;
  pGammaTPC_pion_hNeg = NULL;
  pGammaRDM_pion_hNeg = NULL;
  for(int i=0;i<10;i++) pGammaTPCVsPt_pion_hNeg[i] = NULL;
  // pi- - h-
  pDelta_antiPion_hNeg = NULL;
  pGammaTPC_antiPion_hNeg = NULL;
  pGammaRDM_antiPion_hNeg = NULL;
  for(int i=0;i<10;i++) pGammaTPCVsPt_antiPion_hNeg[i] = NULL;
  //kaon - Inclusive
  // ka+ - h+
  pDelta_kaon_hPos = NULL;
  pGammaTPC_kaon_hPos = NULL;
  pGammaRDM_kaon_hPos = NULL;
  for(int i=0;i<10;i++) pGammaTPCVsPt_kaon_hPos[i] = NULL;
  // ka- - h+
  pDelta_antiKaon_hPos = NULL;
  pGammaTPC_antiKaon_hPos = NULL;
  pGammaRDM_antiKaon_hPos = NULL;
  for(int i=0;i<10;i++) pGammaTPCVsPt_antiKaon_hPos[i] = NULL;
  // ka+ - h-
  pDelta_kaon_hNeg = NULL;
  pGammaTPC_kaon_hNeg = NULL;
  pGammaRDM_kaon_hNeg = NULL;
  for(int i=0;i<10;i++) pGammaTPCVsPt_kaon_hNeg[i] = NULL;
  // ka- - h-
  pDelta_antiKaon_hNeg = NULL;
  pGammaTPC_antiKaon_hNeg = NULL;
  pGammaRDM_antiKaon_hNeg = NULL;
  for(int i=0;i<10;i++) pGammaTPCVsPt_antiKaon_hNeg[i] = NULL;
  //proton - Inclusive
  // p+ - h+
  pDelta_proton_hPos = NULL;
  pGammaTPC_proton_hPos = NULL;
  pGammaRDM_proton_hPos = NULL;
  for(int i=0;i<10;i++) pGammaTPCVsPt_proton_hPos[i] = NULL;
  // p- - h+
  pDelta_antiProton_hPos = NULL;
  pGammaTPC_antiProton_hPos = NULL;
  pGammaRDM_antiProton_hPos = NULL;
  for(int i=0;i<10;i++) pGammaTPCVsPt_antiProton_hPos[i] = NULL;
  // p+ - h-
  pDelta_proton_hNeg = NULL;
  pGammaTPC_proton_hNeg = NULL;
  pGammaRDM_proton_hNeg = NULL;
  for(int i=0;i<10;i++) pGammaTPCVsPt_proton_hNeg[i] = NULL;
  // p- - h-
  pDelta_antiProton_hNeg = NULL;
  pGammaTPC_antiProton_hNeg = NULL;
  pGammaRDM_antiProton_hNeg = NULL;
  for(int i=0;i<10;i++) pGammaTPCVsPt_antiProton_hNeg[i] = NULL;


  //pion - pion
  //pi+ - pi+
  pDelta_pion_pion = NULL;
  pGammaTPC_pion_pion = NULL;
  pGammaRDM_pion_pion = NULL;
  //pi- - pi-
  pDelta_antiPion_antiPion = NULL;
  pGammaTPC_antiPion_antiPion = NULL;
  pGammaRDM_antiPion_antiPion = NULL;
  //pi+ - pi-
  pDelta_pion_antiPion = NULL;
  pGammaTPC_pion_antiPion = NULL;
  pGammaRDM_pion_antiPion = NULL;
  //kaon - kaon
  //ka+ - ka+
  pDelta_kaon_kaon = NULL;
  pGammaTPC_kaon_kaon = NULL;
  pGammaRDM_kaon_kaon = NULL;
  //ka- - ka-
  pDelta_antiKaon_antiKaon = NULL;
  pGammaTPC_antiKaon_antiKaon = NULL;
  pGammaRDM_antiKaon_antiKaon = NULL;
  //ka+ - ka-
  pDelta_kaon_antiKaon = NULL;
  pGammaTPC_kaon_antiKaon = NULL;
  pGammaRDM_kaon_antiKaon = NULL;
  //proton - proton
  //p+ - p+
  pDelta_proton_proton = NULL;
  pGammaTPC_proton_proton = NULL;
  pGammaRDM_proton_proton = NULL;
  //p- - p-
  pDelta_antiProton_antiProton = NULL;
  pGammaTPC_antiProton_antiProton = NULL;
  pGammaRDM_antiProton_antiProton = NULL;
  //p+ - p-
  pDelta_proton_antiProton = NULL;
  pGammaTPC_proton_antiProton = NULL;
  pGammaRDM_proton_antiProton = NULL;


  //pion - kaon
  // pi+ - k+
  pDelta_pion_kaon = NULL;
  pGammaTPC_pion_kaon = NULL;
  pGammaRDM_pion_kaon = NULL;
  // pi- - k+
  pDelta_antiPion_kaon = NULL;
  pGammaTPC_antiPion_kaon = NULL;
  pGammaRDM_antiPion_kaon = NULL;
  // pi+ - k-
  pDelta_pion_antiKaon = NULL;
  pGammaTPC_pion_antiKaon = NULL;
  pGammaRDM_pion_antiKaon = NULL;
  // pi- - k-
  pDelta_antiPion_antiKaon = NULL;
  pGammaTPC_antiPion_antiKaon = NULL;
  pGammaRDM_antiPion_antiKaon = NULL;
  //kaon - proton
  // k+ - p+
  pDelta_kaon_proton = NULL;
  pGammaTPC_kaon_proton = NULL;
  pGammaRDM_kaon_proton = NULL;
  // k- - p+
  pDelta_antiKaon_proton = NULL;
  pGammaTPC_antiKaon_proton = NULL;
  pGammaRDM_antiKaon_proton = NULL;
  // k+ - p-
  pDelta_kaon_antiProton = NULL;
  pGammaTPC_kaon_antiProton = NULL;
  pGammaRDM_kaon_antiProton = NULL;
  // k- - p-
  pDelta_antiKaon_antiProton = NULL;
  pGammaTPC_antiKaon_antiProton = NULL;
  pGammaRDM_antiKaon_antiProton = NULL;
  //proton - pi
  // p+ - pi+
  pDelta_proton_pion = NULL;
  pGammaTPC_proton_pion = NULL;
  pGammaRDM_proton_pion = NULL;
  // p- - pi+
  pDelta_antiProton_pion = NULL;
  pGammaTPC_antiProton_pion = NULL;
  pGammaRDM_antiProton_pion = NULL;
  // p+ - pi-
  pDelta_proton_antiPion = NULL;
  pGammaTPC_proton_antiPion = NULL;
  pGammaRDM_proton_antiPion = NULL;
  // p- - pi-
  pDelta_antiProton_antiPion = NULL;
  pGammaTPC_antiProton_antiPion = NULL;
  pGammaRDM_antiProton_antiPion = NULL;

  //lambda - h
  //lambda - h+
  pDelta_lambda_hPos = NULL;
  pGammaTPC_lambda_hPos = NULL;
  pGammaRDM_lambda_hPos = NULL;
  //antiLambda - h+
  pDelta_antiLambda_hPos = NULL;
  pGammaTPC_antiLambda_hPos = NULL;
  pGammaRDM_antiLambda_hPos = NULL;
  //lambda - h-
  pDelta_lambda_hNeg = NULL;
  pGammaTPC_lambda_hNeg = NULL;
  pGammaRDM_lambda_hNeg = NULL;
  //antiLambda - h-
  pDelta_antiLambda_hNeg = NULL;
  pGammaTPC_antiLambda_hNeg = NULL;
  pGammaRDM_antiLambda_hNeg = NULL;

  //lambda - p
  //lambda - p+
  pDelta_lambda_proton = NULL;
  pGammaTPC_lambda_proton = NULL;
  pGammaRDM_lambda_proton = NULL;
  //antiLambda - p+
  pDelta_antiLambda_proton = NULL;
  pGammaTPC_antiLambda_proton = NULL;
  pGammaRDM_antiLambda_proton = NULL;
  //lambda - p-
  pDelta_lambda_antiProton = NULL;
  pGammaTPC_lambda_antiProton = NULL;
  pGammaRDM_lambda_antiProton = NULL;
  //antiLambda - p-
  pDelta_antiLambda_antiProton = NULL;
  pGammaTPC_antiLambda_antiProton = NULL;
  pGammaRDM_antiLambda_antiProton = NULL;

  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
}

//---------------------------------------------------
AliAnalysisTaskChiralVorticalEffect::~AliAnalysisTaskChiralVorticalEffect()
{
  // destructor
  if(fOutputList) {
    delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
  }
}



//---------------------------------------------------
void AliAnalysisTaskChiralVorticalEffect::UserCreateOutputObjects()
{
  fOutputList = new TList();
  fOutputList->SetName(GetName());
  fOutputList->SetOwner(kTRUE);

  // event-wise
  hEvtCount = new TH1I("evtCount","",20, 1, 21);
  hEvtCount->GetXaxis()->SetBinLabel(1,"All");
  hEvtCount->GetXaxis()->SetBinLabel(2,"Info");
  hEvtCount->GetXaxis()->SetBinLabel(3,"Evt");
  hEvtCount->GetXaxis()->SetBinLabel(4,"Cent");

  hEvtCount->GetXaxis()->SetBinLabel(10,"Manager");
  hEvtCount->GetXaxis()->SetBinLabel(11,"Handler");
  hEvtCount->GetXaxis()->SetBinLabel(12,"AOD");
  hEvtCount->GetXaxis()->SetBinLabel(13,"PID");
  hEvtCount->GetXaxis()->SetBinLabel(14,"Utils");
  // hEvtCount->GetXaxis()->SetBinLabel(15,"MultSel");

  hEvtCount->GetXaxis()->SetBinLabel(16,"Random plane");
  hEvtCount->GetXaxis()->SetBinLabel(17,"TPC plane");
  hEvtCount->GetXaxis()->SetBinLabel(18,"VZERO plane");
  hEvtCount->GetXaxis()->SetBinLabel(19,"ZDC plane");
  hEvtCount->GetXaxis()->SetBinLabel(20,"loops end");
  fOutputList->Add(hEvtCount);

  // 10h
  TString runNumList[91]={
    "139510","139507","139505","139503","139465","139438","139437","139360","139329","139328","139314","139310",
    "139309","139173","139107","139105","139038","139037","139036","139029","139028","138872","138871","138870",
    "138837","138732","138730","138666","138662","138653","138652","138638","138624","138621","138583","138582",
    "138579","138578","138534","138469","138442","138439","138438","138396","138364","138275","138225","138201",
    "138197","138192","138190","137848","137844","137752","137751","137724","137722","137718","137704","137693",
    "137692","137691","137686","137685","137639","137638","137608","137595","137549","137546","137544","137541",
    "137539","137531","137530","137443","137441","137440","137439","137434","137432","137431","137430","137243",
    "137236","137235","137232","137231","137230","137162","137161"
  };
  // 11h
  // TString runNum[39]={"170387","170040","170268","170228","170207","169838","170159","170204","170311","170084",
  //                     "169835","170088","170593","170203","170270","169846","170163","170388","170155","170083",
  //                     "170572","169837","169855","170306","170269","170089","170309","170091","170081","170230",
  //                     "170085","170315","170027","170193","170312","170313","170308","169858","169859"};
  hRunNumBin = new TH1I("runNumBin","",100,0,100);
  for (int i=0; i<91; ++i) {    
    hRunNumBin->GetXaxis()->SetBinLabel(i+1,runNumList[i].Data());
  }
  fOutputList->Add(hRunNumBin);

  hCent = new TH1D("centrality","",100,0,100);
  fOutputList->Add(hCent);
  hCentCorr[0] = new TH2D("centcorr0","",100,0,100,100,0,100);
  hCentCorr[1] = new TH2D("centcorr1","",100,0,100,100,0,100);
  for (int i=0; i<2; ++i) fOutputList->Add(hCentCorr[i]);

  hVxy[0] = new TH2D("vxy0","",100,-0.5,0.5,100,-0.5,0.5);
  hVxy[1] = new TH2D("vxy1","",100,-0.5,0.5,100,-0.5,0.5);
  hVz[0]  = new TH1D("vz0","",200,-50,50);
  hVz[1]  = new TH1D("vz1","",200,-50,50);
  for (int i=0; i<2; ++i) fOutputList->Add(hVxy[i]);
  for (int i=0; i<2; ++i) fOutputList->Add(hVz[i]);


  //Random Event Plane
  hPsi2RDM = new TH3D("hPsiRDM","",100, 0, 100, 10, 0, 10, 180, 0, TMath::Pi());
  hPsi3RDM = new TH3D("hPsiRDM","",100, 0, 100, 10, 0, 10, 180, 0, 2./3.*TMath::Pi());
  fOutputList->Add(hPsi2RDM);
  fOutputList->Add(hPsi3RDM);

  //TPC Event Plane
  hPsi2TPC = new TH3D("hPsiTPC","",100, 0, 100, 10, 0, 10, 180, 0, TMath::Pi());
  hPsi3TPC = new TH3D("hPsiTPC","",100, 0, 100, 10, 0, 10, 180, 0, 2./3.*TMath::Pi());
  fOutputList->Add(hPsi2TPC);
  fOutputList->Add(hPsi3TPC);

  // track-wise
  hPt[0] = new TH1D("hPtBeforeCut", "", 200, 0., 20.);
  hPt[1] = new TH1D("hPtAfterCut", "", 200, 0., 20.);
  for (int i=0; i<2; ++i) fOutputList->Add(hPt[i]);
  hEta[0] = new TH1D("hEtaBeforeCut", "", 200, -10., 10.);
  hEta[1] = new TH1D("hEtaAfterCut",  "", 200, -10., 10.);
  for (int i=0; i<2; ++i) fOutputList->Add(hEta[i]);
  hPhi[0] = new TH1D("hPhiBeforeCut", "", 400, 0, 2*TMath::Pi());
  hPhi[1] = new TH1D("hPhiAfterCut", "", 400, 0, 2*TMath::Pi());
  for (int i=0; i<2; ++i) fOutputList->Add(hPhi[i]);
  hNhits[0] = new TH1D("hNhitsBeforeCut", "", 200, 0., 200.);
  hNhits[1] = new TH1D("hNhitsAfterCut",  "", 200, 0., 200.);
  for (int i=0; i<2; ++i) fOutputList->Add(hNhits[i]);
  hDcaXy[0] = new TH1D("hDcaXyBeforeCut", "", 100, 0., 10.);
  hDcaXy[1] = new TH1D("hDcaXyAfterCut",  "", 100, 0., 10.);
  for (int i=0; i<2; ++i) fOutputList->Add(hDcaXy[i]);
  hDcaZ[0] = new TH1D("hDcaZBeforeCut", "", 100, 0., 10.);
  hDcaZ[1] = new TH1D("hDcaZAfterCut",  "", 100, 0., 10.);
  for (int i=0; i<2; ++i) fOutputList->Add(hDcaZ[i]);
  hPDedx = new TH2D("hPDedx", "", 400, -10., 10., 400, 0, 1000);
  fOutputList->Add(hPDedx);


  //V0-wise
  hV0Pt = new TH1D("hV0Pt","", 200, 0., 20.);
  hV0Eta = new TH1D("hV0Eta","", 200, -10., 10.);
  hV0DcaToPrimVertex = new TH1D("hV0DcaToPrimVertex","",200, 0., 20.);
  hV0CPA = new TH1D("hV0CPA","", 1000, 0.9, 1.);
  hV0DecayLength = new TH1D("hV0DecayLength","",500,0,500.);

  fOutputList->Add(hV0Pt);
  fOutputList->Add(hV0Eta);
  fOutputList->Add(hV0DcaToPrimVertex);
  fOutputList->Add(hV0CPA);
  fOutputList->Add(hV0DecayLength);

  hLambdaPt[0] = new TH1D("hLambdaPt_bfMassCut","", 200, 0., 20.);
  hLambdaPt[1] = new TH1D("hLambdaPt_afMassCut","", 200, 0., 20.);
  hLambdaEta[0] = new TH1D("hLambdaEta_bfMassCut","",200, -10., 10.);
  hLambdaEta[1] = new TH1D("hLambdaEta_afMassCut","",200, -10., 10.);
  hLambdaDcaToPrimVertex[0] = new TH1D("hLambdaDcaToPrimVertex_bfMassCut","",200, 0., 20.);
  hLambdaDcaToPrimVertex[1] = new TH1D("hLambdaDcaToPrimVertex_afMassCut","",200, 0., 20.);
  hLambdaCPA[0] = new TH1D("hLambdaCPA_bfMassCut","", 1000, 0.9, 1.);
  hLambdaCPA[1] = new TH1D("hLambdaCPA_afMassCut","", 1000, 0.9, 1.);
  hLambdaDecayLength[0] = new TH1D("hLambdaDecayLength_bfMassCut","", 500, 0., 500.);
  hLambdaDecayLength[1] = new TH1D("hLambdaDecayLength_afMassCut","", 500, 0., 500.);
  hLambdaMass[0] = new TH1D("hLambdaMass_bfMassCut","",1000,1.,1.25);
  hLambdaMass[1] = new TH1D("hLambdaMass_afMassCut","",1000,1.,1.25);
  for (int i=0; i<10; ++i) hLambdaMassCent[i] = new TH1D(Form("hLambdaMassCent%i",i),"",1000,1.,1.25);

  for( int i=0; i<2; i++) fOutputList->Add(hLambdaPt[i]);
  for( int i=0; i<2; i++) fOutputList->Add(hLambdaEta[i]);
  for( int i=0; i<2; i++) fOutputList->Add(hLambdaDcaToPrimVertex[i]);
  for( int i=0; i<2; i++) fOutputList->Add(hLambdaCPA[i]);
  for( int i=0; i<2; i++) fOutputList->Add(hLambdaDecayLength[i]);
  for( int i=0; i<2; i++) fOutputList->Add(hLambdaMass[i]);
  for (int i=0; i<10; ++i) fOutputList->Add(hLambdaMassCent[i]);

  hAntiLambdaPt[0] = new TH1D("hAntiLambdaPt_bfMassCut","", 200, 0., 20.);
  hAntiLambdaPt[1] = new TH1D("hAntiLambdaPt_afMassCut","", 200, 0., 20.);
  hAntiLambdaEta[0] = new TH1D("hAntiLambdaEta_bfMassCut","",200, -10., 10.);
  hAntiLambdaEta[1] = new TH1D("hAntiLambdaEta_afMassCut","",200, -10., 10.);
  hAntiLambdaDcaToPrimVertex[0] = new TH1D("hAntiLambdaDcaToPrimVertex_bfMassCut","",200, 0., 20.);
  hAntiLambdaDcaToPrimVertex[1] = new TH1D("hAntiLambdaDcaToPrimVertex_afMassCut","",200, 0., 20.);
  hAntiLambdaCPA[0] = new TH1D("hAntiLambdaCPA_bfMassCut","", 1000, 0.9, 1.);
  hAntiLambdaCPA[1] = new TH1D("hAntiLambdaCPA_afMassCut","", 1000, 0.9, 1.);
  hAntiLambdaDecayLength[0] = new TH1D("hAntiLambdaDecayLength_bfMassCut","", 500, 0., 500.);
  hAntiLambdaDecayLength[1] = new TH1D("hAntiLambdaDecayLength_afMassCut","", 500, 0., 500.);
  hAntiLambdaMass[0] = new TH1D("hAntiLambdaMass_bfMassCut","",1000,1.,1.25);
  hAntiLambdaMass[1] = new TH1D("hAntiLambdaMass_afMassCut","",1000,1.,1.25);
  for (int i=0; i<10; ++i) hAntiLambdaMassCent[i] = new TH1D(Form("hAntiLambdaMassCent%i",i),"",1000,1.,1.25);

  for( int i=0; i<2; i++) fOutputList->Add(hAntiLambdaPt[i]);
  for( int i=0; i<2; i++) fOutputList->Add(hAntiLambdaEta[i]);
  for( int i=0; i<2; i++) fOutputList->Add(hAntiLambdaDcaToPrimVertex[i]);
  for( int i=0; i<2; i++) fOutputList->Add(hAntiLambdaCPA[i]);
  for( int i=0; i<2; i++) fOutputList->Add(hAntiLambdaDecayLength[i]);
  for( int i=0; i<2; i++) fOutputList->Add(hAntiLambdaMass[i]);
  for (int i=0; i<10; ++i) fOutputList->Add(hAntiLambdaMassCent[i]);

  //Inclusive
  //Delta
  pDelta_hPos_hPos = new TProfile("pDelta_hPos_hPos","",20,0.,100.);
  pDelta_hNeg_hNeg = new TProfile("pDelta_hNeg_hNeg","",20,0.,100.);
  pDelta_hPos_hNeg = new TProfile("pDelta_hPos_hNeg","",20,0.,100.);
  //Gamma
  //TPC Plane
  pGammaTPC_hPos_hPos = new TProfile("pGammaTPC_hPos_hPos","",20,0.,100.);
  pGammaTPC_hNeg_hNeg = new TProfile("pGammaTPC_hNeg_hNeg","",20,0.,100.);
  pGammaTPC_hPos_hNeg = new TProfile("pGammaTPC_hPos_hNeg","",20,0.,100.);
  pCosCosTPC_hPos_hPos = new TProfile("pCosCosTPC_hPos_hPos","",20,0.,100.);
  pCosCosTPC_hNeg_hNeg = new TProfile("pCosCosTPC_hNeg_hNeg","",20,0.,100.);
  pCosCosTPC_hPos_hNeg = new TProfile("pCosCosTPC_hPos_hNeg","",20,0.,100.);
  pSinSinTPC_hPos_hPos = new TProfile("pSinSinTPC_hPos_hPos","",20,0.,100.);
  pSinSinTPC_hNeg_hNeg = new TProfile("pSinSinTPC_hNeg_hNeg","",20,0.,100.);
  pSinSinTPC_hPos_hNeg = new TProfile("pSinSinTPC_hPos_hNeg","",20,0.,100.);
  //Random Plane 
  pGammaRDM_hPos_hPos = new TProfile("pGammaRDM_hPos_hPos","",20,0.,100.);
  pGammaRDM_hNeg_hNeg = new TProfile("pGammaRDM_hNeg_hNeg","",20,0.,100.);
  pGammaRDM_hPos_hNeg = new TProfile("pGammaRDM_hPos_hNeg","",20,0.,100.);
  pCosCosRDM_hPos_hPos = new TProfile("pCosCosRDM_hPos_hPos","",20,0.,100.);
  pCosCosRDM_hNeg_hNeg = new TProfile("pCosCosRDM_hNeg_hNeg","",20,0.,100.);
  pCosCosRDM_hPos_hNeg = new TProfile("pCosCosRDM_hPos_hNeg","",20,0.,100.);
  pSinSinRDM_hPos_hPos = new TProfile("pSinSinRDM_hPos_hPos","",20,0.,100.);
  pSinSinRDM_hNeg_hNeg = new TProfile("pSinSinRDM_hNeg_hNeg","",20,0.,100.);
  pSinSinRDM_hPos_hNeg = new TProfile("pSinSinRDM_hPos_hNeg","",20,0.,100.);
  //Diff
  for(int i=0;i<10;i++) pGammaTPCVsDeltaPt_hPos_hPos[i] = new TProfile(Form("pGammaTPCVsDeltaPt_hPos_hPos_cent%i",i),"",15,0.,6.);
  for(int i=0;i<10;i++) pGammaTPCVsDeltaPt_hPos_hNeg[i] = new TProfile(Form("pGammaTPCVsDeltaPt_hPos_hNeg_cent%i",i),"",15,0.,6.);
  for(int i=0;i<10;i++) pGammaTPCVsDeltaPt_hNeg_hNeg[i] = new TProfile(Form("pGammaTPCVsDeltaPt_hNeg_hNeg_cent%i",i),"",15,0.,6.);

  for(int i=0;i<10;i++) pGammaTPVVsMeanPt_hPos_hPos[i] = new TProfile(Form("pGammaTPVVsMeanPt_hPos_hPos_cent%i",i),"",15,0.,6.);
  for(int i=0;i<10;i++) pGammaTPVVsMeanPt_hPos_hNeg[i] = new TProfile(Form("pGammaTPVVsMeanPt_hPos_hNeg_cent%i",i),"",15,0.,6.);
  for(int i=0;i<10;i++) pGammaTPVVsMeanPt_hNeg_hNeg[i] = new TProfile(Form("pGammaTPVVsMeanPt_hNeg_hNeg_cent%i",i),"",15,0.,6.);

  for(int i=0;i<10;i++) pGammaTPCVsDeltaEta_hPos_hPos[i] = new TProfile(Form("pGammaTPCVsDeltaEta_hPos_hPos_cent%i",i),"",10,0.,2.);
  for(int i=0;i<10;i++) pGammaTPCVsDeltaEta_hPos_hNeg[i] = new TProfile(Form("pGammaTPCVsDeltaEta_hPos_hNeg_cent%i",i),"",10,0.,2.);
  for(int i=0;i<10;i++) pGammaTPCVsDeltaEta_hNeg_hNeg[i] = new TProfile(Form("pGammaTPCVsDeltaEta_hNeg_hNeg_cent%i",i),"",10,0.,2.);


  //pion - Inclusive
  // pi+ - h+
  pDelta_pion_hPos = new TProfile("pDelta_pion_hPos","",20,0.,100.);
  pGammaTPC_pion_hPos = new TProfile("pGammaTPC_pion_hPos","",20,0.,100.);
  pGammaRDM_pion_hPos = new TProfile("pGammaRDM_pion_hPos","",20,0.,100.);
  for(int i=0;i<10;i++) pGammaTPCVsPt_pion_hPos[i] = new TProfile(Form("pGammaTPCVsPt_pion_hPos_cent%i",i),"",6,0.,3.);
  // pi- - h+
  pDelta_antiPion_hPos = new TProfile("pDelta_antiPion_hPos","",20,0.,100.);
  pGammaTPC_antiPion_hPos = new TProfile("pGammaTPC_antiPion_hPos","",20,0.,100.);
  pGammaRDM_antiPion_hPos = new TProfile("pGammaRDM_antiPion_hPos","",20,0.,100.);
  for(int i=0;i<10;i++) pGammaTPCVsPt_antiPion_hPos[i] = new TProfile(Form("pGammaTPCVsPt_antiPion_hPos_cent%i",i),"",6,0.,3.);
  // pi+ - h-
  pDelta_pion_hNeg = new TProfile("pDelta_pion_hNeg","",20,0.,100.);
  pGammaTPC_pion_hNeg = new TProfile("pGammaTPC_pion_hNeg","",20,0.,100.);
  pGammaRDM_pion_hNeg = new TProfile("pGammaRDM_pion_hNeg","",20,0.,100.);
  for(int i=0;i<10;i++) pGammaTPCVsPt_pion_hNeg[i] = new TProfile(Form("pGammaTPCVsPt_pion_hNeg_cent%i",i),"",6,0.,3.);
  // pi- - h-
  pDelta_antiPion_hNeg = new TProfile("pDelta_antiPion_hNeg","",20,0.,100.);
  pGammaTPC_antiPion_hNeg = new TProfile("pGammaTPC_antiPion_hNeg","",20,0.,100.);
  pGammaRDM_antiPion_hNeg = new TProfile("pGammaRDM_antiPion_hNeg","",20,0.,100.);
  for(int i=0;i<10;i++) pGammaTPCVsPt_antiPion_hNeg[i] = new TProfile(Form("pGammaTPCVsPt_antiPion_hNeg_cent%i",i),"",6,0.,3.);
  //kaon - Inclusive
  // ka+ - h+
  pDelta_kaon_hPos = new TProfile("pDelta_kaon_hPos","",20,0.,100.);
  pGammaTPC_kaon_hPos = new TProfile("pGammaTPC_kaon_hPos","",20,0.,100.);
  pGammaRDM_kaon_hPos = new TProfile("pGammaRDM_kaon_hPos","",20,0.,100.);
  for(int i=0;i<10;i++) pGammaTPCVsPt_kaon_hPos[i] = new TProfile(Form("pGammaTPCVsPt_kaon_hPos_cent%i",i),"",6,0.,3.);
  // ka- - h+
  pDelta_antiKaon_hPos = new TProfile("pDelta_antiKaon_hPos","",20,0.,100.);
  pGammaTPC_antiKaon_hPos = new TProfile("pGammaTPC_antiKaon_hPos","",20,0.,100.);
  pGammaRDM_antiKaon_hPos = new TProfile("pGammaRDM_antiKaon_hPos","",20,0.,100.);
  for(int i=0;i<10;i++) pGammaTPCVsPt_antiKaon_hPos[i] = new TProfile(Form("pGammaTPCVsPt_antiKaon_hPos_cent%i",i),"",6,0.,3.);
  // ka+ - h-
  pDelta_kaon_hNeg = new TProfile("pDelta_kaon_hNeg","",20,0.,100.);
  pGammaTPC_kaon_hNeg = new TProfile("pGammaTPC_kaon_hNeg","",20,0.,100.);
  pGammaRDM_kaon_hNeg = new TProfile("pGammaRDM_kaon_hNeg","",20,0.,100.);
  for(int i=0;i<10;i++) pGammaTPCVsPt_kaon_hNeg[i] = new TProfile(Form("pGammaTPCVsPt_kaon_hNeg_cent%i",i),"",6,0.,3.);
  // ka- - h-
  pDelta_antiKaon_hNeg = new TProfile("pDelta_antiKaon_hNeg","",20,0.,100.);
  pGammaTPC_antiKaon_hNeg = new TProfile("pGammaTPC_antiKaon_hNeg","",20,0.,100.);
  pGammaRDM_antiKaon_hNeg = new TProfile("pGammaRDM_antiKaon_hNeg","",20,0.,100.);
  for(int i=0;i<10;i++) pGammaTPCVsPt_antiKaon_hNeg[i] = new TProfile(Form("pGammaTPCVsPt_antiKaon_hNeg_cent%i",i),"",6,0.,3.);
  //proton - Inclusive
  // p+ - h+
  pDelta_proton_hPos = new TProfile("pDelta_proton_hPos","",20,0.,100.);
  pGammaTPC_proton_hPos = new TProfile("pGammaTPC_proton_hPos","",20,0.,100.);
  pGammaRDM_proton_hPos = new TProfile("pGammaRDM_proton_hPos","",20,0.,100.);
  for(int i=0;i<10;i++) pGammaTPCVsPt_proton_hPos[i] = new TProfile(Form("pGammaTPCVsPt_proton_hPos_cent%i",i),"",6,0.,3.);
  // p- - h+
  pDelta_antiProton_hPos = new TProfile("pDelta_antiProton_hPos","",20,0.,100.);
  pGammaTPC_antiProton_hPos = new TProfile("pGammaTPC_antiProton_hPos","",20,0.,100.);
  pGammaRDM_antiProton_hPos = new TProfile("pGammaRDM_antiProton_hPos","",20,0.,100.);
  for(int i=0;i<10;i++) pGammaTPCVsPt_antiProton_hPos[i] = new TProfile(Form("pGammaTPCVsPt_antiProton_hPos_cent%i",i),"",6,0.,3.);
  // p+ - h-
  pDelta_proton_hNeg = new TProfile("pDelta_proton_hNeg","",20,0.,100.);
  pGammaTPC_proton_hNeg = new TProfile("pGammaTPC_proton_hNeg","",20,0.,100.);
  pGammaRDM_proton_hNeg = new TProfile("pGammaRDM_proton_hNeg","",20,0.,100.);
  for(int i=0;i<10;i++) pGammaTPCVsPt_proton_hNeg[i] = new TProfile(Form("pGammaTPCVsPt_proton_hNeg_cent%i",i),"",6,0.,3.);
  // p- - h-
  pDelta_antiProton_hNeg = new TProfile("pDelta_antiProton_hNeg","",20,0.,100.);
  pGammaTPC_antiProton_hNeg = new TProfile("pGammaTPC_antiProton_hNeg","",20,0.,100.);
  pGammaRDM_antiProton_hNeg = new TProfile("pGammaRDM_antiProton_hNeg","",20,0.,100.);
  for(int i=0;i<10;i++) pGammaTPCVsPt_antiProton_hNeg[i] = new TProfile(Form("pGammaTPCVsPt_antiProton_hNeg_cent%i",i),"",6,0.,3.);


  //pion - pion
  //pi+ - pi+
  pDelta_pion_pion = new TProfile("pDelta_pion_pion","",20,0.,100.);
  pGammaTPC_pion_pion = new TProfile("pGammaTPC_pion_pion","",20,0.,100.);
  pGammaRDM_pion_pion = new TProfile("pGammaRDM_pion_pion","",20,0.,100.);
  //pi- - pi-
  pDelta_antiPion_antiPion = new TProfile("pDelta_antiPion_antiPion","",20,0.,100.);
  pGammaTPC_antiPion_antiPion = new TProfile("pGammaTPC_antiPion_antiPion","",20,0.,100.);
  pGammaRDM_antiPion_antiPion = new TProfile("pGammaRDM_antiPion_antiPion","",20,0.,100.);
  //pi+ - pi-
  pDelta_pion_antiPion = new TProfile("pDelta_pion_antiPion","",20,0.,100.);
  pGammaTPC_pion_antiPion = new TProfile("pGammaTPC_pion_antiPion","",20,0.,100.);
  pGammaRDM_pion_antiPion = new TProfile("pGammaRDM_pion_antiPion","",20,0.,100.);
  //kaon - kaon
  //ka+ - ka+
  pDelta_kaon_kaon = new TProfile("pDelta_kaon_kaon","",20,0.,100.);
  pGammaTPC_kaon_kaon = new TProfile("pGammaTPC_kaon_kaon","",20,0.,100.);
  pGammaRDM_kaon_kaon = new TProfile("pGammaRDM_kaon_kaon","",20,0.,100.);
  //ka- - ka-
  pDelta_antiKaon_antiKaon = new TProfile("pDelta_antiKaon_antiKaon","",20,0.,100.);
  pGammaTPC_antiKaon_antiKaon = new TProfile("pGammaTPC_antiKaon_antiKaon","",20,0.,100.);
  pGammaRDM_antiKaon_antiKaon = new TProfile("pGammaRDM_antiKaon_antiKaon","",20,0.,100.);
  //ka+ - ka-
  pDelta_kaon_antiKaon = new TProfile("pDelta_kaon_antiKaon","",20,0.,100.);
  pGammaTPC_kaon_antiKaon = new TProfile("pGammaTPC_kaon_antiKaon","",20,0.,100.);
  pGammaRDM_kaon_antiKaon = new TProfile("pGammaRDM_kaon_antiKaon","",20,0.,100.);
  //proton - proton
  //p+ - p+
  pDelta_proton_proton = new TProfile("pDelta_proton_proton","",20,0.,100.);
  pGammaTPC_proton_proton = new TProfile("pGammaTPC_proton_proton","",20,0.,100.);
  pGammaRDM_proton_proton = new TProfile("pGammaRDM_proton_proton","",20,0.,100.);
  //p- - p-
  pDelta_antiProton_antiProton = new TProfile("pDelta_antiProton_antiProton","",20,0.,100.);
  pGammaTPC_antiProton_antiProton = new TProfile("pGammaTPC_antiProton_antiProton","",20,0.,100.);
  pGammaRDM_antiProton_antiProton = new TProfile("pGammaRDM_antiProton_antiProton","",20,0.,100.);
  //p+ - p-
  pDelta_proton_antiProton = new TProfile("pDelta_proton_antiProton","",20,0.,100.);
  pGammaTPC_proton_antiProton = new TProfile("pGammaTPC_proton_antiProton","",20,0.,100.);
  pGammaRDM_proton_antiProton = new TProfile("pGammaRDM_proton_antiProton","",20,0.,100.);


  //pion - kaon
  // pi+ - k+
  pDelta_pion_kaon = new TProfile("pDelta_pion_kaon","",20,0.,100.);
  pGammaTPC_pion_kaon = new TProfile("pGammaTPC_pion_kaon","",20,0.,100.);
  pGammaRDM_pion_kaon = new TProfile("pGammaRDM_pion_kaon","",20,0.,100.);
  // pi- - k+
  pDelta_antiPion_kaon = new TProfile("pDelta_antiPion_kaon","",20,0.,100.);
  pGammaTPC_antiPion_kaon = new TProfile("pGammaTPC_antiPion_kaon","",20,0.,100.);
  pGammaRDM_antiPion_kaon = new TProfile("pGammaRDM_antiPion_kaon","",20,0.,100.);
  // pi+ - k-
  pDelta_pion_antiKaon = new TProfile("pDelta_pion_antiKaon","",20,0.,100.);
  pGammaTPC_pion_antiKaon = new TProfile("pGammaTPC_pion_antiKaon","",20,0.,100.);
  pGammaRDM_pion_antiKaon = new TProfile("pGammaRDM_pion_antiKaon","",20,0.,100.);
  // pi- - k-
  pDelta_antiPion_antiKaon = new TProfile("pDelta_antiPion_antiKaon","",20,0.,100.);
  pGammaTPC_antiPion_antiKaon = new TProfile("pGammaTPC_antiPion_antiKaon","",20,0.,100.);
  pGammaRDM_antiPion_antiKaon = new TProfile("pGammaRDM_antiPion_antiKaon","",20,0.,100.);
  //kaon - proton
  // k+ - p+
  pDelta_kaon_proton = new TProfile("pDelta_kaon_proton","",20,0.,100.);
  pGammaTPC_kaon_proton = new TProfile("pGammaTPC_kaon_proton","",20,0.,100.);
  pGammaRDM_kaon_proton = new TProfile("pGammaRDM_kaon_proton","",20,0.,100.);
  // k- - p+
  pDelta_antiKaon_proton = new TProfile("pDelta_antiKaon_proton","",20,0.,100.);
  pGammaTPC_antiKaon_proton = new TProfile("pGammaTPC_antiKaon_proton","",20,0.,100.);
  pGammaRDM_antiKaon_proton = new TProfile("pGammaRDM_antiKaon_proton","",20,0.,100.);
  // k+ - p-
  pDelta_kaon_antiProton = new TProfile("pDelta_kaon_antiProton","",20,0.,100.);
  pGammaTPC_kaon_antiProton = new TProfile("pGammaTPC_kaon_antiProton","",20,0.,100.);
  pGammaRDM_kaon_antiProton = new TProfile("pGammaRDM_kaon_antiProton","",20,0.,100.);
  // k- - p-
  pDelta_antiKaon_antiProton = new TProfile("pDelta_antiKaon_antiProton","",20,0.,100.);
  pGammaTPC_antiKaon_antiProton = new TProfile("pGammaTPC_antiKaon_antiProton","",20,0.,100.);
  pGammaRDM_antiKaon_antiProton = new TProfile("pGammaRDM_antiKaon_antiProton","",20,0.,100.);
  //proton - pi
  // p+ - pi+
  pDelta_proton_pion = new TProfile("pDelta_proton_pion","",20,0.,100.);
  pGammaTPC_proton_pion = new TProfile("pGammaTPC_proton_pion","",20,0.,100.);
  pGammaRDM_proton_pion = new TProfile("pGammaRDM_proton_pion","",20,0.,100.);
  // p- - pi+
  pDelta_antiProton_pion = new TProfile("pDelta_antiProton_pion","",20,0.,100.);
  pGammaTPC_antiProton_pion = new TProfile("pGammaTPC_antiProton_pion","",20,0.,100.);
  pGammaRDM_antiProton_pion = new TProfile("pGammaRDM_antiProton_pion","",20,0.,100.);
  // p+ - pi-
  pDelta_proton_antiPion = new TProfile("pDelta_proton_antiPion","",20,0.,100.);
  pGammaTPC_proton_antiPion = new TProfile("pGammaTPC_proton_antiPion","",20,0.,100.);
  pGammaRDM_proton_antiPion = new TProfile("pGammaRDM_proton_antiPion","",20,0.,100.);
  // p- - pi-
  pDelta_antiProton_antiPion = new TProfile("pDelta_antiProton_antiPion","",20,0.,100.);
  pGammaTPC_antiProton_antiPion = new TProfile("pGammaTPC_antiProton_antiPion","",20,0.,100.);
  pGammaRDM_antiProton_antiPion = new TProfile("pGammaRDM_antiProton_antiPion","",20,0.,100.);



  //Inclusive
  //Delta
  fOutputList->Add(pDelta_hPos_hPos);
  fOutputList->Add(pDelta_hNeg_hNeg);
  fOutputList->Add(pDelta_hPos_hNeg);
  //Gamma
  //TPC Plane
  fOutputList->Add(pGammaTPC_hPos_hPos);
  fOutputList->Add(pGammaTPC_hNeg_hNeg);
  fOutputList->Add(pGammaTPC_hPos_hNeg);
  fOutputList->Add(pCosCosTPC_hPos_hPos);
  fOutputList->Add(pCosCosTPC_hNeg_hNeg);
  fOutputList->Add(pCosCosTPC_hPos_hNeg);
  fOutputList->Add(pSinSinTPC_hPos_hPos);
  fOutputList->Add(pSinSinTPC_hNeg_hNeg);
  fOutputList->Add(pSinSinTPC_hPos_hNeg);
  //Random Plane 
  fOutputList->Add(pGammaRDM_hPos_hPos);
  fOutputList->Add(pGammaRDM_hNeg_hNeg);
  fOutputList->Add(pGammaRDM_hPos_hNeg);
  fOutputList->Add(pCosCosRDM_hPos_hPos);
  fOutputList->Add(pCosCosRDM_hNeg_hNeg);
  fOutputList->Add(pCosCosRDM_hPos_hNeg);
  fOutputList->Add(pSinSinRDM_hPos_hPos);
  fOutputList->Add(pSinSinRDM_hNeg_hNeg);
  fOutputList->Add(pSinSinRDM_hPos_hNeg);
  //Diff
  for(int i=0;i<10;i++) fOutputList->Add(pGammaTPCVsDeltaPt_hPos_hPos[i]);
  for(int i=0;i<10;i++) fOutputList->Add(pGammaTPCVsDeltaPt_hPos_hNeg[i]);
  for(int i=0;i<10;i++) fOutputList->Add(pGammaTPCVsDeltaPt_hNeg_hNeg[i]);

  for(int i=0;i<10;i++) fOutputList->Add(pGammaTPVVsMeanPt_hPos_hPos[i]);
  for(int i=0;i<10;i++) fOutputList->Add(pGammaTPVVsMeanPt_hPos_hNeg[i]);
  for(int i=0;i<10;i++) fOutputList->Add(pGammaTPVVsMeanPt_hNeg_hNeg[i]);

  for(int i=0;i<10;i++) fOutputList->Add(pGammaTPCVsDeltaEta_hPos_hPos[i]);
  for(int i=0;i<10;i++) fOutputList->Add(pGammaTPCVsDeltaEta_hPos_hNeg[i]);
  for(int i=0;i<10;i++) fOutputList->Add(pGammaTPCVsDeltaEta_hNeg_hNeg[i]);


  //pion - Inclusive
  // pi+ - h+
  fOutputList->Add(pDelta_pion_hPos);
  fOutputList->Add(pGammaTPC_pion_hPos);
  fOutputList->Add(pGammaRDM_pion_hPos);
  for(int i=0;i<10;i++) fOutputList->Add(pGammaTPCVsPt_pion_hPos[i]);
  // pi- - h+
  fOutputList->Add(pDelta_antiPion_hPos);
  fOutputList->Add(pGammaTPC_antiPion_hPos);
  fOutputList->Add(pGammaRDM_antiPion_hPos);
  for(int i=0;i<10;i++) fOutputList->Add(pGammaTPCVsPt_antiPion_hPos[i]);
  // pi+ - h-
  fOutputList->Add(pDelta_pion_hNeg);
  fOutputList->Add(pGammaTPC_pion_hNeg);
  fOutputList->Add(pGammaRDM_pion_hNeg);
  for(int i=0;i<10;i++) fOutputList->Add(pGammaTPCVsPt_pion_hNeg[i]);
  // pi- - h-
  fOutputList->Add(pDelta_antiPion_hNeg);
  fOutputList->Add(pGammaTPC_antiPion_hNeg);
  fOutputList->Add(pGammaRDM_antiPion_hNeg);
  for(int i=0;i<10;i++) fOutputList->Add(pGammaTPCVsPt_antiPion_hNeg[i]);
  //kaon - Inclusive
  // ka+ - h+
  fOutputList->Add(pDelta_kaon_hPos);
  fOutputList->Add(pGammaTPC_kaon_hPos);
  fOutputList->Add(pGammaRDM_kaon_hPos);
  for(int i=0;i<10;i++) fOutputList->Add(pGammaTPCVsPt_kaon_hPos[i]);
  // ka- - h+
  fOutputList->Add(pDelta_antiKaon_hPos);
  fOutputList->Add(pGammaTPC_antiKaon_hPos);
  fOutputList->Add(pGammaRDM_antiKaon_hPos);
  for(int i=0;i<10;i++) fOutputList->Add(pGammaTPCVsPt_antiKaon_hPos[i]);
  // ka+ - h-
  fOutputList->Add(pDelta_kaon_hNeg);
  fOutputList->Add(pGammaTPC_kaon_hNeg);
  fOutputList->Add(pGammaRDM_kaon_hNeg);
  for(int i=0;i<10;i++) fOutputList->Add(pGammaTPCVsPt_kaon_hNeg[i]);
  // ka- - h-
  fOutputList->Add(pDelta_antiKaon_hNeg);
  fOutputList->Add(pGammaTPC_antiKaon_hNeg);
  fOutputList->Add(pGammaRDM_antiKaon_hNeg);
  for(int i=0;i<10;i++) fOutputList->Add(pGammaTPCVsPt_antiKaon_hNeg[i]);
  //proton - Inclusive
  // p+ - h+
  fOutputList->Add(pDelta_proton_hPos);
  fOutputList->Add(pGammaTPC_proton_hPos);
  fOutputList->Add(pGammaRDM_proton_hPos);
  for(int i=0;i<10;i++) fOutputList->Add(pGammaTPCVsPt_proton_hPos[i]);
  // p- - h+
  fOutputList->Add(pDelta_antiProton_hPos);
  fOutputList->Add(pGammaTPC_antiProton_hPos);
  fOutputList->Add(pGammaRDM_antiProton_hPos);
  for(int i=0;i<10;i++) fOutputList->Add(pGammaTPCVsPt_antiProton_hPos[i]);
  // p+ - h-
  fOutputList->Add(pDelta_proton_hNeg);
  fOutputList->Add(pGammaTPC_proton_hNeg);
  fOutputList->Add(pGammaRDM_proton_hNeg);
  for(int i=0;i<10;i++) fOutputList->Add(pGammaTPCVsPt_proton_hNeg[i]);
  // p- - h-
  fOutputList->Add(pDelta_antiProton_hNeg);
  fOutputList->Add(pGammaTPC_antiProton_hNeg);
  fOutputList->Add(pGammaRDM_antiProton_hNeg);
  for(int i=0;i<10;i++) fOutputList->Add(pGammaTPCVsPt_antiProton_hNeg[i]);


  //pion - pion
  //pi+ - pi+
  fOutputList->Add(pDelta_pion_pion);
  fOutputList->Add(pGammaTPC_pion_pion);
  fOutputList->Add(pGammaRDM_pion_pion);
  //pi- - pi-
  fOutputList->Add(pDelta_antiPion_antiPion);
  fOutputList->Add(pGammaTPC_antiPion_antiPion);
  fOutputList->Add(pGammaRDM_antiPion_antiPion);
  //pi+ - pi-
  fOutputList->Add(pDelta_pion_antiPion);
  fOutputList->Add(pGammaTPC_pion_antiPion);
  fOutputList->Add(pGammaRDM_pion_antiPion);
  //kaon - kaon
  //ka+ - ka+
  fOutputList->Add(pDelta_kaon_kaon);
  fOutputList->Add(pGammaTPC_kaon_kaon);
  fOutputList->Add(pGammaRDM_kaon_kaon);
  //ka- - ka-
  fOutputList->Add(pDelta_antiKaon_antiKaon);
  fOutputList->Add(pGammaTPC_antiKaon_antiKaon);
  fOutputList->Add(pGammaRDM_antiKaon_antiKaon);
  //ka+ - ka-
  fOutputList->Add(pDelta_kaon_antiKaon);
  fOutputList->Add(pGammaTPC_kaon_antiKaon);
  fOutputList->Add(pGammaRDM_kaon_antiKaon);
  //proton - proton
  //p+ - p+
  fOutputList->Add(pDelta_proton_proton);
  fOutputList->Add(pGammaTPC_proton_proton);
  fOutputList->Add(pGammaRDM_proton_proton);
  //p- - p-
  fOutputList->Add(pDelta_antiProton_antiProton);
  fOutputList->Add(pGammaTPC_antiProton_antiProton);
  fOutputList->Add(pGammaRDM_antiProton_antiProton);
  //p+ - p-
  fOutputList->Add(pDelta_proton_antiProton);
  fOutputList->Add(pGammaTPC_proton_antiProton);
  fOutputList->Add(pGammaRDM_proton_antiProton);


  //pion - kaon
  // pi+ - k+
  fOutputList->Add(pDelta_pion_kaon);
  fOutputList->Add(pGammaTPC_pion_kaon);
  fOutputList->Add(pGammaRDM_pion_kaon);
  // pi- - k+
  fOutputList->Add(pDelta_antiPion_kaon);
  fOutputList->Add(pGammaTPC_antiPion_kaon);
  fOutputList->Add(pGammaRDM_antiPion_kaon);
  // pi+ - k-
  fOutputList->Add(pDelta_pion_antiKaon);
  fOutputList->Add(pGammaTPC_pion_antiKaon);
  fOutputList->Add(pGammaRDM_pion_antiKaon);
  // pi- - k-
  fOutputList->Add(pDelta_antiPion_antiKaon);
  fOutputList->Add(pGammaTPC_antiPion_antiKaon);
  fOutputList->Add(pGammaRDM_antiPion_antiKaon);
  //kaon - proton
  // k+ - p+
  fOutputList->Add(pDelta_kaon_proton);
  fOutputList->Add(pGammaTPC_kaon_proton);
  fOutputList->Add(pGammaRDM_kaon_proton);
  // k- - p+
  fOutputList->Add(pDelta_antiKaon_proton);
  fOutputList->Add(pGammaTPC_antiKaon_proton);
  fOutputList->Add(pGammaRDM_antiKaon_proton);
  // k+ - p-
  fOutputList->Add(pDelta_kaon_antiProton);
  fOutputList->Add(pGammaTPC_kaon_antiProton);
  fOutputList->Add(pGammaRDM_kaon_antiProton);
  // k- - p-
  fOutputList->Add(pDelta_antiKaon_antiProton);
  fOutputList->Add(pGammaTPC_antiKaon_antiProton);
  fOutputList->Add(pGammaRDM_antiKaon_antiProton);
  //proton - pi
  // p+ - pi+
  fOutputList->Add(pDelta_proton_pion);
  fOutputList->Add(pGammaTPC_proton_pion);
  fOutputList->Add(pGammaRDM_proton_pion);
  // p- - pi+
  fOutputList->Add(pDelta_antiProton_pion);
  fOutputList->Add(pGammaTPC_antiProton_pion);
  fOutputList->Add(pGammaRDM_antiProton_pion);
  // p+ - pi-
  fOutputList->Add(pDelta_proton_antiPion);
  fOutputList->Add(pGammaTPC_proton_antiPion);
  fOutputList->Add(pGammaRDM_proton_antiPion);
  // p- - pi-
  fOutputList->Add(pDelta_antiProton_antiPion);
  fOutputList->Add(pGammaTPC_antiProton_antiPion);
  fOutputList->Add(pGammaRDM_antiProton_antiPion);


  //lambda - h
  //lambda - h+
  pDelta_lambda_hPos = new TProfile("pDelta_lambda_hPos","",20,0.,100.);
  pGammaTPC_lambda_hPos = new TProfile("pGammaTPC_lambda_hPos","",20,0.,100.);
  pGammaRDM_lambda_hPos = new TProfile("pGammaRDM_lambda_hPos","",20,0.,100.);
  //antiLambda - h+
  pDelta_antiLambda_hPos = new TProfile("pDelta_antiLambda_hPos","",20,0.,100.);
  pGammaTPC_antiLambda_hPos = new TProfile("pGammaTPC_antiLambda_hPos","",20,0.,100.);
  pGammaRDM_antiLambda_hPos = new TProfile("pGammaRDM_antiLambda_hPos","",20,0.,100.);
  //lambda - h-
  pDelta_lambda_hNeg = new TProfile("pDelta_lambda_hNeg","",20,0.,100.);
  pGammaTPC_lambda_hNeg = new TProfile("pGammaTPC_lambda_hNeg","",20,0.,100.);
  pGammaRDM_lambda_hNeg = new TProfile("pGammaRDM_lambda_hNeg","",20,0.,100.);
  //antiLambda - h-
  pDelta_antiLambda_hNeg = new TProfile("pDelta_antiLambda_hNeg","",20,0.,100.);
  pGammaTPC_antiLambda_hNeg = new TProfile("pGammaTPC_antiLambda_hNeg","",20,0.,100.);
  pGammaRDM_antiLambda_hNeg = new TProfile("pGammaRDM_antiLambda_hNeg","",20,0.,100.);

  //lambda - p
  //lambda - p+
  pDelta_lambda_proton = new TProfile("pDelta_lambda_proton","",20,0.,100.);
  pGammaTPC_lambda_proton = new TProfile("pGammaTPC_lambda_proton","",20,0.,100.);
  pGammaRDM_lambda_proton = new TProfile("pGammaRDM_lambda_proton","",20,0.,100.);
  //antiLambda - p+
  pDelta_antiLambda_proton = new TProfile("pDelta_antiLambda_proton","",20,0.,100.);
  pGammaTPC_antiLambda_proton = new TProfile("pGammaTPC_antiLambda_proton","",20,0.,100.);
  pGammaRDM_antiLambda_proton = new TProfile("pGammaRDM_antiLambda_proton","",20,0.,100.);
  //lambda - p-
  pDelta_lambda_antiProton = new TProfile("pDelta_lambda_antiProton","",20,0.,100.);
  pGammaTPC_lambda_antiProton = new TProfile("pGammaTPC_lambda_antiProton","",20,0.,100.);
  pGammaRDM_lambda_antiProton = new TProfile("pGammaRDM_lambda_antiProton","",20,0.,100.);
  //antiLambda - p-
  pDelta_antiLambda_antiProton = new TProfile("pDelta_antiLambda_antiProton","",20,0.,100.);
  pGammaTPC_antiLambda_antiProton = new TProfile("pGammaTPC_antiLambda_antiProton","",20,0.,100.);
  pGammaRDM_antiLambda_antiProton = new TProfile("pGammaRDM_antiLambda_antiProton","",20,0.,100.);

  fOutputList->Add(pDelta_lambda_hPos);
  fOutputList->Add(pGammaTPC_lambda_hPos);
  fOutputList->Add(pGammaRDM_lambda_hPos);
  fOutputList->Add(pDelta_antiLambda_hPos);
  fOutputList->Add(pGammaTPC_antiLambda_hPos);
  fOutputList->Add(pGammaRDM_antiLambda_hPos);
  fOutputList->Add(pDelta_lambda_hNeg);
  fOutputList->Add(pGammaTPC_lambda_hNeg);
  fOutputList->Add(pGammaRDM_lambda_hNeg);
  fOutputList->Add(pDelta_antiLambda_hNeg);
  fOutputList->Add(pGammaTPC_antiLambda_hNeg);
  fOutputList->Add(pGammaRDM_antiLambda_hNeg);
  fOutputList->Add(pDelta_lambda_proton);
  fOutputList->Add(pGammaTPC_lambda_proton);
  fOutputList->Add(pGammaRDM_lambda_proton);
  fOutputList->Add(pDelta_antiLambda_proton);
  fOutputList->Add(pGammaTPC_antiLambda_proton);
  fOutputList->Add(pGammaRDM_antiLambda_proton);
  fOutputList->Add(pDelta_lambda_antiProton);
  fOutputList->Add(pGammaTPC_lambda_antiProton);
  fOutputList->Add(pGammaRDM_lambda_antiProton);
  fOutputList->Add(pDelta_antiLambda_antiProton);
  fOutputList->Add(pGammaTPC_antiLambda_antiProton);
  fOutputList->Add(pGammaRDM_antiLambda_antiProton);

  PostData(1,fOutputList);
}

//---------------------------------------------------
void AliAnalysisTaskChiralVorticalEffect::UserExec(Option_t *)
{
  hEvtCount->Fill(1);

  AliAnalysisManager* manager = AliAnalysisManager::GetAnalysisManager();
  if (!manager) {
    AliError(Form("%s: Could not get Analysis Manager", GetName()));
  } else hEvtCount->Fill(10);
  AliAODInputHandler* handler = (AliAODInputHandler*)manager->GetInputEventHandler();
  if (!handler) {
    AliError(Form("%s: Could not get Input Handler", GetName()));
  } else hEvtCount->Fill(11);
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fAOD) {
    AliError(Form("%s: Could not get AOD event", GetName()));
    return;
  } else hEvtCount->Fill(12);
  AliPIDResponse* fPID = handler->GetPIDResponse();
  if (!fPID) {
    AliError(Form("%s: Could not get PIDResponse", GetName()));
  } else hEvtCount->Fill(13);
  AliAnalysisUtils* fUtils = new AliAnalysisUtils();
  if (!fUtils) {
    AliError(Form("%s: Could not get AliAnalysisUtils", GetName()));
  } else hEvtCount->Fill(14);
  // AliMultSelection* fMultSel = (AliMultSelection*)InputEvent()->FindListObject("MultSelection");
  // if (!fMultSel) {
  //   AliError(Form("%s: Could not get AliMultSelection", GetName()));
  // } else hEvtCount->Fill(15);
  if (!manager || !handler || !fAOD || !fUtils) return;
  hEvtCount->Fill(2);

  //------------------
  // event-wise
  //------------------

  // runNumber
  runNum = fAOD->GetRunNumber();
  runNumBin = GetRunNumBin(runNum);
  if (runNumBin<0) return;
  hRunNumBin->Fill(runNumBin);

  // vertex
  AliAODVertex* fVtx = fAOD->GetPrimaryVertex();
  vtx[0] = (double)fVtx->GetX();
  vtx[1] = (double)fVtx->GetY();
  vtx[2] = (double)fVtx->GetZ();
  double vzSPD  = fAOD->GetPrimaryVertexSPD()->GetZ();
  hVxy[0]->Fill(vtx[0], vtx[1]);
  hVz[0]->Fill(vtx[2]);
  if (fabs(vtx[0])<1e-6 || fabs(vtx[1])<1e-6 || fabs(vtx[2])<1e-6) return;


  if (fabs(vtx[2])>10) return;
  if (fabs(vtx[2]-vzSPD)>0.5) return;
  hVxy[1]->Fill(vtx[0], vtx[1]);
  hVz[1]->Fill(vtx[2]);
  for (int i = 0; i < 20; ++i) {
    if (vtx[2] > -10+i*1 && vtx[2] < -10+(i+1)*1) {vzBin = i; break;}
  }
  if (vzBin<0) return;

  // pileup
  fUtils->SetUseOutOfBunchPileUp(true);
  fUtils->SetUseMVPlpSelection(true);
  // fUtils->SetMinPlpContribMV(5);
  bool isPileup = fUtils->IsPileUpEvent(fAOD);
  // bool isPileup = fUtils->IsPileUpMV(fAOD); // pp, p-Pb
  if (isPileup) return;
  hEvtCount->Fill(3); 

  // centrality
  // run1
  cent = fAOD->GetCentrality()->GetCentralityPercentile("V0M");
  double cent2 = fAOD->GetCentrality()->GetCentralityPercentile("CL1");
  double cent3 = fAOD->GetCentrality()->GetCentralityPercentile("TRK");
  hCentCorr[0]->Fill(cent,cent2);
  if (fabs(cent-cent2)>7.5) return;
  if (fabs(cent-cent3)>5) return; // ANA-280
  if (cent<0 || cent>=100) return;
  hCentCorr[1]->Fill(cent,cent2);
  centBin = (int)cent/10;
  hCent->Fill(cent);
  hEvtCount->Fill(4);
  
  //------------------
  //* loop trk
  //------------------

  vector<double> vecPt;
  vector<double> vecEta;
  vector<double> vecPhi;
  vector<short>  vecID;
  vector<int>    vecType;

  //TPC QxQy
  double sumCos2Phi = 0.;
  double sumSin2Phi = 0.;
  double sumCos3Phi = 0.;
  double sumSin3Phi = 0.;

  int nTrk = fAOD->GetNumberOfTracks();
  for (int iTrk = 0; iTrk < nTrk; ++iTrk) {
    AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(iTrk));
    if (!track) {
      AliError(Form("%s: Could not get Track", GetName()));
      continue;
    }
    if (!track->TestFilterBit(fFilterBit)) continue;
    double pt     = track->Pt();
    double eta    = track->Eta();
    double phi    = track->Phi();
    int    charge = track->Charge();
    int    nhits  = track->GetTPCNcls();
    double dedx   = track->GetTPCsignal();
    double chi2   = track->Chi2perNDF();
    hPt[0]->Fill(pt);
    hEta[0]->Fill(eta);
    hPhi[0]->Fill(phi);
    hNhits[0]->Fill(nhits);

    if (pt<fPtMin || pt>fPtMax) continue;
    if (fabs(eta)>fEtaMax) continue;
    if (fabs(nhits)<fNhitsMin) continue;
    if (chi2>fChi2Max) continue;
    if (dedx<fDeDxMin) continue;

    hPt[1]->Fill(pt);
    hEta[1]->Fill(eta);
    hPhi[1]->Fill(phi);
    hNhits[1]->Fill(nhits);
    hPDedx->Fill(track->P()*charge, dedx);

    //------------------
    // dca
    //------------------
    double mag = fAOD->GetMagneticField();
    double dcaxy  = 999.;
    double dcaz   = 999.;
    double r[3];
    double dca[2];
    double cov[3];
    bool proptodca = track->PropagateToDCA(fVtx, mag, 100., dca, cov);
    if (track->GetXYZ(r)) {
      dcaxy = r[0];
      dcaz  = r[1];
    } else {
      double dcax = r[0] - vtx[0];
      double dcay = r[1] - vtx[1];
      dcaz  = r[2] - vtx[2];
      dcaxy = sqrt(dcax*dcax + dcay*dcay);
      // dcaxy = dca[0];
    }
    hDcaXy[0]->Fill(dcaxy);
    if (fabs(dcaxy)>fDcaXyMax) continue;//DCAxy cut
    hDcaXy[1]->Fill(dcaxy);
    hDcaZ[0]->Fill(dcaz);
    if (fabs(dcaz)>fDcaZMax) continue;//DCAz cut
    hDcaZ[1]->Fill(dcaz);

    //PID
    //TPC
    float nSigmaTPCPi = fPID->NumberOfSigmasTPC(track, (AliPID::EParticleType)AliPID::kPion);
    float nSigmaTPCK  = fPID->NumberOfSigmasTPC(track, (AliPID::EParticleType)AliPID::kKaon);
    float nSigmaTPCP  = fPID->NumberOfSigmasTPC(track, (AliPID::EParticleType)AliPID::kProton);
    //TOF
    float nSigmaTOFPi = fPID->NumberOfSigmasTOF(track, (AliPID::EParticleType)AliPID::kPion);
    float nSigmaTOFK  = fPID->NumberOfSigmasTOF(track, (AliPID::EParticleType)AliPID::kKaon);
    float nSigmaTOFP  = fPID->NumberOfSigmasTOF(track, (AliPID::EParticleType)AliPID::kProton);

    bool isPion   = ((TMath::Abs(nSigmaTPCPi) <= 2.) && (TMath::Abs(nSigmaTOFPi) <= 2.)); 
    bool isKaon   = ((TMath::Abs(nSigmaTPCK)  <= 2.) && (TMath::Abs(nSigmaTOFK)  <= 2.)); 
    bool isProton = ((TMath::Abs(nSigmaTPCP)  <= 2.) && (TMath::Abs(nSigmaTOFP)  <= 2.)); 

    int type = 0;
    if(charge > 0) {
      type = 999;
      if(isPion) type = 211;
      if(isKaon) type = 321;
      if(isProton) type = 2212;
    }
    if(charge < 0) {
      type = -999;
      if(isPion) type = -211;
      if(isKaon) type = -321;
      if(isProton) type = -2212;
    }
    if (type == 0) {
      AliError("Wrong Track!");
      return;
    }

    //ID
    int id = track->GetID();

    vecPhi.push_back(phi);
    vecEta.push_back(eta);
    vecPt.push_back(pt);
    vecType.push_back(type);
    vecID.push_back(id);

    //Qn
    double cos2phi = TMath::Cos(2 * phi);
    double sin2phi = TMath::Sin(2 * phi);
    double cos3phi = TMath::Cos(3 * phi);
    double sin3phi = TMath::Sin(3 * phi);
    sumCos2Phi += cos2phi;
    sumSin2Phi += sin2phi;
    sumCos3Phi += cos3phi;
    sumSin3Phi += sin3phi;
  }


  vector<double> vecLambdaPhi;
  vector<short>  vecLambdaPosID;
  vector<short>  vecLambdaNegID;

  vector<double> vecAntiLambdaPhi;
  vector<short>  vecAntiLambdaPosID;
  vector<short>  vecAntiLambdaNegID;

  int nV0sTot = fAOD->GetNumberOfV0s();
    for (int iV0 = 0; iV0 < nV0sTot; iV0++) {
      AliAODv0 *v0 = fAOD->GetV0(iV0);
      if(!v0) continue;
      //Basic kinematic variable
      double pt  = v0->Pt();
      double eta = v0->PseudoRapV0();
      double dcaToPV = v0->DcaV0ToPrimVertex();//DCA to Primary Vertex
      double CPA = v0->CosPointingAngle(vtx);//cosine pointing angle
      double dl  = v0->DecayLengthV0(vtx);
      hV0Pt->Fill(pt);
      hV0Eta->Fill(eta);
      hV0DcaToPrimVertex->Fill(dcaToPV);
      hV0CPA->Fill(CPA);
      hV0DecayLength->Fill(dl);

      //V0 cut
      if(!IsGoodV0(v0)) continue;
      //V0 daughters cut
      AliAODTrack *nTrack = dynamic_cast<AliAODTrack *>(v0->GetDaughter(1));
      AliAODTrack *pTrack = dynamic_cast<AliAODTrack *>(v0->GetDaughter(0));
      if(!(IsGoodDaughterTrack(nTrack)) || !(IsGoodDaughterTrack(pTrack))) continue;
      float nDcaPV = v0->DcaNegToPrimVertex();
      float pDcaPV = v0->DcaPosToPrimVertex();
      if( nDcaPV<fDaughtersDCAToPrimVtxMin || pDcaPV<fDaughtersDCAToPrimVtxMin) continue;

      //daughers Nsigma
      float nSigTPCPosProton = TMath::Abs(fPID->NumberOfSigmasTPC(pTrack,AliPID::kProton));//TPC p+
      float nSigTOFPosProton = TMath::Abs(fPID->NumberOfSigmasTPC(pTrack,AliPID::kProton));//TOF p+
      float nSigTPCNegPion   = TMath::Abs(fPID->NumberOfSigmasTPC(nTrack,AliPID::kPion));//TPC π-
      float nSigTOFNegPion   = TMath::Abs(fPID->NumberOfSigmasTPC(nTrack,AliPID::kPion));//TOF π-

      float nSigTPCPosPion   = TMath::Abs(fPID->NumberOfSigmasTPC(pTrack,AliPID::kPion));//TPC π+
      float nSigTOFPosPion   = TMath::Abs(fPID->NumberOfSigmasTPC(pTrack,AliPID::kPion));//TOF π+
      float nSigTPCNegProton = TMath::Abs(fPID->NumberOfSigmasTPC(nTrack,AliPID::kProton));//TPC p-
      float nSigTOFNegProton = TMath::Abs(fPID->NumberOfSigmasTPC(nTrack,AliPID::kProton));//TOF p-

      TVector2 Vt(v0->MomV0X(), v0->MomV0Y());
      double phi = Vt.Phi();
      short id_posDaughter = v0->GetPosID();
      short id_negDaughter = v0->GetNegID();

      double massLambda     = v0->MassLambda();
      double massAntiLambda = v0->MassAntiLambda();

      //Lambda PID 
      //Λ-->(p+)+(π-)
      if (nSigTPCPosProton < fDaughtersNsigma && //p+
          nSigTOFPosProton < fDaughtersNsigma && //p+
          nSigTPCNegPion   < fDaughtersNsigma && //π-
          nSigTOFNegPion   < fDaughtersNsigma)   //π-
      {
            hLambdaPt[0]->Fill(pt);
            hLambdaEta[0]->Fill(eta);
            hLambdaDcaToPrimVertex[0]->Fill(dcaToPV);
            hLambdaCPA[0]->Fill(CPA);
            hLambdaDecayLength[0]->Fill(dl);
            hLambdaMass[0]->Fill(massLambda);

          if (TMath::Abs(massLambda - fMassMean) < fLambdaMassCut)
          {
            hLambdaPt[1]->Fill(pt);
            hLambdaEta[1]->Fill(eta);
            hLambdaDcaToPrimVertex[1]->Fill(dcaToPV);
            hLambdaCPA[1]->Fill(CPA);
            hLambdaDecayLength[1]->Fill(dl);
            hLambdaMass[1]->Fill(massLambda);
            hLambdaMassCent[centBin]->Fill(massLambda);

            
            vecLambdaPhi.push_back(phi);
            vecLambdaPosID.push_back(id_posDaughter);
            vecLambdaNegID.push_back(id_negDaughter);
          }
        
      }

      //AntiLambda PID
      //(-Λ)-->(p-)+(π+)
      if (nSigTPCNegProton < fDaughtersNsigma &&
          nSigTOFNegProton < fDaughtersNsigma &&
          nSigTPCPosPion   < fDaughtersNsigma &&
          nSigTOFPosPion   < fDaughtersNsigma) 
      {
            hAntiLambdaPt[0]->Fill(pt);
            hAntiLambdaEta[0]->Fill(eta);
            hAntiLambdaDcaToPrimVertex[0]->Fill(dcaToPV);
            hAntiLambdaCPA[0]->Fill(CPA);
            hAntiLambdaDecayLength[0]->Fill(dl);
            hAntiLambdaMass[0]->Fill(massAntiLambda);

          if (TMath::Abs(massAntiLambda - fMassMean) < fLambdaMassCut)
          {

            hAntiLambdaEta[1]->Fill(eta);
            hAntiLambdaPt[1]->Fill(pt);
            hAntiLambdaDcaToPrimVertex[1]->Fill(dcaToPV);
            hAntiLambdaCPA[1]->Fill(CPA);
            hAntiLambdaDecayLength[1]->Fill(dl);
            hAntiLambdaMass[1]->Fill(massAntiLambda);
            hAntiLambdaMassCent[centBin]->Fill(massAntiLambda);

            vecAntiLambdaPhi.push_back(phi);
            vecAntiLambdaPosID.push_back(id_posDaughter);
            vecAntiLambdaNegID.push_back(id_negDaughter);
          }
        }
      
    }

  //RDM Plane
  double psi2RDM = gRandom->Uniform(0,TMath::Pi());
  hPsi2RDM -> Fill(runNumBin, centBin, psi2RDM);
  double psi3RDM = gRandom->Uniform(0,2/3*TMath::Pi());
  hPsi2RDM -> Fill(runNumBin, centBin, psi3RDM);
  hEvtCount->Fill(16);

  //TPC Plane
  TVector2 Q2TPC;
  Q2TPC.Set(sumCos2Phi,sumSin2Phi);
  double psi2TPC = Q2TPC.Phi()/2.;
  hPsi2TPC -> Fill(runNumBin, centBin, psi2TPC);
  TVector2 Q3TPC;
  Q3TPC.Set(sumCos3Phi,sumSin3Phi);
  double psi3TPC = Q3TPC.Phi()/3.;
  hPsi3TPC -> Fill(runNumBin, centBin, psi3TPC);
  hEvtCount->Fill(17);

  //VZERO Plane
  AliAODVZERO* fVZERO = fAOD->GetVZEROData();
  if (!fVZERO) return;
  double psi2V0A = -999, psi2V0C = -999, psi2V0M = -999; 
  if (!GetVZEROPlane(psi2V0A, psi2V0C, psi2V0M, fVZERO)) return;
  hEvtCount->Fill(18);

  //ZDC Plane
  AliAODZDC* fZDC = fAOD->GetZDCData();
  if (!fZDC) return;
  double psi1ZNA = -999, psi1ZNC = -999;
  if(!GetZDCPlane(psi1ZNA, psi1ZNC, fZDC)) return;
  hEvtCount->Fill(19);

  for (vector<double>::size_type iTrk = 0; iTrk < vecPt.size(); iTrk++) {
    short  id_1  = vecID[iTrk];
    double pt_1  = vecPt[iTrk];
    double eta_1 = vecEta[iTrk];
    double phi_1 = vecPhi[iTrk];
    int  type_1  = vecType[iTrk];
    for (vector<double>::size_type jTrk = 0; jTrk < vecPt.size(); jTrk++) {
      short  id_2  = vecID[jTrk];
      if(id_2 == id_1) continue;
      double pt_2  = vecPt[jTrk];
      double eta_2 = vecEta[jTrk];
      double phi_2 = vecPhi[jTrk];
      int   type_2 = vecType[jTrk];
      
      //TPC Event Plane deducted autocorrelation
      TVector2 q;
      double qx = TMath::Cos(2 * phi_1) + TMath::Cos(2 * phi_2);
      double qy = TMath::Sin(2 * phi_1) + TMath::Sin(2 * phi_2);
      q.Set(qx,qy);
      TVector2 QTmp = Q2TPC - q;
      double psiTPC_deAutoCor = QTmp.Phi()/mHarmonic;

      //delta
      double delta = TMath::Cos(phi_1 - phi_2);
      //gamma
      double gammaRDM  = TMath::Cos(phi_1 + phi_2 - 2 * psi2RDM);
      double gammaTPC  = TMath::Cos(phi_1 + phi_2 - 2 * psiTPC_deAutoCor);
      //coscos sinsin
      double coscosRDM = TMath::Cos(phi_1 - psi2RDM) * TMath::Cos(phi_2 - psi2RDM);
      double sinsinRDM = TMath::Sin(phi_1 - psi2RDM) * TMath::Sin(phi_2 - psi2RDM);
      double coscosTPC = TMath::Cos(phi_1 - psi2TPC) * TMath::Cos(phi_2 - psi2TPC);
      double sinsinTPC = TMath::Sin(phi_1 - psi2TPC) * TMath::Sin(phi_2 - psi2TPC);

      //Inclusive
      if(type_1 > 0 && type_2 > 0)
      {
        pDelta_hPos_hPos->Fill(cent,delta);
        pGammaRDM_hPos_hPos->Fill(cent,gammaRDM);
        pGammaTPC_hPos_hPos->Fill(cent,gammaTPC);
        pCosCosRDM_hPos_hPos->Fill(cent,coscosRDM);
        pSinSinRDM_hPos_hPos->Fill(cent,sinsinRDM);
        pCosCosTPC_hPos_hPos->Fill(cent,coscosTPC);
        pSinSinTPC_hPos_hPos->Fill(cent,sinsinTPC);
        pGammaTPCVsDeltaPt_hPos_hPos[centBin]  ->Fill(TMath::Abs(pt_1 - pt_2), gammaTPC);
        pGammaTPVVsMeanPt_hPos_hPos[centBin]   ->Fill((pt_1 + pt_2)/2., gammaTPC);
        pGammaTPCVsDeltaEta_hPos_hPos[centBin] ->Fill(TMath::Abs(eta_1 - eta_2), gammaTPC);
      }
      if(type_1 < 0 && type_2 < 0)
      {
        pDelta_hNeg_hNeg->Fill(cent,delta);
        pGammaRDM_hNeg_hNeg->Fill(cent,gammaRDM);
        pGammaTPC_hNeg_hNeg->Fill(cent,gammaTPC);
        pCosCosRDM_hNeg_hNeg->Fill(cent,coscosRDM);
        pSinSinRDM_hNeg_hNeg->Fill(cent,sinsinRDM);
        pCosCosTPC_hNeg_hNeg->Fill(cent,coscosTPC);
        pSinSinTPC_hNeg_hNeg->Fill(cent,sinsinTPC);
        pGammaTPCVsDeltaPt_hNeg_hNeg[centBin]  ->Fill(TMath::Abs(pt_1 - pt_2), gammaTPC);
        pGammaTPVVsMeanPt_hNeg_hNeg[centBin]   ->Fill((pt_1 + pt_2)/2., gammaTPC);
        pGammaTPCVsDeltaEta_hNeg_hNeg[centBin] ->Fill(TMath::Abs(eta_1 - eta_2), gammaTPC);
      }
      if(type_1 > 0 && type_2 < 0)
      {
        pDelta_hPos_hNeg->Fill(cent,delta);
        pGammaRDM_hPos_hNeg->Fill(cent,gammaRDM);
        pGammaTPC_hPos_hNeg->Fill(cent,gammaTPC);
        pCosCosRDM_hPos_hNeg->Fill(cent,coscosRDM);
        pSinSinRDM_hPos_hNeg->Fill(cent,sinsinRDM);
        pCosCosTPC_hPos_hNeg->Fill(cent,coscosTPC);
        pSinSinTPC_hPos_hNeg->Fill(cent,sinsinTPC);
        pGammaTPCVsDeltaPt_hPos_hNeg[centBin]  ->Fill(TMath::Abs(pt_1 - pt_2), gammaTPC);
        pGammaTPVVsMeanPt_hPos_hNeg[centBin]   ->Fill((pt_1 + pt_2)/2., gammaTPC);
        pGammaTPCVsDeltaEta_hPos_hNeg[centBin] ->Fill(TMath::Abs(eta_1 - eta_2), gammaTPC);
      }

      // //pion - Inclusive
      // // pi+ - h+
      // if(type_1 ==  211 && type_2 > 0)
      // {
      //   pDelta_pion_hPos->Fill(cent,delta);
      //   pGammaTPC_pion_hPos->Fill(cent,gammaTPC);
      //   pGammaRDM_pion_hPos->Fill(cent,gammaRDM);
      //   pGammaTPCVsPt_pion_hPos[centBin]->Fill(pt_1,gammaTPC);
      // }
      // // pi- - h+
      // if(type_1 == -211 && type_2 > 0)
      // {
      //   pDelta_antiPion_hPos->Fill(cent,delta);
      //   pGammaTPC_antiPion_hPos->Fill(cent,gammaTPC);
      //   pGammaRDM_antiPion_hPos->Fill(cent,gammaRDM);
      //   pGammaTPCVsPt_antiPion_hPos[centBin]->Fill(pt_1,gammaTPC);
      // }
      // // pi+ - h-
      // if(type_1 ==  211 && type_2 < 0)
      // {
      //   pDelta_pion_hNeg->Fill(cent,delta);
      //   pGammaTPC_pion_hNeg->Fill(cent,gammaTPC);
      //   pGammaRDM_pion_hNeg->Fill(cent,gammaRDM);
      //   pGammaTPCVsPt_pion_hNeg[centBin]->Fill(pt_1,gammaTPC);
      // }
      // // pi- - h-
      // if(type_1 == -211 && type_2 < 0)
      // {
      //   pDelta_antiPion_hNeg->Fill(cent,delta);
      //   pGammaTPC_antiPion_hNeg->Fill(cent,gammaTPC);
      //   pGammaRDM_antiPion_hNeg->Fill(cent,gammaRDM);
      //   pGammaTPCVsPt_antiPion_hNeg[centBin]->Fill(pt_1,gammaTPC);
      // }
      // //kaon - Inclusive
      // // k+ - h+
      // if(type_1 ==  321 && type_2 > 0)
      // {
      //   pDelta_kaon_hPos->Fill(cent,delta);
      //   pGammaTPC_kaon_hPos->Fill(cent,gammaTPC);
      //   pGammaRDM_kaon_hPos->Fill(cent,gammaRDM);
      //   pGammaTPCVsPt_kaon_hPos[centBin]->Fill(pt_1,gammaTPC);
      // }
      // // k- - h+
      // if(type_1 == -321 && type_2 > 0)
      // {
      //   pDelta_antiKaon_hPos->Fill(cent,delta);
      //   pGammaTPC_antiKaon_hPos->Fill(cent,gammaTPC);
      //   pGammaRDM_antiKaon_hPos->Fill(cent,gammaRDM);
      //   pGammaTPCVsPt_antiKaon_hPos[centBin]->Fill(pt_1,gammaTPC);
      // }
      // // k+ - h-
      // if(type_1 ==  321 && type_2 < 0)
      // {
      //   pDelta_kaon_hNeg->Fill(cent,delta);
      //   pGammaTPC_kaon_hNeg->Fill(cent,gammaTPC);
      //   pGammaRDM_kaon_hNeg->Fill(cent,gammaRDM);
      //   pGammaTPCVsPt_kaon_hNeg[centBin]->Fill(pt_1,gammaTPC);
      // }
      // // k- - h-
      // if(type_1 == -321 && type_2 < 0)
      // {
      //   pDelta_antiKaon_hNeg->Fill(cent,delta);
      //   pGammaTPC_antiKaon_hNeg->Fill(cent,gammaTPC);
      //   pGammaRDM_antiKaon_hNeg->Fill(cent,gammaRDM);
      //   pGammaTPCVsPt_antiKaon_hNeg[centBin]->Fill(pt_1,gammaTPC);
      // }
      //proton - Inclusive
      // p+ - h+
      if(type_1 ==  2212 && type_2 > 0)
      {
        pDelta_proton_hPos->Fill(cent,delta);
        pGammaTPC_proton_hPos->Fill(cent,gammaTPC);
        pGammaRDM_proton_hPos->Fill(cent,gammaRDM);
        pGammaTPCVsPt_proton_hPos[centBin]->Fill(pt_1,gammaTPC);
      }// p- - h+
      if(type_1 == -2212 && type_2 > 0)
      {
        pDelta_antiProton_hPos->Fill(cent,delta);
        pGammaTPC_antiProton_hPos->Fill(cent,gammaTPC);
        pGammaRDM_antiProton_hPos->Fill(cent,gammaRDM);
        pGammaTPCVsPt_antiProton_hPos[centBin]->Fill(pt_1,gammaTPC);
      }// p+ - h-
      if(type_1 ==  2212 && type_2 < 0) 
      {
        pDelta_proton_hNeg->Fill(cent,delta);
        pGammaTPC_proton_hNeg->Fill(cent,gammaTPC);
        pGammaRDM_proton_hNeg->Fill(cent,gammaRDM);
        pGammaTPCVsPt_proton_hNeg[centBin]->Fill(pt_1,gammaTPC);
      }// p- - h-
      if(type_1 == -2212 && type_2 < 0) 
      {
        pDelta_antiProton_hNeg->Fill(cent,delta);
        pGammaTPC_antiProton_hNeg->Fill(cent,gammaTPC);
        pGammaRDM_antiProton_hNeg->Fill(cent,gammaRDM);
        pGammaTPCVsPt_antiProton_hNeg[centBin]->Fill(pt_1,gammaTPC);
      }

      // //pion - pion
      // //pi+ - pi+
      // if(type_1 == 211 && type_2 == 211)
      // {
      //   pDelta_pion_pion->Fill(cent,delta);
      //   pGammaTPC_pion_pion->Fill(cent,gammaTPC);
      //   pGammaRDM_pion_pion->Fill(cent,gammaRDM);
      // }
      // //pi- - pi-
      // if(type_1 == -211 && type_2 == -211)
      // {
      //   pDelta_antiPion_antiPion->Fill(cent,delta);
      //   pGammaTPC_antiPion_antiPion->Fill(cent,gammaTPC);
      //   pGammaRDM_antiPion_antiPion->Fill(cent,gammaRDM);
      // }
      // //pi+ - pi-
      // if(type_1 == 211 && type_2 == -211)
      // {
      //   pDelta_pion_antiPion->Fill(cent,delta);
      //   pGammaTPC_pion_antiPion->Fill(cent,gammaTPC);
      //   pGammaRDM_pion_antiPion->Fill(cent,gammaRDM);
      // }
      // //kaon - kaon
      // //k+ - k+
      // if(type_1 == 321 && type_2 == 321)
      // {
      //   pDelta_kaon_kaon->Fill(cent,delta);
      //   pGammaTPC_kaon_kaon->Fill(cent,gammaTPC);
      //   pGammaRDM_kaon_kaon->Fill(cent,gammaRDM);
      // }
      // //k- - k-
      // if(type_1 == -321 && type_2 == -321)      
      // {
      //   pDelta_antiKaon_antiKaon->Fill(cent,delta);
      //   pGammaTPC_antiKaon_antiKaon->Fill(cent,gammaTPC);
      //   pGammaRDM_antiKaon_antiKaon->Fill(cent,gammaRDM);
      // }
      // //k+ - k-
      // if(type_1 == 321 && type_2 == -321)      
      // {
      //   pDelta_kaon_antiKaon->Fill(cent,delta);
      //   pGammaTPC_kaon_antiKaon->Fill(cent,gammaTPC);
      //   pGammaRDM_kaon_antiKaon->Fill(cent,gammaRDM);
      // }
      //proton - proton
      //p+ - p+
      if(type_1 == 2212 && type_2 == 2212)      
      {
        pDelta_proton_proton->Fill(cent,delta);
        pGammaTPC_proton_proton->Fill(cent,gammaTPC);
        pGammaRDM_proton_proton->Fill(cent,gammaRDM);
      }
      //p- - p-
      if(type_1 == -2212 && type_2 == -2212)      
      {
        pDelta_antiProton_antiProton->Fill(cent,delta);
        pGammaTPC_antiProton_antiProton->Fill(cent,gammaTPC);
        pGammaRDM_antiProton_antiProton->Fill(cent,gammaRDM);
      }
      //p+ - p-
      if(type_1 == 2212 && type_2 == -2212)      
      {
        pDelta_proton_antiProton->Fill(cent,delta);
        pGammaTPC_proton_antiProton->Fill(cent,gammaTPC);
        pGammaRDM_proton_antiProton->Fill(cent,gammaRDM);
      }


      // //pion - kaon
      // // pi+ - k+
      // if(type_1 ==  211 && type_2 == 321)
      // {
      //   pDelta_pion_kaon->Fill(cent,delta);
      //   pGammaTPC_pion_kaon->Fill(cent,gammaTPC);
      //   pGammaRDM_pion_kaon->Fill(cent,gammaRDM);
      // }
      // // pi- - k+
      // if(type_1 == -211 && type_2 == 321)      
      // {
      //   pDelta_antiPion_kaon->Fill(cent,delta);
      //   pGammaTPC_antiPion_kaon->Fill(cent,gammaTPC);
      //   pGammaRDM_antiPion_kaon->Fill(cent,gammaRDM);
      // }
      // // pi+ - k-
      // if(type_1 ==  211 && type_2 == -321)       
      // {
      //   pDelta_pion_antiKaon->Fill(cent,delta);
      //   pGammaTPC_pion_antiKaon->Fill(cent,gammaTPC);
      //   pGammaRDM_pion_antiKaon->Fill(cent,gammaRDM);
      // }
      // // pi- - k-
      // if(type_1 == -211 && type_2 == -321)       
      // {
      //   pDelta_antiPion_antiKaon->Fill(cent,delta);
      //   pGammaTPC_antiPion_antiKaon->Fill(cent,gammaTPC);
      //   pGammaRDM_antiPion_antiKaon->Fill(cent,gammaRDM);
      // }
      // //kaon - proton
      // // k+ - p+
      // if(type_1 ==  321 && type_2 == 2212)
      // {
      //   pDelta_kaon_proton->Fill(cent,delta);
      //   pGammaTPC_kaon_proton->Fill(cent,gammaTPC);
      //   pGammaRDM_kaon_proton->Fill(cent,gammaRDM);
      // } 
      // // k- - p+
      // if(type_1 == -321 && type_2 == 2212)
      // {
      //   pDelta_antiKaon_proton->Fill(cent,delta);
      //   pGammaTPC_antiKaon_proton->Fill(cent,gammaTPC);
      //   pGammaRDM_antiKaon_proton->Fill(cent,gammaRDM);
      // }    
      // // k+ - p-
      // if(type_1 ==  321 && type_2 == -2212)
      // {
      //   pDelta_kaon_antiProton->Fill(cent,delta);
      //   pGammaTPC_kaon_antiProton->Fill(cent,gammaTPC);
      //   pGammaRDM_kaon_antiProton->Fill(cent,gammaRDM);
      // } 
      // // k- - p-
      // if(type_1 == -321 && type_2 == -2212)
      // {
      //   pDelta_antiKaon_antiProton->Fill(cent,delta);
      //   pGammaTPC_antiKaon_antiProton->Fill(cent,gammaTPC);
      //   pGammaRDM_antiKaon_antiProton->Fill(cent,gammaRDM);
      // } 
      // //proton - pi
      // // p+ - pi+
      // if(type_1 ==  2212 && type_2 > 211)
      // {
      //   pDelta_proton_pion->Fill(cent,delta);
      //   pGammaTPC_proton_pion->Fill(cent,gammaTPC);
      //   pGammaRDM_proton_pion->Fill(cent,gammaRDM);
      // }
      // // p- - pi+
      // if(type_1 == -2212 && type_2 > 211)
      // {
      //   pDelta_antiProton_pion->Fill(cent,delta);
      //   pGammaTPC_antiProton_pion->Fill(cent,gammaTPC);
      //   pGammaRDM_antiProton_pion->Fill(cent,gammaRDM);
      // }
      // // p+ - pi-
      // if(type_1 ==  2212 && type_2 == -211)
      // {
      //   pDelta_proton_antiPion->Fill(cent,delta);
      //   pGammaTPC_proton_antiPion->Fill(cent,gammaTPC);
      //   pGammaRDM_proton_antiPion->Fill(cent,gammaRDM);
      // }
      // // p- - pi-
      // if(type_1 == -2212 && type_2 == -211)
      // {
      //   pDelta_antiProton_antiPion->Fill(cent,delta);
      //   pGammaTPC_antiProton_antiPion->Fill(cent,gammaTPC);
      //   pGammaRDM_antiProton_antiPion->Fill(cent,gammaRDM);
      // }
    }

    //Lambda
    for (vector<double>::size_type jLambda = 0; jLambda < vecLambdaPhi.size(); jLambda++)
    {
      double phi_lambda =  vecLambdaPhi[jLambda];
      short  id_posDaughter = vecLambdaPosID[jLambda];
      short  id_negDaughter = vecLambdaNegID[jLambda];
      if(id_1 == id_posDaughter || id_1 == id_negDaughter) continue;

      //TPC Event Plane deducted autocorrelation
      double qx = TMath::Cos(mHarmonic * phi_1);
      double qy = TMath::Sin(mHarmonic * phi_1);
      if(find(vecID.begin(),vecID.end(), id_posDaughter) != vecID.end()|| 
         find(vecID.begin(),vecID.end(), id_negDaughter) != vecID.end()) 
      {
        qx += TMath::Cos(mHarmonic * vecPhi[iTrk]);
        qy += TMath::Cos(mHarmonic * vecPhi[iTrk]);
      }
      TVector2 q;
      q.Set(qx,qy);
      TVector2 QTmp = Q2TPC - q;
      double psiTPC_deAutoCor = QTmp.Phi()/mHarmonic;

      //delta
      double delta = TMath::Cos(phi_lambda - phi_1);
      //gamma
      double gammaRDM  = TMath::Cos(phi_lambda + phi_1 - 2 * psi2RDM);
      double gammaTPC  = TMath::Cos(phi_lambda + phi_1 - 2 * psiTPC_deAutoCor);

      //lambda - h+
      if(type_1 > 0)
      {
        pDelta_lambda_hPos->Fill(cent,delta);
        pGammaTPC_lambda_hPos->Fill(cent,gammaTPC);
        pGammaRDM_lambda_hPos->Fill(cent,gammaRDM);
      }
      //lambda - h-
      if(type_1 < 0)
      {
        pDelta_lambda_hNeg->Fill(cent,delta);
        pGammaTPC_lambda_hNeg->Fill(cent,gammaTPC);
        pGammaRDM_lambda_hNeg->Fill(cent,gammaRDM);
      }
      
      //lambda - p+
      if(type_1 == 2212)
      {
        pDelta_lambda_proton->Fill(cent,delta);
        pGammaTPC_lambda_proton->Fill(cent,gammaTPC);
        pGammaRDM_lambda_proton->Fill(cent,gammaRDM);
      }
      //lambda - p-
      if(type_1 == -2212)
      {
        pDelta_lambda_antiProton->Fill(cent,delta);
        pGammaTPC_lambda_antiProton->Fill(cent,gammaTPC);
        pGammaRDM_lambda_antiProton->Fill(cent,gammaRDM);
      }
    }

    //Anti-Lambda
    for (vector<double>::size_type jAntiLambda = 0; jAntiLambda < vecAntiLambdaPhi.size(); jAntiLambda++)
    {
      double phi_antiLambda =  vecAntiLambdaPhi[jAntiLambda];
      short  id_posDaughter = vecAntiLambdaPosID[jAntiLambda];
      short  id_negDaughter = vecAntiLambdaNegID[jAntiLambda];
      if(id_1 == id_posDaughter || id_1 == id_negDaughter) continue;

      //TPC Event Plane deducted autocorrelation
      double qx = TMath::Cos(mHarmonic * phi_1);
      double qy = TMath::Sin(mHarmonic * phi_1);
      if(find(vecID.begin(),vecID.end(), id_posDaughter) != vecID.end()|| 
         find(vecID.begin(),vecID.end(), id_negDaughter) != vecID.end()) 
      {
        qx += TMath::Cos(mHarmonic * vecPhi[iTrk]);
        qy += TMath::Cos(mHarmonic * vecPhi[iTrk]);
      }
      TVector2 q;
      q.Set(qx,qy);
      TVector2 QTmp = Q2TPC - q;
      double psiTPC_deAutoCor = QTmp.Phi()/mHarmonic;

      //delta
      double delta = TMath::Cos(phi_antiLambda - phi_1);
      //gamma
      double gammaRDM  = TMath::Cos(phi_antiLambda + phi_1 - 2 * psi2RDM);
      double gammaTPC  = TMath::Cos(phi_antiLambda + phi_1 - 2 * psiTPC_deAutoCor);

      //antiLambda - h+
      if(type_1 > 0)
      {
        pDelta_antiLambda_hPos->Fill(cent,delta);
        pGammaTPC_antiLambda_hPos->Fill(cent,gammaTPC);
        pGammaRDM_antiLambda_hPos->Fill(cent,gammaRDM);
      }
      //antiLambda - h-
      if(type_1 < 0)
      {
        pDelta_antiLambda_hNeg->Fill(cent,delta);
        pGammaTPC_antiLambda_hNeg->Fill(cent,gammaTPC);
        pGammaRDM_antiLambda_hNeg->Fill(cent,gammaRDM);
      }
      
      //antiLambda - p+
      if(type_1 == 2212)
      {
        pDelta_antiLambda_proton->Fill(cent,delta);
        pGammaTPC_antiLambda_proton->Fill(cent,gammaTPC);
        pGammaRDM_antiLambda_proton->Fill(cent,gammaRDM);
      }
      //antiLambda - p-
      if(type_1 == -2212)
      {
        pDelta_antiLambda_antiProton->Fill(cent,delta);
        pGammaTPC_antiLambda_antiProton->Fill(cent,gammaTPC);
        pGammaRDM_antiLambda_antiProton->Fill(cent,gammaRDM);
      }
    }
  }

  hEvtCount->Fill(20);
  PostData(1,fOutputList);
}

//---------------------------------------------------
int AliAnalysisTaskChiralVorticalEffect::GetRunNumBin(int runNum)
{
  int runNumBin=-1;
  // 10h
  int runNumList[91]={
    139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328, 139314, 139310,
    139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029, 139028, 138872, 138871, 138870,
    138837, 138732, 138730, 138666, 138662, 138653, 138652, 138638, 138624, 138621, 138583, 138582,
    138579, 138578, 138534, 138469, 138442, 138439, 138438, 138396, 138364, 138275, 138225, 138201,
    138197, 138192, 138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693,
    137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595, 137549, 137546, 137544, 137541,
    137539, 137531, 137530, 137443, 137441, 137440, 137439, 137434, 137432, 137431, 137430, 137243,
    137236, 137235, 137232, 137231, 137230, 137162, 137161
  };
  // 11h
  // int runNumList[39]={170387,170040,170268,170228,170207,169838,170159,170204,170311,170084,
  //                     169835,170088,170593,170203,170270,169846,170163,170388,170155,170083,
  //                     170572,169837,169855,170306,170269,170089,170309,170091,170081,170230,
  //                     170085,170315,170027,170193,170312,170313,170308,169858,169859};
  for (int i = 0; i < 91; ++i) {
    if (runNum==runNumList[i]) {runNumBin=i; break;}
    else continue;
  }
  return runNumBin;
}


bool AliAnalysisTaskChiralVorticalEffect::IsGoodV0(AliAODv0 *aodV0) {
  if (!aodV0) {
    AliError(Form("ERROR: Could not retrieve aodV0"));
    return kFALSE;
  }
  //* Offline reconstructed V0 only
  if ( aodV0->GetOnFlyStatus() ) return kFALSE;
  //* Get daughters and check them     
  AliAODTrack *myTrackNegTest = dynamic_cast<AliAODTrack *>(aodV0->GetDaughter(1));
  AliAODTrack *myTrackPosTest = dynamic_cast<AliAODTrack *>(aodV0->GetDaughter(0));
  if (!myTrackPosTest || !myTrackNegTest) {
    Printf("strange analysis::UserExec:: Error:Could not retreive one of the daughter track\n");
    return kFALSE;
  }
  //* Unlike signs of daughters
  if ( myTrackNegTest->Charge() == myTrackPosTest->Charge() ) return kFALSE;
  //* Cosinus of pointing angle      
  double dCPA = aodV0->CosPointingAngle(vtx);
  //* cut on Cosinus of pointing angle
  if ( dCPA < fV0CPAMin ) return kFALSE;
  //* DCA of V0
  double dV0Dca = aodV0->DcaV0ToPrimVertex();
  if ( TMath::Abs(dV0Dca) > fV0DCAToPrimVtxMax ) return kFALSE;
  //* V0 parh length before decay
  double dDecayLength = aodV0->DecayLengthV0(vtx);
  if ( dDecayLength > fV0DecayLengthMax ) return kFALSE;
  if ( dDecayLength < fV0DecayLengthMin ) return kFALSE;
  //* DCA between daughters
  double dDCA = aodV0->DcaV0Daughters();
  if ( dDCA > fV0DcaBetweenDaughtersMax ) return kFALSE;
  double dPt = aodV0->Pt();
  if( dPt < fV0PtMin ) return kFALSE;
  double dRapidity = aodV0->RapLambda();
  if( TMath::Abs(dRapidity) > fV0RapidityMax ) return kFALSE;
  return kTRUE;
}

bool AliAnalysisTaskChiralVorticalEffect::IsGoodDaughterTrack(const AliAODTrack *t)
{
  //* TPC refit
  if ( !t->IsOn(AliAODTrack::kTPCrefit) ) return kFALSE;
  //* Maximum value of transverse momentum 
  double dPt = t->Pt();
  if (dPt > fDaughtersPtMax) return kFALSE;
  //* Maximum value of pseudorapidity
  double dEta = t->Eta();
  if (TMath::Abs(dEta) > fDaughtersEtaMax) return kFALSE;
  //* Minimum number of clusters
  float nCrossedRowsTPC = t->GetTPCClusterInfo(2,1);
  if (nCrossedRowsTPC < fDaughtersTPCNclsMin) return kFALSE;
  //* Findable clusters > 0
  int findable = t->GetTPCNclsF();
  if (findable <= 0) return kFALSE;
  //* [number of crossed rows]>0.8 * [number of findable clusters].
  if (nCrossedRowsTPC/findable < 0.8) return kFALSE;
  return kTRUE;
}


bool AliAnalysisTaskChiralVorticalEffect::GetZDCPlaneLucas(double &psiA , double &psiC, AliAODZDC* fZDC)
{
  const auto min_energy = 1.e-6;
  const std::array<double, 5> phiZNA{0, 7.*TMath::PiOver4(), 5.*TMath::PiOver4(), TMath::PiOver4(), 3.*TMath::PiOver4()};
  const std::array<double, 5> phiZNC{0, 5.*TMath::PiOver4(), 7.*TMath::PiOver4(), 3.*TMath::PiOver4(), TMath::PiOver4()};

  auto energyZNA = fZDC->GetZNATowerEnergy();
  auto energyZNC = fZDC->GetZNCTowerEnergy();

  double sumEnergyZNA = 0;
  double sumEnergyZNC = 0;
  double QxZNA = 0., QyZNA = 0., QxZNC = 0., QyZNC = 0.;

  for (int iTower = 0; iTower < 4; ++iTower) {
    if (energyZNA[iTower] < min_energy || energyZNC[iTower] < min_energy) return;
    QxZNA += energyZNA[iTower+1] * TMath::Cos(phiZNA[iTower]);
    QyZNA += energyZNA[iTower+1] * TMath::Cos(phiZNA[iTower]);
    QxZNC += energyZNC[iTower+1] * TMath::Cos(phiZNC[iTower]);
    QyZNC += energyZNC[iTower+1] * TMath::Cos(phiZNC[iTower]);
    sumEnergyZNA += energyZNA[iTower+1];
    sumEnergyZNC += energyZNC[iTower+1];
  }
  QxZNA /= sumEnergyZNA;
  QyZNA /= sumEnergyZNA;
  QxZNC /= sumEnergyZNC;
  QyZNC /= sumEnergyZNC;

  TVector2 QZNA(QxZNA,QyZNA);
  TVector2 QZNC(QxZNC,QyZNC);

  psiA = QZNA.Phi();
  psiC = QZNC.Phi();

  hPsi1ZNA[0]->Fill(runNumBin,centBin,psiA);
  hPsi1ZNC[0]->Fill(runNumBin,centBin,psiC);

  double sumEnergyZNAEqualize = -999;
  double sumEnergyZNCEqualize = -999;
  double sumEnergyZNARecenter = -999;
  double sumEnergyZNCRecenter = -999;
  
  //TODO
  TH1D* hVtxMeanCentRun[2];
  double vtxMean[2];
  double vtxSigma[2];
  for(int i=0;i<2;i++) {
    vtxMean[i] = hVtxMeanCentRun[i]->GetBinContent(runNumBin);
    vtxSigma[i] = hVtxMeanCentRun[i]->GetBinError(runNumBin);
  }
  double vtxNSigma[2];
  for(int i=0;i<2;i++) {
    vtxNSigma[i] = (vtx[i]-vtxMean[i])/vtxSigma[i];
  }

  if(bZDCEqualize) {
    TVector2 QZNAEqualize = GetEqualizeQVector(sumEnergyZNAEqualize, energyZNA + 1, phiZNA.data() + 1);
    TVector2 QZNCEqualize = GetEqualizeQVector(sumEnergyZNCEqualize, energyZNC + 1, phiZNC.data() + 1);
    std::vector<double> variables = {cent, vtxNSigma[0], vtxNSigma[1], vtx[2]};
    if(bZDCEqualize && bZDCRecenter) {
      TVector2 QZNARecenter = GetRecenterQVector(sumEnergyZNARecenter, QZNAEqualize, variables.data());
      TVector2 QZNCRecenter = GetRecenterQVector(sumEnergyZNCRecenter, QZNAEqualize, variables.data());
    }
  }
    //TODO
    // if(bZDCShift) {
    //   // shifting
    //   p_shiftCoeffSinA_ = (TProfile2D*)listEPShifting_->FindObject(Form("shiftCoeffSinA_%i", runNum));
    //   p_shiftCoeffCosA_ = (TProfile2D*)listEPShifting_->FindObject(Form("shiftCoeffCosA_%i", runNum));
    //   p_shiftCoeffSinC_ = (TProfile2D*)listEPShifting_->FindObject(Form("shiftCoeffSinC_%i", runNum));
    //   p_shiftCoeffCosC_ = (TProfile2D*)listEPShifting_->FindObject(Form("shiftCoeffCosC_%i", runNum));
    //   p_shiftCoeffSinFull_ = (TProfile2D*)listEPShifting_->FindObject(Form("shiftCoeffSinFull_%i", runNum));
    //   p_shiftCoeffCosFull_ = (TProfile2D*)listEPShifting_->FindObject(Form("shiftCoeffCosFull_%i", runNum));

    //   for (int i = 1; i <= 20; ++i) {
    //     double shiftSinA = p_shiftCoeffSinA_->GetBinContent(centBin+1,i);
    //     double shiftCosA = p_shiftCoeffCosA_->GetBinContent(centBin+1,i);
    //     double shiftSinC = p_shiftCoeffSinC_->GetBinContent(centBin+1,i);
    //     double shiftCosC = p_shiftCoeffCosC_->GetBinContent(centBin+1,i);
    //     double shiftSinFull = p_shiftCoeffSinFull_->GetBinContent(centBin+1,i);
    //     double shiftCosFull = p_shiftCoeffCosFull_->GetBinContent(centBin+1,i);
    //     psiA += (1./i)*(-shiftSinA*cos(i*psiA)+shiftCosA*sin(i*psiA));
    //     psiC += (1./i)*(-shiftSinC*cos(i*psiC)+shiftCosC*sin(i*psiC));
    //     psiFull += (1./i)*(-shiftSinFull*cos(i*psiFull)+shiftCosFull*sin(i*psiFull));
    //   }

    //   hPsi1ZNA[2]->Fill(runNumBin,centBin,psiA);
    //   hPsi1ZNC[2]->Fill(runNumBin,centBin,psiC);
    //   hPsi1ZDC[2]->Fill(runNumBin,centBin,psiFull);
    // }
}


TVector2 AliAnalysisTaskChiralVorticalEffect::GetEqualizeQVector(double &energySum, const double *rawEnergies, const double *phis) 
{
    //TODO
    TH1D*fAverageGainIn;
    const double min_energy = 1.e-6;
    unsigned int i_group = 0;
    double energyEq[4];
    for (int i = 0; i < 4; ++i) {
      if (i == i_group) i_group = i; 
      energyEq[i] = rawEnergies[i] / fAverageGainIn->GetBinContent(i_group+1) * fAverageGainIn->GetBinContent(i+1);
    }
    fAverageGainIn = NULL;
    delete fAverageGainIn;
  
    double QxEq = 0., QyEq = 0.;
    for (int iTower = 0; iTower < 4; ++iTower) {
      if (energyEq[iTower] < min_energy) return;
      QxEq += energyEq[iTower] * TMath::Cos(phis[iTower]);
      QyEq += energyEq[iTower] * TMath::Cos(phis[iTower]);
      energySum += energyEq[iTower];
    }
    QxEq /= energySum;
    QyEq /= energySum;

    TVector2 QEq;
    QEq.Set(QxEq,QyEq);
    return QEq;
}


// void AliAnalysisTaskChiralVorticalEffect::MakeZDCCalibFiles(AliAODZDC* fZDC) 
// {
//   for (auto i = 0u; i < 4; ++i) {
//     fAverageGainOut->Fill(i, channel_energies[i]);
//     fGainQAbefore->Fill(i, channel_energies[i]);
//   }

// }


TVector2 AliAnalysisTaskChiralVorticalEffect::GetRecenterQVector(double &energySum, const TVector2 QVector, double* variables) 
{
  //TODO
  THnD*fSumWin;
  THnD*fSumXin;
  THnD*fSumYin;
  int fMinEntries;
  std::vector<Int_t> bins(fSumWin->GetNdimensions());
  for (auto i = 0; i < fSumWin->GetNdimensions(); ++i) {
    bins[i] = fSumWin->GetAxis(i)->FindBin(variables[i]);
  }
  auto bin  = fSumWin->GetBin(bins.data());
  auto n    = fSumWin->GetBinContent(bin);
  auto xsum = fSumXin->GetBinContent(bin);
  auto ysum = fSumYin->GetBinContent(bin);
      
  auto xmean = xsum / n;
  auto ymean = ysum / n;

  //auto xwidth = 1.;
  //auto ywidth = 1.;
  // if (fEqualization) {
  //   auto xerror = fSumXin->GetBinError2(bin);
  //   auto yerror = fSumYin->GetBinError2(bin);
  //   xwidth = sqrt(std::fabs(xerror/n - xmean*xmean));
  //   ywidth = sqrt(std::fabs(yerror/n - ymean*ymean));
  // }
  TVector2 QMean;
  QMean.Set(xmean,ymean);
  return QVector - QMean;
}


//---------------------------------------------------
void AliAnalysisTaskChiralVorticalEffect::Terminate(Option_t *) 
{

}


bool AliAnalysisTaskChiralVorticalEffect::GetVZEROPlane(double &psiA, double &psiC, AliAODVZERO* fVZERO);
{
  AliEventplane* V0EP = fAOD->GetEventplane();
  double psiAraw = TVector2::Phi_0_2pi(V0EP->GetEventplane("V0A", fAOD)); if (psiA>TMath::Pi()) psiA-= TMath::Pi();
  double psiCraw = TVector2::Phi_0_2pi(V0EP->GetEventplane("V0C", fAOD)); if (psiC>TMath::Pi()) psiC-= TMath::Pi();

  if(bVZEROEqualize) {
    //V0 info
    double Qxa2 = 0., Qya2 = 0.;
    double Qxc2 = 0., Qyc2 = 0.; 
    double sumMa = 0., sumMc = 0., sumMm = 0.;

    for (Int_t iV0 = 0; iV0 < 64; iV0++) {
      double phiV0 = TMath::PiOver4()*(0.5 + iV0 % 8);
      double multv0 = aodV0->GetMultiplicity(iV0);
      if (iV0 < 32) {
        double multCorC = -10;
        if (iV0 < 8) multCorC = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(1);
        else if (iV0 >= 8  && iV0 < 16) multCorC = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(9);
        else if (iV0 >= 16 && iV0 < 24) multCorC = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(17);
        else if (iV0 >= 24 && iV0 < 32) multCorC = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(25);
        if (multCorC < 0) continue;
        Qxc2 += TMath::Cos(2.*phiV0) * multCorC;
        Qyc2 += TMath::Sin(2.*phiV0) * multCorC;
        sumMc = sumMc + multCorC;
      } else {
        double multCorA = -10;
        if (iV0 >= 32 && iV0 < 40) multCorA = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(33);
        else if (iV0 >= 40 && iV0 < 48) multCorA = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(41);
        else if (iV0 >= 48 && iV0 < 56) multCorA = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(49);
        else if (iV0 >= 56 && iV0 < 64) multCorA = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(57);
        if (multCorA < 0) continue;
        Qxa2 += TMath::Cos(2.*phiV0) * multCorA;
        Qya2 += TMath::Sin(2.*phiV0) * multCorA;
        sumMa = sumMa + multCorA;
      }
    }
    if (sumMa <= 0 || sumMc <= 0) return;
    psiA = TMath::ATan2(Qya2, Qxa2)/2.;
    psiC = TMath::ATan2(Qyc2, Qxc2)/2.;

    if(bVZERORecenter) {
      double Qxa2Cor = (Qxa2 - fQx2mV0A->GetBinContent(iCentSPD+1))/fQx2sV0A->GetBinContent(iCentSPD+1);
      double Qya2Cor = (Qya2 - fQy2mV0A->GetBinContent(iCentSPD+1))/fQy2sV0A->GetBinContent(iCentSPD+1);
      double Qxc2Cor = (Qxc2 - fQx2mV0C->GetBinContent(iCentSPD+1))/fQx2sV0C->GetBinContent(iCentSPD+1);
      double Qyc2Cor = (Qyc2 - fQy2mV0C->GetBinContent(iCentSPD+1))/fQy2sV0C->GetBinContent(iCentSPD+1);

      psiA = TMath::ATan2(Qya2Cor, Qxa2Cor)/2.;
      psiC = TMath::ATan2(Qyc2Cor, Qxc2Cor)/2.;
    }
  }
  // if (bQAVZERO) {
  //     fPsiA[centrCode]->Fill(evPlAngV0A);
  //     fPsiC[centrCode]->Fill(evPlAngV0C);
  //     fPsiAvsPsiC[centrCode]->Fill(evPlAngV0A, evPlAngV0C);

  //     fQxavsV0Bef->Fill(centV0, Qxa2);
  //     fQyavsV0Bef->Fill(centV0, Qya2);
  //     fQxcvsV0Bef->Fill(centV0, Qxc2);
  //     fQycvsV0Bef->Fill(centV0, Qyc2);

  //     fQxavsVtxZBef->Fill(vtxZ, Qxa2);
  //     fQyavsVtxZBef->Fill(vtxZ, Qya2);
  //     fQxcvsVtxZBef->Fill(vtxZ, Qxc2);
  //     fQycvsVtxZBef->Fill(vtxZ, Qyc2);

  //     fQxavsV0Aft->Fill(centV0, Qxa2Cor);
  //     fQyavsV0Aft->Fill(centV0, Qya2Cor);
  //     fQxcvsV0Aft->Fill(centV0, Qxc2Cor);
  //     fQycvsV0Aft->Fill(centV0, Qyc2Cor);

  //     fQxavsVtxZAft->Fill(vtxZ, Qxa2Cor);
  //     fQyavsVtxZAft->Fill(vtxZ, Qya2Cor);
  //     fQxcvsVtxZAft->Fill(vtxZ, Qxc2Cor);
  //     fQycvsVtxZAft->Fill(vtxZ, Qyc2Cor);
  //   }
  return kTRUE;
}


bool AliAnalysisTaskChiralVorticalEffect::GetZDCPlaneJacopo(double &psiA , double &psiC, AliAODZDC* fZDC) {

  if(RunNum!=fCachedRunNum) {
    fbFlagIsPosMagField = kFALSE;
    Int_t dRun15hPos[] = {246390, 246391, 246392, 246994, 246991, 246989, 246984, 246982, 246980, 246948, 246945, 246928, 246851, 246847, 246846, 246845, 246844, 246810, 246809, 246808, 246807, 246805, 246804, 246766, 246765, 246763, 246760, 246759, 246758, 246757, 246751, 246750, 246495, 246493, 246488, 246487, 246434, 246431, 246428, 246424};
    for (Int_t i=0; i<40; i++) {
        if(RunNum==dRun15hPos[i]) fbFlagIsPosMagField = kTRUE;
    }
  }

    const Double_t * towZNCraw = aodZDC->GetZNCTowerEnergy();
    const Double_t * towZNAraw = aodZDC->GetZNATowerEnergy();
    
    // Get centroid from ZDCs *******************************************************
    
    Double_t Enucl = (RunNum < 209122 ? 1380. : 2511.);
    Double_t xyZNC[2]={0.,0.}, xyZNA[2]={0.,0.};
    Double_t towZNC[5]={0.}, towZNA[5]={0.};
    
    Double_t ZNCcalib=1., ZNAcalib=1.;
    
    // equalize gain of all towers
    if(RunNum!=fCachedRunNum) {
        for(Int_t i=0; i<5; i++) {
            fTowerGainEq[0][i] = (TH1D*)(fTowerEqList->FindObject(Form("fZNCTower[%d][%d]",RunNum,i)));
            fTowerGainEq[1][i] = (TH1D*)(fTowerEqList->FindObject(Form("fZNATower[%d][%d]",RunNum,i)));
        }
    }
    for(Int_t i=0; i<5; i++) {
        if(fTowerGainEq[0][i]) towZNC[i] = towZNCraw[i]*fTowerGainEq[0][i]->GetBinContent(fTowerGainEq[0][i]->FindBin(Centrality));
        if(fTowerGainEq[1][i]) towZNA[i] = towZNAraw[i]*fTowerGainEq[1][i]->GetBinContent(fTowerGainEq[1][i]->FindBin(Centrality));
    }
    
    if(RunNum>=245829) towZNA[2] = 0.;
    Double_t zncEnergy=0., znaEnergy=0.;
    for(Int_t i=0; i<5; i++){
        zncEnergy += towZNC[i];
        znaEnergy += towZNA[i];
    }
    if(RunNum>=245829) znaEnergy *= 8./7.;
    Double_t fZNCen = towZNC[0]/Enucl;
    Double_t fZNAen = towZNA[0]/Enucl;
    
    const Double_t x[4] = {-1.75, 1.75, -1.75, 1.75};
    const Double_t y[4] = {-1.75, -1.75, 1.75, 1.75};
    Double_t numXZNC=0., numYZNC=0., denZNC=0., cZNC, wZNC, EZNC, SumEZNC=0.;
    Double_t numXZNA=0., numYZNA=0., denZNA=0., cZNA, wZNA, EZNA, SumEZNA=0., BadChOr;
    Bool_t fAllChONZNC=kTRUE, fAllChONZNA=kTRUE;
    
    for(Int_t i=0; i<4; i++){
        // get energy
        EZNC = towZNC[i+1];
        SumEZNC += EZNC;
        
        // build centroid
        wZNC = TMath::Power(EZNC, fZDCGainAlpha);
        numXZNC += x[i]*wZNC;
        numYZNC += y[i]*wZNC;
        denZNC += wZNC;
        
        // get energy
        if(i==1) {
            EZNA = towZNA[0]-towZNA[1]-towZNA[3]-towZNA[4];
        } else {
            EZNA = towZNA[i+1];
        }
        SumEZNA += EZNA;
        
        // build centroid
        wZNA = TMath::Power(EZNA, fZDCGainAlpha);
        numXZNA += x[i]*wZNA;
        numYZNA += y[i]*wZNA;
        denZNA += wZNA;
    }
    if(denZNC>0.){
        Double_t nSpecnC = SumEZNC/Enucl;
        cZNC = 1.89358-0.71262/(nSpecnC+0.71789);
        xyZNC[0] = cZNC*numXZNC/denZNC;
        xyZNC[1] = cZNC*numYZNC/denZNC;
        denZNC *= cZNC;
    }
    else{
        xyZNC[0] = xyZNC[1] = 0.;
    }
    if(denZNA>0.){
        Double_t nSpecnA = SumEZNA/Enucl;
        cZNA = 1.89358-0.71262/(nSpecnA+0.71789);
        xyZNA[0] = cZNA*numXZNA/denZNA;
        xyZNA[1] = cZNA*numYZNA/denZNA;
        denZNA *= cZNA;
    }
    else{
        xyZNA[0] = xyZNA[1] = 0.;
    }
    
    fZDCFlowVect[0]->Set(xyZNC[0],xyZNC[1]);
    fZDCFlowVect[1]->Set(xyZNA[0],xyZNA[1]);
    fZDCFlowVect[0]->SetMult(denZNC);
    fZDCFlowVect[1]->SetMult(denZNA);
    
    // RE-CENTER ZDC Q-VECTORS ***************************************************
    
    Int_t qb[4] = {0};
    if(fbFlagIsPosMagField) { qb[0]=0; qb[1]=1; qb[2]=4; qb[3]=5; }
    else                    { qb[0]=2; qb[1]=3; qb[2]=6; qb[3]=7; }
    
    // get re-centered QM*
    //  Double_t QMCrec = denZNC;
    //  Double_t QMArec = denZNA;
    //  if(fAvEZDCCRbRPro && fAvEZDCARbRPro) {
    //    Int_t runbin = fAvEZDCCRbRPro->GetXaxis()->FindBin(Form("%d",RunNum));
    //    Int_t cenbin = fAvEZDCCRbRPro->GetYaxis()->FindBin(Centrality);
    //    QMCrec -= fAvEZDCCRbRPro->GetBinContent(runbin,cenbin);
    //    QMArec -= fAvEZDCARbRPro->GetBinContent(runbin,cenbin);
    //  }
    
    Bool_t IsGoodEvent = kTRUE;
    
    if(fZDCCalibList) {
        if(RunNum!=fCachedRunNum) {
            // get histos of run
            fZDCQHist[0] = (TProfile*)(fZDCCalibList->FindObject(Form("Run %d",RunNum))->FindObject(Form("fCRCZDCQVecC[%d][%d]",RunNum,0)));
            fZDCQHist[1] = (TProfile*)(fZDCCalibList->FindObject(Form("Run %d",RunNum))->FindObject(Form("fCRCZDCQVecC[%d][%d]",RunNum,1)));
            fZDCQHist[2] = (TProfile*)(fZDCCalibList->FindObject(Form("Run %d",RunNum))->FindObject(Form("fCRCZDCQVecA[%d][%d]",RunNum,0)));
            fZDCQHist[3] = (TProfile*)(fZDCCalibList->FindObject(Form("Run %d",RunNum))->FindObject(Form("fCRCZDCQVecA[%d][%d]",RunNum,1)));
            
            for(Int_t k=0; k<4; k++) {
                fZDCVtxHist[k] = (TProfile3D*)(fZDCCalibList->FindObject(Form("Run %d",RunNum))->FindObject(Form("fCRCZDCQVecVtxPos540[%d][%d]",RunNum,k)));
                
                fZDCEcomTotHist[k] = (TProfile2D*)(fZDCCalibList->FindObject(Form("fCRCZDCQVecEComTot[%d]",k)));
            }
            
            for(Int_t c=0; c<10; c++) {
                for(Int_t k=0; k<4; k++) {
                    fZDCVtxCenHist[c][k] = (TProfile3D*)(fZDCCalibList->FindObject(Form("fCRCZDCQVecVtxPosCen[%d][%d]",c,k)));
                }
            }
            
            for(Int_t k=0; k<4; k++) {
                fZDCVtxFitHist[k] = (TH3D*)fZDCCalibList->FindObject(Form("TH3SlopeRunCenVtx[%d]",k));
                if(fZDCVtxFitHist[k]) {
                    fZDCVtxFitHist[k]->Sumw2(kFALSE);
                    Int_t runbin = fZDCVtxFitHist[k]->GetXaxis()->FindBin(Form("%d",RunNum));
                    fZDCVtxFitHist[k]->GetXaxis()->SetRange(runbin,runbin);
                    for(Int_t i=0; i<3; i++) {
                        fZDCVtxFitHist[k]->GetZaxis()->SetRange(i+1,i+1);
                        fZDCVtxFitCenProjHist[k][i] = (TH1D*)fZDCVtxFitHist[k]->Project3D("y")->Clone(Form("proj[%d][%d]",k,i));
                    }
                }
            }
            
            for(Int_t c=0; c<10; c++) {
                for(Int_t k=0; k<8; k++) {
                    fZDCVtxCenHistMagPol[c][k] = (TProfile3D*)(fZDCCalibList->FindObject(Form("fZDCVtxCenHistMagPol[%d][%d]",c,k)));
                }
            }
            
            //      for(Int_t i=0; i<10; i++) {
            //        for(Int_t z=0; z<10; z++) {
            //          for(Int_t k=0; k<4; k++) {
            //            fZDCQVecVtxCenEZDC3D[i][z][k] = (TH3D*)(fZDCCalibList->FindObject(Form("ZDCQVecVtxCenEZDC3D[%d][%d][%d]",i,z,qb[k])));
            //          }
            //        }
            //      }
            
        }
        
        // ZDCN-C
        Double_t QCRe = fZDCFlowVect[0]->X();
        Double_t QCIm = fZDCFlowVect[0]->Y();
        Double_t QMC  = fZDCFlowVect[0]->GetMult();
        // ZDCN-A
        Double_t QARe = fZDCFlowVect[1]->X();
        Double_t QAIm = fZDCFlowVect[1]->Y();
        Double_t QMA  = fZDCFlowVect[1]->GetMult();
        
        Double_t QCReR=QCRe, QCImR=QCIm, QAReR=QARe, QAImR=QAIm;
        
        // STEP #1: re-center vs centrality (1%) vs run number
        
        if (fZDCQHist[0]) {
            Double_t AvQCRe = fZDCQHist[0]->GetBinContent(fZDCQHist[0]->FindBin(Centrality));
            Double_t SDQCRe = fZDCQHist[0]->GetBinError(fZDCQHist[0]->FindBin(Centrality));
            Double_t AvQCIm = fZDCQHist[1]->GetBinContent(fZDCQHist[1]->FindBin(Centrality));
            Double_t SDQCIm = fZDCQHist[1]->GetBinError(fZDCQHist[1]->FindBin(Centrality));
            
            Double_t AvQARe = fZDCQHist[2]->GetBinContent(fZDCQHist[2]->FindBin(Centrality));
            Double_t SDQARe = fZDCQHist[2]->GetBinError(fZDCQHist[2]->FindBin(Centrality));
            Double_t AvQAIm = fZDCQHist[3]->GetBinContent(fZDCQHist[3]->FindBin(Centrality));
            Double_t SDQAIm = fZDCQHist[3]->GetBinError(fZDCQHist[3]->FindBin(Centrality));
            
            if(AvQCRe && AvQCIm && QMC>0. && sqrt(QCRe*QCRe+QCIm*QCIm)>1.E-6) {
                QCReR = QCRe-AvQCRe;
                QCImR = QCIm-AvQCIm;
                fZDCFlowVect[0]->Set(QCReR,QCImR);
            }
            
            if(AvQARe && AvQAIm && QMA>0. && sqrt(QARe*QARe+QAIm*QAIm)>1.E-6) {
                QAReR = QARe-AvQARe;
                QAImR = QAIm-AvQAIm;
                fZDCFlowVect[1]->Set(QAReR,QAImR);
            }
        }
        
        // STEP #2: re-center vs primary vtx vs centrality (10%)
        
        if (fZDCVtxCenHist[fCenBin][0]) {
            Bool_t pass = kTRUE;
            if(fVtxPosCor[0] < fZDCVtxCenHist[fCenBin][0]->GetXaxis()->GetXmin() || fVtxPosCor[0] > fZDCVtxCenHist[fCenBin][0]->GetXaxis()->GetXmax()) pass = kFALSE;
            if(fVtxPosCor[1] < fZDCVtxCenHist[fCenBin][0]->GetYaxis()->GetXmin() || fVtxPosCor[1] > fZDCVtxCenHist[fCenBin][0]->GetYaxis()->GetXmax()) pass = kFALSE;
            if(fVtxPosCor[2] < fZDCVtxCenHist[fCenBin][0]->GetZaxis()->GetXmin() || fVtxPosCor[2] > fZDCVtxCenHist[fCenBin][0]->GetZaxis()->GetXmax()) pass = kFALSE;
            if(!pass) {
                IsGoodEvent = kFALSE;
            } else {
                QCReR -= fZDCVtxCenHist[fCenBin][0]->GetBinContent(fZDCVtxCenHist[fCenBin][0]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
                QCImR -= fZDCVtxCenHist[fCenBin][1]->GetBinContent(fZDCVtxCenHist[fCenBin][1]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
                QAReR -= fZDCVtxCenHist[fCenBin][2]->GetBinContent(fZDCVtxCenHist[fCenBin][2]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
                QAImR -= fZDCVtxCenHist[fCenBin][3]->GetBinContent(fZDCVtxCenHist[fCenBin][3]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
                fZDCFlowVect[0]->Set(QCReR,QCImR);
                fZDCFlowVect[1]->Set(QAReR,QAImR);
            }
        }
        
        // STEP #3: re-center vs primary vtx vs run number
        
        if (fZDCVtxHist[0]) {
            QCReR -= fZDCVtxHist[0]->GetBinContent(fZDCVtxHist[0]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
            QCImR -= fZDCVtxHist[1]->GetBinContent(fZDCVtxHist[1]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
            fZDCFlowVect[0]->Set(QCReR,QCImR);
            QAReR -= fZDCVtxHist[2]->GetBinContent(fZDCVtxHist[2]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
            QAImR -= fZDCVtxHist[3]->GetBinContent(fZDCVtxHist[3]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
            fZDCFlowVect[1]->Set(QAReR,QAImR);
        }
        
        // STEP #4: re-center vs centrality vs energy in the common tower
        
        if(fZDCEcomTotHist[0]) {
            QCReR -= fZDCEcomTotHist[0]->GetBinContent(fZDCEcomTotHist[0]->FindBin(Centrality,fZNCen));
            QCImR -= fZDCEcomTotHist[1]->GetBinContent(fZDCEcomTotHist[1]->FindBin(Centrality,fZNCen));
            fZDCFlowVect[0]->Set(QCReR,QCImR);
            QAReR -= fZDCEcomTotHist[2]->GetBinContent(fZDCEcomTotHist[2]->FindBin(Centrality,fZNAen));
            QAImR -= fZDCEcomTotHist[3]->GetBinContent(fZDCEcomTotHist[3]->FindBin(Centrality,fZNAen));
            fZDCFlowVect[1]->Set(QAReR,QAImR);
        }
        
        // STEP #5: re-center vs vtx vs cen vs run number (through fits)
        
        if(fZDCVtxFitHist[0]) {
            for (Int_t i=0; i<3; i++) {
                QCReR -= fVtxPosCor[i]*fZDCVtxFitCenProjHist[0][i]->Interpolate(Centrality);
                QCImR -= fVtxPosCor[i]*fZDCVtxFitCenProjHist[1][i]->Interpolate(Centrality);
                QAReR -= fVtxPosCor[i]*fZDCVtxFitCenProjHist[2][i]->Interpolate(Centrality);
                QAImR -= fVtxPosCor[i]*fZDCVtxFitCenProjHist[3][i]->Interpolate(Centrality);
            }
            fZDCFlowVect[0]->Set(QCReR,QCImR);
            fZDCFlowVect[1]->Set(QAReR,QAImR);
        }
        
        // second iteration (2D)
        
        if (fZDCVtxCenHistMagPol[fCenBin][0]) {
            if(fbFlagIsPosMagField) {
                QCReR -= fZDCVtxCenHistMagPol[fCenBin][0]->GetBinContent(fZDCVtxCenHistMagPol[fCenBin][0]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
                QCImR -= fZDCVtxCenHistMagPol[fCenBin][1]->GetBinContent(fZDCVtxCenHistMagPol[fCenBin][1]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
                QAReR -= fZDCVtxCenHistMagPol[fCenBin][4]->GetBinContent(fZDCVtxCenHistMagPol[fCenBin][4]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
                QAImR -= fZDCVtxCenHistMagPol[fCenBin][5]->GetBinContent(fZDCVtxCenHistMagPol[fCenBin][5]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
            } else {
                QCReR -= fZDCVtxCenHistMagPol[fCenBin][2]->GetBinContent(fZDCVtxCenHistMagPol[fCenBin][2]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
                QCImR -= fZDCVtxCenHistMagPol[fCenBin][3]->GetBinContent(fZDCVtxCenHistMagPol[fCenBin][3]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
                QAReR -= fZDCVtxCenHistMagPol[fCenBin][6]->GetBinContent(fZDCVtxCenHistMagPol[fCenBin][6]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
                QAImR -= fZDCVtxCenHistMagPol[fCenBin][7]->GetBinContent(fZDCVtxCenHistMagPol[fCenBin][7]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
            }
            fZDCFlowVect[0]->Set(QCReR,QCImR);
            fZDCFlowVect[1]->Set(QAReR,QAImR);
        }
        
        // STEP #6: re-center vs centrality vs total energy vs vtx
        
        //    Int_t EZDCCBin = fCRCZDCQVecDummyEZDCBins[fCenBin]->GetXaxis()->FindBin(QMCrec)-1;
        //    Int_t EZDCABin = fCRCZDCQVecDummyEZDCBins[fCenBin]->GetXaxis()->FindBin(QMArec)-1;
        //
        //    if(fZDCQVecVtxCenEZDC3D[0][0][0]) {
        //      printf("doing step 6 \n");
        //      Bool_t pass2=kTRUE;
        //      // exclude events with vtx outside of range
        //      if(fVtxPosCor[0] < fZDCQVecVtxCenEZDC3D[0][0][0]->GetXaxis()->GetXmin() || fVtxPosCor[0] > fZDCQVecVtxCenEZDC3D[0][0][0]->GetXaxis()->GetXmax()) pass2 = kFALSE;
        //      if(fVtxPosCor[1] < fZDCQVecVtxCenEZDC3D[0][0][0]->GetYaxis()->GetXmin() || fVtxPosCor[1] > fZDCQVecVtxCenEZDC3D[0][0][0]->GetYaxis()->GetXmax()) pass2 = kFALSE;
        //      if(fVtxPosCor[2] < fZDCQVecVtxCenEZDC3D[0][0][0]->GetZaxis()->GetXmin() || fVtxPosCor[2] > fZDCQVecVtxCenEZDC3D[0][0][0]->GetZaxis()->GetXmax()) pass2 = kFALSE;
        //      // exclude events with very low or very high total energy
        //      if(fabs(QMCrec)>100. || fabs(QMArec)>100.) pass2 = kFALSE;
        //      // get EZDC bin
        //      Int_t EZDCCBin = fCRCZDCQVecDummyEZDCBins[fCenBin]->GetXaxis()->FindBin(QMCrec)-1;
        //      Int_t EZDCABin = fCRCZDCQVecDummyEZDCBins[fCenBin]->GetXaxis()->FindBin(QMArec)-1;
        //      if(EZDCCBin<0) EZDCCBin=0;
        //      if(EZDCCBin>9) EZDCCBin=9;
        //      if(EZDCABin<0) EZDCABin=0;
        //      if(EZDCABin>9) EZDCABin=9;
        //      if(pass2) {
        //        // check if possible to interpolate
        //        Bool_t bInterp = kTRUE;
        //        Int_t bx = fZDCQVecVtxCenEZDC3D[0][0][0]->GetXaxis()->FindBin(fVtxPosCor[0]);
        //        Int_t by = fZDCQVecVtxCenEZDC3D[0][0][0]->GetYaxis()->FindBin(fVtxPosCor[1]);
        //        Int_t bz = fZDCQVecVtxCenEZDC3D[0][0][0]->GetZaxis()->FindBin(fVtxPosCor[2]);
        //        if(bx==1 || bx==fZDCQVecVtxCenEZDC3D[0][0][0]->GetXaxis()->GetNbins()) bInterp = kFALSE;
        //        if(by==1 || by==fZDCQVecVtxCenEZDC3D[0][0][0]->GetYaxis()->GetNbins()) bInterp = kFALSE;
        //        if(bz==1 || bz==fZDCQVecVtxCenEZDC3D[0][0][0]->GetZaxis()->GetNbins()) bInterp = kFALSE;
        //        if(bInterp) {
        //          QCReR -= fZDCQVecVtxCenEZDC3D[fCenBin][EZDCCBin][0]->Interpolate(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]);
        //          QCImR -= fZDCQVecVtxCenEZDC3D[fCenBin][EZDCCBin][1]->Interpolate(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]);
        //          QAReR -= fZDCQVecVtxCenEZDC3D[fCenBin][EZDCABin][2]->Interpolate(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]);
        //          QAImR -= fZDCQVecVtxCenEZDC3D[fCenBin][EZDCABin][3]->Interpolate(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]);
        //        } else {
        //          QCReR -= fZDCQVecVtxCenEZDC3D[fCenBin][EZDCCBin][0]->GetBinContent(bx,by,bz);
        //          QCImR -= fZDCQVecVtxCenEZDC3D[fCenBin][EZDCCBin][1]->GetBinContent(bx,by,bz);
        //          QAReR -= fZDCQVecVtxCenEZDC3D[fCenBin][EZDCABin][2]->GetBinContent(bx,by,bz);
        //          QAImR -= fZDCQVecVtxCenEZDC3D[fCenBin][EZDCABin][3]->GetBinContent(bx,by,bz);
        //        }
        //      } else {
        //        IsGoodEvent = kFALSE;
        //      }
        //    }
        
    } else {
        printf("WARNING: no list provided for ZDC Q-vector re-centering ! \n");
    }
    
    Double_t xyZNCfinal[2]={fZDCFlowVect[0]->X(),fZDCFlowVect[0]->Y()};
    Double_t xyZNAfinal[2]={-fZDCFlowVect[1]->X(),fZDCFlowVect[1]->Y()}; // this is not a bug: QAReR --> -QAReR
    if(!IsGoodEvent) {
        xyZNCfinal[0]=0.; xyZNCfinal[1]=0.;
        xyZNAfinal[0]=0.; xyZNAfinal[1]=0.;
    }
    fFlowEvent->SetZDC2Qsub(xyZNCfinal,denZNC,xyZNAfinal,denZNA);
    
    // save run number
    fCachedRunNum = RunNum;
}





