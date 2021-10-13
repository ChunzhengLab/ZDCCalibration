#include <iostream>
#include <cstdlib>
#include <sys/time.h>
#include <vector>
// ROOT classes
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TList.h"
#include "TH1I.h"
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
#include "TLorentzVector.h"
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
// Alice classes
#include "AliCentrality.h"
#include "AliEventplane.h"
#include "AliEventCuts.h"
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
// #include "AliMultSelection.h"

#include "AliAnalysisTaskEE.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskEE);

//---------------------------------------------------

AliAnalysisTaskEE::AliAnalysisTaskEE() : 
  AliAnalysisTaskSE()
{
  fileNameEPMean_="alien:///alice/cern.ch/user/q/qshou/EPCalib/LHC11h_pass2/recentering/Merged_AnalysisResults.root";
  fileEPMean_=NULL;
  listEPMean_=NULL;
  fileNameEPShifting_="alien:///alice/cern.ch/user/q/qshou/EPCalib/LHC11h_pass2/shifting/Merged_AnalysisResults.root";
  fileEPShifting_=NULL;
  listEPShifting_=NULL;
  p_XAMean_=NULL;
  p_XCMean_=NULL;
  p_YAMean_=NULL;
  p_YCMean_=NULL;
  XAMean_=0;
  XCMean_=0;
  YAMean_=0;
  YCMean_=0;
  p_shiftCoeffSinA_=NULL;
  p_shiftCoeffCosA_=NULL;
  p_shiftCoeffSinC_=NULL;
  p_shiftCoeffCosC_=NULL;
  p_shiftCoeffSinFull_=NULL;
  p_shiftCoeffCosFull_=NULL;

  fOutputList=NULL;
  h_evtCount=NULL;
  h_cent=NULL;
  h_mult=NULL;
  for (int i=0; i<4; ++i) {
    h_ZNA[i]=NULL;
    h_ZNC[i]=NULL;
  }
  h_XYA=NULL;
  h_XYC=NULL;
  for (int i=0; i<3; ++i) {
    for (int j=0; j<10; ++j) {
      h_psiZDC_A[i][j]=NULL;
      h_psiZDC_C[i][j]=NULL;
      h_psiZDC_full[i][j]=NULL;
    }
  }
  p_psiZDC_reso=NULL;
  h_pt=NULL;
  h_eta[0]=NULL;
  h_eta[1]=NULL;
  h_phi=NULL;
  h_dcaXy[0]=NULL;
  h_dcaXy[1]=NULL;
  h_dcaZ[0]=NULL;
  h_dcaZ[1]=NULL;
  for (int i=0; i<4; ++i) {
    h_nHits[i]=NULL;
  }
  for (int i=0; i<6; ++i) {
    h_nSigmaE_p[i]=NULL;
  }
  h_mass_epem=NULL;
}

//---------------------------------------------------

AliAnalysisTaskEE::AliAnalysisTaskEE(const char *name) : 
  AliAnalysisTaskSE(name)
{
  fileNameEPMean_="alien:///alice/cern.ch/user/q/qshou/EPCalib/LHC11h_pass2/recentering/Merged_AnalysisResults.root";
  fileEPMean_=NULL;
  listEPMean_=NULL;
  fileNameEPShifting_="alien:///alice/cern.ch/user/q/qshou/EPCalib/LHC11h_pass2/shifting/Merged_AnalysisResults.root";
  fileEPShifting_=NULL;
  listEPShifting_=NULL;
  p_XAMean_=NULL;
  p_XCMean_=NULL;
  p_YAMean_=NULL;
  p_YCMean_=NULL;
  XAMean_=0;
  XCMean_=0;
  YAMean_=0;
  YCMean_=0;
  p_shiftCoeffSinA_=NULL;
  p_shiftCoeffCosA_=NULL;
  p_shiftCoeffSinC_=NULL;
  p_shiftCoeffCosC_=NULL;
  p_shiftCoeffSinFull_=NULL;
  p_shiftCoeffCosFull_=NULL;

  fOutputList=NULL;
  h_evtCount=NULL;
  h_cent=NULL;
  h_mult=NULL;
  for (int i=0; i<4; ++i) {
    h_ZNA[i]=NULL;
    h_ZNC[i]=NULL;
  }
  h_XYA=NULL;
  h_XYC=NULL;
  for (int i=0; i<3; ++i) {
    for (int j=0; j<10; ++j) {
      h_psiZDC_A[i][j]=NULL;
      h_psiZDC_C[i][j]=NULL;
      h_psiZDC_full[i][j]=NULL;
    }
  }
  p_psiZDC_reso=NULL;
  h_pt=NULL;
  h_eta[0]=NULL;
  h_eta[1]=NULL;
  h_phi=NULL;
  h_dcaXy[0]=NULL;
  h_dcaXy[1]=NULL;
  h_dcaZ[0]=NULL;
  h_dcaZ[1]=NULL;
  for (int i=0; i<4; ++i) {
    h_nHits[i]=NULL;
  }
  for (int i=0; i<6; ++i) {
    h_nSigmaE_p[i]=NULL;
  }
  h_mass_epem=NULL;

  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
}

//---------------------------------------------------

AliAnalysisTaskEE::~AliAnalysisTaskEE()
{}

//---------------------------------------------------

void AliAnalysisTaskEE::Terminate(Option_t *) 
{}

//---------------------------------------------------

void AliAnalysisTaskEE::UserCreateOutputObjects()
{
  fOutputList = new TList();
  fOutputList->SetName(GetName());
  fOutputList->SetOwner(kTRUE);

  fileEPMean_ = TFile::Open(fileNameEPMean_.Data(),"READ");
  if (!fileEPMean_) {
    AliError(Form("%s: Could not retrieve EP-Recentering file.", GetName()));
    return;
  }
  listEPMean_ = (TList*)fileEPMean_->Get("output");
  if (!listEPMean_) {
    AliError(Form("%s: Could not retrieve EP-Recentering list.", GetName()));
    return;
  }
  fileEPShifting_ = TFile::Open(fileNameEPShifting_.Data(),"READ");
  if (!fileEPShifting_) {
    AliError(Form("%s: Could not retrieve EP-Shifting file.", GetName()));
    return;
  }
  listEPShifting_ = (TList*)fileEPShifting_->Get("output");
  if (!listEPShifting_) {
    AliError(Form("%s: Could not retrieve EP-Shifting list.", GetName()));
    return;
  }

  // event-wise
  h_evtCount = new TH1I("evtCount", "", 20, 1, 21);
  h_evtCount->GetXaxis()->SetBinLabel(1,"All");
  h_evtCount->GetXaxis()->SetBinLabel(2,"Info.");
  h_evtCount->GetXaxis()->SetBinLabel(3,"Evt");
  h_evtCount->GetXaxis()->SetBinLabel(10,"Manager");
  h_evtCount->GetXaxis()->SetBinLabel(11,"Handler");
  h_evtCount->GetXaxis()->SetBinLabel(12,"AOD");
  h_evtCount->GetXaxis()->SetBinLabel(13,"PID");
  h_evtCount->GetXaxis()->SetBinLabel(14,"Utils");
  // h_evtCount->GetXaxis()->SetBinLabel(15,"MultSel");
  fOutputList->Add(h_evtCount);

  h_cent = new TH1D("centrality", "", 100, 0, 100);
  fOutputList->Add(h_cent);

  h_mult = new TH1D("mult", "", 3000, 0, 3000);
  fOutputList->Add(h_mult);

  h_ZNA[0] = new TH1D("zna_1", "", 40000, 0, 40000);
  h_ZNA[1] = new TH1D("zna_2", "", 40000, 0, 40000);
  h_ZNA[2] = new TH1D("zna_3", "", 40000, 0, 40000);
  h_ZNA[3] = new TH1D("zna_4", "", 40000, 0, 40000);
  h_ZNC[0] = new TH1D("znc_1", "", 40000, 0, 40000);
  h_ZNC[1] = new TH1D("znc_2", "", 40000, 0, 40000);
  h_ZNC[2] = new TH1D("znc_3", "", 40000, 0, 40000);
  h_ZNC[3] = new TH1D("znc_4", "", 40000, 0, 40000);
  for (int i=0; i<4; ++i) {
    fOutputList->Add(h_ZNA[i]); 
    fOutputList->Add(h_ZNC[i]);
  }
  
  for (int i=0; i<10; ++i) {
    h_psiZDC_A[0][i] = new TH1D(Form("psiZDC_A_raw_cent%i",i), "", 100, -2*TMath::Pi(), 2*TMath::Pi());
    h_psiZDC_C[0][i] = new TH1D(Form("psiZDC_C_raw_cent%i",i), "", 100, -2*TMath::Pi(), 2*TMath::Pi());
    h_psiZDC_full[0][i] = new TH1D(Form("psiZDC_full_raw_cent%i",i), "", 100, -2*TMath::Pi(), 2*TMath::Pi());
    h_psiZDC_A[1][i] = new TH1D(Form("psiZDC_A_afRecent_cent%i",i), "", 100, -2*TMath::Pi(), 2*TMath::Pi());
    h_psiZDC_C[1][i] = new TH1D(Form("psiZDC_C_afRecent_cent%i",i), "", 100, -2*TMath::Pi(), 2*TMath::Pi());
    h_psiZDC_full[1][i] = new TH1D(Form("psiZDC_full_afRecent_cent%i",i), "", 100, -2*TMath::Pi(), 2*TMath::Pi());
    h_psiZDC_A[2][i] = new TH1D(Form("psiZDC_A_afShift_cent%i",i), "", 100, -2*TMath::Pi(), 2*TMath::Pi());
    h_psiZDC_C[2][i] = new TH1D(Form("psiZDC_C_afShift_cent%i",i), "", 100, -2*TMath::Pi(), 2*TMath::Pi());
    h_psiZDC_full[2][i] = new TH1D(Form("psiZDC_full_afShift_cent%i",i), "", 100, -2*TMath::Pi(), 2*TMath::Pi());
    fOutputList->Add(h_psiZDC_A[0][i]);
    fOutputList->Add(h_psiZDC_C[0][i]);
    fOutputList->Add(h_psiZDC_full[0][i]);
    fOutputList->Add(h_psiZDC_A[1][i]);
    fOutputList->Add(h_psiZDC_C[1][i]);
    fOutputList->Add(h_psiZDC_full[1][i]);
    fOutputList->Add(h_psiZDC_A[2][i]);
    fOutputList->Add(h_psiZDC_C[2][i]);
    fOutputList->Add(h_psiZDC_full[2][i]);
  }

  p_psiZDC_reso = new TProfile("psiZDC_reso","",10,0,10);
  fOutputList->Add(p_psiZDC_reso);

  h_XYA = new TH2D("XYA","",200,-1,1,200,-1,1);
  h_XYC = new TH2D("XYC","",200,-1,1,200,-1,1);
  fOutputList->Add(h_XYA);
  fOutputList->Add(h_XYC);

  // track-wise
  h_pt = new TH1D("pt", "", 200, 0, 20);
  fOutputList->Add(h_pt);

  h_eta[0] = new TH1D("eta_bfCut", "", 400, -10, 10);
  h_eta[1] = new TH1D("eta_afCut",  "", 400, -10, 10);
  fOutputList->Add(h_eta[0]);
  fOutputList->Add(h_eta[1]);

  h_phi = new TH1D("phi", "", 1000, 0, 2*TMath::Pi());
  fOutputList->Add(h_phi);

  h_dcaXy[0] = new TH1D("dcaXy_bfCut", "", 200, 0, 10);
  h_dcaXy[1] = new TH1D("dcaXy_afCut",  "", 200, 0, 10);
  fOutputList->Add(h_dcaXy[0]);
  fOutputList->Add(h_dcaXy[1]);

  h_dcaZ[0] = new TH1D("dcaZ_bfCut", "", 200, 0, 10);
  h_dcaZ[1] = new TH1D("dcaZ_afCut",  "", 200, 0, 10);
  fOutputList->Add(h_dcaZ[0]);
  fOutputList->Add(h_dcaZ[1]);

  h_nHits[0] = new TH1D("nHitsTPC_bfCut", "", 200, 0, 200);
  h_nHits[1] = new TH1D("nHitsTPC_afCut",  "", 200, 0, 200);
  h_nHits[2] = new TH1D("nHitsITS_bfCut", "", 200, 0, 200);
  h_nHits[3] = new TH1D("nHitsITS_afCut",  "", 200, 0, 200);
  for (int i=0; i<4; ++i) {    
    fOutputList->Add(h_nHits[i]);
  }

  h_nSigmaE_p[0] = new TH2D("nSigmaETPC_p_bfCut", "", 500, -10, 10, 800, -40, 40);
  h_nSigmaE_p[1] = new TH2D("nSigmaETPC_p_afCut", "", 500, -10, 10, 800, -40, 40);
  h_nSigmaE_p[2] = new TH2D("nSigmaEITS_p_bfCut", "", 500, -10, 10, 800, -40, 40);
  h_nSigmaE_p[3] = new TH2D("nSigmaEITS_p_afCut", "", 500, -10, 10, 800, -40, 40);
  h_nSigmaE_p[4] = new TH2D("nSigmaETOF_p_bfCut", "", 500, -10, 10, 800, -40, 40);
  h_nSigmaE_p[5] = new TH2D("nSigmaETOF_p_afCut", "", 500, -10, 10, 800, -40, 40);
  for (int i=0; i<6; ++i) {    
    fOutputList->Add(h_nSigmaE_p[i]);
  }

  h_mass_epem = new TH1D("mass_epem", "", 200, 0, 4);
  fOutputList->Add(h_mass_epem);

  PostData(1,fOutputList);
}

//---------------------------------------------------

void AliAnalysisTaskEE::UserExec(Option_t *)
{
  h_evtCount->Fill(1);

  AliAnalysisManager* manager = AliAnalysisManager::GetAnalysisManager();
  if (!manager) {
    AliError(Form("%s: Could not get Analysis Manager", GetName()));
  } else h_evtCount->Fill(10);
 
  AliAODInputHandler* handler = (AliAODInputHandler*)manager->GetInputEventHandler();
  if (!handler) {
    AliError(Form("%s: Could not get Input Handler", GetName()));
  } else h_evtCount->Fill(11);
 
  AliAODEvent* fAOD = (AliAODEvent*)InputEvent();
  if (!fAOD) {
    AliError(Form("%s: Could not get AOD event", GetName()));
  } else h_evtCount->Fill(12);
 
  AliPIDResponse* fPID = handler->GetPIDResponse();
  if (!fPID) {
    AliError(Form("%s: Could not get PIDResponse", GetName()));
  } else h_evtCount->Fill(13);
 
  // AliAnalysisUtils* fUtils = new AliAnalysisUtils();
  // if (!fUtils) {
  //   AliError(Form("%s: Could not get AliAnalysisUtils", GetName()));
  // } else h_evtCount->Fill(14);
 
  // AliMultSelection* fMultSel = (AliMultSelection*)InputEvent()->FindListObject("MultSelection");
  // if (!fMultSel) {
  //   AliError(Form("%s: Could not get AliMultSelection", GetName()));
  // } else h_evtCount->Fill(15);
 
  // if (!manager || !handler || !fAOD || !fPID || !fUtils || !fMultSel) return;
  if (!manager || !handler || !fAOD || !fPID) return;
  h_evtCount->Fill(2);

  //------------------
  // event
  //------------------

  int runNum = fAOD->GetRunNumber();
  int runNumBin = GetRunNumBin(runNum);
  if (runNumBin<0) return;

  // vertex
  AliAODVertex* fVtx = fAOD->GetPrimaryVertex();
  double vx    = fVtx->GetX();
  double vy    = fVtx->GetY();
  double vz    = fVtx->GetZ();
  double vzSPD = fAOD->GetPrimaryVertexSPD()->GetZ();
  if (fabs(vz)>10) return;

  h_evtCount->Fill(3);

  // pileup
  // fUtils->SetUseOutOfBunchPileUp(true);
  // bool isPileup = fUtils->IsPileUpEvent(fAOD);

  // bool isPileup = fUtils->IsPileUpMV(fAOD);
  // if (isPileup) return;

  // centrality
  double cent = -999;
  // run1
  cent = fAOD->GetCentrality()->GetCentralityPercentile("V0M");
  // cent = event->GetCentrality()->GetCentralityPercentile("CL1");
  // cent = event->GetCentrality()->GetCentralityPercentile("TRK");
  // run2
  // if (collision=="pPb") {
  //   double centV0A = fAOD->GetCentrality()->GetCentralityPercentile("V0A");
  //   double centV0M = fAOD->GetCentrality()->GetCentralityPercentile("V0M");
  //   double centTRK = fAOD->GetCentrality()->GetCentralityPercentile("TRK");
  //   double centSPD = fAOD->GetCentrality()->GetCentralityPercentile("CL1");
  //   double centZNA = fAOD->GetCentrality()->GetCentralityPercentile("ZNA");
  //   double centZNC = fAOD->GetCentrality()->GetCentralityPercentile("ZNC");
  //   double centCND = fAOD->GetCentrality()->GetCentralityPercentile("CND");
  //   cent = centV0A;
  // } else if (collision=="PbPb") {
  //   // https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PACentStudiesRun2
  //   // double lPercentile = fMultSel->GetMultiplicityPercentile("V0M");
  //   // int lEvSelCode    = fMultSel->GetEvSelCode();
  //   // if (lEvSelCode>1) return;
  //   // cent = lPercentile;

  //   // https://twiki.cern.ch/twiki/bin/viewauth/ALICE/AliDPGtoolsEventProp
  //   cent = fEventCuts.GetCentrality();
  // }
  if (cent<0) return;
  int centBin=(int)cent/10;
  h_cent->Fill(cent);

  //------------------
  // event plane
  //------------------

  // fEP = fAOD->GetEventplane();
  // double fEPV0=-999, fEPV0A=-999, fEPV0C=-999;
  // if (fEP) {
  //   fEPV0  = TVector2::Phi_0_2pi(fEP->GetEventplane("V0",  fAOD)); if (fEPV0>TMath::Pi()) fEPV0-=TMath::Pi();
  //   fEPV0A = TVector2::Phi_0_2pi(fEP->GetEventplane("V0A", fAOD)); if (fEPV0A>TMath::Pi()) fEPV0A-=TMath::Pi();
  //   fEPV0C = TVector2::Phi_0_2pi(fEP->GetEventplane("V0C", fAOD)); if (fEPV0C>TMath::Pi()) fEPV0C-=TMath::Pi();
  // } 
  //cout<<fEPV0<<" "<<fEPV0A<<" "<<fEPV0C<<endl;

  AliAODZDC* aodZDC = fAOD->GetZDCData();
  if (!aodZDC) return;
  const double* energyZNA  = aodZDC->GetZNATowerEnergy();
  const double* energyZNC  = aodZDC->GetZNCTowerEnergy();
  if (energyZNA[1]<0 ||
      energyZNA[2]<0 ||
      energyZNA[3]<0 ||
      energyZNA[4]<0 ||
      energyZNC[1]<0 ||
      energyZNC[2]<0 ||
      energyZNC[3]<0 ||
      energyZNC[4]<0) return;
  const double x[4] = {-1.75, 1.75, -1.75, 1.75};
  const double y[4] = {-1.75, -1.75, 1.75, 1.75};
  double xEA=0, xEC=0, yEA=0, yEC=0, EA=0, EC=0;
  for (int i=0; i<4; ++i) {
    h_ZNA[i]->Fill(energyZNA[i+1]);
    h_ZNC[i]->Fill(energyZNC[i+1]);

    xEA += x[i]*energyZNA[i+1];
    yEA += y[i]*energyZNA[i+1];
    xEC += x[i]*energyZNC[i+1];
    yEC += y[i]*energyZNC[i+1];

    EA += energyZNA[i+1];
    EC += energyZNC[i+1];
  }
  double XA = xEA/EA;
  double XC = xEC/EC;
  double YA = yEA/EA;
  double YC = yEC/EC;
  double psiA = atan2(YA,XA);
  double psiC = atan2(YC,XC);
  double psiFull = atan2(YA+YC,XA+XC);
  h_psiZDC_A[0][centBin]->Fill(psiA);
  h_psiZDC_C[0][centBin]->Fill(psiC);
  h_psiZDC_full[0][centBin]->Fill(psiFull);

  // recentering
  double XA_mean=0,XC_mean=0,YA_mean=0,YC_mean=0;
  GetEPMean(runNum,centBin,vx,vy);
  XA -= XAMean_;
  XC -= XCMean_;
  YA -= YAMean_;
  YC -= YCMean_;
  h_XYA->Fill(XA,YA);
  h_XYC->Fill(XC,YC);
  psiA = atan2(YA,XA);
  psiC = atan2(YC,XC);
  psiFull = atan2(YA+YC,XA+XC);
  h_psiZDC_A[1][centBin]->Fill(psiA);
  h_psiZDC_C[1][centBin]->Fill(psiC);
  h_psiZDC_full[1][centBin]->Fill(psiFull);

  // shifting
  p_shiftCoeffSinA_ = (TProfile2D*)listEPShifting_->FindObject(Form("shiftCoeffSinA_%i", runNum));
  p_shiftCoeffCosA_ = (TProfile2D*)listEPShifting_->FindObject(Form("shiftCoeffCosA_%i", runNum));
  p_shiftCoeffSinC_ = (TProfile2D*)listEPShifting_->FindObject(Form("shiftCoeffSinC_%i", runNum));
  p_shiftCoeffCosC_ = (TProfile2D*)listEPShifting_->FindObject(Form("shiftCoeffCosC_%i", runNum));
  p_shiftCoeffSinFull_ = (TProfile2D*)listEPShifting_->FindObject(Form("shiftCoeffSinFull_%i", runNum));
  p_shiftCoeffCosFull_ = (TProfile2D*)listEPShifting_->FindObject(Form("shiftCoeffCosFull_%i", runNum));

  for (int i = 1; i <= 20; ++i) {
    double shiftSinA = p_shiftCoeffSinA_->GetBinContent(centBin+1,i);
    double shiftCosA = p_shiftCoeffCosA_->GetBinContent(centBin+1,i);
    double shiftSinC = p_shiftCoeffSinC_->GetBinContent(centBin+1,i);
    double shiftCosC = p_shiftCoeffCosC_->GetBinContent(centBin+1,i);
    double shiftSinFull = p_shiftCoeffSinFull_->GetBinContent(centBin+1,i);
    double shiftCosFull = p_shiftCoeffCosFull_->GetBinContent(centBin+1,i);
    psiA += (1./i)*(-shiftSinA*cos(i*psiA)+shiftCosA*sin(i*psiA));
    psiC += (1./i)*(-shiftSinC*cos(i*psiC)+shiftCosC*sin(i*psiC));
    psiFull += (1./i)*(-shiftSinFull*cos(i*psiFull)+shiftCosFull*sin(i*psiFull));
  }
  h_psiZDC_A[2][centBin]->Fill(psiA);
  h_psiZDC_C[2][centBin]->Fill(psiC);
  h_psiZDC_full[2][centBin]->Fill(psiFull);

  p_psiZDC_reso->Fill(centBin,cos(psiA-psiC));

  //------------------
  // track
  //------------------

  int nTrk = fAOD->GetNumberOfTracks();
  h_mult->Fill(nTrk);

  std::vector<TLorentzVector> v_eminus;
  std::vector<TLorentzVector> v_eplus;

  for (int iTrk=0; iTrk<nTrk; ++iTrk) {
    AliAODTrack *track = (AliAODTrack*)fAOD->GetTrack(iTrk);
    if (!track) {
      AliError(Form("%s: Could not get Track", GetName()));
      continue;
    }

    // if (!track->TestFilterBit(1)) continue;

    double p     = track->P();
    double pt     = track->Pt();
    double eta    = track->Eta();
    double phi    = track->Phi();
    int   charge = track->Charge();
    double dedx   = track->GetTPCsignal();
    int   nhits_tpc = track->GetTPCNcls();
    int   nhits_its = track->GetITSNcls();
    
    if (pt<0.4) continue;
    h_pt->Fill(pt);
    h_phi->Fill(phi);
    h_eta[0]->Fill(eta);
    if (fabs(eta)>0.8) continue;
    h_eta[1]->Fill(eta);
    h_nHits[0]->Fill(nhits_tpc);
    h_nHits[2]->Fill(nhits_its);
    if (fabs(nhits_tpc)<50) continue;
    if (fabs(nhits_its)<4) continue;
    h_nHits[1]->Fill(nhits_tpc);
    h_nHits[3]->Fill(nhits_its);

    //------------------
    // dca
    //------------------

    double dcaxy  = 999.;
    double dcaz   = 999.;
    double r[3];
    double dca[2];
    double cov[3];
    double mag = fAOD->GetMagneticField();
    bool proptodca = track->PropagateToDCA(fVtx, mag, 100., dca, cov);
    if (track->GetXYZ(r)) {
      dcaxy = r[0];
      dcaz  = r[1];
    } else {
      double dcax = r[0] - vx;
      double dcay = r[1] - vy;
      dcaz  = r[2] - vz;
      dcaxy = sqrt(dcax*dcax + dcay*dcay);
      // dcaxy = dca[0];
    }
    h_dcaXy[0]->Fill(dcaxy);
    if (fabs(dcaxy)>1) continue;
    h_dcaXy[1]->Fill(dcaxy);
    h_dcaZ[0]->Fill(dcaz);
    if (fabs(dcaz)>3) continue;
    h_dcaZ[1]->Fill(dcaz);

    //------------------
    // PID
    //------------------

    // double nSigmaTPCPi = fPID->NumberOfSigmasTPC(track, (AliPID::EParticleType)AliPID::kPion);
    // double nSigmaTPCK  = fPID->NumberOfSigmasTPC(track, (AliPID::EParticleType)AliPID::kKaon);
    // double nSigmaTPCP  = fPID->NumberOfSigmasTPC(track, (AliPID::EParticleType)AliPID::kProton);
    // double nSigmaTOFPi = fPID->NumberOfSigmasTOF(track, (AliPID::EParticleType)AliPID::kPion);
    // double nSigmaTOFK  = fPID->NumberOfSigmasTOF(track, (AliPID::EParticleType)AliPID::kKaon);
    // double nSigmaTOFP  = fPID->NumberOfSigmasTOF(track, (AliPID::EParticleType)AliPID::kProton);

    if (p<0.4 || p>3.5) continue;
    double nSigmaTPC_e = fPID->NumberOfSigmasTPC(track, (AliPID::EParticleType)AliPID::kElectron);
    double nSigmaITS_e = fPID->NumberOfSigmasITS(track, (AliPID::EParticleType)AliPID::kElectron);
    double nSigmaTOF_e = fPID->NumberOfSigmasTOF(track, (AliPID::EParticleType)AliPID::kElectron);
    h_nSigmaE_p[0]->Fill(p*charge, nSigmaTPC_e);
    h_nSigmaE_p[2]->Fill(p*charge, nSigmaITS_e);
    h_nSigmaE_p[4]->Fill(p*charge, nSigmaTOF_e);
    bool isElectron = false;
    if (nSigmaITS_e>-4 && nSigmaITS_e<1 &&
        fabs(nSigmaTOF_e)<3 &&
        nSigmaTPC_e>-1.5 && nSigmaTPC_e<3) isElectron = true;
    if (isElectron) {
      h_nSigmaE_p[1]->Fill(p*charge, nSigmaTPC_e);
      h_nSigmaE_p[3]->Fill(p*charge, nSigmaITS_e);
      h_nSigmaE_p[5]->Fill(p*charge, nSigmaTOF_e);

      TLorentzVector e;
      e.SetXYZM(track->Px(), track->Py(), track->Pz(), 0.511e-3);
      if (charge<0) v_eminus.push_back(e);
      else v_eplus.push_back(e);
    }

    //------------------
    // pair
    //------------------
    for (unsigned i=0; i<v_eminus.size(); ++i) {
      for (unsigned j=0; j<v_eplus.size(); ++j) {
        TLorentzVector e = v_eminus.at(i);
        TLorentzVector p = v_eplus.at(j);
        TLorentzVector pair = e+p;
        h_mass_epem->Fill(pair.M());
      }
    }
  }

  PostData(1,fOutputList);
}

//---------------------------------------------------

int AliAnalysisTaskEE::GetRunNumBin(int runNum)
{
  int runNumBin=-1;
  int runNumList[39]={170387,170040,170268,170228,170207,169838,170159,170204,170311,170084,169835,170088,170593,170203,170270,169846,170163,170388,170155,170083,170572,169837,169855,170306,170269,170089,170309,170091,170081,170230,170085,170315,170027,170193,170312,170313,170308,169858,169859};    
  for (int i=0; i<39; ++i) {
    if (runNum==runNumList[i]) {runNumBin=i; break;}
    else continue;
  }

  return runNumBin;
}

void AliAnalysisTaskEE::GetEPMean(int runNum, int centBin, double vx, double vy)
{
  p_XAMean_ = (TProfile3D*)listEPMean_->FindObject(Form("XA_%i", runNum));
  p_XCMean_ = (TProfile3D*)listEPMean_->FindObject(Form("XC_%i", runNum));
  p_YAMean_ = (TProfile3D*)listEPMean_->FindObject(Form("YA_%i", runNum));
  p_YCMean_ = (TProfile3D*)listEPMean_->FindObject(Form("YC_%i", runNum));

  int vxBin = p_XAMean_->GetYaxis()->FindBin(vx);
  int vyBin = p_XAMean_->GetZaxis()->FindBin(vy);
  XAMean_ = p_XAMean_->GetBinContent(centBin+1,vxBin,vyBin);
  XCMean_ = p_XCMean_->GetBinContent(centBin+1,vxBin,vyBin);
  YAMean_ = p_YAMean_->GetBinContent(centBin+1,vxBin,vyBin);
  YCMean_ = p_YCMean_->GetBinContent(centBin+1,vxBin,vyBin);
}
