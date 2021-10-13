#include <iostream>
#include <cstdlib>
#include <sys/time.h>
// ROOT classes
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TList.h"
#include "TH1I.h"
#include "TH3I.h"
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

#include "AliAnalysisTaskEPCalib.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskEPCalib);

//---------------------------------------------------
AliAnalysisTaskEPCalib::AliAnalysisTaskEPCalib() : 
  AliAnalysisTaskSE()
{
  FillNUAWeight_=false;
  NUAWeightOn_  =false;
  FillRecenter_ =false;
  FillShift_    =false;
  TPCCalibOn_   =false;

  fileNameNUAEntry_="alien:///alice/cern.ch/user/q/qshou/EPCalib/LHC11h_pass2/recentering/Merged_AnalysisResults.root";
  fileNUAEntry_=NULL;
  listNUAEntry_=NULL;
  runNum_      =-999;
  vzBin_       =-999;
  centBin_     =-999;
  harmonic_    =2;
  // p_XAMean_=NULL;
  // p_XCMean_=NULL;
  // p_YAMean_=NULL;
  // p_YCMean_=NULL;
  // XAMean_=0;
  // XCMean_=0;
  // YAMean_=0;
  // YCMean_=0;

  fOutputList=NULL;
  h_evtCount =NULL;
  h_runNumBin=NULL;
  h_cent     =NULL;
  for (int i=0; i<2; ++i) h_vxy[i]=NULL;
  for (int i=0; i<2; ++i) h_vz[i] =NULL;
  for (int i=0; i<2; ++i) h_phi[i]=NULL;
  for (int iCent=0; iCent<10; ++iCent) {
    for (int i = 0; i < 2; ++i) {
      h_NUAEntryWrite[iCent][i] = NULL;
    }
  }
  for (int iCent=0; iCent<10; ++iCent) {
    for (int i=0; i<6; ++i) {
      h_psiTPC[iCent][i] = NULL;
    }
  }
}

//---------------------------------------------------
AliAnalysisTaskEPCalib::AliAnalysisTaskEPCalib(const char *name) : 
  AliAnalysisTaskSE(name)
{
  FillNUAWeight_=false;
  NUAWeightOn_  =false;
  FillRecenter_ =false;
  FillShift_    =false;
  TPCCalibOn_   =false;

  fileNameNUAEntry_="alien:///alice/cern.ch/user/q/qshou/EPCalib/LHC11h_pass2/recentering/Merged_AnalysisResults.root";
  fileNUAEntry_=NULL;
  listNUAEntry_=NULL;
  runNum_      =-999;
  vzBin_       =-999;
  centBin_     =-999;
  harmonic_    =2;
  // p_XAMean_=NULL;
  // p_XCMean_=NULL;
  // p_YAMean_=NULL;
  // p_YCMean_=NULL;
  // XAMean_=0;
  // XCMean_=0;
  // YAMean_=0;
  // YCMean_=0;

  fOutputList=NULL;
  h_evtCount =NULL;
  h_runNumBin=NULL;
  h_cent     =NULL;
  for (int i=0; i<2; ++i) h_vxy[i]=NULL;
  for (int i=0; i<2; ++i) h_vz[i] =NULL;
  for (int i=0; i<2; ++i) h_phi[i]=NULL;
  for (int iCent=0; iCent<10; ++iCent) {
    for (int i = 0; i < 2; ++i) {
      h_NUAEntryWrite[iCent][i] = NULL;
    }
  }
  for (int iCent=0; iCent<10; ++iCent) {
    for (int i=0; i<6; ++i) {
      h_psiTPC[iCent][i] = NULL;
    }
  }

  // for (int i=0; i<4; ++i) {    
  //   h_ZNA[i]=NULL; 
  //   h_ZNC[i]=NULL;
  // }
  // for (int i=0; i<39; ++i) {    
  //   p_XA[i]=NULL;
  //   p_XC[i]=NULL;
  //   p_YA[i]=NULL;
  //   p_YC[i]=NULL;
  // }
  // for (int i=0; i<2; ++i) {
  //   h_XYA[i]=NULL;
  //   h_XYC[i]=NULL;
  //   h_psiZDC_A[i]=NULL;
  //   h_psiZDC_C[i]=NULL;
  //   h_psiZDC_full[i]=NULL;
  // }
  // for (int i=0; i<39; ++i) {
  //   p_shiftCoeffSinA[i]=NULL;
  //   p_shiftCoeffCosA[i]=NULL;
  //   p_shiftCoeffSinC[i]=NULL;
  //   p_shiftCoeffCosC[i]=NULL;
  //   p_shiftCoeffSinFull[i]=NULL;
  //   p_shiftCoeffCosFull[i]=NULL;
  // }

  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
}

//---------------------------------------------------
AliAnalysisTaskEPCalib::~AliAnalysisTaskEPCalib()
{}

//---------------------------------------------------
void AliAnalysisTaskEPCalib::Terminate(Option_t *) 
{}

//---------------------------------------------------
void AliAnalysisTaskEPCalib::UserCreateOutputObjects()
{
  fOutputList = new TList();
  fOutputList->SetName(GetName());
  fOutputList->SetOwner(kTRUE);

  if (NUAWeightOn_) {
    fileNUAEntry_ = TFile::Open(fileNameNUAEntry_.Data(),"READ");
    if (!fileNUAEntry_) {
      AliError(Form("%s: Wrong NUAEntry file.", GetName()));
      return;
    }
    listNUAEntry_ = (TList*)fileNUAEntry_->Get("output");
    if (!listNUAEntry_) {
      AliError(Form("%s: Wrong NUAEntry file.", GetName()));
      return;
    }
  }

  // event-wise
  h_evtCount = new TH1I("evtCount","",20, 1, 21);
  h_evtCount->GetXaxis()->SetBinLabel(1,"All");
  h_evtCount->GetXaxis()->SetBinLabel(2,"Info");
  h_evtCount->GetXaxis()->SetBinLabel(3,"Evt");
  h_evtCount->GetXaxis()->SetBinLabel(4,"Cent");
  h_evtCount->GetXaxis()->SetBinLabel(5,"TPCEP");
  h_evtCount->GetXaxis()->SetBinLabel(10,"Manager");
  h_evtCount->GetXaxis()->SetBinLabel(11,"Handler");
  h_evtCount->GetXaxis()->SetBinLabel(12,"AOD");
  // h_evtCount->GetXaxis()->SetBinLabel(13,"Utils");
  // h_evtCount->GetXaxis()->SetBinLabel(15,"MultSel");
  fOutputList->Add(h_evtCount);

  TString runNum[39]={"170387","170040","170268","170228","170207","169838","170159","170204","170311","170084",
                      "169835","170088","170593","170203","170270","169846","170163","170388","170155","170083",
                      "170572","169837","169855","170306","170269","170089","170309","170091","170081","170230",
                      "170085","170315","170027","170193","170312","170313","170308","169858","169859"};
  h_runNumBin = new TH1I("runNumBin","",39,0,39);
  for (int i=0; i<39; ++i) {    
    h_runNumBin->GetXaxis()->SetBinLabel(i,runNum[i].Data());
  }
  fOutputList->Add(h_runNumBin);

  h_cent = new TH1D("centrality","",100,0,100);
  fOutputList->Add(h_cent);
  h_vxy[0] = new TH2D("vxy0","",100,-0.5,0.5,100,-0.5,0.5);
  h_vxy[1] = new TH2D("vxy1","",100,-0.5,0.5,100,-0.5,0.5);
  h_vz[0]  = new TH1D("vz0","",200,-50,50);
  h_vz[1]  = new TH1D("vz1","",200,-50,50);
  for (int i=0; i<2; ++i) {
    fOutputList->Add(h_vxy[i]);
    fOutputList->Add(h_vz[i]);
  }

  h_phi[0] = new TH1D("phi0", "",360,0,2*TMath::Pi());
  h_phi[1] = new TH1D("phi1", "",360,0,2*TMath::Pi());
  fOutputList->Add(h_phi[0]);
  fOutputList->Add(h_phi[1]);

  for (int iCent=0; iCent<10; ++iCent) {
    h_NUAEntryWrite[iCent][0] = new TH3I(Form("NUAEntryPos_cent%i_run",iCent),"",20,-10,10,360,0,2*TMath::Pi(),20,-0.8,0.8);
    h_NUAEntryWrite[iCent][1] = new TH3I(Form("NUAEntryNeg_cent%i_run",iCent),"",20,-10,10,360,0,2*TMath::Pi(),20,-0.8,0.8);
    if (FillNUAWeight_) {
      fOutputList->Add(h_NUAEntryWrite[iCent][0]);
      fOutputList->Add(h_NUAEntryWrite[iCent][1]);
    }
  }

  for (int iCent=0; iCent<10; ++iCent) {
    for (int i=0; i<6; ++i) {
      h_psiTPC[iCent][i] = new TH1D(Form("psiTPC_cent%i_%i",iCent,i), "", 360, -2*TMath::Pi(), 2*TMath::Pi());
      fOutputList->Add(h_psiTPC[iCent][i]);
    }
  }

  // h_ZNA[0] = new TH1D("zna_1", "", 40000, 0, 40000);
  // h_ZNA[1] = new TH1D("zna_2", "", 40000, 0, 40000);
  // h_ZNA[2] = new TH1D("zna_3", "", 40000, 0, 40000);
  // h_ZNA[3] = new TH1D("zna_4", "", 40000, 0, 40000);
  // h_ZNC[0] = new TH1D("znc_1", "", 40000, 0, 40000);
  // h_ZNC[1] = new TH1D("znc_2", "", 40000, 0, 40000);
  // h_ZNC[2] = new TH1D("znc_3", "", 40000, 0, 40000);
  // h_ZNC[3] = new TH1D("znc_4", "", 40000, 0, 40000);
  // for (int i=0; i<39; ++i) {
  //   p_XA[i] = new TProfile3D((TString("XA_")+=runNum[i]).Data(),"",10,0,10,40,-2,2,40,-2,2);
  //   p_XC[i] = new TProfile3D((TString("XC_")+=runNum[i]).Data(),"",10,0,10,40,-2,2,40,-2,2);
  //   p_YA[i] = new TProfile3D((TString("YA_")+=runNum[i]).Data(),"",10,0,10,40,-2,2,40,-2,2);
  //   p_YC[i] = new TProfile3D((TString("YC_")+=runNum[i]).Data(),"",10,0,10,40,-2,2,40,-2,2);
  // }
  // if (!doShift_) {
  //   for (int i=0; i<4; ++i) {
  //     fOutputList->Add(h_ZNA[i]); 
  //     fOutputList->Add(h_ZNC[i]);
  //   }

  //   for (int i=0; i<39; ++i) {
  //     fOutputList->Add(p_XA[i]);
  //     fOutputList->Add(p_XC[i]);
  //     fOutputList->Add(p_YA[i]);
  //     fOutputList->Add(p_YC[i]);
  //   }
  // }
  // h_XYA[0] = new TH2D("XYA_raw","",200,-1,1,200,-1,1);
  // h_XYC[0] = new TH2D("XYC_raw","",200,-1,1,200,-1,1);
  // h_XYA[1] = new TH2D("XYA_afR","",200,-1,1,200,-1,1);
  // h_XYC[1] = new TH2D("XYC_afR","",200,-1,1,200,-1,1);
  // fOutputList->Add(h_XYA[0]);
  // fOutputList->Add(h_XYC[0]);
  // h_psiZDC_A[0] = new TH1D("psiZDC_A_raw", "", 100, -2*TMath::Pi(), 2*TMath::Pi());
  // h_psiZDC_C[0] = new TH1D("psiZDC_C_raw", "", 100, -2*TMath::Pi(), 2*TMath::Pi());
  // h_psiZDC_full[0] = new TH1D("psiZDC_full_raw", "", 100, -2*TMath::Pi(), 2*TMath::Pi()); 
  // h_psiZDC_A[1] = new TH1D("psiZDC_A_afR", "", 100, -2*TMath::Pi(), 2*TMath::Pi());
  // h_psiZDC_C[1] = new TH1D("psiZDC_C_afR", "", 100, -2*TMath::Pi(), 2*TMath::Pi());
  // h_psiZDC_full[1] = new TH1D("psiZDC_full_afR", "", 100, -2*TMath::Pi(), 2*TMath::Pi()); 
  // fOutputList->Add(h_psiZDC_A[0]);
  // fOutputList->Add(h_psiZDC_C[0]);
  // fOutputList->Add(h_psiZDC_full[0]);
  // for (int i=0; i<39; ++i) {
  //   p_shiftCoeffSinA[i] = new TProfile2D((TString("shiftCoeffSinA_")+=runNum[i]).Data(),"",10,0,10,20,0,20);
  //   p_shiftCoeffCosA[i] = new TProfile2D((TString("shiftCoeffCosA_")+=runNum[i]).Data(),"",10,0,10,20,0,20);
  //   p_shiftCoeffSinC[i] = new TProfile2D((TString("shiftCoeffSinC_")+=runNum[i]).Data(),"",10,0,10,20,0,20);
  //   p_shiftCoeffCosC[i] = new TProfile2D((TString("shiftCoeffCosC_")+=runNum[i]).Data(),"",10,0,10,20,0,20);
  //   p_shiftCoeffSinFull[i] = new TProfile2D((TString("shiftCoeffSinFull_")+=runNum[i]).Data(),"",10,0,10,20,0,20);
  //   p_shiftCoeffCosFull[i] = new TProfile2D((TString("shiftCoeffCosFull_")+=runNum[i]).Data(),"",10,0,10,20,0,20);
  // }
  // if (doShift_) {
  //   fOutputList->Add(h_XYA[1]);
  //   fOutputList->Add(h_XYC[1]);
  //   fOutputList->Add(h_psiZDC_A[1]);
  //   fOutputList->Add(h_psiZDC_C[1]);
  //   fOutputList->Add(h_psiZDC_full[1]);
  //   for (int i=0; i<39; ++i) {
  //     fOutputList->Add(p_shiftCoeffSinA[i]);
  //     fOutputList->Add(p_shiftCoeffCosA[i]);
  //     fOutputList->Add(p_shiftCoeffSinC[i]);
  //     fOutputList->Add(p_shiftCoeffCosC[i]);
  //     fOutputList->Add(p_shiftCoeffSinFull[i]);
  //     fOutputList->Add(p_shiftCoeffCosFull[i]);
  //   }
  // }

  PostData(1,fOutputList);
}

//---------------------------------------------------
void AliAnalysisTaskEPCalib::UserExec(Option_t *)
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
  AliAnalysisUtils* fUtils = new AliAnalysisUtils();
  if (!fUtils) {
    AliError(Form("%s: Could not get AliAnalysisUtils", GetName()));
  } else h_evtCount->Fill(13);
  // AliMultSelection* fMultSel = (AliMultSelection*)InputEvent()->FindListObject("MultSelection");
  // if (!fMultSel) {
  //   AliError(Form("%s: Could not get AliMultSelection", GetName()));
  // } else h_evtCount->Fill(15);
  if (!manager || !handler || !fAOD) return;
  h_evtCount->Fill(2);

  //------------------
  // event
  //------------------
  runNum_    = fAOD->GetRunNumber();
  int runNumBin = GetRunNumBin(runNum_);
  if (runNumBin<0) return;
  h_runNumBin->Fill(runNumBin);
  // run by run
  for (int iCent=0; iCent<10; ++iCent) {
    TString name = h_NUAEntryWrite[iCent][0]->GetName();
    name += runNum_;
    h_NUAEntryWrite[iCent][0]->SetName(name);
    name = h_NUAEntryWrite[iCent][1]->GetName();
    name += runNum_;
    h_NUAEntryWrite[iCent][1]->SetName(name);
  }

  // vertex
  AliAODVertex* fVtx = fAOD->GetPrimaryVertex();
  float vx    = fVtx->GetX();
  float vy    = fVtx->GetY();
  float vz    = fVtx->GetZ();
  // float vzSPD = fAOD->GetPrimaryVertexSPD()->GetZ();
  h_vxy[0]->Fill(vx, vy);
  h_vz[0]->Fill(vz);
  if (fabs(vx)<1e-6 || fabs(vy)<1e-6 || fabs(vz)<1e-6) return;
  // https://twiki.cern.ch/twiki/bin/viewauth/ALICE/AliDPGtoolsEventProp
  // AliEventCuts fEventCuts;
  // if (!fEventCuts.AcceptEvent(fAOD)) {
  //   // PostData(1, fOutputList);
  //   return;
  // }
  if (fabs(vz)>10) return;
  h_vxy[1]->Fill(vx, vy);
  h_vz[1]->Fill(vz);
  for (int i = 0; i < 20; ++i) {
    if (vz > -10+i*1 && vz < -10+(i+1)*1) {vzBin_ = i; break;}
  }
  if (vzBin_<0) return;

  // pileup
  fUtils->SetUseOutOfBunchPileUp(true);
  bool isPileup = fUtils->IsPileUpEvent(fAOD);
  // bool isPileup = fUtils->IsPileUpMV(fAOD);
  if (isPileup) return;
  h_evtCount->Fill(3);

  // centrality
  double cent = -999;
  // run1
  cent = fAOD->GetCentrality()->GetCentralityPercentile("V0M");
  // cent = event->GetCentrality()->GetCentralityPercentile("CL1");
  // cent = event->GetCentrality()->GetCentralityPercentile("TRK");
  if (cent<0) return;
  centBin_ = (int)cent/10;
  h_cent->Fill(cent);
  h_evtCount->Fill(4);

  // V0 plane
  // fEP = fAOD->GetEventplane();
  // double fEPV0=-999, fEPV0A=-999, fEPV0C=-999;
  // if (fEP) {
  //   fEPV0  = TVector2::Phi_0_2pi(fEP->GetEventplane("V0",  fAOD)); if (fEPV0>TMath::Pi()) fEPV0-=TMath::Pi();
  //   fEPV0A = TVector2::Phi_0_2pi(fEP->GetEventplane("V0A", fAOD)); if (fEPV0A>TMath::Pi()) fEPV0A-=TMath::Pi();
  //   fEPV0C = TVector2::Phi_0_2pi(fEP->GetEventplane("V0C", fAOD)); if (fEPV0C>TMath::Pi()) fEPV0C-=TMath::Pi();
  // } //cout<<fEPV0<<" "<<fEPV0A<<" "<<fEPV0C<<endl;

  if (!TPCPlane(fAOD)) return;
  h_evtCount->Fill(5);

  PostData(1,fOutputList);
}

//---------------------------------------------------
int AliAnalysisTaskEPCalib::GetRunNumBin(int runNum)
{
  int runNumBin=-1;
  int runNumList[39]={170387,170040,170268,170228,170207,169838,170159,170204,170311,170084,
                      169835,170088,170593,170203,170270,169846,170163,170388,170155,170083,
                      170572,169837,169855,170306,170269,170089,170309,170091,170081,170230,
                      170085,170315,170027,170193,170312,170313,170308,169858,169859};    
  for (int i = 0; i < 39; ++i) {
    if (runNum==runNumList[i]) {runNumBin=i; break;}
    else continue;
  }
  return runNumBin;
}

//---------------------------------------------------
bool AliAnalysisTaskEPCalib::TPCPlane(AliAODEvent* fAOD)
{
  TH3I* h_NUAEntryPosRead = NULL;
  TH3I* h_NUAEntryNegRead = NULL;
  if (NUAWeightOn_) {
    h_NUAEntryPosRead = (TH3I*)listNUAEntry_->FindObject(Form("NUAEntryPos_cent%i_run%i",centBin_,runNum_));
    h_NUAEntryNegRead = (TH3I*)listNUAEntry_->FindObject(Form("NUAEntryNeg_cent%i_run%i",centBin_,runNum_));
  }

  AliAODVertex* fVtx = fAOD->GetPrimaryVertex();
  double vtx[3] = {(double)fVtx->GetX(), (double)fVtx->GetY(), (double)fVtx->GetZ()};
  double mag    = fAOD->GetMagneticField();

  double sumCos[3]  = {0};
  double sumSin[3]  = {0};
  double sumCosA[3] = {0};
  double sumSinA[3] = {0};
  double sumCosC[3] = {0};
  double sumSinC[3] = {0};

  int nTrk = fAOD->GetNumberOfTracks();
  for (int iTrk = 0; iTrk < nTrk; ++iTrk) {
    AliAODTrack* track = (AliAODTrack*)fAOD->GetTrack(iTrk);
    if (!track) {
      AliError(Form("%s: Could not get Track", GetName()));
      continue;
    }
    if (!track->TestFilterBit(768)) continue;
    double pt     = track->Pt();
    double eta    = track->Eta();
    double phi    = track->Phi();
    int    charge = track->Charge();
    int    nhits  = track->GetTPCNcls();
    double dedx   = track->GetTPCsignal();
    double chi2   = track->Chi2perNDF();
    if (pt<0.2 || pt>5) continue;
    if (fabs(eta)>0.8) continue;
    if (fabs(nhits)<70) continue;
    if (chi2<0.1 || chi2>4.) continue;
    if (dedx<10) continue;
    double dcaxy = 999.;
    double dcaz  = 999.;
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
    if (fabs(dcaxy)>2.4) continue;
    if (fabs(dcaz)>3.2) continue;
    h_phi[0]->Fill(phi);

    if (FillNUAWeight_) {
      if (charge>0)      h_NUAEntryWrite[centBin_][0]->Fill(vtx[2],phi,eta);
      else if (charge<0) h_NUAEntryWrite[centBin_][1]->Fill(vtx[2],phi,eta);
    }

    double NUAWeight = 1;
    double cosnphi = NUAWeight*cos(harmonic_*phi);
    double sinnphi = NUAWeight*sin(harmonic_*phi);
    sumCos[0] += cosnphi;
    sumSin[0] += sinnphi;
    if (eta>0) {
      sumCosA[0] += cosnphi;
      sumSinA[0] += sinnphi;
    } else {
      sumCosC[0] += cosnphi;
      sumSinC[0] += sinnphi;
    }

    if (NUAWeightOn_) {
      int NUAEntry = 0;
      int phiBin = h_NUAEntryPosRead->GetYaxis()->FindBin(phi);
      int etaBin = h_NUAEntryPosRead->GetZaxis()->FindBin(eta);
      if (charge>0)      NUAEntry = h_NUAEntryPosRead->GetBinContent(vzBin_,phiBin,etaBin);
      else if (charge<0) NUAEntry = h_NUAEntryNegRead->GetBinContent(vzBin_,phiBin,etaBin);
      if (NUAEntry==0) continue;
      NUAWeight = 1./NUAEntry;

      h_phi[1]->Fill(phi, NUAEntry);
      cosnphi = NUAWeight*cos(harmonic_*phi);
      sinnphi = NUAWeight*sin(harmonic_*phi);
      sumCos[1] += cosnphi;
      sumSin[1] += sinnphi;
      if (eta>0) {
        sumCosA[1] += cosnphi;
        sumSinA[1] += sinnphi;
      } else {
        sumCosC[1] += cosnphi;
        sumSinC[1] += sinnphi;
      }
    }


    // recentering
    if (FillRecenter_) {
      // if (eta>0 && vtx[2]>0) {
      //   cosRecenterFW->Fill(runNumberPointer,centBin_,cosnphi);
      //   sinRecenterFW->Fill(runNumberPointer,centBin_,sinnphi);
      // }
      // else if (eta>0 && vtx[2]<0) {
      //   cosRecenterW->Fill(runNumberPointer,centBin_,cosnphi);
      //   sinRecenterW->Fill(runNumberPointer,centBin_,sinnphi);
      // }
      // else if(eta<0 && vtx[2]>0) {
      //   cosRecenterE->Fill(runNumberPointer,centBin_,cosnphi);
      //   sinRecenterE->Fill(runNumberPointer,centBin_,sinnphi);
      // }
      // else if(eta<0 && vtx[2]<0) {
      //   cosRecenterFE->Fill(runNumberPointer,centBin_,cosnphi);
      //   sinRecenterFE->Fill(runNumberPointer,centBin_,sinnphi);
      // }
    }

    if (TPCCalibOn_) {
      // double cosMean = cosRecenterW->GetBinContent(centBin_+1,vxBin,vyBin);
      // double sinMean = sinRecenterE->GetBinContent(centBin_+1,vxBin,vyBin);
      // sumCos[2] += (cosnphi-cosMean);
      // sumSin[2] += (sinnphi-sinMean);
      // if (eta>0) {
      //   sumCosA[2] += (cosnphi-cosMean);
      //   sumSinA[2] += (sinnphi-sinMean);
      // } else {
      //   sumCosC[2] += (cosnphi-cosMean);
      //   sumSinC[2] += (sinnphi-sinMean);
      // }
    }
  }

  if (sumCos[0] <1e-6 && sumSin[0] <1e-6) return false;
  if (sumCosA[0]<1e-6 && sumSinA[0]<1e-6) return false;
  if (sumCosC[0]<1e-6 && sumSinC[0]<1e-6) return false;
  TVector2 Q0, Q0A, Q0C;
  Q0.Set(sumCos[0], sumSin[0]);
  Q0A.Set(sumCosA[0], sumSinA[0]);
  Q0C.Set(sumCosC[0], sumSinC[0]);
  double psi0  = Q0.Phi()/harmonic_; 
  double psi0A = Q0A.Phi()/harmonic_; 
  double psi0C = Q0C.Phi()/harmonic_; 
  h_psiTPC[centBin_][0]->Fill(psi0);
  h_psiTPC[centBin_][1]->Fill(psi0A);
  h_psiTPC[centBin_][2]->Fill(psi0C);

  if (NUAWeightOn_) {
    if (sumCos[1] <1e-6 && sumSin[1] <1e-6) return false;
    if (sumCosA[1]<1e-6 && sumSinA[1]<1e-6) return false;
    if (sumCosC[1]<1e-6 && sumSinC[1]<1e-6) return false;
    TVector2 Q1, Q1A, Q1C;
    Q1.Set(sumCos[1], sumSin[1]);
    Q1A.Set(sumCosA[1], sumSinA[1]);
    Q1C.Set(sumCosC[1], sumSinC[1]);
    double psi1  = Q1.Phi()/harmonic_; 
    double psi1A = Q1A.Phi()/harmonic_; 
    double psi1C = Q1C.Phi()/harmonic_; 
    h_psiTPC[centBin_][3]->Fill(psi1);
    h_psiTPC[centBin_][4]->Fill(psi1A);
    h_psiTPC[centBin_][5]->Fill(psi1C);
  }

  if (TPCCalibOn_) {
    // if (sumCos[2] <1e-6 && sumSin[2] <1e-6) return false;
    // if (sumCosA[2]<1e-6 && sumSinA[2]<1e-6) return false;
    // if (sumCosC[2]<1e-6 && sumSinC[2]<1e-6) return false;
    // TVector2 Q1, Q1A, Q1C;
    // Q1.Set(sumCos[2], sumSin[2]);
    // Q1A.Set(sumCosA[2], sumSinA[2]);
    // Q1C.Set(sumCosC[2], sumSinC[2]);
    // double psi1  = Q1.Phi()/harmonic_; 
    // double psi1A = Q1A.Phi()/harmonic_; 
    // double psi1C = Q1C.Phi()/harmonic_; 
    // h_psiTPC[6]->Fill(psi1);
    // h_psiTPC[7]->Fill(psi1A);
    // h_psiTPC[8]->Fill(psi1C);

    // p_shiftCoeffSinA_ = (TProfile2D*)listEPShifting_->FindObject(Form("shiftCoeffSinA_%i", runNum));
    // p_shiftCoeffCosA_ = (TProfile2D*)listEPShifting_->FindObject(Form("shiftCoeffCosA_%i", runNum));
    // p_shiftCoeffSinC_ = (TProfile2D*)listEPShifting_->FindObject(Form("shiftCoeffSinC_%i", runNum));
    // p_shiftCoeffCosC_ = (TProfile2D*)listEPShifting_->FindObject(Form("shiftCoeffCosC_%i", runNum));
    // p_shiftCoeffSinFull_ = (TProfile2D*)listEPShifting_->FindObject(Form("shiftCoeffSinFull_%i", runNum));
    // p_shiftCoeffCosFull_ = (TProfile2D*)listEPShifting_->FindObject(Form("shiftCoeffCosFull_%i", runNum));
    // for (int i = 1; i <= 20; ++i) {\
    //   double shiftSin = p_shiftCoeffSin_->GetBinContent(centBin_+1,i);
    //   double shiftCos = p_shiftCoeffCos_->GetBinContent(centBin_+1,i);
    //   double shiftSinA = p_shiftCoeffSinA_->GetBinContent(centBin_+1,i);
    //   double shiftCosA = p_shiftCoeffCosA_->GetBinContent(centBin_+1,i);
    //   double shiftSinC = p_shiftCoeffSinC_->GetBinContent(centBin_+1,i);
    //   double shiftCosC = p_shiftCoeffCosC_->GetBinContent(centBin_+1,i);
    //   psi1  += (2/i/harmonic_)*(-shiftSin*cos(i*psi1)+shiftCos*sin(i*psi1));
    //   psi1A += (2/i/harmonic_)*(-shiftSinA*cos(i*psi1A)+shiftCosA*sin(i*psi1A));
    //   psi1C += (2/i/harmonic_)*(-shiftSinC*cos(i*psi1C)+shiftCosC*sin(i*psi1C));
    // }
    // h_psiTPC[9]->Fill(psi1);
    // h_psiTPC[10]->Fill(psi1A);
    // h_psiTPC[11]->Fill(psi1C);
  }

  if (FillShift_) {
    // for (int i = 1; i <= 20; ++i) {
    //   p_shiftCoeffCos[runNumBin] ->Fill(centBin_, i-1, cos(i*harmonic_*psiFull));
    //   p_shiftCoeffSin[runNumBin] ->Fill(centBin_, i-1, sin(i*harmonic_*psiFull));
    //   p_shiftCoeffCosA[runNumBin]->Fill(centBin_, i-1, cos(i*harmonic_*psiA));
    //   p_shiftCoeffSinA[runNumBin]->Fill(centBin_, i-1, sin(i*harmonic_*psiA));
    //   p_shiftCoeffCosC[runNumBin]->Fill(centBin_, i-1, cos(i*harmonic_*psiC));
    //   p_shiftCoeffSinC[runNumBin]->Fill(centBin_, i-1, sin(i*harmonic_*psiC));
    // }
  }

  return true;
}

// void AliAnalysisTaskEPCalib::GetEPMean(int runNum, int centBin, double vx, double vy)
// {
//   p_XAMean_ = (TProfile3D*)listEPMean_->FindObject(Form("XA_%i", runNum));
//   p_XCMean_ = (TProfile3D*)listEPMean_->FindObject(Form("XC_%i", runNum));
//   p_YAMean_ = (TProfile3D*)listEPMean_->FindObject(Form("YA_%i", runNum));
//   p_YCMean_ = (TProfile3D*)listEPMean_->FindObject(Form("YC_%i", runNum));
//   int vxBin = p_XAMean_->GetYaxis()->FindBin(vx);
//   int vyBin = p_XAMean_->GetZaxis()->FindBin(vy);
//   XAMean_ = p_XAMean_->GetBinContent(centBin+1,vxBin,vyBin);
//   XCMean_ = p_XCMean_->GetBinContent(centBin+1,vxBin,vyBin);
//   YAMean_ = p_YAMean_->GetBinContent(centBin+1,vxBin,vyBin);
//   YCMean_ = p_YCMean_->GetBinContent(centBin+1,vxBin,vyBin);

//   // ZDC plane
//   AliAODZDC* aodZDC = fAOD->GetZDCData();
//   if (!aodZDC) return;
//   const double* energyZNA  = aodZDC->GetZNATowerEnergy();
//   const double* energyZNC  = aodZDC->GetZNCTowerEnergy();
//   if (energyZNA[1]<0 ||
//       energyZNA[2]<0 ||
//       energyZNA[3]<0 ||
//       energyZNA[4]<0 ||
//       energyZNC[1]<0 ||
//       energyZNC[2]<0 ||
//       energyZNC[3]<0 ||
//       energyZNC[4]<0) return;
//   const double x[4] = {-1.75, 1.75, -1.75, 1.75};
//   const double y[4] = {-1.75, -1.75, 1.75, 1.75};
//   double xEA=0, xEC=0, yEA=0, yEC=0, EA=0, EC=0;
//   for (int i = 0; i < 4; ++i) {
//     if (!doShift_) {
//       h_ZNA[i]->Fill(energyZNA[i+1]);
//       h_ZNC[i]->Fill(energyZNC[i+1]);
//     }
//     xEA += x[i]*energyZNA[i+1];
//     yEA += y[i]*energyZNA[i+1];
//     xEC += x[i]*energyZNC[i+1];
//     yEC += y[i]*energyZNC[i+1];
//     EA += energyZNA[i+1];
//     EC += energyZNC[i+1];
//   }

//   double XA = xEA/EA;
//   double XC = xEC/EC;
//   double YA = yEA/EA;
//   double YC = yEC/EC;
//   if (!doShift_) {
//     p_XA[runNumBin]->Fill(centBin,vx,vy,XA);
//     p_XC[runNumBin]->Fill(centBin,vx,vy,XC);
//     p_YA[runNumBin]->Fill(centBin,vx,vy,YA);
//     p_YC[runNumBin]->Fill(centBin,vx,vy,YC);
//   }
//   h_XYA[0]->Fill(XA,YA);
//   h_XYC[0]->Fill(XC,YC);
//   double psiA = atan2(YA,XA);
//   double psiC = atan2(YC,XC);
//   double psiFull = atan2(YA+YC,XA+XC);
//   h_psiZDC_A[0]->Fill(psiA);
//   h_psiZDC_C[0]->Fill(psiC);
//   h_psiZDC_full[0]->Fill(psiFull);

//   if (doShift_) {
//     double XA_mean=0,XC_mean=0,YA_mean=0,YC_mean=0;
//     GetEPMean(runNum,centBin,vx,vy);
//     XA -= XAMean_;
//     XC -= XCMean_;
//     YA -= YAMean_;
//     YC -= YCMean_;
//     h_XYA[1]->Fill(XA,YA);
//     h_XYC[1]->Fill(XC,YC);
//     psiA = atan2(YA,XA);
//     psiC = atan2(YC,XC);
//     psiFull = atan2(YA+YC,XA+XC);
//     h_psiZDC_A[1]->Fill(psiA);
//     h_psiZDC_C[1]->Fill(psiC);
//     h_psiZDC_full[1]->Fill(psiFull);

//     for (int i = 1; i <= 20; ++i) {
//       // shiftCoeffSinA[centBin][i-1] += sin(i*psiA);
//       // shiftCoeffCosA[centBin][i-1] += cos(i*psiA);
//       // shiftCoeffSinC[centBin][i-1] += sin(i*psiC);
//       // shiftCoeffCosC[centBin][i-1] += cos(i*psiC);
//       // shiftCoeffSinFull[centBin][i-1] += sin(i*psiFull);
//       // shiftCoeffCosFull[centBin][i-1] += cos(i*psiFull);
//       p_shiftCoeffSinA[runNumBin]->Fill(centBin, i-1, sin(i*psiA));
//       p_shiftCoeffCosA[runNumBin]->Fill(centBin, i-1, cos(i*psiA));
//       p_shiftCoeffSinC[runNumBin]->Fill(centBin, i-1, sin(i*psiC));
//       p_shiftCoeffCosC[runNumBin]->Fill(centBin, i-1, cos(i*psiC));
//       p_shiftCoeffSinFull[runNumBin]->Fill(centBin, i-1, sin(i*psiFull));
//       p_shiftCoeffCosFull[runNumBin]->Fill(centBin, i-1, cos(i*psiFull));
//     }
//   }
// }
