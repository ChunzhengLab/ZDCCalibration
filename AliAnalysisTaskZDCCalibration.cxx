/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* AliAnaysisTaskMyTask
 *
 * empty task which can serve as a starting point for building an analysis
 * as an example, one histogram is filled
 */

#include <iostream>
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "THnSparse.h"
#include "TList.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisUtils.h"
#include "AliAnalysisTaskZDCCalibration.h"

class AliAnalysisTaskZDCCalibration;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskZDCCalibration) // classimp: necessary for root

AliAnalysisTaskZDCCalibration::AliAnalysisTaskZDCCalibration() : AliAnalysisTaskSE(), 
    fDataSet("10h"),
    bFirstFillHistVetex(kTRUE),
    bFillHistForGE(kTRUE),
    bFillHistForRC(kFALSE),
    bFillHistForSF(kFALSE),
    bCalculateV2(kFALSE),
    fFilterBit(1),
    bApplyGE(kFALSE),
    bApplyRC(kFALSE),
    bApplySF(kFALSE),
    bFillQAHist(kFALSE),

    fZDCCalibrationList(nullptr),

    fCentrality(-999),
    fAOD(nullptr),
    fZDC(nullptr),
    fUtils(nullptr),
    fOutputList(nullptr),
    fQAList(nullptr),
    fHistCent(nullptr),
    fHist2DVxVyTot(nullptr),
    fHistVzTot(nullptr),
    
    fHist2DForMeanVxVy(nullptr),
    fHistForMeanVz(nullptr),

    fProfileForZNCGE(nullptr),
    fProfileForZNAGE(nullptr),

    fHn4DForZNCQxRC(nullptr),
    fHn4DForZNCQyRC(nullptr),
    fHn4DForZNCMtRC(nullptr),
    fHn4DForZNAQxRC(nullptr),
    fHn4DForZNAQyRC(nullptr),
    fHn4DForZNAMtRC(nullptr),

    fProfile2DForCosC(nullptr),
    fProfile2DForSinC(nullptr),
    fProfile2DForCosA(nullptr),
    fProfile2DForSinA(nullptr)
{
    for (size_t i = 0; i < 2; i++) {fHistPt[i] = nullptr;}
    for (size_t i = 0; i < 2; i++) {fHistEta[i] = nullptr;}
    for (size_t i = 0; i < 2; i++) {fHistPhi[i] = nullptr;}
    
    for (size_t i = 0; i < 3; i++) {fVtx[i] = -999;}
    for (size_t i = 0; i < 2; i++) {fHist2DCentCorr[i] = nullptr;}
    for (size_t i = 0; i < 2; i++)
    {
      fProfileQxCCentTot[i] = nullptr;
      fProfileQyCCentTot[i] = nullptr;
      fProfileQxACentTot[i] = nullptr;
      fProfileQyACentTot[i] = nullptr;
      fProfileQxAQxCCentTot[i] = nullptr;
      fProfileQxAQyCCentTot[i] = nullptr;
      fProfileQyAQxCCentTot[i] = nullptr;
      fProfileQyAQyCCentTot[i] = nullptr;
    }

    for (size_t i = 0; i < 3; i++)
    {
      fHist2DPsiCCentTot[i] = nullptr;
      fHist2DPsiACentTot[i] = nullptr;
    }

    for (size_t iRun = 0; iRun < fnRunMax; iRun++)
    {
      fRunList[iRun]                  = nullptr;
      fQAListThisRun[iRun]            = nullptr;
      fHist2DVxVy[iRun]               = nullptr;
      fHistVz[iRun]                   = nullptr;
      fProfileZNATowerEnergy[iRun]    = nullptr;
      fProfileZNCTowerEnergy[iRun]    = nullptr;
      fProfileZNATowerEnergyBfGE[iRun]= nullptr;
      fProfileZNCTowerEnergyBfGE[iRun]= nullptr;

      //vx vy sigma vz
      fHn4DQxZNCCentVxVySigmaVz[iRun] = nullptr;
      fHn4DQyZNCCentVxVySigmaVz[iRun] = nullptr;
      fHn4DMtZNCCentVxVySigmaVz[iRun] = nullptr;
      fHn4DQxZNACentVxVySigmaVz[iRun] = nullptr;
      fHn4DQyZNACentVxVySigmaVz[iRun] = nullptr;
      fHn4DMtZNACentVxVySigmaVz[iRun] = nullptr;

      fProfile2DShiftCentiCosC[iRun] = nullptr;
      fProfile2DShiftCentiSinC[iRun] = nullptr;
      fProfile2DShiftCentiCosA[iRun] = nullptr;
      fProfile2DShiftCentiSinA[iRun] = nullptr;

      for (size_t i = 0; i < 2; i++)
      {
        fProfileQxCCent[iRun][i] = nullptr;
        fProfileQyCCent[iRun][i] = nullptr;
        fProfileQxACent[iRun][i] = nullptr;
        fProfileQyACent[iRun][i] = nullptr;

        fProfileQxAQxCCent[iRun][i] = nullptr;
        fProfileQxAQyCCent[iRun][i] = nullptr;
        fProfileQyAQxCCent[iRun][i] = nullptr;
        fProfileQyAQyCCent[iRun][i] = nullptr;
      }

      for (size_t i = 0; i < 3; i++)
      {
        fHist2DPsiCCent[iRun][i] = nullptr;
        fHist2DPsiACent[iRun][i] = nullptr;
      }
    }
    
    for (size_t i = 0; i < 3; i++)
    {
      fProfileV2PsiZNCCent[i];
      fProfileV2PsiZNACent[i];      
    }
    
    for (size_t i = 0; i < 10; i++)
    {
      for (size_t j = 0; j < 3; i++)
      {
        fProfileV2PsiZNCVsPt[i][j];
        fProfilePhiPsiZNCCent[i][j];
        fProfileV2PsiZNAVsPt[i][j];
        fProfilePhiPsiZNACent[i][j];
      }
    }
    
}
//_____________________________________________________________________________
AliAnalysisTaskZDCCalibration::AliAnalysisTaskZDCCalibration(const char* name) : AliAnalysisTaskSE(name),
    fDataSet("10h"),
    bFirstFillHistVetex(kTRUE),
    bFillHistForGE(kTRUE),
    bFillHistForRC(kFALSE),
    bFillHistForSF(kFALSE),
    bCalculateV2(kFALSE),
    fFilterBit(1),
    bApplyGE(kFALSE),
    bApplyRC(kFALSE),
    bApplySF(kFALSE),
    bFillQAHist(kFALSE),

    fZDCCalibrationList(nullptr),

    fCentrality(-999),
    fAOD(nullptr),
    fZDC(nullptr),
    fUtils(nullptr),
    fOutputList(nullptr),
    fQAList(nullptr),
    fHistCent(nullptr),
    fHist2DVxVyTot(nullptr),
    fHistVzTot(nullptr),
    
    fHist2DForMeanVxVy(nullptr),
    fHistForMeanVz(nullptr),

    fProfileForZNCGE(nullptr),
    fProfileForZNAGE(nullptr),

    fHn4DForZNCQxRC(nullptr),
    fHn4DForZNCQyRC(nullptr),
    fHn4DForZNCMtRC(nullptr),
    fHn4DForZNAQxRC(nullptr),
    fHn4DForZNAQyRC(nullptr),
    fHn4DForZNAMtRC(nullptr),

    fProfile2DForCosC(nullptr),
    fProfile2DForSinC(nullptr),
    fProfile2DForCosA(nullptr),
    fProfile2DForSinA(nullptr)
{
    for (size_t i = 0; i < 2; i++) {fHistPt[i] = nullptr;}
    for (size_t i = 0; i < 2; i++) {fHistEta[i] = nullptr;}
    for (size_t i = 0; i < 2; i++) {fHistPhi[i] = nullptr;}
    
    for (size_t i = 0; i < 3; i++) {fVtx[i] = -999;}
    for (size_t i = 0; i < 2; i++) {fHist2DCentCorr[i] = nullptr;}
    for (size_t i = 0; i < 2; i++)
    {
      fProfileQxCCentTot[i] = nullptr;
      fProfileQyCCentTot[i] = nullptr;
      fProfileQxACentTot[i] = nullptr;
      fProfileQyACentTot[i] = nullptr;
      fProfileQxAQxCCentTot[i] = nullptr;
      fProfileQxAQyCCentTot[i] = nullptr;
      fProfileQyAQxCCentTot[i] = nullptr;
      fProfileQyAQyCCentTot[i] = nullptr;
    }

    for (size_t i = 0; i < 3; i++)
    {
      fHist2DPsiCCentTot[i] = nullptr;
      fHist2DPsiACentTot[i] = nullptr;
    }

    for (size_t iRun = 0; iRun < fnRunMax; iRun++)
    {
      fRunList[iRun]                  = nullptr;
      fQAListThisRun[iRun]            = nullptr;
      fHist2DVxVy[iRun]               = nullptr;
      fHistVz[iRun]                   = nullptr;
      fProfileZNATowerEnergy[iRun]    = nullptr;
      fProfileZNCTowerEnergy[iRun]    = nullptr;
      fProfileZNATowerEnergyBfGE[iRun]= nullptr;
      fProfileZNCTowerEnergyBfGE[iRun]= nullptr;

      //vx vy sigma vz
      fHn4DQxZNCCentVxVySigmaVz[iRun] = nullptr;
      fHn4DQyZNCCentVxVySigmaVz[iRun] = nullptr;
      fHn4DMtZNCCentVxVySigmaVz[iRun] = nullptr;
      fHn4DQxZNACentVxVySigmaVz[iRun] = nullptr;
      fHn4DQyZNACentVxVySigmaVz[iRun] = nullptr;
      fHn4DMtZNACentVxVySigmaVz[iRun] = nullptr;

      fProfile2DShiftCentiCosC[iRun] = nullptr;
      fProfile2DShiftCentiSinC[iRun] = nullptr;
      fProfile2DShiftCentiCosA[iRun] = nullptr;
      fProfile2DShiftCentiSinA[iRun] = nullptr;

      for (size_t i = 0; i < 2; i++)
      {
        fProfileQxCCent[iRun][i] = nullptr;
        fProfileQyCCent[iRun][i] = nullptr;
        fProfileQxACent[iRun][i] = nullptr;
        fProfileQyACent[iRun][i] = nullptr;

        fProfileQxAQxCCent[iRun][i] = nullptr;
        fProfileQxAQyCCent[iRun][i] = nullptr;
        fProfileQyAQxCCent[iRun][i] = nullptr;
        fProfileQyAQyCCent[iRun][i] = nullptr;
      }

      for (size_t i = 0; i < 3; i++)
      {
        fHist2DPsiCCent[iRun][i] = nullptr;
        fHist2DPsiACent[iRun][i] = nullptr;
      }
    }
    
    for (size_t i = 0; i < 3; i++)
    {
      fProfileV2PsiZNCCent[i];
      fProfileV2PsiZNACent[i];      
    }
    
    for (size_t i = 0; i < 10; i++)
    {
      for (size_t j = 0; j < 3; i++)
      {
        fProfileV2PsiZNCVsPt[i][j];
        fProfilePhiPsiZNCCent[i][j];
        fProfileV2PsiZNAVsPt[i][j];
        fProfilePhiPsiZNACent[i][j];
      }
    }
    
    // constructor
    DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
                                        // this chain is created by the analysis manager, so no need to worry about it, 
                                        // it does its work automatically
    DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms 
                                        // you can add more output objects by calling DefineOutput(2, classname::Class())
                                        // if you add more output objects, make sure to call PostData for all of them, and to
                                        // make changes to your AddTask macro!
}
//_____________________________________________________________________________
AliAnalysisTaskZDCCalibration::~AliAnalysisTaskZDCCalibration()
{
    // destructor
    if(fOutputList) delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    if(fUtils) delete fUtils;  
}
//_____________________________________________________________________________
void AliAnalysisTaskZDCCalibration::UserCreateOutputObjects()
{
    // create output objects
    //
    // this function is called ONCE at the start of your analysis (RUNTIME)
    // here you ceate the histograms that you want to use 
    //
    // the histograms are in this case added to a tlist, this list is in the end saved
    // to an output file
    //
    fOutputList = new TList();          // this is a list which will contain all of your histograms
                                        // at the end of the analysis, the contents of this list are written
                                        // to the output file
    fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them
                                        // if requested (dont worry about this now)
    fHistCent = new TH1D("fHistCent","Centrality", 100, 0., 100.);
    fHist2DCentCorr[0] = new TH2D("fHist2DCentCorrBfCut","Centrality(before Cut) V0M vs. TRK", 100, 0., 100., 100, 0., 100.);
    fHist2DCentCorr[1] = new TH2D("fHist2DCentCorrAfCut","Centrality(after Cut) V0M vs. TRK", 100, 0., 100., 100, 0., 100.);
    fHist2DVxVyTot = new TH2D("fHist2DVxVyTot","Vx vs.Vy (all run) ",1000,-0.5,0.5,1000,-0.5,0.5);
    fHistVzTot = new TH1D("fHistVzTot","Vz (all run)",200,-20,20);

    fOutputList -> Add(fHistCent);
    fOutputList -> Add(fHist2DCentCorr[0]);
    fOutputList -> Add(fHist2DCentCorr[1]);
    fOutputList -> Add(fHist2DVxVyTot);
    fOutputList -> Add(fHistVzTot);

    if (bCalculateV2)
    {
      fHistPt[0]  = new TH1D("fHistPtBfCut","pT before Cut",500,0.,50.);
      fHistEta[0] = new TH1D("fHistEtaBfCut","Eta before Cut",300,-3,3);
      fHistPhi[0] = new TH1D("fHistPhiBfCut","Phi before Cut",360,-TMath::TwoPi(),TMath::TwoPi());
      fHistPt[1]  = new TH1D("fHistPtAffCut","pT after Cut",500,0.,50.);
      fHistEta[1] = new TH1D("fHistEtaAffCut","Eta after Cut",300,-3,3);
      fHistPhi[1] = new TH1D("fHistPhiAffCut","Phi after Cut",360,-TMath::TwoPi(),TMath::TwoPi());
      fOutputList->Add(fHistPt[0]);
      fOutputList->Add(fHistEta[0]);
      fOutputList->Add(fHistPhi[0]);
      fOutputList->Add(fHistPt[1]);
      fOutputList->Add(fHistEta[1]);
      fOutputList->Add(fHistPhi[1]);

      fProfileV2PsiZNCCent[0] = new TProfile("fProfileV2PsiZNCCentAfGE","v2 vs. centrality ZNC After GE",10,0.,100.);
      fProfileV2PsiZNACent[0] = new TProfile("fProfileV2PsiZNACentAfGE","v2 vs. centrality ZNA After GE",10,0.,100.);
      fProfileV2PsiZNCCent[1] = new TProfile("fProfileV2PsiZNCCentAfRC","v2 vs. centrality ZNC After RC",10,0.,100.);
      fProfileV2PsiZNACent[1] = new TProfile("fProfileV2PsiZNACentAfRC","v2 vs. centrality ZNA After RC",10,0.,100.);
      fProfileV2PsiZNCCent[2] = new TProfile("fProfileV2PsiZNCCentAfSF","v2 vs. centrality ZNC After SF",10,0.,100.);
      fProfileV2PsiZNACent[2] = new TProfile("fProfileV2PsiZNACentAfSF","v2 vs. centrality ZNA After SF",10,0.,100.);
      for (size_t iCentBin = 0; iCentBin < 10; iCentBin++)
      {
        fProfileV2PsiZNCVsPt[iCentBin][0] = new TProfile(Form("fProfileV2PsiZNCVsPtAfGECent%d",iCentBin),Form("V2 vs. pT ZNC After GE Cent%d",iCentBin),120,0.2,6.2);
        fProfileV2PsiZNAVsPt[iCentBin][0] = new TProfile(Form("fProfileV2PsiZNAVsPtAfGECent%d",iCentBin),Form("V2 vs. pT ZNA After GE Cent%d",iCentBin),120,0.2,6.2);
        fProfilePhiPsiZNCCent[iCentBin][0] = new TProfile(Form("fProfilePhiPsiZNCAfGECent%d",iCentBin),Form("phi ZNC After GE Cent%d",iCentBin),360,-TMath::TwoPi(),TMath::TwoPi());
        fProfilePhiPsiZNACent[iCentBin][0] = new TProfile(Form("fProfilePhiPsiZNAAfGECent%d",iCentBin),Form("phi ZNA After GE Cent%d",iCentBin),360,-TMath::TwoPi(),TMath::TwoPi());
        fProfileV2PsiZNCVsPt[iCentBin][1] = new TProfile(Form("fProfileV2PsiZNCVsPtAfRCCent%d",iCentBin),Form("V2 vs. pT ZNC After RC Cent%d",iCentBin),120,0.2,6.2);
        fProfileV2PsiZNAVsPt[iCentBin][1] = new TProfile(Form("fProfileV2PsiZNAVsPtAfRCCent%d",iCentBin),Form("V2 vs. pT ZNA After RC Cent%d",iCentBin),120,0.2,6.2);
        fProfilePhiPsiZNCCent[iCentBin][1] = new TProfile(Form("fProfilePhiPsiZNCAfRCCent%d",iCentBin),Form("phi ZNC After RC Cent%d",iCentBin),360,-TMath::TwoPi(),TMath::TwoPi());
        fProfilePhiPsiZNACent[iCentBin][1] = new TProfile(Form("fProfilePhiPsiZNAAfRCCent%d",iCentBin),Form("phi ZNA After RC Cent%d",iCentBin),360,-TMath::TwoPi(),TMath::TwoPi());
        fProfileV2PsiZNCVsPt[iCentBin][2] = new TProfile(Form("fProfileV2PsiZNCVsPtAfSFCent%d",iCentBin),Form("V2 vs. pT ZNC After SF Cent%d",iCentBin),120,0.2,6.2);
        fProfileV2PsiZNAVsPt[iCentBin][2] = new TProfile(Form("fProfileV2PsiZNAVsPtAfSFCent%d",iCentBin),Form("V2 vs. pT ZNA After SF Cent%d",iCentBin),120,0.2,6.2);
        fProfilePhiPsiZNCCent[iCentBin][2] = new TProfile(Form("fProfilePhiPsiZNCAfSFCent%d",iCentBin),Form("phi ZNC After SF Cent%d",iCentBin),360,-TMath::TwoPi(),TMath::TwoPi());
        fProfilePhiPsiZNACent[iCentBin][2] = new TProfile(Form("fProfilePhiPsiZNAAfSFCent%d",iCentBin),Form("phi ZNA After SF Cent%d",iCentBin),360,-TMath::TwoPi(),TMath::TwoPi());
      }
      for (size_t i = 0; i < 3; i++)
      {
        fOutputList->Add(fProfileV2PsiZNCCent[i]);
        fOutputList->Add(fProfileV2PsiZNACent[i]);
        for (size_t iCentBin = 0; iCentBin < 10; iCentBin++)
        {
          fOutputList->Add(fProfileV2PsiZNCVsPt[iCentBin][i]);
          fOutputList->Add(fProfilePhiPsiZNCCent[iCentBin][i]);
          fOutputList->Add(fProfileV2PsiZNAVsPt[iCentBin][i]);
          fOutputList->Add(fProfilePhiPsiZNACent[iCentBin][i]);
        }
      }  
    }
    

    if(bFillQAHist)
    {
      fQAList = new TList();
      fQAList ->SetOwner(kTRUE);
      fQAList ->SetName("QA Total");

      fProfileQxACentTot[0] = new TProfile("fProfileQxACentBfRC","<XA> before Recentering",100,0.,100.);
      fProfileQyACentTot[0] = new TProfile("fProfileQyACentBfRC","<YA> before Recentering",100,0.,100.);
      fProfileQxCCentTot[0] = new TProfile("fProfileQxCCentBfRC","<XC> before Recentering",100,0.,100.);
      fProfileQyCCentTot[0] = new TProfile("fProfileQyCCentBfRC","<YC> before Recentering",100,0.,100.);
      fProfileQxACentTot[1] = new TProfile("fProfileQxACentAfRC","<XA> after Recentering",100,0.,100.);
      fProfileQyACentTot[1] = new TProfile("fProfileQyACentAfRC","<YA> after Recentering",100,0.,100.);
      fProfileQxCCentTot[1] = new TProfile("fProfileQxCCentAfRC","<XC> after Recentering",100,0.,100.);
      fProfileQyCCentTot[1] = new TProfile("fProfileQyCCentAfRC","<YC> after Recentering",100,0.,100.);

      fQAList -> Add(fProfileQxACentTot[0]);
      fQAList -> Add(fProfileQyACentTot[0]);
      fQAList -> Add(fProfileQxCCentTot[0]);
      fQAList -> Add(fProfileQyCCentTot[0]);
      fQAList -> Add(fProfileQxACentTot[1]);
      fQAList -> Add(fProfileQyACentTot[1]);
      fQAList -> Add(fProfileQxCCentTot[1]);
      fQAList -> Add(fProfileQyCCentTot[1]);

      fProfileQxAQxCCentTot[0] = new TProfile("fProfileQxAQxCCentTotBfRC","<XAXC> before Recentering", 100, 0., 100.);
      fProfileQxAQyCCentTot[0] = new TProfile("fProfileQxAQyCCentTotBfRC","<XAYC> before Recentering", 100, 0., 100.);
      fProfileQyAQxCCentTot[0] = new TProfile("fProfileQyAQxCCentTotBfRC","<YAXC> before Recentering", 100, 0., 100.);
      fProfileQyAQyCCentTot[0] = new TProfile("fProfileQyAQyCCentTotBfRC","<YAYC> before Recentering", 100, 0., 100.);
      fProfileQxAQxCCentTot[1] = new TProfile("fProfileQxAQxCCentTotAfRC","<XAXC> after Recentering", 100, 0., 100.);        
      fProfileQxAQyCCentTot[1] = new TProfile("fProfileQxAQyCCentTotAfRC","<XAYC> after Recentering", 100, 0., 100.);        
      fProfileQyAQxCCentTot[1] = new TProfile("fProfileQyAQxCCentTotAfRC","<YAXC> after Recentering", 100, 0., 100.);        
      fProfileQyAQyCCentTot[1] = new TProfile("fProfileQyAQyCCentTotAfRC","<YAYC> after Recentering", 100, 0., 100.);    

      fQAList -> Add(fProfileQxAQxCCentTot[0]);
      fQAList -> Add(fProfileQxAQyCCentTot[0]);
      fQAList -> Add(fProfileQyAQxCCentTot[0]);
      fQAList -> Add(fProfileQyAQyCCentTot[0]);
      fQAList -> Add(fProfileQxAQxCCentTot[1]);
      fQAList -> Add(fProfileQxAQyCCentTot[1]);
      fQAList -> Add(fProfileQyAQxCCentTot[1]);
      fQAList -> Add(fProfileQyAQyCCentTot[1]);

      fHist2DPsiCCentTot[0] = new TH2D("fHist2DPsiCCentGE","Hist2D Cent vs. PsiC after Gain Equalization",10,0.,100.,180,0.,TMath::TwoPi());
      fHist2DPsiACentTot[0] = new TH2D("fHist2DPsiACentGE","Hist2D Cent vs. PsiA after Gain Equalization",10,0.,100.,180,0.,TMath::TwoPi());
      fHist2DPsiCCentTot[1] = new TH2D("fHist2DPsiCCentRC","Hist2D Cent vs. PsiC after Recentering",10,0.,100.,180,0.,TMath::TwoPi());
      fHist2DPsiACentTot[1] = new TH2D("fHist2DPsiACentRC","Hist2D Cent vs. PsiA after Recentering",10,0.,100.,180,0.,TMath::TwoPi());
      fHist2DPsiCCentTot[2] = new TH2D("fHist2DPsiCCentSF","Hist2D Cent vs. PsiC after Shifting",10,0.,100.,180,0.,TMath::TwoPi());
      fHist2DPsiACentTot[2] = new TH2D("fHist2DPsiACentSF","Hist2D Cent vs. PsiA after Shifting",10,0.,100.,180,0.,TMath::TwoPi());

      fQAList -> Add(fHist2DPsiCCentTot[0]);
      fQAList -> Add(fHist2DPsiACentTot[0]);
      fQAList -> Add(fHist2DPsiCCentTot[1]);
      fQAList -> Add(fHist2DPsiACentTot[1]);
      fQAList -> Add(fHist2DPsiCCentTot[2]);
      fQAList -> Add(fHist2DPsiACentTot[2]);

      fOutputList->Add(fQAList);
    }

    int runNumList10h[91]={
      139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328, 139314, 139310,
      139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029, 139028, 138872, 138871, 138870,
      138837, 138732, 138730, 138666, 138662, 138653, 138652, 138638, 138624, 138621, 138583, 138582,
      138579, 138578, 138534, 138469, 138442, 138439, 138438, 138396, 138364, 138275, 138225, 138201,
      138197, 138192, 138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693,
      137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595, 137549, 137546, 137544, 137541,
      137539, 137531, 137530, 137443, 137441, 137440, 137439, 137434, 137432, 137431, 137430, 137243,
      137236, 137235, 137232, 137231, 137230, 137162, 137161
    };  
    //11h
    int runNumList11h[39]={
      170387 ,170040, 170268, 170228, 170207, 169838, 170159, 170204, 170311, 170084,
      169835 ,170088, 170593, 170203, 170270, 169846, 170163, 170388, 170155, 170083,
      170572 ,169837, 169855, 170306, 170269, 170089, 170309, 170091, 170081, 170230,
      170085 ,170315, 170027, 170193, 170312, 170313, 170308, 169858, 169859
    };
     
    int fAllRunsNum = -999;
    if (fDataSet.Contains("10h")) fAllRunsNum = 91;
    if (fDataSet.Contains("11h")) fAllRunsNum = 39;
    
    for (int iRun = 0; iRun < fAllRunsNum; iRun++)
    {
      fRunList[iRun] = new TList();
      fRunList[iRun] -> SetOwner(kTRUE);
      if (fDataSet.Contains("10h")) fRunList[iRun] -> SetName(Form("%i", runNumList10h[iRun]));
      if (fDataSet.Contains("11h")) fRunList[iRun] -> SetName(Form("%i", runNumList11h[iRun]));

      fHist2DVxVy[iRun] = new TH2D("hist2DVxVy","hist2DVxVy",1000,-0.5,0.5,1000,-0.5,0.5);
      fHistVz[iRun] = new TH1D("histVz","histVz",200,-20,20);
      fRunList[iRun] -> Add(fHist2DVxVy[iRun]);
      fRunList[iRun] -> Add(fHistVz[iRun]);

      if (bFillHistForGE)
      {
        fProfileZNCTowerEnergy[iRun] = new TProfile("profileZNCTowerEnergy","ProfileZNCTowerEnergy",5,0,5);
        fProfileZNATowerEnergy[iRun] = new TProfile("profileZNATowerEnergy","ProfileZNATowerEnergy",5,0,5);
        fRunList[iRun] -> Add(fProfileZNCTowerEnergy[iRun]);
        fRunList[iRun] -> Add(fProfileZNATowerEnergy[iRun]);
      }

      if (bApplyGE)
      {
        fProfileForZNAGE = new TProfile();
        fProfileForZNAGE = new TProfile();
        fProfileZNCTowerEnergyBfGE[iRun] = new TProfile("profileZNCTowerEnergyAfGE","ProfileZNCTowerEnergyAfGE",5,0,5);
        fProfileZNATowerEnergyBfGE[iRun] = new TProfile("profileZNATowerEnergyAfGE","ProfileZNATowerEnergyAfGE",5,0,5);
        fRunList[iRun] -> Add(fProfileZNCTowerEnergyBfGE[iRun]);
        fRunList[iRun] -> Add(fProfileZNATowerEnergyBfGE[iRun]);
      }
      
      if (bFillHistForRC)
      {
        int   nbins[4] = {100,  3,  3,  3};
        double xmin[4] = {  0,  0,  0,  0};
        double xmax[4] = {100,  3,  3,  3};

        fHn4DQxZNCCentVxVySigmaVz[iRun] = new THnSparseD("fHn4DQxZNCCentVxVySigmaVz","fHn4DQxZNCCentVxVySigmaVz",4,nbins,xmin,xmax);
        fHn4DQyZNCCentVxVySigmaVz[iRun] = new THnSparseD("fHn4DQyZNCCentVxVySigmaVz","fHn4DQyZNCCentVxVySigmaVz",4,nbins,xmin,xmax);
        fHn4DMtZNCCentVxVySigmaVz[iRun] = new THnSparseD("fHn4DMtZNCCentVxVySigmaVz","fHn4DMtZNCCentVxVySigmaVz",4,nbins,xmin,xmax);
        fHn4DQxZNACentVxVySigmaVz[iRun] = new THnSparseD("fHn4DQxZNACentVxVySigmaVz","fHn4DQxZNACentVxVySigmaVz",4,nbins,xmin,xmax);
        fHn4DQyZNACentVxVySigmaVz[iRun] = new THnSparseD("fHn4DQyZNACentVxVySigmaVz","fHn4DQyZNACentVxVySigmaVz",4,nbins,xmin,xmax);
        fHn4DMtZNACentVxVySigmaVz[iRun] = new THnSparseD("fHn4DMtZNACentVxVySigmaVz","fHn4DMtZNACentVxVySigmaVz",4,nbins,xmin,xmax);
        
        fRunList[iRun] -> Add(fHn4DQxZNCCentVxVySigmaVz[iRun]);
        fRunList[iRun] -> Add(fHn4DQyZNCCentVxVySigmaVz[iRun]);
        fRunList[iRun] -> Add(fHn4DMtZNCCentVxVySigmaVz[iRun]);
        fRunList[iRun] -> Add(fHn4DQxZNACentVxVySigmaVz[iRun]);
        fRunList[iRun] -> Add(fHn4DQyZNACentVxVySigmaVz[iRun]);
        fRunList[iRun] -> Add(fHn4DMtZNACentVxVySigmaVz[iRun]);
      }
      
      if (bApplyRC)
      {
        fHn4DForZNAQxRC = new THnSparseD();
        fHn4DForZNAQyRC = new THnSparseD();
        fHn4DForZNAMtRC = new THnSparseD();
        fHn4DForZNCQxRC = new THnSparseD();
        fHn4DForZNCQyRC = new THnSparseD();
        fHn4DForZNCMtRC = new THnSparseD();
      }

      if (bFillHistForSF)
      {
        fProfile2DShiftCentiCosC[iRun] = new TProfile2D("fProfile2DShiftCentiCosC","fProfile2DShiftCentiCosC",10,0.,100.,20,0.,20.);
        fProfile2DShiftCentiSinC[iRun] = new TProfile2D("fProfile2DShiftCentiSinC","fProfile2DShiftCentiSinC",10,0.,100.,20,0.,20.);
        fProfile2DShiftCentiCosA[iRun] = new TProfile2D("fProfile2DShiftCentiCosA","fProfile2DShiftCentiCosA",10,0.,100.,20,0.,20.);
        fProfile2DShiftCentiSinA[iRun] = new TProfile2D("fProfile2DShiftCentiSinA","fProfile2DShiftCentiSinA",10,0.,100.,20,0.,20.);
        fRunList[iRun] -> Add(fProfile2DShiftCentiCosC[iRun]);
        fRunList[iRun] -> Add(fProfile2DShiftCentiSinC[iRun]);
        fRunList[iRun] -> Add(fProfile2DShiftCentiCosA[iRun]);
        fRunList[iRun] -> Add(fProfile2DShiftCentiSinA[iRun]);
      }

      if (bApplySF)
      {
        fProfile2DForCosC = new TProfile2D();
        fProfile2DForSinC = new TProfile2D();
        fProfile2DForCosA = new TProfile2D();
        fProfile2DForSinA = new TProfile2D();
      }
      
      
      if (bFillQAHist)
      {
        fQAListThisRun[iRun] = new TList();
        fQAListThisRun[iRun] ->SetOwner();
        fQAListThisRun[iRun] ->SetName("QA");

        if (fDataSet.Contains("10h")) 
        {
          fProfileQxACent[iRun][0] = new TProfile("fProfileQxACentBfRC",Form("<XA> before Recentering Run%d",runNumList10h[iRun]),100,0.,100.);
          fProfileQyACent[iRun][0] = new TProfile("fProfileQyACentBfRC",Form("<YA> before Recentering Run%d",runNumList10h[iRun]),100,0.,100.);
          fProfileQxCCent[iRun][0] = new TProfile("fProfileQxCCentBfRC",Form("<XC> before Recentering Run%d",runNumList10h[iRun]),100,0.,100.);
          fProfileQyCCent[iRun][0] = new TProfile("fProfileQyCCentBfRC",Form("<YC> before Recentering Run%d",runNumList10h[iRun]),100,0.,100.);
          fProfileQxACent[iRun][1] = new TProfile("fProfileQxACentAfRC",Form("<XA> after Recentering Run%d",runNumList10h[iRun]),100,0.,100.);
          fProfileQyACent[iRun][1] = new TProfile("fProfileQyACentAfRC",Form("<YA> after Recentering Run%d",runNumList10h[iRun]),100,0.,100.);
          fProfileQxCCent[iRun][1] = new TProfile("fProfileQxCCentAfRC",Form("<XC> after Recentering Run%d",runNumList10h[iRun]),100,0.,100.);
          fProfileQyCCent[iRun][1] = new TProfile("fProfileQyCCentAfRC",Form("<YC> after Recentering Run%d",runNumList10h[iRun]),100,0.,100.);
        }

        if (fDataSet.Contains("11h")) 
        {
          fProfileQxACent[iRun][0] = new TProfile("fProfileQxACentBfRC",Form("<XA> before Recentering Run%d",runNumList11h[iRun]),100,0.,100.);
          fProfileQyACent[iRun][0] = new TProfile("fProfileQyACentBfRC",Form("<YA> before Recentering Run%d",runNumList11h[iRun]),100,0.,100.);
          fProfileQxCCent[iRun][0] = new TProfile("fProfileQxCCentBfRC",Form("<XC> before Recentering Run%d",runNumList11h[iRun]),100,0.,100.);
          fProfileQyCCent[iRun][0] = new TProfile("fProfileQyCCentBfRC",Form("<YC> before Recentering Run%d",runNumList11h[iRun]),100,0.,100.);
          fProfileQxACent[iRun][1] = new TProfile("fProfileQxACentAfRC",Form("<XA> after Recentering Run%d",runNumList11h[iRun]),100,0.,100.);
          fProfileQyACent[iRun][1] = new TProfile("fProfileQyACentAfRC",Form("<YA> after Recentering Run%d",runNumList11h[iRun]),100,0.,100.);
          fProfileQxCCent[iRun][1] = new TProfile("fProfileQxCCentAfRC",Form("<XC> after Recentering Run%d",runNumList11h[iRun]),100,0.,100.);
          fProfileQyCCent[iRun][1] = new TProfile("fProfileQyCCentAfRC",Form("<YC> after Recentering Run%d",runNumList11h[iRun]),100,0.,100.);
        }

        fQAListThisRun[iRun]->Add(fProfileQxACent[iRun][0]);
        fQAListThisRun[iRun]->Add(fProfileQyACent[iRun][0]);
        fQAListThisRun[iRun]->Add(fProfileQxCCent[iRun][0]);
        fQAListThisRun[iRun]->Add(fProfileQyCCent[iRun][0]);
        fQAListThisRun[iRun]->Add(fProfileQxACent[iRun][1]);
        fQAListThisRun[iRun]->Add(fProfileQyACent[iRun][1]);
        fQAListThisRun[iRun]->Add(fProfileQxCCent[iRun][1]);
        fQAListThisRun[iRun]->Add(fProfileQyCCent[iRun][1]);


        if (fDataSet.Contains("10h")) 
        {
          fProfileQxAQxCCent[iRun][0] = new TProfile("fProfileQxAQxCCentBfRCThisRun",Form("<XAXC> before Recentering Run%d", runNumList10h[iRun]), 100, 0., 100.);
          fProfileQxAQyCCent[iRun][0] = new TProfile("fProfileQxAQyCCentBfRCThisRun",Form("<XAYC> before Recentering Run%d", runNumList10h[iRun]), 100, 0., 100.); 
          fProfileQyAQxCCent[iRun][0] = new TProfile("fProfileQyAQxCCentBfRCThisRun",Form("<YAXC> before Recentering Run%d", runNumList10h[iRun]), 100, 0., 100.); 
          fProfileQyAQyCCent[iRun][0] = new TProfile("fProfileQyAQyCCentBfRCThisRun",Form("<YAYC> before Recentering Run%d", runNumList10h[iRun]), 100, 0., 100.); 
          fProfileQxAQxCCent[iRun][1] = new TProfile("fProfileQxAQxCCentAfRCThisRun",Form("<XAXC> after Recentering Run%d", runNumList10h[iRun]), 100, 0., 100.);
          fProfileQxAQyCCent[iRun][1] = new TProfile("fProfileQxAQyCCentAfRCThisRun",Form("<XAYC> after Recentering Run%d", runNumList10h[iRun]), 100, 0., 100.); 
          fProfileQyAQxCCent[iRun][1] = new TProfile("fProfileQyAQxCCentAfRCThisRun",Form("<YAXC> after Recentering Run%d", runNumList10h[iRun]), 100, 0., 100.); 
          fProfileQyAQyCCent[iRun][1] = new TProfile("fProfileQyAQyCCentAfRCThisRun",Form("<YAYC> after Recentering Run%d", runNumList10h[iRun]), 100, 0., 100.); 
        }

        if (fDataSet.Contains("11h")) 
        {
          fProfileQxAQxCCent[iRun][0] = new TProfile("fProfileQxAQxCCentBfRCThisRun",Form("<XAXC> before Recentering Run%d", runNumList11h[iRun]), 100, 0., 100.);
          fProfileQxAQyCCent[iRun][0] = new TProfile("fProfileQxAQyCCentBfRCThisRun",Form("<XAYC> before Recentering Run%d", runNumList11h[iRun]), 100, 0., 100.); 
          fProfileQyAQxCCent[iRun][0] = new TProfile("fProfileQyAQxCCentBfRCThisRun",Form("<YAXC> before Recentering Run%d", runNumList11h[iRun]), 100, 0., 100.); 
          fProfileQyAQyCCent[iRun][0] = new TProfile("fProfileQyAQyCCentBfRCThisRun",Form("<YAYC> before Recentering Run%d", runNumList11h[iRun]), 100, 0., 100.); 
          fProfileQxAQxCCent[iRun][1] = new TProfile("fProfileQxAQxCCentAfRCThisRun",Form("<XAXC> after Recentering Run%d", runNumList11h[iRun]), 100, 0., 100.);
          fProfileQxAQyCCent[iRun][1] = new TProfile("fProfileQxAQyCCentAfRCThisRun",Form("<XAYC> after Recentering Run%d", runNumList11h[iRun]), 100, 0., 100.); 
          fProfileQyAQxCCent[iRun][1] = new TProfile("fProfileQyAQxCCentAfRCThisRun",Form("<YAXC> after Recentering Run%d", runNumList11h[iRun]), 100, 0., 100.); 
          fProfileQyAQyCCent[iRun][1] = new TProfile("fProfileQyAQyCCentAfRCThisRun",Form("<YAYC> after Recentering Run%d", runNumList11h[iRun]), 100, 0., 100.); 
        }

        fQAListThisRun[iRun]->Add(fProfileQxAQxCCent[iRun][0]);
        fQAListThisRun[iRun]->Add(fProfileQxAQyCCent[iRun][0]);
        fQAListThisRun[iRun]->Add(fProfileQyAQxCCent[iRun][0]);
        fQAListThisRun[iRun]->Add(fProfileQyAQyCCent[iRun][0]);
        fQAListThisRun[iRun]->Add(fProfileQxAQxCCent[iRun][1]);
        fQAListThisRun[iRun]->Add(fProfileQxAQyCCent[iRun][1]);
        fQAListThisRun[iRun]->Add(fProfileQyAQxCCent[iRun][1]);
        fQAListThisRun[iRun]->Add(fProfileQyAQyCCent[iRun][1]);

        if (fDataSet.Contains("10h")) 
        {
          fHist2DPsiCCent[iRun][0] = new TH2D("fHist2DPsiCCentGE",Form("Hist2D Cent vs. PsiC after Gain Equalization Run%d",runNumList10h[iRun]),10,0.,100.,50,0.,TMath::TwoPi());
          fHist2DPsiACent[iRun][0] = new TH2D("fHist2DPsiACentGE",Form("Hist2D Cent vs. PsiA after Gain Equalization Run%d",runNumList10h[iRun]),10,0.,100.,50,0.,TMath::TwoPi());
          fHist2DPsiCCent[iRun][1] = new TH2D("fHist2DPsiCCentRC",Form("Hist2D Cent vs. PsiC after Recentering Run%d",runNumList10h[iRun]),10,0.,100.,50,0.,TMath::TwoPi());
          fHist2DPsiACent[iRun][1] = new TH2D("fHist2DPsiACentRC",Form("Hist2D Cent vs. PsiA after Recentering Run%d",runNumList10h[iRun]),10,0.,100.,50,0.,TMath::TwoPi());
          fHist2DPsiCCent[iRun][2] = new TH2D("fHist2DPsiCCentSF",Form("Hist2D Cent vs. PsiC after Shifting Run%d",runNumList10h[iRun]),10,0.,100.,50,0.,TMath::TwoPi());
          fHist2DPsiACent[iRun][2] = new TH2D("fHist2DPsiACentSF",Form("Hist2D Cent vs. PsiA after Shifting Run%d",runNumList10h[iRun]),10,0.,100.,50,0.,TMath::TwoPi());
        }
        
        if (fDataSet.Contains("11h")) 
        {
          fHist2DPsiCCent[iRun][0] = new TH2D("fHist2DPsiCCentGE",Form("Hist2D Cent vs. PsiC after Gain Equalization Run%d",runNumList11h[iRun]),10,0.,100.,50,0.,TMath::TwoPi());
          fHist2DPsiACent[iRun][0] = new TH2D("fHist2DPsiACentGE",Form("Hist2D Cent vs. PsiA after Gain Equalization Run%d",runNumList11h[iRun]),10,0.,100.,50,0.,TMath::TwoPi());
          fHist2DPsiCCent[iRun][1] = new TH2D("fHist2DPsiCCentRC",Form("Hist2D Cent vs. PsiC after Recentering Run%d",runNumList11h[iRun]),10,0.,100.,50,0.,TMath::TwoPi());
          fHist2DPsiACent[iRun][1] = new TH2D("fHist2DPsiACentRC",Form("Hist2D Cent vs. PsiA after Recentering Run%d",runNumList11h[iRun]),10,0.,100.,50,0.,TMath::TwoPi());
          fHist2DPsiCCent[iRun][2] = new TH2D("fHist2DPsiCCentSF",Form("Hist2D Cent vs. PsiC after Shifting Run%d",runNumList11h[iRun]),10,0.,100.,50,0.,TMath::TwoPi());
          fHist2DPsiACent[iRun][2] = new TH2D("fHist2DPsiACentSF",Form("Hist2D Cent vs. PsiA after Shifting Run%d",runNumList11h[iRun]),10,0.,100.,50,0.,TMath::TwoPi());
        }

        fQAListThisRun[iRun]->Add(fHist2DPsiCCent[iRun][0]);
        fQAListThisRun[iRun]->Add(fHist2DPsiACent[iRun][0]);
        fQAListThisRun[iRun]->Add(fHist2DPsiCCent[iRun][1]);
        fQAListThisRun[iRun]->Add(fHist2DPsiACent[iRun][1]);
        fQAListThisRun[iRun]->Add(fHist2DPsiCCent[iRun][2]);
        fQAListThisRun[iRun]->Add(fHist2DPsiACent[iRun][2]);

        fRunList[iRun]->Add(fQAListThisRun[iRun]);
      }

      fOutputList -> Add(fRunList[iRun]);
    }
    

    PostData(1, fOutputList);
}
//_____________________________________________________________________________
void AliAnalysisTaskZDCCalibration::UserExec(Option_t *)
{
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());    
    if(!fAOD) return;
    fZDC = fAOD -> GetZDCData();
    if(!fZDC) return;
    fUtils = new AliAnalysisUtils();
    if(!fUtils) return;

    // run Number
    int runNum = fAOD->GetRunNumber();
    int runNumBin = GetRunNumBin(runNum);
    if (runNumBin == -1) return; 

    //pile Up
    if (fDataSet.Contains("10h"))
    {
      fUtils->SetUseOutOfBunchPileUp(true);
      fUtils->SetUseMVPlpSelection(true);
      bool isPileup = fUtils->IsPileUpEvent(fAOD);
      if (isPileup) return;
    }
    if (fDataSet.Contains("11h"))
    {
      bool isPileup = fAOD->IsPileupFromSPD();
      if (isPileup) return;
    }
    double centV0M = fAOD->GetCentrality()->GetCentralityPercentile("V0M");
    double centCL1 = fAOD->GetCentrality()->GetCentralityPercentile("CL1");
    double centTRK = fAOD->GetCentrality()->GetCentralityPercentile("TRK");
    fHist2DCentCorr[0]->Fill(centV0M,centCL1);
    if (fabs(centV0M-centCL1)>7.5) return;
    if (centV0M<0 || centV0M>=100) return;
    fHist2DCentCorr[1]->Fill(centV0M,centCL1);
    fCentrality = centV0M;
    fHistCent->Fill(centV0M);

    // vertex
    fAOD -> GetPrimaryVertex() -> GetXYZ(fVtx);
    double vzSPD  = fAOD->GetPrimaryVertexSPD()->GetZ();
    if (fabs(fVtx[2])>10) return;
    if (fabs(fVtx[2]-vzSPD)>0.5) return;
    if (fabs(fVtx[0])<1e-6 || fabs(fVtx[1])<1e-6 || fabs(fVtx[2])<1e-6) return;
    fHist2DVxVy[runNumBin] ->Fill(fVtx[0],fVtx[1]);
    fHistVz[runNumBin] ->Fill(fVtx[2]);
    fHist2DVxVyTot->Fill(fVtx[0],fVtx[1]);
    fHistVzTot->Fill(fVtx[2]);

    TList* fZDCCalibrationListThisRun = nullptr;
    if(!bFirstFillHistVetex) fZDCCalibrationListThisRun = (TList*)fZDCCalibrationList->FindObject(Form("%d",runNum));

    //Plane
    Double_t PsiZNAGE = -999;
    Double_t PsiZNCGE = -999;
    Double_t PsiZNARC = -999;
    Double_t PsiZNCRC = -999;
    Double_t PsiZNASF = -999;
    Double_t PsiZNCSF = -999;

    const double x[4] = {-1.75, 1.75, -1.75, 1.75};
    const double y[4] = {-1.75, -1.75, 1.75, 1.75};
    const double* EZNCRaw = fZDC->GetZNCTowerEnergy();
    const double* EZNARaw = fZDC->GetZNATowerEnergy();

    if(bFillHistForGE) {
      for (int iTower = 0; iTower < 5; ++iTower) {
        fProfileZNCTowerEnergy[runNumBin] -> Fill(iTower+0.5, EZNCRaw[iTower]);
        fProfileZNATowerEnergy[runNumBin] -> Fill(iTower+0.5, EZNARaw[iTower]);
      }
    }

    double EZNC[5] = {0.,0.,0.,0.,0.}; 
    double EZNA[5] = {0.,0.,0.,0.,0.};
    double QxC = 0.;
    double QyC = 0.;
    double MC  = 0.;
    double QxA = 0.;
    double QyA = 0.;
    double MA  = 0.;

    if(bApplyGE) {
      fProfileForZNCGE = (TProfile*)fZDCCalibrationListThisRun->FindObject("profileZNCTowerEnergy"); 
      fProfileForZNAGE = (TProfile*)fZDCCalibrationListThisRun->FindObject("profileZNATowerEnergy"); 
      for (int iTower = 0; iTower < 5; ++iTower) {
        EZNC[iTower] = EZNCRaw[iTower] / (fProfileForZNCGE->GetBinContent(iTower + 1)) * (fProfileForZNCGE->GetBinContent(2));//here to select the ref Tower
        EZNA[iTower] = EZNARaw[iTower] / (fProfileForZNAGE->GetBinContent(iTower + 1)) * (fProfileForZNAGE->GetBinContent(2));
      }

      for (int iTower = 0; iTower < 4; ++iTower) {
        QxC += EZNC[iTower + 1] * x[iTower];
        QyC += EZNC[iTower + 1] * y[iTower];
        MC  += EZNC[iTower + 1];
        QxA += EZNA[iTower + 1] * x[iTower];
        QyA += EZNA[iTower + 1] * y[iTower];
        MA  += EZNA[iTower + 1];
      }
      if(MC < 1.e-6) return;
      if(MA < 1.e-6) return;
      //QA
      for (int iTower = 0; iTower < 5; ++iTower) {
        fProfileZNCTowerEnergyBfGE[runNumBin] -> Fill(iTower+0.5, EZNC[iTower]);
        fProfileZNATowerEnergyBfGE[runNumBin] -> Fill(iTower+0.5, EZNA[iTower]);
      }
    }

    //GetVexBin
    int vxBin = -1;
    int vyBin = -1;
    int vzBin = -1;
    double fillPosition[4] = {fCentrality,-999.,-999.,-999.};

    if (bFillHistForRC || bApplyRC)
    {
      fHist2DForMeanVxVy = (TH2D*)fZDCCalibrationListThisRun->FindObject("hist2DVxVy"); 
      fHistForMeanVz     = (TH1D*)fZDCCalibrationListThisRun->FindObject("histVz");
      double vxMean = fHist2DForMeanVxVy->ProjectionX()->GetMean();
      double vyMean = fHist2DForMeanVxVy->ProjectionY()->GetMean();
      double vxSigmaMean = fHist2DForMeanVxVy->ProjectionX()->GetRMS();
      double vySigmaMean = fHist2DForMeanVxVy->ProjectionY()->GetRMS();
      if(vxSigmaMean < 1.e-6 || vySigmaMean < 1.e-6) return;

      double vxNSigma = (fVtx[0] - vxMean)/vxSigmaMean;
      double vyNSigma = (fVtx[1] - vyMean)/vySigmaMean;

      if     (vxNSigma>-3.   && vxNSigma<=-0.42)  {vxBin=1;}
      else if(vxNSigma>-0.42 && vxNSigma<= 0.42)  {vxBin=2;}
      else if(vxNSigma> 0.42 && vxNSigma<= 3.  )  {vxBin=3;}

      if     (vyNSigma>-3.   && vyNSigma<=-0.42)  {vyBin=1;}
      else if(vyNSigma>-0.42 && vyNSigma<= 0.42)  {vyBin=2;}
      else if(vyNSigma> 0.42 && vyNSigma<= 3.  )  {vyBin=3;}
     
      if     (fVtx[2]> -10.  && fVtx[2] <= -2.5)  {vzBin=1;}
      else if(fVtx[2]> -2.5  && fVtx[2] <=  2.5)  {vzBin=2;}
      else if(fVtx[2]>  2.5  && fVtx[2] <   10.)  {vzBin=3;}

      if(vxBin <= 0 && vyBin <= 0 && vzBin <= 0) return;
      
      fillPosition[1] = vxBin-0.5;
      fillPosition[2] = vyBin-0.5;
      fillPosition[3] = vzBin-0.5;
    }

    if(bFillHistForRC) {
      fHn4DQxZNCCentVxVySigmaVz[runNumBin] -> Fill(fillPosition,QxC);
      fHn4DQyZNCCentVxVySigmaVz[runNumBin] -> Fill(fillPosition,QyC);
      fHn4DMtZNCCentVxVySigmaVz[runNumBin] -> Fill(fillPosition,MC);

      fHn4DQxZNACentVxVySigmaVz[runNumBin] -> Fill(fillPosition,QxA);
      fHn4DQyZNACentVxVySigmaVz[runNumBin] -> Fill(fillPosition,QyA);
      fHn4DMtZNACentVxVySigmaVz[runNumBin] -> Fill(fillPosition,MA);
    }

    QxC /= MC;
    QyC /= MC;
    QxA /= MA;
    QyA /= MA;

    if(bFillQAHist && bApplyGE) {
      fProfileQxCCent[runNumBin][0]->Fill(fCentrality,QxC);
      fProfileQyCCent[runNumBin][0]->Fill(fCentrality,QyC);
      fProfileQxACent[runNumBin][0]->Fill(fCentrality,-QxA);
      fProfileQyACent[runNumBin][0]->Fill(fCentrality,QyA);
      fProfileQxAQxCCent[runNumBin][0]->Fill(fCentrality,-QxA*QxC);
      fProfileQxAQyCCent[runNumBin][0]->Fill(fCentrality,-QxA*QyC);
      fProfileQyAQxCCent[runNumBin][0]->Fill(fCentrality,QyA*QxC);
      fProfileQyAQyCCent[runNumBin][0]->Fill(fCentrality,QyA*QyC);

      fProfileQxCCentTot[0]->Fill(fCentrality,QxC);
      fProfileQyCCentTot[0]->Fill(fCentrality,QyC);
      fProfileQxACentTot[0]->Fill(fCentrality,-QxA);
      fProfileQyACentTot[0]->Fill(fCentrality,QyA);
      fProfileQxAQxCCentTot[0]->Fill(fCentrality,-QxA*QxC);
      fProfileQxAQyCCentTot[0]->Fill(fCentrality,-QxA*QyC);
      fProfileQyAQxCCentTot[0]->Fill(fCentrality,QyA*QxC);
      fProfileQyAQyCCentTot[0]->Fill(fCentrality,QyA*QyC);

      double psiC = atan2(QyC,QxC)>0. ? atan2(QyC,QxC) : atan2(QyC,QxC)+TMath::TwoPi();
      double psiA = atan2(QyA,-QxA)>0. ? atan2(QyA,-QxA) : atan2(QyA,-QxA)+TMath::TwoPi();
      PsiZNCGE = psiC;
      PsiZNAGE = psiA;

      fHist2DPsiCCent[runNumBin][0]->Fill(fCentrality,psiC);
      fHist2DPsiACent[runNumBin][0]->Fill(fCentrality,psiA);

      fHist2DPsiCCentTot[0]->Fill(fCentrality,psiC);
      fHist2DPsiACentTot[0]->Fill(fCentrality,psiA);
    }

    double QxCMean = 0.;
    double QyCMean = 0.;
    double MCMean  = 0.;
    double QxAMean = 0.;
    double QyAMean = 0.;
    double MAMean  = 0.;

    if(bApplyRC) {
      fHn4DForZNCQxRC = (THnSparseD*)fZDCCalibrationListThisRun->FindObject("fHn4DQxZNCCentVxVySigmaVz");
      fHn4DForZNCQyRC = (THnSparseD*)fZDCCalibrationListThisRun->FindObject("fHn4DQyZNCCentVxVySigmaVz");
      fHn4DForZNCMtRC = (THnSparseD*)fZDCCalibrationListThisRun->FindObject("fHn4DMtZNCCentVxVySigmaVz");
      fHn4DForZNAQxRC = (THnSparseD*)fZDCCalibrationListThisRun->FindObject("fHn4DQxZNACentVxVySigmaVz");
      fHn4DForZNAQyRC = (THnSparseD*)fZDCCalibrationListThisRun->FindObject("fHn4DQyZNACentVxVySigmaVz");
      fHn4DForZNAMtRC = (THnSparseD*)fZDCCalibrationListThisRun->FindObject("fHn4DMtZNACentVxVySigmaVz");

      QxCMean = fHn4DForZNCQxRC -> GetBinContent(fHn4DForZNCQxRC->GetBin(fillPosition));
      QyCMean = fHn4DForZNCQyRC -> GetBinContent(fHn4DForZNCQyRC->GetBin(fillPosition));
      MCMean  = fHn4DForZNCMtRC -> GetBinContent(fHn4DForZNCMtRC->GetBin(fillPosition));
      QxAMean = fHn4DForZNAQxRC -> GetBinContent(fHn4DForZNAQxRC->GetBin(fillPosition));
      QyAMean = fHn4DForZNAQyRC -> GetBinContent(fHn4DForZNAQyRC->GetBin(fillPosition));
      MAMean  = fHn4DForZNAMtRC -> GetBinContent(fHn4DForZNAMtRC->GetBin(fillPosition));

      if(MCMean < 1.e-6) return;
      if(MAMean < 1.e-6) return;
      
      QxCMean /= MCMean; 
      QyCMean /= MCMean; 
      QxAMean /= MAMean; 
      QyAMean /= MAMean;

      QxC -= QxCMean; 
      QyC -= QyCMean; 
      QxA -= QxAMean; 
      QyA -= QyAMean; 
    }

    if(bFillQAHist && bApplyRC) {
      fProfileQxCCent[runNumBin][1]->Fill(fCentrality,QxC);
      fProfileQyCCent[runNumBin][1]->Fill(fCentrality,QyC);
      fProfileQxACent[runNumBin][1]->Fill(fCentrality,-QxA);
      fProfileQyACent[runNumBin][1]->Fill(fCentrality,QyA);
      fProfileQxAQxCCent[runNumBin][1]->Fill(fCentrality,-QxA*QxC);
      fProfileQxAQyCCent[runNumBin][1]->Fill(fCentrality,-QxA*QyC);
      fProfileQyAQxCCent[runNumBin][1]->Fill(fCentrality,QyA*QxC);
      fProfileQyAQyCCent[runNumBin][1]->Fill(fCentrality,QyA*QyC);


      fProfileQxCCentTot[1]->Fill(fCentrality,QxC);
      fProfileQyCCentTot[1]->Fill(fCentrality,QyC);
      fProfileQxACentTot[1]->Fill(fCentrality,-QxA);
      fProfileQyACentTot[1]->Fill(fCentrality,QyA);
      fProfileQxAQxCCentTot[1]->Fill(fCentrality,-QxA*QxC);
      fProfileQxAQyCCentTot[1]->Fill(fCentrality,-QxA*QyC);
      fProfileQyAQxCCentTot[1]->Fill(fCentrality,QyA*QxC);
      fProfileQyAQyCCentTot[1]->Fill(fCentrality,QyA*QyC);

      double psiC = atan2(QyC,QxC)>0. ? atan2(QyC,QxC) : atan2(QyC,QxC)+TMath::TwoPi();
      double psiA = atan2(QyA,-QxA)>0. ? atan2(QyA,-QxA) : atan2(QyA,-QxA)+TMath::TwoPi();
      PsiZNCRC = psiC;
      PsiZNARC = psiA;

      fHist2DPsiCCent[runNumBin][1]->Fill(fCentrality,psiC);
      fHist2DPsiACent[runNumBin][1]->Fill(fCentrality,psiA);
      fHist2DPsiCCentTot[1]->Fill(fCentrality,psiC);
      fHist2DPsiACentTot[1]->Fill(fCentrality,psiA);
    }

    if (bFillHistForSF)
    {
      double psiC = atan2(QyC,QxC)>0. ? atan2(QyC,QxC) : atan2(QyC,QxC)+TMath::TwoPi();
      double psiA = atan2(QyA,-QxA)>0. ? atan2(QyA,-QxA) : atan2(QyA,-QxA)+TMath::TwoPi();
      for (int i = 1; i <= 20; ++i) {
        fProfile2DShiftCentiCosC[runNumBin]->Fill(fCentrality, i-0.5, TMath::Cos(i * psiC));
        fProfile2DShiftCentiSinC[runNumBin]->Fill(fCentrality, i-0.5, TMath::Sin(i * psiC));
        fProfile2DShiftCentiCosA[runNumBin]->Fill(fCentrality, i-0.5, TMath::Cos(i * psiA));
        fProfile2DShiftCentiSinA[runNumBin]->Fill(fCentrality, i-0.5, TMath::Sin(i * psiA));
      }
    }
    
    if (bApplySF)
    {
      fProfile2DForCosC = (TProfile2D*)fZDCCalibrationListThisRun->FindObject("fProfile2DShiftCentiCosC");
      fProfile2DForSinC = (TProfile2D*)fZDCCalibrationListThisRun->FindObject("fProfile2DShiftCentiSinC");
      fProfile2DForCosA = (TProfile2D*)fZDCCalibrationListThisRun->FindObject("fProfile2DShiftCentiCosA");
      fProfile2DForSinA = (TProfile2D*)fZDCCalibrationListThisRun->FindObject("fProfile2DShiftCentiSinA");
      double psiC = atan2(QyC,QxC)>0. ? atan2(QyC,QxC) : atan2(QyC,QxC)+TMath::TwoPi();
      double psiA = atan2(QyA,-QxA)>0. ? atan2(QyA,-QxA) : atan2(QyA,-QxA)+TMath::TwoPi();
      for (int i = 1; i <= 20; i++)
      {
        double shiftCosC = fProfile2DForCosC->GetBinContent(fProfile2DForCosC->GetXaxis()->FindBin(fCentrality),i);
        double shiftSinC = fProfile2DForSinC->GetBinContent(fProfile2DForSinC->GetXaxis()->FindBin(fCentrality),i);
        double shiftCosA = fProfile2DForCosA->GetBinContent(fProfile2DForCosA->GetXaxis()->FindBin(fCentrality),i);
        double shiftSinA = fProfile2DForSinA->GetBinContent(fProfile2DForSinA->GetXaxis()->FindBin(fCentrality),i);
        psiC += (2./i) * (-shiftSinC * TMath::Cos(i*psiC) + shiftCosC * TMath::Sin(i*psiC));
        psiA += (2./i) * (-shiftSinA * TMath::Cos(i*psiA) + shiftCosA * TMath::Sin(i*psiA));
      }
      PsiZNCSF = psiC;
      PsiZNASF = psiA;

      if (bFillQAHist)
      {
        fHist2DPsiCCent[runNumBin][2]->Fill(fCentrality, psiC);
        fHist2DPsiACent[runNumBin][2]->Fill(fCentrality, psiA);
        fHist2DPsiCCentTot[2]->Fill(fCentrality,psiC);
        fHist2DPsiACentTot[2]->Fill(fCentrality,psiA);
      }
    }

    int centBin = centV0M/10;
    if (bCalculateV2) {
      Int_t nTracks(fAOD->GetNumberOfTracks());
      for(Int_t iTrack(0); iTrack < nTracks; iTrack++) {
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(iTrack));
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
        fHistPt[0]->Fill(pt);
        fHistEta[0]->Fill(eta);
        fHistPhi[0]->Fill(phi);
    
        if (pt < 0.2 || pt > 30.) continue;
        if (fabs(eta) > 0.8) continue;
        if (fabs(nhits) < 70) continue;
        if (chi2 > 4.) continue;
    
        fHistPt[1]->Fill(pt);
        fHistEta[1]->Fill(eta);
        fHistPhi[1]->Fill(phi);
    
        //------------------
        // dca
        //------------------
        double mag = fAOD->GetMagneticField();
        double dcaxy  = 999.;
        double dcaz   = 999.;
        double r[3];
        double dca[2];
        double cov[3];
        bool proptodca = track->PropagateToDCA(fAOD -> GetPrimaryVertex(), mag, 100., dca, cov);
        if (track->GetXYZ(r)) {
          dcaxy = r[0];
          dcaz  = r[1];
        } else {
          double dcax = r[0] - fVtx[0];
          double dcay = r[1] - fVtx[1];
          dcaz  = r[2] - fVtx[2];
          dcaxy = sqrt(dcax*dcax + dcay*dcay);
          // dcaxy = dca[0];
        }
        if (fabs(dcaxy) > 2.4) continue;//DCAxy cut
        if (fabs(dcaz) > 3.2) continue;//DCAz cut

        //------------------
        // FLOW
        //------------------
        Double_t v2ZNCGETmp = TMath::Cos(2*(phi-PsiZNCGE));
        Double_t v2ZNAGETmp = TMath::Cos(2*(phi-PsiZNAGE));
        Double_t v2ZNCRCTmp = TMath::Cos(2*(phi-PsiZNCRC));
        Double_t v2ZNARCTmp = TMath::Cos(2*(phi-PsiZNARC));
        Double_t v2ZNCSFTmp = TMath::Cos(2*(phi-PsiZNCSF));
        Double_t v2ZNASFTmp = TMath::Cos(2*(phi-PsiZNASF));

        fProfileV2PsiZNCCent[0]->Fill(centV0M,v2ZNCGETmp);
        fProfileV2PsiZNACent[0]->Fill(centV0M,v2ZNAGETmp);
        fProfileV2PsiZNCCent[1]->Fill(centV0M,v2ZNCRCTmp);
        fProfileV2PsiZNACent[1]->Fill(centV0M,v2ZNARCTmp);
        fProfileV2PsiZNCCent[2]->Fill(centV0M,v2ZNCSFTmp);
        fProfileV2PsiZNACent[2]->Fill(centV0M,v2ZNASFTmp);

        fProfilePhiPsiZNCCent[centBin][0]->Fill(phi-PsiZNCGE);
        fProfilePhiPsiZNACent[centBin][0]->Fill(phi-PsiZNAGE);
        fProfilePhiPsiZNCCent[centBin][1]->Fill(phi-PsiZNCRC);
        fProfilePhiPsiZNACent[centBin][1]->Fill(phi-PsiZNARC);
        fProfilePhiPsiZNCCent[centBin][2]->Fill(phi-PsiZNCSF);
        fProfilePhiPsiZNACent[centBin][2]->Fill(phi-PsiZNASF);

        fProfileV2PsiZNCVsPt[centBin][0]->Fill(pt,v2ZNCGETmp);
        fProfileV2PsiZNAVsPt[centBin][0]->Fill(pt,v2ZNAGETmp);
        fProfileV2PsiZNCVsPt[centBin][1]->Fill(pt,v2ZNCRCTmp);
        fProfileV2PsiZNAVsPt[centBin][1]->Fill(pt,v2ZNARCTmp);
        fProfileV2PsiZNCVsPt[centBin][2]->Fill(pt,v2ZNCSFTmp);
        fProfileV2PsiZNAVsPt[centBin][2]->Fill(pt,v2ZNASFTmp);
      }
    }
    
                                                        // continue until all the tracks are processed
    PostData(1, fOutputList);                           // stream the results the analysis of this event to
                                                        // the output manager which will take care of writing
                                                        // it to a file
}
//_____________________________________________________________________________
void AliAnalysisTaskZDCCalibration::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________

int AliAnalysisTaskZDCCalibration::GetRunNumBin(int runNum)
{
  int runNumBin=-1;
  // 10h
  int runNumList10h[91]={
    139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328, 139314, 139310,
    139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029, 139028, 138872, 138871, 138870,
    138837, 138732, 138730, 138666, 138662, 138653, 138652, 138638, 138624, 138621, 138583, 138582,
    138579, 138578, 138534, 138469, 138442, 138439, 138438, 138396, 138364, 138275, 138225, 138201,
    138197, 138192, 138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693,
    137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595, 137549, 137546, 137544, 137541,
    137539, 137531, 137530, 137443, 137441, 137440, 137439, 137434, 137432, 137431, 137430, 137243,
    137236, 137235, 137232, 137231, 137230, 137162, 137161
  };  
  //11h
  int runNumList11h[39]={
    170387,170040,170268,170228,170207,169838,170159,170204,170311,170084,
    169835,170088,170593,170203,170270,169846,170163,170388,170155,170083,
    170572,169837,169855,170306,170269,170089,170309,170091,170081,170230,
    170085,170315,170027,170193,170312,170313,170308,169858,169859
  };
  
  if (fDataSet.Contains("10h")) {
    for (int i = 0; i < 91; ++i) {
      if (runNum==runNumList10h[i]) {runNumBin=i; break;}
      else continue;
    }
  }
  if (fDataSet.Contains("11h")) { 
    for (int i = 0; i < 39; ++i) {
      if (runNum==runNumList11h[i]) {runNumBin=i; break;}
      else continue;
    }
  }

  return runNumBin;
}
