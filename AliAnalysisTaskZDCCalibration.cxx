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

#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "THnSparse.h"
#include "TList.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskZDCCalibration.h"
#include "AliAnalysisUtils.h"

class AliAnalysisTaskZDCCalibration;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskZDCCalibration) // classimp: necessary for root

AliAnalysisTaskZDCCalibration::AliAnalysisTaskZDCCalibration() : AliAnalysisTaskSE(), 
    fAOD(0),fZDC(0),fUtils(0),fHistCent(0),fOutputList(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
    for (size_t i = 0; i < 2; i++)
    {
      fProfileQxAQxCCentTot[i] = NULL;
      fProfileQxAQyCCentTot[i] = NULL;
      fProfileQyAQxCCentTot[i] = NULL;
      fProfileQyAQyCCentTot[i] = NULL;
    }

    for (size_t iRun = 0; iRun < fnRunMax; iRun++)
    {
      fRunList[iRun]                  = NULL;
      fHist2DVxVy[iRun]               = NULL;
      fHistVz[iRun]                   = NULL;
      fProfileZNATowerEnergy[iRun]    = NULL;
      fProfileZNCTowerEnergy[iRun]    = NULL;

      //vx vy sigma vz
      fHn4DQxZNACentVxVySigmaVz[iRun] = NULL;
      fHn4DQyZNACentVxVySigmaVz[iRun] = NULL;
      fHn4DMtZNACentVxVySigmaVz[iRun] = NULL;
      fHn4DQxZNCCentVxVySigmaVz[iRun] = NULL;
      fHn4DQyZNCCentVxVySigmaVz[iRun] = NULL;
      fHn4DMtZNCCentVxVySigmaVz[iRun] = NULL;
        
      for (size_t i = 0; i < 2; i++)
      {
        fProfileQxAQxCCent[iRun][i] = NULL;
        fProfileQxAQyCCent[iRun][i] = NULL;
        fProfileQyAQxCCent[iRun][i] = NULL;
        fProfileQyAQyCCent[iRun][i] = NULL;

        fHist2DPsiACentBin[iRun][i] = NULL;
        fHist2DPsiCCentBin[iRun][i] = NULL;
      }
    }
}
//_____________________________________________________________________________
AliAnalysisTaskZDCCalibration::AliAnalysisTaskZDCCalibration(const char* name) : AliAnalysisTaskSE(name),
    fAOD(0),fZDC(0),fUtils(0),fHistCent(0),fOutputList(0)
{
    for (size_t i = 0; i < 2; i++)
    {
      fProfileQxAQxCCentTot[i] = NULL;
      fProfileQxAQyCCentTot[i] = NULL;
      fProfileQyAQxCCentTot[i] = NULL;
      fProfileQyAQyCCentTot[i] = NULL;
    }

    for (size_t iRun = 0; iRun < fnRunMax; iRun++)
    {
      fRunList[iRun]                  = NULL;
      fHist2DVxVy[iRun]               = NULL;
      fHistVz[iRun]                   = NULL;
      fProfileZNATowerEnergy[iRun]    = NULL;
      fProfileZNCTowerEnergy[iRun]    = NULL;

      //vx vy sigma vz
      fHn4DQxZNACentVxVySigmaVz[iRun] = NULL;
      fHn4DQyZNACentVxVySigmaVz[iRun] = NULL;
      fHn4DMtZNACentVxVySigmaVz[iRun] = NULL;
      fHn4DQxZNCCentVxVySigmaVz[iRun] = NULL;
      fHn4DQyZNCCentVxVySigmaVz[iRun] = NULL;
      fHn4DMtZNCCentVxVySigmaVz[iRun] = NULL;
        
      for (size_t i = 0; i < 2; i++)
      {
        fProfileQxAQxCCent[iRun][i] = NULL;
        fProfileQxAQyCCent[iRun][i] = NULL;
        fProfileQyAQxCCent[iRun][i] = NULL;
        fProfileQyAQyCCent[iRun][i] = NULL;

        fHist2DPsiACentBin[iRun][i] = NULL;
        fHist2DPsiCCentBin[iRun][i] = NULL;
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
    if(fOutputList) {
        delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    }
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
    fHistCent = new TH1D("fHistCent","Hist Centrality", 100, 0., 100.);
    fHist2DCentCorr[0] = new TH2D("fHist2DCentCorrBfCut","Centrality V0M vs. TRK", 100, 0., 100., 100, 0., 100.);
    fHist2DCentCorr[1] = new TH2D("fHist2DCentCorrBfCut","Centrality V0M vs. TRK", 100, 0., 100., 100, 0., 100.);

    fProfileQxAQxCCentTot[0] = new TProfile("fProfileQxAQxCCentTotBfEC","<XAXC> before Recentering", 100, 0., 100.);
    fProfileQxAQyCCentTot[0] = new TProfile("fProfileQxAQyCCentTotBfEC","<XAYC> before Recentering", 100, 0., 100.);
    fProfileQyAQxCCentTot[0] = new TProfile("fProfileQyAQxCCentTotBfEC","<YAXC> before Recentering", 100, 0., 100.);
    fProfileQyAQyCCentTot[0] = new TProfile("fProfileQyAQyCCentTotBfEC","<YAYC> before Recentering", 100, 0., 100.);

    fProfileQxAQxCCentTot[1] = new TProfile("fProfileQxAQxCCentTotAfEC","<XAXC> after Recentering", 100, 0., 100.);        
    fProfileQxAQyCCentTot[1] = new TProfile("fProfileQxAQyCCentTotAfEC","<XAYC> after Recentering", 100, 0., 100.);        
    fProfileQyAQxCCentTot[1] = new TProfile("fProfileQyAQxCCentTotAfEC","<YAXC> after Recentering", 100, 0., 100.);        
    fProfileQyAQyCCentTot[1] = new TProfile("fProfileQyAQyCCentTotAfEC","<YAYC> after Recentering", 100, 0., 100.);    


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


    for (int iRun = 0; iRun < fnRunMax; iRun++)
    {
      fRunList[iRun] = new TList();
      fRunList[iRun] -> SetOwner(kTRUE);
      fRunList[iRun] -> SetName(Form("%i", runNumList[iRun]));


      fHist2DVxVy[iRun] = new TH2D("hist2DVxVy","hist2DVxVy",200,-1,1,200,-1,1);
      fHistVz[iRun] = new TH1D("histVz","histVz",200,-10.,10.); 
      fRunList[iRun] -> Add(fHist2DVxVy[iRun]);
      fRunList[iRun] -> Add(fHistVz[iRun]);

      if (bFillHistForGE)
      {
        fProfileZNCTowerEnergy[iRun] = new TProfile("fProfileZNCTowerEnergy","fProfileZNCTowerEnergy",5,0,5);
        fProfileZNATowerEnergy[iRun] = new TProfile("fProfileZNATowerEnergy","fProfileZNATowerEnergy",5,0,5);
        fRunList[iRun] -> Add(fProfileZNCTowerEnergy[iRun]);
        fRunList[iRun] -> Add(fProfileZNATowerEnergy[iRun]);
      }

      if (bApplyGE)
      {
        fProfileForZNAGE = new TProfile();
        fProfileForZNAGE = new TProfile();
      }
      
      if (bFillHistForRC)
      {
        int   nbins[4] = {100,  3,  3,   3};
        double xmin[4] = {  0,  0,  0, -10};
        double xmax[4] = {100,  3,  3,  10};
        fHn4DQxZNACentVxVySigmaVz[iRun] = new THnSparseD("fHn4DQxZNACentVxVySigmaVz","fHn4DQxZNACentVxVySigmaVz",4,nbins,xmin,xmax);
        fHn4DQyZNACentVxVySigmaVz[iRun] = new THnSparseD("fHn4DQyZNACentVxVySigmaVz","fHn4DQyZNACentVxVySigmaVz",4,nbins,xmin,xmax);
        fHn4DMtZNACentVxVySigmaVz[iRun] = new THnSparseD("fHn4DMtZNACentVxVySigmaVz","fHn4DMtZNACentVxVySigmaVz",4,nbins,xmin,xmax);
        fHn4DQxZNCCentVxVySigmaVz[iRun] = new THnSparseD("fHn4DQxZNCCentVxVySigmaVz","fHn4DQxZNCCentVxVySigmaVz",4,nbins,xmin,xmax);
        fHn4DQyZNCCentVxVySigmaVz[iRun] = new THnSparseD("fHn4DQyZNCCentVxVySigmaVz","fHn4DQyZNCCentVxVySigmaVz",4,nbins,xmin,xmax);
        fHn4DMtZNCCentVxVySigmaVz[iRun] = new THnSparseD("fHn4DMtZNCCentVxVySigmaVz","fHn4DMtZNCCentVxVySigmaVz",4,nbins,xmin,xmax);
        fRunList[iRun] -> Add(fHn4DQxZNACentVxVySigmaVz[iRun]);
        fRunList[iRun] -> Add(fHn4DQyZNACentVxVySigmaVz[iRun]);
        fRunList[iRun] -> Add(fHn4DMtZNACentVxVySigmaVz[iRun]);
        fRunList[iRun] -> Add(fHn4DQxZNCCentVxVySigmaVz[iRun]);
        fRunList[iRun] -> Add(fHn4DQyZNCCentVxVySigmaVz[iRun]);
        fRunList[iRun] -> Add(fHn4DMtZNCCentVxVySigmaVz[iRun]);
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

      fOutputList -> Add(fHistCent);
      fOutputList -> Add(fHist2DCentCorr[0]);
      fOutputList -> Add(fHist2DCentCorr[1]);
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
    fUtils->SetUseOutOfBunchPileUp(true);
    fUtils->SetUseMVPlpSelection(true);
    bool isPileup = fUtils->IsPileUpEvent(fAOD);
    if (isPileup) return;
    
    double centV0M = fAOD->GetCentrality()->GetCentralityPercentile("V0M");
    double centCL1 = fAOD->GetCentrality()->GetCentralityPercentile("CL1");
    double centTRK = fAOD->GetCentrality()->GetCentralityPercentile("TRK");
    fHist2DCentCorr[0]->Fill(centV0M,centCL1);
    if (fabs(centV0M-centCL1)>7.5) return;
    if (fabs(centV0M-centTRK)>5) return; // ANA-280
    if (centV0M<0 || centV0M>=100) return;
    fHist2DCentCorr[1]->Fill(centV0M,centCL1);
    fCentrality = centV0M;

    // vertex
    fAOD -> GetPrimaryVertex() -> GetXYZ(fVtx);
    double vzSPD  = fAOD->GetPrimaryVertexSPD()->GetZ();
    if (fabs(fVtx[2])>10) return;
    if (fabs(fVtx[2]-vzSPD)>0.5) return;
    if (fabs(fVtx[0])<1e-6 || fabs(fVtx[1])<1e-6 || fabs(fVtx[2])<1e-6) return;
    fHist2DVxVy[runNumBin] ->Fill(fVtx[0],fVtx[1]);
    fHistVz[runNumBin] ->Fill(fVtx[2]);

    //GetVexBin
    int vxBin = -1;
    int vyBin = -1;
    if(bGetVetexBin){
      fHist2DForMeanVxVy = (TH2D*)fVetexList->FindObject(Form("%d",runNumBin))->FindObject("fHist2DVxVy"); 
      fHistForMeanVz     = (TH1D*)fVetexList->FindObject(Form("%d",runNumBin))->FindObject("fHistVz");
      double vxMean = fHist2DForMeanVxVy->ProjectionX()->GetMean();
      double vyMean = fHist2DForMeanVxVy->ProjectionY()->GetMean();
      double vxSigmaMean = fHist2DForMeanVxVy->ProjectionX()->GetRMS();
      double vySigmaMean = fHist2DForMeanVxVy->ProjectionY()->GetRMS();

      if(vxSigmaMean < 1.e-6 || vxSigmaMean < 1.e-6) return;

      double vxSigma = (fVtx[0] - vxMean)/vxSigmaMean;
      double vySigma = (fVtx[1] - vyMean)/vxSigmaMean;

      //TODO
      if     (vxSigma>xx && vxSigma<xx){vxBin=1;}
      else if(vxSigma>xx && vxSigma<xx){vxBin=2;}
      else if(vxSigma>xx && vxSigma<xx){vxBin=3;}

      if     (vySigma>xx && vySigma<xx){vxBin=1;}
      else if(vySigma>xx && vySigma<xx){vxBin=2;}
      else if(vySigma>xx && vySigma<xx){vxBin=3;}

      if(vxBin <= 0 && vyBin <= 0) return;
    
    

    const double x[4] = {-1.75, 1.75, -1.75, 1.75};
    const double y[4] = {-1.75, -1.75, 1.75, 1.75};
    const double* EZNARaw = fZDC->GetZNATowerEnergy();
    const double* EZNCRaw = fZDC->GetZNCTowerEnergy();

    double EZNC[4] = {0.,0.,0.,0.}; 
    double EZNA[4] = {0.,0.,0.,0.};

    if(bFillHistForGE) {
      for (int iTower = 0; iTower < 5; ++iTower) {
        fProfileZNCTowerEnergy[runNumBin] -> Fill(iTower+0.5, EZNCRaw[iTower]);
        fProfileZNATowerEnergy[runNumBin] -> Fill(iTower+0.5, EZNARaw[iTower]);
      }
    }

    if(bApplyGE && fZDCGEList && !bFillHistForGE) {
      fProfileForZNCGE = (TProfile*)fZDCGEList->FindObject(Form("%d",runNumBin))->FindObject("fProfileZNCTowerEnergy"); 
      fProfileForZNAGE = (TProfile*)fZDCGEList->FindObject(Form("%d",runNumBin))->FindObject("fProfileZNATowerEnergy"); 
      for (int iTower = 0; iTower < 4; ++iTower) {
        EZNC[iTower] = EZNARaw[iTower + 1] * (fProfileForZNCGE->GetBinContent(iTower + 1.5)) / (fProfileZNCTowerEnergy[runNumBin]->GetBinContent(1.5));
        EZNA[iTower] = EZNARaw[iTower + 1] * (fProfileForZNAGE->GetBinContent(iTower + 1.5)) / (fProfileZNATowerEnergy[runNumBin]->GetBinContent(1.5));
      }
    }

    double QxC = 0.;
    double QyC = 0.;
    double MC  = 0.;
    double QxA = 0.;
    double QyA = 0.;
    double MA  = 0.;

    for (int iTower = 0; iTower < 4; ++iTower) {
      QxC += EZNC[iTower] * x[iTower];
      QyC += EZNC[iTower] * y[iTower];
      MC  += EZNC[iTower];
      QxA += EZNA[iTower] * x[iTower];
      QyA += EZNA[iTower] * y[iTower];
      MA  += EZNA[iTower];
    }

    if(MC < 1.e-6) return;
    if(MA < 1.e-6) return;

    double fillPosition[4] = {fCentrality,vxBin-0.5,vyBin-0.5,fVtx[2]};
    if(bFillHistForRC) {
      fHn4DQxZNCCentVxVySigmaVz[runNum] -> Fill(fillPosition,QxC);
      fHn4DQyZNCCentVxVySigmaVz[runNum] -> Fill(fillPosition,QyC);
      fHn4DMtZNCCentVxVySigmaVz[runNum] -> Fill(fillPosition,MC);

      fHn4DQxZNACentVxVySigmaVz[runNum] -> Fill(fillPosition,QxA);
      fHn4DQyZNACentVxVySigmaVz[runNum] -> Fill(fillPosition,QyA);
      fHn4DQyZNACentVxVySigmaVz[runNum] -> Fill(fillPosition,MA);
    }

    QxC /= MC; 
    QyC /= MC; 
    QxA /= MA; 
    QyA /= MA; 
    
    fProfileQxAQxCCent[runNum][0]->Fill(fCentrality,QxA*QxC);
    fProfileQxAQyCCent[runNum][0]->Fill(fCentrality,QxA*QyC);
    fProfileQyAQxCCent[runNum][0]->Fill(fCentrality,QyA*QxC);
    fProfileQyAQyCCent[runNum][0]->Fill(fCentrality,QyA*QyC);

    fProfileQxAQxCCentTot[0]->Fill(fCentrality,QxA*QxC);
    fProfileQxAQyCCentTot[0]->Fill(fCentrality,QxA*QyC);
    fProfileQyAQxCCentTot[0]->Fill(fCentrality,QyA*QxC);
    fProfileQyAQyCCentTot[0]->Fill(fCentrality,QyA*QyC);

    fHist2DPsiACentBin[fnRunMax][0]->Fill(fCentrality,atan2(QyA,QxA));
    fHist2DPsiCCentBin[fnRunMax][0]->Fill(fCentrality,atan2(QyC,QxC));

    if(bApplyRC) {
      
      double QxCMean = fHn4DQxZNCCentVxVySigmaVz[runNum] -> GetBinContent(fHn4DQxZNCCentVxVySigmaVz[runNum]->GetBin(fillPosition));
      double QyCMean = fHn4DQyZNCCentVxVySigmaVz[runNum] -> GetBinContent(fHn4DQyZNCCentVxVySigmaVz[runNum]->GetBin(fillPosition));
      double MCMean  = fHn4DMtZNCCentVxVySigmaVz[runNum] -> GetBinContent(fHn4DMtZNCCentVxVySigmaVz[runNum]->GetBin(fillPosition));

      double QxAMean = fHn4DQxZNACentVxVySigmaVz[runNum] -> GetBinContent(fHn4DQxZNACentVxVySigmaVz[runNum]->GetBin(fillPosition));
      double QyAMean = fHn4DQyZNACentVxVySigmaVz[runNum] -> GetBinContent(fHn4DQyZNACentVxVySigmaVz[runNum]->GetBin(fillPosition));
      double MAMean  = fHn4DMtZNACentVxVySigmaVz[runNum] -> GetBinContent(fHn4DMtZNACentVxVySigmaVz[runNum]->GetBin(fillPosition));

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
    

      fProfileQxAQxCCent[runNum][1]->Fill(fCentrality,QxA*QxC);
      fProfileQxAQyCCent[runNum][1]->Fill(fCentrality,QxA*QyC);
      fProfileQyAQxCCent[runNum][1]->Fill(fCentrality,QyA*QxC);
      fProfileQyAQyCCent[runNum][1]->Fill(fCentrality,QyA*QyC);

      fProfileQxAQxCCentTot[1]->Fill(fCentrality,QxA*QxC);
      fProfileQxAQyCCentTot[1]->Fill(fCentrality,QxA*QyC);
      fProfileQyAQxCCentTot[1]->Fill(fCentrality,QyA*QxC);
      fProfileQyAQyCCentTot[1]->Fill(fCentrality,QyA*QyC);

      fHist2DPsiACentBin[fnRunMax][1]->Fill(fCentrality,atan2(QyA,QxA));
      fHist2DPsiCCentBin[fnRunMax][1]->Fill(fCentrality,atan2(QyC,QxC));
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
int AliAnalysisTaskZDCCalibration::GetCentBin(double centrality)
{
  int CentBin=-1;
  if (centrality>0.  && centrality<5. ) CentBin=0;
  if (centrality>5.  && centrality<10.) CentBin=1;
  if (centrality>10. && centrality<20.) CentBin=2;
  if (centrality>20. && centrality<30.) CentBin=3;
  if (centrality>30. && centrality<40.) CentBin=4;
  if (centrality>40. && centrality<50.) CentBin=5;
  if (centrality>50. && centrality<60.) CentBin=6;
  if (centrality>60. && centrality<70.) CentBin=7;
  if (centrality>70. && centrality<80.) CentBin=8;
  if (centrality>80. && centrality<90.) CentBin=9;
  if (centrality>90. && centrality<100.) CentBin=10;
  return CentBin;
}

int AliAnalysisTaskZDCCalibration::GetRunNumBin(int runNum)
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