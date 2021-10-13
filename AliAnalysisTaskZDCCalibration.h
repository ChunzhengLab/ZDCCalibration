/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskZDCCalibration_H
#define AliAnalysisTaskZDCCalibration_H

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskZDCCalibration : public AliAnalysisTaskSE  
{
    public:
                                AliAnalysisTaskZDCCalibration();
                                AliAnalysisTaskZDCCalibration(const char *name);
        virtual                 ~AliAnalysisTaskZDCCalibration();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);

        virtual void            FirstFillHistVetex(bool fillHistVetex)          {bFirstFillHistVetex = fillHistVetex;}
        virtual void            FillHistForGainEqualization(bool fillHistForGE) {bFillHistForGE = fillHistForGE;}
        virtual void            FillHistForQVZDCRecentering(bool fillHistForRC) {bFillHistForRC = fillHistForRC;}
        virtual void            FillHistForZDCPlaneShifting(bool fillHitsForSF) {bFillHistForSF = fillHitsForSF;}

        virtual void            ApplyGainEqualization(bool applyGE) {bApplyGE = applyGE;}
        virtual void            ApplyQVZDCRecentering(bool applyRC) {bApplyRC = applyRC;}
        virtual void            ApplyZDCPlaneShifting(bool applySF) {bApplySF = applySF;}

        virtual void            SetZDCCalibrationList(TList* const wlist)   {this->fZDCCalibrationList = wlist;}

        virtual void            FillQAHist(bool fillQAHist)                 {bFillQAHist = fillQAHist;}

    private:
        bool                    bFirstFillHistVetex;
        bool                    bFillHistForGE;
        bool                    bFillHistForRC;
        bool                    bFillHistForSF;

        bool                    bApplyGE;
        bool                    bApplyRC;
        bool                    bApplySF;
        bool                    bFillQAHist;

        TList*                  fZDCCalibrationList;

        double                  fVtx[3];
        double                  fCentrality;
        
        AliAODEvent*            fAOD;           //! input event
        AliAODZDC*              fZDC;
        AliAnalysisUtils*       fUtils;
        TList*                  fOutputList;    //! output list
        TList*                  fQAList;
        TH1D*                   fHistCent;
        TH2D*                   fHist2DCentCorr[2];
        TH2D*                   fHist2DVxVyTot;
        TH1D*                   fHistVzTot;

        TProfile*               fProfileQxCCentTot[2];
        TProfile*               fProfileQyCCentTot[2];
        TProfile*               fProfileQxACentTot[2];
        TProfile*               fProfileQyACentTot[2];
        TProfile*               fProfileQxAQxCCentTot[2];
        TProfile*               fProfileQxAQyCCentTot[2];
        TProfile*               fProfileQyAQxCCentTot[2];
        TProfile*               fProfileQyAQyCCentTot[2];

        TH2D*                   fHist2DPsiCCentTot[3];
        TH2D*                   fHist2DPsiACentTot[3];

        const static int        fnRunMax = 91;
        TList*                  fRunList[fnRunMax];
        TList*                  fQAListThisRun[fnRunMax];

        //Get Vtx
        //Write
        TH2D*                   fHist2DVxVy[fnRunMax];
        TH1D*                   fHistVz[fnRunMax];
        //Read
        TH2D*                   fHist2DForMeanVxVy;
        TH1D*                   fHistForMeanVz;    

        //For GE
        //Write
        TProfile*               fProfileZNCTowerEnergy[fnRunMax];
        TProfile*               fProfileZNATowerEnergy[fnRunMax];
        TProfile*               fProfileZNCTowerEnergyBfGE[fnRunMax];
        TProfile*               fProfileZNATowerEnergyBfGE[fnRunMax];
        //Read
        TProfile*               fProfileForZNCGE;
        TProfile*               fProfileForZNAGE;

        //For RC
        //vxsigma vysigma vz
        //Write
        THnSparseD*              fHn4DQxZNCCentVxVySigmaVz[fnRunMax];
        THnSparseD*              fHn4DQyZNCCentVxVySigmaVz[fnRunMax];
        THnSparseD*              fHn4DMtZNCCentVxVySigmaVz[fnRunMax];
        THnSparseD*              fHn4DQxZNACentVxVySigmaVz[fnRunMax];
        THnSparseD*              fHn4DQyZNACentVxVySigmaVz[fnRunMax];
        THnSparseD*              fHn4DMtZNACentVxVySigmaVz[fnRunMax];
        //Read
        THnSparseD*              fHn4DForZNCQxRC;
        THnSparseD*              fHn4DForZNCQyRC;
        THnSparseD*              fHn4DForZNCMtRC;
        THnSparseD*              fHn4DForZNAQxRC;
        THnSparseD*              fHn4DForZNAQyRC;
        THnSparseD*              fHn4DForZNAMtRC;

        //For SF
        //Write
        TProfile2D*              fProfile2DShiftCentiCosC[fnRunMax];
        TProfile2D*              fProfile2DShiftCentiSinC[fnRunMax];
        TProfile2D*              fProfile2DShiftCentiCosA[fnRunMax];
        TProfile2D*              fProfile2DShiftCentiSinA[fnRunMax];
        //Read
        TProfile2D*              fProfile2DForCosC;
        TProfile2D*              fProfile2DForSinC;
        TProfile2D*              fProfile2DForCosA;
        TProfile2D*              fProfile2DForSinA;

        //Corr
        TProfile*                fProfileQxCCent[fnRunMax][2];
        TProfile*                fProfileQyCCent[fnRunMax][2];
        TProfile*                fProfileQxACent[fnRunMax][2];
        TProfile*                fProfileQyACent[fnRunMax][2];
        TProfile*                fProfileQxAQxCCent[fnRunMax][2];
        TProfile*                fProfileQxAQyCCent[fnRunMax][2];
        TProfile*                fProfileQyAQxCCent[fnRunMax][2];
        TProfile*                fProfileQyAQyCCent[fnRunMax][2];

        //Psi
        TH2D*                    fHist2DPsiCCent[fnRunMax][3];
        TH2D*                    fHist2DPsiACent[fnRunMax][3];

        int                     GetRunNumBin(int runNum);

        AliAnalysisTaskZDCCalibration(const AliAnalysisTaskZDCCalibration&); // not implemented
        AliAnalysisTaskZDCCalibration& operator=(const AliAnalysisTaskZDCCalibration&); // not implemented        

        ClassDef(AliAnalysisTaskZDCCalibration, 1);
};

#endif