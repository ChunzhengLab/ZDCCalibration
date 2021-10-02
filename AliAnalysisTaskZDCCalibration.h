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

        virtual void            FillHistForGainEqualization(bool fillHistForGE) {bFillHistForGE = fillHistForGE;}
        virtual void            FillHistForQVZDCRecentering(bool fillHistForRC) {bFillHistForRC = fillHistForRC;}

        virtual void            ApplyGainEqualization(bool applyGE) {bApplyGE = applyGE;}
        virtual void            ApplyQVZDCRecentering(bool applyRC) {bApplyRC = applyRC;}

        virtual void            SetGainEqualizationList(TList* const wlist) {this->fZDCGEList = wlist;}
        virtual void            SetQVZDCRecenteringList(TList* const wlist) {this->fZDCRCList = wlist;}

    private:
        bool                    bFillHistForGE;
        bool                    bFillHistForRC;

        bool                    bApplyGE;
        bool                    bApplyRC;
        TList*                  fZDCGEList;
        TList*                  fZDCRCList;

        double                  fVtx[3];
        double                  fCentrality;
        double                  fCentBin;
        
        AliAODEvent*            fAOD;           //! input event
        AliAODZDC*              fZDC;
        AliAnalysisUtils*       fUtils;
        TList*                  fOutputList;    //! output list
        TH1D*                   fHistCent;
        TH2D*                   fHist2DCentCorr[2];

        TProfile*               fProfileQxAQxCCentTot[2];
        TProfile*               fProfileQxAQyCCentTot[2];
        TProfile*               fProfileQyAQxCCentTot[2];
        TProfile*               fProfileQyAQyCCentTot[2];

        const static int fnRunMax = 91;
        TList*                  fRunList[fnRunMax];
        TH2D*                   fHist2DVxVy[fnRunMax];
        TH1D*                   fHistVz[fnRunMax];
        TProfile*               fProfileZNATowerEnergy[fnRunMax];
        TProfile*               fProfileZNCTowerEnergy[fnRunMax];

        //vxsigma vysigma vz
        THnSparse*                    fHn4DQxZNACentVxVySigmaVz[fnRunMax];
        THnSparse*                    fHn4DQyZNACentVxVySigmaVz[fnRunMax];
        THnSparse*                    fHn4DMtZNACentVxVySigmaVz[fnRunMax];
        THnSparse*                    fHn4DQxZNCCentVxVySigmaVz[fnRunMax];
        THnSparse*                    fHn4DQyZNCCentVxVySigmaVz[fnRunMax];
        THnSparse*                    fHn4DMtZNCCentVxVySigmaVz[fnRunMax];
        
        TProfile*               fProfileQxAQxCCent[fnRunMax][2];
        TProfile*               fProfileQxAQyCCent[fnRunMax][2];
        TProfile*               fProfileQyAQxCCent[fnRunMax][2];
        TProfile*               fProfileQyAQyCCent[fnRunMax][2];

        TH2D*                   fHist2DPsiACentBin[fnRunMax][2];
        TH2D*                   fHist2DPsiCCentBin[fnRunMax][2];

        int                     GetCentBin(double centrality);
        int                     GetRunNumBin(int runNum);
        int                     GetVxySigmaBin(double v);

        AliAnalysisTaskZDCCalibration(const AliAnalysisTaskZDCCalibration&); // not implemented
        AliAnalysisTaskZDCCalibration& operator=(const AliAnalysisTaskZDCCalibration&); // not implemented        

        ClassDef(AliAnalysisTaskZDCCalibration, 1);
};

#endif
