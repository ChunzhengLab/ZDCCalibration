#ifndef AliAnalysisTaskEE_cxx
#define AliAnalysisTaskEE_cxx

//class TList;
//class TH1F;
//class TH2F;
//class TProfile;
//class AliAnalysisUtils;

//#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskEE : public AliAnalysisTaskSE
{
 public:
  AliAnalysisTaskEE();
  AliAnalysisTaskEE(const char *name);
  virtual ~AliAnalysisTaskEE();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);

  void SetFileNameEPMean(TString name) {fileNameEPMean_ = name;}
  void SetFileNameEPShifting(TString name) {fileNameEPShifting_ = name;}

 private:
  int GetRunNumBin(int runNum);
  void GetEPMean(int runNum, int centBin, double vx, double vy);

  // AliEventplane* fEP;
  
  TString fileNameEPMean_;
  TFile* fileEPMean_;
  TList* listEPMean_;
  TString fileNameEPShifting_;
  TFile* fileEPShifting_;
  TList* listEPShifting_;
  TProfile3D* p_XAMean_;
  TProfile3D* p_XCMean_;
  TProfile3D* p_YAMean_;
  TProfile3D* p_YCMean_;
  double XAMean_;
  double XCMean_;
  double YAMean_;
  double YCMean_;
  TProfile2D* p_shiftCoeffSinA_;
  TProfile2D* p_shiftCoeffCosA_;
  TProfile2D* p_shiftCoeffSinC_;
  TProfile2D* p_shiftCoeffCosC_;
  TProfile2D* p_shiftCoeffSinFull_;
  TProfile2D* p_shiftCoeffCosFull_;

  TList*     fOutputList;
  TH1I*      h_evtCount;
  TH1D*      h_cent;
  TH1D*      h_mult;
  TH1D*      h_ZNA[4];
  TH1D*      h_ZNC[4];
  TH2D*      h_XYA;
  TH2D*      h_XYC;
  TH1D*      h_psiZDC_A[3][10];
  TH1D*      h_psiZDC_C[3][10];
  TH1D*      h_psiZDC_full[3][10];
  TProfile*  p_psiZDC_reso;
  TH1D*      h_pt;
  TH1D*      h_eta[2];
  TH1D*      h_phi;
  TH1D*      h_dcaXy[2];
  TH1D*      h_dcaZ[2];
  TH1D*      h_nHits[4];
  TH2D*      h_nSigmaE_p[6];
  TH1D*      h_mass_epem;

  AliAnalysisTaskEE(const AliAnalysisTaskEE&);
  AliAnalysisTaskEE& operator=(const AliAnalysisTaskEE&);

  ClassDef(AliAnalysisTaskEE, 1);
};
#endif
