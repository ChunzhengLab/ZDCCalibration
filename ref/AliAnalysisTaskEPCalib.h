#ifndef AliAnalysisTaskEPCalib_cxx
#define AliAnalysisTaskEPCalib_cxx

//class TList;
//class TH1F;
//class TH2F;
//class TProfile;
//class AliAnalysisUtils;

//#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskEPCalib : public AliAnalysisTaskSE
{
 public:

  AliAnalysisTaskEPCalib();
  AliAnalysisTaskEPCalib(const char *name);
  virtual ~AliAnalysisTaskEPCalib();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);

  void SetFileNameNUAEntry(TString name) {fileNameNUAEntry_ = name;} 
  bool FillNUAWeight(bool b)             {FillNUAWeight_    = b;} 
  bool NUAWeightOn(bool b)               {NUAWeightOn_      = b;} 
  void FillRecenter(bool b)              {FillRecenter_     = b;} 
  void FillShift(bool b)                 {FillShift_        = b;} 
  void TPCCalibOn(bool b)                {TPCCalibOn_       = b;} 
  void SetHarmonic(double n)             {harmonic_         = n;}

 private:
  int GetRunNumBin(int runNum);
  bool TPCPlane(AliAODEvent* fAOD);
  // void GetEPMean(int runNum, int centBin, double vx, double vy);

  bool FillNUAWeight_;
  bool NUAWeightOn_  ;
  bool FillRecenter_ ;
  bool FillShift_    ;
  bool TPCCalibOn_   ;

  TString fileNameNUAEntry_;
  TFile*  fileNUAEntry_;
  TList*  listNUAEntry_;
  int     runNum_  ;
  double  vzBin_   ;
  int     centBin_ ;
  double  harmonic_;

  // bool doShift_;
  // TString fileNameEPMean_;
  // TFile* fileEPMean_;
  // TList* listEPMean_;
  // TProfile3D* p_XAMean_;
  // TProfile3D* p_XCMean_;
  // TProfile3D* p_YAMean_;
  // TProfile3D* p_YCMean_;
  // double XAMean_;
  // double XCMean_;
  // double YAMean_;
  // double YCMean_;

  TList*      fOutputList;
  TH1I*       h_evtCount;
  TH1I*       h_runNumBin;
  TH1D*       h_cent;
  TH2D*       h_vxy[2];
  TH1D*       h_vz[2];
  TH1D*       h_phi[2];
  TH3I*       h_NUAEntryWrite[10][2];
  TH1D*       h_psiTPC[10][6];

  // TH1D*       h_ZNA[4];
  // TH1D*       h_ZNC[4];
  // TProfile3D* p_XA[39];
  // TProfile3D* p_XC[39];
  // TProfile3D* p_YA[39];
  // TProfile3D* p_YC[39];
  // TH2D*       h_XYA[2];
  // TH2D*       h_XYC[2];
  // TH1D*       h_psiZDC_A[2];
  // TH1D*       h_psiZDC_C[2];
  // TH1D*       h_psiZDC_full[2];
  // TProfile2D* p_shiftCoeffSinA[39];
  // TProfile2D* p_shiftCoeffCosA[39];
  // TProfile2D* p_shiftCoeffSinC[39];
  // TProfile2D* p_shiftCoeffCosC[39];
  // TProfile2D* p_shiftCoeffSinFull[39];
  // TProfile2D* p_shiftCoeffCosFull[39];

  AliAnalysisTaskEPCalib(const AliAnalysisTaskEPCalib&);
  AliAnalysisTaskEPCalib& operator=(const AliAnalysisTaskEPCalib&);

  ClassDef(AliAnalysisTaskEPCalib, 1);
};
#endif
