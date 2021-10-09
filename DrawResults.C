#include "TFile.h"
#include "TProfile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"

void DrawResults() {
  TFile* inputFile = TFile::Open("AnalysisResults.root");
  TDirectory* inputDir = inputFile->GetDirectory("myTask");
  TList* listTot = nullptr;
  inputDir->GetObject("Output",listTot);

  TH1D* hCent = (TH1D*)listTot->FindObject("fHistCent");
  TH2D* hCentCorrBfCut = (TH2D*)listTot->FindObject("fHist2DCentCorrBfCut");
  TH2D* hCentCorrAfCut = (TH2D*)listTot->FindObject("fHist2DCentCorrAfCut");
  TH2D* hVxVy = (TH2D*)listTot->FindObject("fHist2DVxVyTot");
  TH1D* hVz = (TH1D*)listTot->FindObject("fHistVzTot");
  TList* listThisRun = (TList*)listTot->FindObject("139510");

  TProfile* pZNCTowerEnergyThisRun[2];
  TProfile* pZNATowerEnergyThisRun[2];

  for(int i = 0; i<2; i++) pZNCTowerEnergyThisRun[i] = nullptr;
  for(int i = 0; i<2; i++) pZNATowerEnergyThisRun[i] = nullptr;

  pZNCTowerEnergyThisRun[0] = (TProfile*)listThisRun->FindObject("profileZNCTowerEnergy");
  pZNATowerEnergyThisRun[0] = (TProfile*)listThisRun->FindObject("profileZNATowerEnergy");
  pZNCTowerEnergyThisRun[1] = (TProfile*)listThisRun->FindObject("profileZNCTowerEnergyBfGE");
  pZNATowerEnergyThisRun[1] = (TProfile*)listThisRun->FindObject("profileZNATowerEnergyBfGE");

  
  TList* QAThisRun = (TList*)listThisRun->FindObject("QA");
  TProfile* pQxCCentThisRun[2];
  TProfile* pQyCCentThisRun[2];
  TProfile* pQxACentThisRun[2];
  TProfile* pQyACentThisRun[2];
  TProfile* pQxAQxCCentThisRun[2];
  TProfile* pQxAQyCCentThisRun[2];
  TProfile* pQyAQxCCentThisRun[2];
  TProfile* pQyAQyCCentThisRun[2];
  TH2D*     hPsiACentThisRun[2];
  TH2D*     hPsiCCentThisRun[2];

  for(int i = 0; i<2; i++) pQxCCentThisRun[i] = nullptr;
  for(int i = 0; i<2; i++) pQyCCentThisRun[i] = nullptr;
  for(int i = 0; i<2; i++) pQxACentThisRun[i] = nullptr;
  for(int i = 0; i<2; i++) pQyACentThisRun[i] = nullptr;
  for(int i = 0; i<2; i++) pQxAQxCCentThisRun[i] = nullptr;
  for(int i = 0; i<2; i++) pQxAQyCCentThisRun[i] = nullptr;
  for(int i = 0; i<2; i++) pQyAQxCCentThisRun[i] = nullptr;
  for(int i = 0; i<2; i++) pQyAQyCCentThisRun[i] = nullptr;
  for(int i = 0; i<2; i++) hPsiACentThisRun[i] = nullptr;
  for(int i = 0; i<2; i++) hPsiCCentThisRun[i] = nullptr;

  pQxCCentThisRun[0] = (TProfile*)QAThisRun->FindObject("fProfileQxCCentBfRC");
  pQyCCentThisRun[0] = (TProfile*)QAThisRun->FindObject("fProfileQyCCentBfRC");
  pQxACentThisRun[0] = (TProfile*)QAThisRun->FindObject("fProfileQxACentBfRC");
  pQyACentThisRun[0] = (TProfile*)QAThisRun->FindObject("fProfileQyACentBfRC");

  pQxAQxCCentThisRun[0] = (TProfile*)QAThisRun->FindObject("fProfileQxAQxCCentBfRCThisRun");
  pQxAQyCCentThisRun[0] = (TProfile*)QAThisRun->FindObject("fProfileQxAQyCCentBfRCThisRun");
  pQyAQxCCentThisRun[0] = (TProfile*)QAThisRun->FindObject("fProfileQyAQxCCentBfRCThisRun");
  pQyAQyCCentThisRun[0] = (TProfile*)QAThisRun->FindObject("fProfileQyAQyCCentBfRCThisRun");
  hPsiACentThisRun[0] = (TH2D*)QAThisRun->FindObject("fHist2DPsiACentBfRC");
  hPsiCCentThisRun[0] = (TH2D*)QAThisRun->FindObject("fHist2DPsiCCentBfRC");

  pQxCCentThisRun[1] = (TProfile*)QAThisRun->FindObject("fProfileQxCCentAfRC");
  pQyCCentThisRun[1] = (TProfile*)QAThisRun->FindObject("fProfileQyCCentAfRC");
  pQxACentThisRun[1] = (TProfile*)QAThisRun->FindObject("fProfileQxACentAfRC");
  pQyACentThisRun[1] = (TProfile*)QAThisRun->FindObject("fProfileQyACentAfRC");

  pQxAQxCCentThisRun[1] = (TProfile*)QAThisRun->FindObject("fProfileQxAQxCCentAfRCThisRun");
  pQxAQyCCentThisRun[1] = (TProfile*)QAThisRun->FindObject("fProfileQxAQyCCentAfRCThisRun");
  pQyAQxCCentThisRun[1] = (TProfile*)QAThisRun->FindObject("fProfileQyAQxCCentAfRCThisRun");
  pQyAQyCCentThisRun[1] = (TProfile*)QAThisRun->FindObject("fProfileQyAQyCCentAfRCThisRun");
  hPsiACentThisRun[1] = (TH2D*)QAThisRun->FindObject("fHist2DPsiACentAfRC");
  hPsiCCentThisRun[1] = (TH2D*)QAThisRun->FindObject("fHist2DPsiCCentAfRC");

  int ci[4];
  TColor* color[4];
  ci[0] = TColor::GetFreeColorIndex();
  color[0] = new TColor(ci[0],   0/255.,  24/255., 113/255.);
  ci[1] = TColor::GetFreeColorIndex();
  color[1] = new TColor(ci[1], 255/255.,  88/255.,  93/255.);
  ci[2] = TColor::GetFreeColorIndex();
  color[2] = new TColor(ci[2], 255/255., 181/255.,  73/255.);
  ci[3] = TColor::GetFreeColorIndex();
  color[3] = new TColor(ci[3], 65/255.,  182/255., 230/255.);

  TCanvas* cGlobalQA = new TCanvas("","",1000,800);
  cGlobalQA->Divide(2,2);
  cGlobalQA->cd(1);
  hCent->Draw();
  cGlobalQA->cd(2);
  hCentCorrAfCut->Draw("cont0");
  cGlobalQA->cd(3);
  hVxVy->Draw("cont0");
  cGlobalQA->cd(4);
  hVz->Draw();
  cGlobalQA->SaveAs("plots/GlobalQA.png");

  TCanvas* cGE = new TCanvas("","",1000,800);
  cGE->Divide(2,2);
  cGE->cd(1);
  pZNCTowerEnergyThisRun[0]->Draw();
  pZNCTowerEnergyThisRun[0]->SetMarkerColor(ci[0]);
  pZNCTowerEnergyThisRun[0]->SetLineColor(ci[0]);
  pZNCTowerEnergyThisRun[0]->SetLineStyle(9);
  pZNCTowerEnergyThisRun[0]->SetLineWidth(2);
  cGE->cd(2);
  pZNATowerEnergyThisRun[0]->Draw();
  pZNATowerEnergyThisRun[0]->SetMarkerColor(ci[1]);
  pZNATowerEnergyThisRun[0]->SetLineColor(ci[1]);
  pZNATowerEnergyThisRun[0]->SetLineStyle(9);
  pZNATowerEnergyThisRun[0]->SetLineWidth(2);

  cGE->cd(3);
  pZNCTowerEnergyThisRun[1]->Draw();
  pZNCTowerEnergyThisRun[1]->SetMarkerColor(ci[0]);
  pZNCTowerEnergyThisRun[1]->SetLineColor(ci[0]);
  pZNCTowerEnergyThisRun[1]->SetLineWidth(2);
  cGE->cd(4);
  pZNATowerEnergyThisRun[1]->Draw();
  pZNATowerEnergyThisRun[1]->SetMarkerColor(ci[1]);
  pZNATowerEnergyThisRun[1]->SetLineColor(ci[1]);
  pZNATowerEnergyThisRun[1]->SetLineWidth(2);
  cGE->SaveAs("plots/GE.png");



  TCanvas* cQ = new TCanvas("","",1000,800);
  cQ->Divide(2,2);
  cQ->cd(1);
  pQxCCentThisRun[0]->SetMarkerColor(ci[0]);
  pQxCCentThisRun[0]->SetLineColor(ci[0]);
  pQxCCentThisRun[0]->Draw("SAME");
  pQxCCentThisRun[1]->SetMarkerColor(ci[1]);
  pQxCCentThisRun[1]->SetLineColor(ci[1]);
  pQxCCentThisRun[1]->Draw("SAME");

  cQ->cd(2);
  pQyCCentThisRun[0]->SetMarkerColor(ci[0]);
  pQyCCentThisRun[0]->SetLineColor(ci[0]);
  pQyCCentThisRun[0]->Draw("SAME");
  pQyCCentThisRun[1]->SetMarkerColor(ci[1]);
  pQyCCentThisRun[1]->SetLineColor(ci[1]);
  pQyCCentThisRun[1]->Draw("SAME");
  cQ->cd(3);
  pQxACentThisRun[0]->SetMarkerColor(ci[0]);
  pQxACentThisRun[0]->SetLineColor(ci[0]);
  pQxACentThisRun[0]->Draw("SAME");
  pQxACentThisRun[1]->SetMarkerColor(ci[1]);
  pQxACentThisRun[1]->SetLineColor(ci[1]);
  pQxACentThisRun[1]->Draw("SAME");
  cQ->cd(4);
  pQyACentThisRun[0]->Draw("SAME");
  pQyACentThisRun[0]->SetMarkerColor(ci[0]);
  pQyACentThisRun[0]->SetLineColor(ci[0]);
  pQyACentThisRun[1]->Draw("SAME");
  pQyACentThisRun[1]->SetMarkerColor(ci[1]);
  pQyACentThisRun[1]->SetLineColor(ci[1]);
  cQ->SaveAs("plots/Q.png");


  TCanvas* cQQ = new TCanvas("","",1000,800);
  cQQ->Divide(2,2);
  cQQ->cd(1);
  pQxAQxCCentThisRun[0]->Draw("SAME");
  pQxAQxCCentThisRun[1]->Draw("SAME");
  cQQ->cd(2);
  pQxAQyCCentThisRun[0]->Draw("SAME");
  pQxAQyCCentThisRun[1]->Draw("SAME");
  cQQ->cd(3);
  pQyAQxCCentThisRun[0]->Draw("SAME");
  pQyAQxCCentThisRun[1]->Draw("SAME");
  cQQ->cd(4);
  pQyAQyCCentThisRun[0]->Draw("SAME");
  pQyAQyCCentThisRun[1]->Draw("SAME");
  cQQ->SaveAs("plots/QQ.png");

  TCanvas* cQQTot = new TCanvas("","",1000,800);
  pQxAQxCCentThisRun[0]->Rebin(2);
  pQxAQxCCentThisRun[1]->Rebin(2);
  pQxAQyCCentThisRun[0]->Rebin(2);
  pQxAQyCCentThisRun[1]->Rebin(2);
  pQyAQxCCentThisRun[0]->Rebin(2);
  pQyAQxCCentThisRun[1]->Rebin(2);
  pQyAQyCCentThisRun[0]->Rebin(2);
  pQyAQyCCentThisRun[1]->Rebin(2);

  pQxAQxCCentThisRun[0]->SetMarkerColor(ci[0]);
  pQxAQxCCentThisRun[0]->SetLineColor(ci[0]);
  pQxAQxCCentThisRun[1]->SetMarkerColor(ci[1]);
  pQxAQxCCentThisRun[1]->SetLineColor(ci[1]);

  pQxAQyCCentThisRun[0]->SetMarkerColor(ci[0]);
  pQxAQyCCentThisRun[0]->SetLineColor(ci[0]);
  pQxAQyCCentThisRun[1]->SetMarkerColor(ci[1]);
  pQxAQyCCentThisRun[1]->SetLineColor(ci[1]);

  pQyAQxCCentThisRun[0]->SetMarkerColor(ci[0]);
  pQyAQxCCentThisRun[0]->SetLineColor(ci[0]);
  pQyAQxCCentThisRun[1]->SetMarkerColor(ci[1]);
  pQyAQxCCentThisRun[1]->SetLineColor(ci[1]);

  pQyAQyCCentThisRun[0]->SetMarkerColor(ci[0]);
  pQyAQyCCentThisRun[0]->SetLineColor(ci[0]);
  pQyAQyCCentThisRun[1]->SetMarkerColor(ci[1]);
  pQyAQyCCentThisRun[1]->SetLineColor(ci[1]);


  pQxAQxCCentThisRun[0]->Draw("SAME");
  pQxAQxCCentThisRun[1]->Draw("SAME");

  pQxAQyCCentThisRun[0]->Draw("SAME");
  pQxAQyCCentThisRun[1]->Draw("SAME");

  pQyAQxCCentThisRun[0]->Draw("SAME");
  pQyAQxCCentThisRun[1]->Draw("SAME");

  pQyAQyCCentThisRun[0]->Draw("SAME");
  pQyAQyCCentThisRun[1]->Draw("SAME");
  cQQTot->SaveAs("plots/QQTot.png");

  TCanvas* cPlanelego = new TCanvas("","",1000,800);
  cPlanelego->Divide(2,2);
  cPlanelego->cd(1);
  hPsiCCentThisRun[0]->Draw("lego");
  cPlanelego->cd(2);
  hPsiACentThisRun[0]->Draw("lego");
  cPlanelego->cd(3);
  hPsiCCentThisRun[1]->Draw("lego");
  cPlanelego->cd(4);
  hPsiACentThisRun[1]->Draw("lego");
  cPlanelego->SaveAs("plots/Planelego.png");

  TCanvas* cPlane = new TCanvas("","",1000,800);
  cPlane->Divide(2,2);
  cPlane->cd(1);
  hPsiCCentThisRun[0]->Draw("cont4");
  cPlane->cd(2);
  hPsiACentThisRun[0]->Draw("cont4");
  cPlane->cd(3);
  hPsiCCentThisRun[1]->Draw("cont4");
  cPlane->cd(4);
  hPsiACentThisRun[1]->Draw("cont4");
  cPlane->SaveAs("plots/Plane.png");
}