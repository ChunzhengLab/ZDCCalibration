#include "TFile.h"
#include "TProfile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TColor.h"

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
  TList* listQATot = (TList*)listTot->FindObject("QA Total");

  TProfile* pQxAQxCCentTot[2];
  TProfile* pQxAQyCCentTot[2];
  TProfile* pQyAQxCCentTot[2];
  TProfile* pQyAQyCCentTot[2];
  for(int i = 0; i<2; i++) pQxAQxCCentTot[i] = nullptr;
  for(int i = 0; i<2; i++) pQxAQyCCentTot[i] = nullptr;
  for(int i = 0; i<2; i++) pQyAQxCCentTot[i] = nullptr;
  for(int i = 0; i<2; i++) pQyAQyCCentTot[i] = nullptr;
  pQxAQxCCentTot[0] = (TProfile*)listQATot->FindObject("fProfileQxAQxCCentTotBfRC");
  pQxAQyCCentTot[0] = (TProfile*)listQATot->FindObject("fProfileQxAQyCCentTotBfRC");
  pQyAQxCCentTot[0] = (TProfile*)listQATot->FindObject("fProfileQyAQxCCentTotBfRC");
  pQyAQyCCentTot[0] = (TProfile*)listQATot->FindObject("fProfileQyAQyCCentTotBfRC");
  pQxAQxCCentTot[1] = (TProfile*)listQATot->FindObject("fProfileQxAQxCCentTotAfRC");
  pQxAQyCCentTot[1] = (TProfile*)listQATot->FindObject("fProfileQxAQyCCentTotAfRC");
  pQyAQxCCentTot[1] = (TProfile*)listQATot->FindObject("fProfileQyAQxCCentTotAfRC");
  pQyAQyCCentTot[1] = (TProfile*)listQATot->FindObject("fProfileQyAQyCCentTotAfRC");

  TList* listThisRun = (TList*)listTot->FindObject("137231");
  TProfile* pZNCTowerEnergyThisRun[2];
  TProfile* pZNATowerEnergyThisRun[2];
  for(int i = 0; i<2; i++) pZNCTowerEnergyThisRun[i] = nullptr;
  for(int i = 0; i<2; i++) pZNATowerEnergyThisRun[i] = nullptr;
  pZNCTowerEnergyThisRun[0] = (TProfile*)listThisRun->FindObject("profileZNCTowerEnergy");
  pZNATowerEnergyThisRun[0] = (TProfile*)listThisRun->FindObject("profileZNATowerEnergy");
  pZNCTowerEnergyThisRun[1] = (TProfile*)listThisRun->FindObject("profileZNCTowerEnergyAfGE");
  pZNATowerEnergyThisRun[1] = (TProfile*)listThisRun->FindObject("profileZNATowerEnergyAfGE");


  TList* QAThisRun = (TList*)listThisRun->FindObject("QA");
  TProfile* pQxCCentThisRun[2];
  TProfile* pQyCCentThisRun[2];
  TProfile* pQxACentThisRun[2];
  TProfile* pQyACentThisRun[2];
  TH2D*     hPsiACentThisRun[2];
  TH2D*     hPsiCCentThisRun[2];

  for(int i = 0; i<2; i++) pQxCCentThisRun[i] = nullptr;
  for(int i = 0; i<2; i++) pQyCCentThisRun[i] = nullptr;
  for(int i = 0; i<2; i++) pQxACentThisRun[i] = nullptr;
  for(int i = 0; i<2; i++) pQyACentThisRun[i] = nullptr;
  for(int i = 0; i<2; i++) hPsiACentThisRun[i] = nullptr;
  for(int i = 0; i<2; i++) hPsiCCentThisRun[i] = nullptr;

  pQxCCentThisRun[0] = (TProfile*)QAThisRun->FindObject("fProfileQxCCentBfRC");
  pQyCCentThisRun[0] = (TProfile*)QAThisRun->FindObject("fProfileQyCCentBfRC");
  pQxACentThisRun[0] = (TProfile*)QAThisRun->FindObject("fProfileQxACentBfRC");
  pQyACentThisRun[0] = (TProfile*)QAThisRun->FindObject("fProfileQyACentBfRC");

  hPsiACentThisRun[0] = (TH2D*)QAThisRun->FindObject("fHist2DPsiACentBfRC");
  hPsiCCentThisRun[0] = (TH2D*)QAThisRun->FindObject("fHist2DPsiCCentBfRC");

  pQxCCentThisRun[1] = (TProfile*)QAThisRun->FindObject("fProfileQxCCentAfRC");
  pQyCCentThisRun[1] = (TProfile*)QAThisRun->FindObject("fProfileQyCCentAfRC");
  pQxACentThisRun[1] = (TProfile*)QAThisRun->FindObject("fProfileQxACentAfRC");
  pQyACentThisRun[1] = (TProfile*)QAThisRun->FindObject("fProfileQyACentAfRC");

  hPsiACentThisRun[1] = (TH2D*)QAThisRun->FindObject("fHist2DPsiACentAfRC");
  hPsiCCentThisRun[1] = (TH2D*)QAThisRun->FindObject("fHist2DPsiCCentAfRC");

  TH1D* hPsiACent[2][10];
  TH1D* hPsiCCent[2][10];
  for(int i = 0; i<2; i++) for (int iCent = 0; iCent < 10; iCent++) {hPsiACent[i][iCent] = nullptr; hPsiCCent[i][iCent] = nullptr;}
  for (int i = 0; i < 2; i++) {
  for (int iCent = 0; iCent < 10; iCent++) {      
      hPsiACent[i][iCent] = (TH1D*)hPsiACentThisRun[i]->ProjectionY(Form("PsiAcent%d_%d",iCent,i),iCent+1,iCent+1);
      hPsiCCent[i][iCent] = (TH1D*)hPsiCCentThisRun[i]->ProjectionY(Form("PsiCcent%d_%d",iCent,i),iCent+1,iCent+1);
      hPsiACent[i][iCent]->GetYaxis()->SetRangeUser(0,450);
      hPsiCCent[i][iCent]->GetYaxis()->SetRangeUser(0,450);
    }
  }

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
  pZNCTowerEnergyThisRun[0]->SetLineStyle(6);
  pZNCTowerEnergyThisRun[0]->SetLineWidth(2);
  cGE->cd(2);
  pZNATowerEnergyThisRun[0]->Draw();
  pZNATowerEnergyThisRun[0]->SetMarkerColor(ci[1]);
  pZNATowerEnergyThisRun[0]->SetLineColor(ci[1]);
  pZNATowerEnergyThisRun[0]->SetLineStyle(6);
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


  TCanvas* cQQTot = new TCanvas("","",1000,800);
  pQxAQxCCentTot[0]->Rebin(4);
  pQxAQxCCentTot[1]->Rebin(4);
  pQxAQyCCentTot[0]->Rebin(4);
  pQxAQyCCentTot[1]->Rebin(4);
  pQyAQxCCentTot[0]->Rebin(4);
  pQyAQxCCentTot[1]->Rebin(4);
  pQyAQyCCentTot[0]->Rebin(4);
  pQyAQyCCentTot[1]->Rebin(4);
  pQxAQxCCentTot[0]->SetMarkerColor(ci[0]);
  pQxAQxCCentTot[0]->SetLineColor(ci[0]);
  pQxAQxCCentTot[0]->SetMarkerStyle(kOpenSquare);
  pQxAQxCCentTot[1]->SetMarkerColor(ci[0]);
  pQxAQxCCentTot[1]->SetLineColor(ci[0]);
  pQxAQxCCentTot[1]->SetMarkerStyle(kFullSquare);

  pQxAQyCCentTot[0]->SetMarkerColor(ci[1]);
  pQxAQyCCentTot[0]->SetLineColor(ci[1]);
  pQxAQyCCentTot[0]->SetMarkerStyle(kOpenSquare);
  pQxAQyCCentTot[1]->SetMarkerColor(ci[1]);
  pQxAQyCCentTot[1]->SetLineColor(ci[1]);
  pQxAQyCCentTot[1]->SetMarkerStyle(kFullSquare);


  pQyAQxCCentTot[0]->SetMarkerColor(ci[2]);
  pQyAQxCCentTot[0]->SetLineColor(ci[2]);
  pQyAQxCCentTot[0]->SetMarkerStyle(kOpenSquare);
  pQyAQxCCentTot[1]->SetMarkerColor(ci[2]);
  pQyAQxCCentTot[1]->SetLineColor(ci[2]);
  pQyAQxCCentTot[1]->SetMarkerStyle(kFullSquare);


  pQyAQyCCentTot[0]->SetMarkerColor(ci[3]);
  pQyAQyCCentTot[0]->SetLineColor(ci[3]);
  pQyAQyCCentTot[0]->SetMarkerStyle(kOpenSquare);
  pQyAQyCCentTot[1]->SetMarkerColor(ci[3]);
  pQyAQyCCentTot[1]->SetLineColor(ci[3]);
  pQyAQyCCentTot[1]->SetMarkerStyle(kFullSquare);
  for(int i = 0; i<2; i++) pQxAQxCCentTot[i] ->GetYaxis()->SetRangeUser(-0.025,0.025);
  for(int i = 0; i<2; i++) pQxAQyCCentTot[i] ->GetYaxis()->SetRangeUser(-0.025,0.025);
  for(int i = 0; i<2; i++) pQyAQxCCentTot[i] ->GetYaxis()->SetRangeUser(-0.025,0.025);
  for(int i = 0; i<2; i++) pQyAQyCCentTot[i] ->GetYaxis()->SetRangeUser(-0.025,0.025);

  pQxAQxCCentTot[0]->Draw("SAME");
  pQxAQxCCentTot[1]->Draw("SAME");

  pQxAQyCCentTot[0]->Draw("SAME");
  pQxAQyCCentTot[1]->Draw("SAME");

  pQyAQxCCentTot[0]->Draw("SAME");
  pQyAQxCCentTot[1]->Draw("SAME");

  pQyAQyCCentTot[0]->Draw("SAME");
  pQyAQyCCentTot[1]->Draw("SAME");

  cQQTot->SaveAs("plots/QQTot.png");

  TCanvas* cQQ = new TCanvas("","",1000,800);
  cQQ->Divide(2,2);
  cQQ->cd(1);
  pQxAQxCCentTot[0]->Draw("SAME");
  pQxAQxCCentTot[1]->Draw("SAME");
  cQQ->cd(2);
  pQxAQyCCentTot[0]->Draw("SAME");
  pQxAQyCCentTot[1]->Draw("SAME");
  cQQ->cd(3);
  pQyAQxCCentTot[0]->Draw("SAME");
  pQyAQxCCentTot[1]->Draw("SAME");
  cQQ->cd(4);
  pQyAQyCCentTot[0]->Draw("SAME");
  pQyAQyCCentTot[1]->Draw("SAME");
  cQQ->SaveAs("plots/QQ.png");

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

  TCanvas* cPlaneCCentBf = new TCanvas("","",1000,800);
  cPlaneCCentBf->Divide(3,3);
  for (int i = 0; i < 9; i++)
  {
    cPlaneCCentBf->cd(i+1);
    hPsiCCent[0][i]->Draw();
  }
  cPlaneCCentBf->SaveAs("plots/PlaneCBf.png");

  TCanvas* cPlaneCCentAf = new TCanvas("","",1000,800);
  cPlaneCCentAf->Divide(3,3);
  for (int i = 0; i < 9; i++)
  {
    cPlaneCCentAf->cd(i+1);
    hPsiCCent[1][i]->Draw();
  }
  cPlaneCCentAf->SaveAs("plots/PlaneCAf.png");

  TCanvas* cPlaneACentBf = new TCanvas("","",1000,800);
  cPlaneACentBf->Divide(3,3);
  for (int i = 0; i < 9; i++)
  {
    cPlaneACentBf->cd(i+1);
    hPsiACent[0][i]->Draw();
  }
  cPlaneACentBf->SaveAs("plots/PlaneABf.png");

  TCanvas* cPlaneACentAf = new TCanvas("","",1000,800);
  cPlaneACentAf->Divide(3,3);
  for (int i = 0; i < 9; i++)
  {
    cPlaneACentAf->cd(i+1);
    hPsiACent[1][i]->Draw();
  }
  cPlaneACentAf->SaveAs("plots/PlaneAAf.png");
}