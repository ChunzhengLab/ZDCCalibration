#include <iostream>
#include "TROOT.h"
#include "TGrid.h"
#include "TFile.h"
#include "AliAnalysisTaskZDCCalibration.h"
#include "AliAnalysisManager.h"
using namespace std;

AliAnalysisTaskZDCCalibration* AddTaskZDCCalibration(TString name = "name")
{
    // get the manager via the static access member. since it's static, you don't need
    // to create an instance of the class here to call the function
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        return 0x0;
    }
    // get the input event handler, again via a static method. 
    // this handler is part of the managing system and feeds events
    // to your task
    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }
    // by default, a file is open for writing. here, we get the filename
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":myTask";      // create a subfolder in the file
    // now we create an instance of your task
    AliAnalysisTaskZDCCalibration* task = new AliAnalysisTaskZDCCalibration(name.Data());   
    if(!task) return 0x0;
    TString dataSet("10h");
    if(dataSet.Contains("10h")) task->SelectCollisionCandidates(AliVEvent::kMB);
    if(dataSet.Contains("11h")) task->SelectCollisionCandidates(AliVEvent::kSemiCentral | AliVEvent::kCentral | AliVEvent::kMB);
    bool isFirstFillVetex = kTRUE;
    task->FirstFillHistVetex(isFirstFillVetex);
    task->SetDataSet(dataSet);

    // add your task to the manager
    mgr->AddTask(task);

    if (!isFirstFillVetex)
    {
      TGrid::Connect("alien://");
      // add list for ZDC Calibration
      TString ZDCCalibrationFileName = "alien:///alice/cern.ch/user/c/chunzhen/ZDCCalibration.root";
      TFile* ZDCCalibrationFile = TFile::Open(ZDCCalibrationFileName,"READ");
      if(!ZDCCalibrationFile) {
          cout << "ERROR: ZDC Calibration File is not found!" << endl;
          exit(1);
      }
      gROOT->cd();
      TDirectory* inputDir = ZDCCalibrationFile->GetDirectory("myTask");
      TList* ZDCCalibrationList = nullptr;
      inputDir->GetObject("Output",ZDCCalibrationList);
      if(ZDCCalibrationList) {
          task->SetZDCCalibrationList(ZDCCalibrationList);
          cout << "ZDC Calibration file: set! (from " <<  ZDCCalibrationFileName.Data() << ")" << endl;
      } else {
          cout << "ERROR: ZDC Calibration file: TList not found!" << endl;
          exit(1);
      }
      delete ZDCCalibrationFile;
    }

    // your task needs input: here we connect the manager to your task
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    // same for the output
    mgr->ConnectOutput(task,1,mgr->CreateContainer("Output", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    // in the end, this macro returns a pointer to your task. this will be convenient later on
    // when you will run your analysis in an analysis train on grid
    return task;
}
