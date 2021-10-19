// include the header of your analysis task here! for classes already compiled by aliBuild,
// precompiled header files (with extension pcm) are available, so that you do not need to
// specify includes for those. for your own task however, you (probably) have not generated a
// pcm file, so we need to include it explicitly
#include "AliAnalysisTaskZDCCalibration.h"

void runAnalysis()
{
    // set if you want to run the analysis locally (kTRUE), or on grid (kFALSE)
    Bool_t local = kFALSE;
    // if you run on grid, specify test mode (kTRUE) or full grid model (kFALSE)
    Bool_t gridTest = kFALSE;
    
    // since we will compile a class, tell root where to look for headers  
#if !defined (__CINT__) || defined (__CLING__)
    gInterpreter->ProcessLine(".include $ROOTSYS/include");
    gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
#else
    gROOT->ProcessLine(".include $ROOTSYS/include");
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
#endif
     
    // create the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisTaskExample");
    AliAODInputHandler *aodH = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodH);

    // load the macro and add the task
    TMacro PIDadd(gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C"));
    AliAnalysisTaskPIDResponse *PIDresponseTask = reinterpret_cast<AliAnalysisTaskPIDResponse *>(PIDadd.Exec());

    TMacro multSelection(gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"));
    AliMultSelectionTask *multSelectionTask = reinterpret_cast<AliMultSelectionTask *>(multSelection.Exec());

    // compile the class and load the add task macro
    // here we have to differentiate between using the just-in-time compiler
    // from root6, or the interpreter of root5
#if !defined (__CINT__) || defined (__CLING__)
    gInterpreter->LoadMacro("AliAnalysisTaskZDCCalibration.cxx++g");
    AliAnalysisTaskZDCCalibration *task = reinterpret_cast<AliAnalysisTaskZDCCalibration*>(gInterpreter->ExecuteMacro("AddTaskZDCCalibration.C"));
#else
    gROOT->LoadMacro("AliAnalysisTaskZDCCalibration.cxx++g");
    gROOT->LoadMacro("AddTaskZDCCalibration.C");
    AliAnalysisTaskZDCCalibration *task = AddTaskZDCCalibration();
#endif


    if(!mgr->InitAnalysis()) return;
    mgr->SetDebugLevel(2);
    mgr->PrintStatus();
    mgr->SetUseProgressBar(1, 25);

    if(local) {
        // if you want to run locally, we need to define some input
        TChain* chain = new TChain("aodTree");
        // add a few files to the chain (change this so that your local files are added)
        chain->Add("AliAOD.root");
        // start the analysis locally, reading the events from the tchain
        mgr->StartAnalysis("local", chain);
    } else {
        // if we want to run on grid, we create and configure the plugin
        AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
        // also specify the include (header) paths on grid
        alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
        // make sure your source files get copied to grid
        alienHandler->SetAdditionalLibs("AliAnalysisTaskZDCCalibration.cxx AliAnalysisTaskZDCCalibration.h");
        alienHandler->SetAnalysisSource("AliAnalysisTaskZDCCalibration.cxx");
        // select the aliphysics version. all other packages
        // are LOADED AUTOMATICALLY!
        alienHandler->SetAliPhysicsVersion("vAN-20181028_ROOT6-1");
        // set the Alien API version
        alienHandler->SetAPIVersion("V1.1x");

        // //10h
        // // select the input data
        // alienHandler->SetGridDataDir("/alice/data/2010/LHC10h");
        // alienHandler->SetDataPattern("ESDs/pass2/AOD160/*/AliAOD.root");
        // // MC has no prefix, data has prefix 000
        // alienHandler->SetRunPrefix("000");
        // // runnumber
        // alienHandler->AddRunNumber(139510);

        //11h
        alienHandler->SetGridDataDir("/alice/data/2011/LHC11h_2");
        alienHandler->SetDataPattern("*ESDs/pass2/AOD145/*AOD.root");
        alienHandler->SetRunPrefix("000");
        alienHandler->AddRunNumber(170387);
        alienHandler->AddRunNumber(170040);
        alienHandler->AddRunNumber(170268);
        alienHandler->AddRunNumber(170228);
        alienHandler->AddRunNumber(170207);
        alienHandler->AddRunNumber(169838);
        alienHandler->AddRunNumber(170159);
        alienHandler->AddRunNumber(170204);
        alienHandler->AddRunNumber(170311);
        alienHandler->AddRunNumber(170084);
    
        // number of files per subjob
        alienHandler->SetSplitMaxInputFileNumber(40);
        alienHandler->SetExecutable("myTask.sh");
        // specify how many seconds your job may take
        alienHandler->SetTTL(10000);
        alienHandler->SetJDLName("myTask.jdl");

        alienHandler->SetOutputToRunNo(kTRUE);
        alienHandler->SetKeepLogs(kTRUE);
        // merging: run with kTRUE to merge on grid
        // after re-running the jobs in SetRunMode("terminate") 
        // (see below) mode, set SetMergeViaJDL(kFALSE) 
        // to collect final results
        alienHandler->SetMaxMergeStages(1);
        alienHandler->SetMergeViaJDL(kTRUE);

        // define the output folders
        alienHandler->SetGridWorkingDir("ZDCCalibration_fillHistForGE");
        //alienHandler->SetGridWorkingDir("ZDCCalibration_fillHistForRC");
        //alienHandler->SetGridWorkingDir("ZDCCalibration_finishZDCCali");
        alienHandler->SetGridOutputDir("Output");

        // connect the alien plugin to the manager
        mgr->SetGridHandler(alienHandler);
        if(gridTest) {
            // speficy on how many files you want to run
            alienHandler->SetNtestFiles(1);
            // and launch the analysis
            alienHandler->SetRunMode("test");
            mgr->StartAnalysis("grid");
        } else {
            // else launch the full grid analysis
            //alienHandler->SetRunMode("terminate");
            alienHandler->SetRunMode("full");
            mgr->StartAnalysis("grid");
        }
    }
}
