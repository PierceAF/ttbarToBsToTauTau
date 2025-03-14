// SKIM - 1: save skimmed ntuple, 0: run on already skimmed ntuple
#define SYST  0
#define SKIM  0

#ifndef PLOTTER
#include "include/treestream.h"
#endif

#if SKIM == 0
#include "Plotting.h" // Specify here the implementations for your Analysis
#else 
#include "Analysis_Skim.h" // Specify here the implementations for your Analysis
#endif

struct settings {

  //-----------------------------------------------------------------------------
  // -- Constants
  //-----------------------------------------------------------------------------
  settings() :
    runOnSkim                ( 1-SKIM ),
    saveSkimmedNtuple        ( SKIM   ),
    doTopPtReweighting       ( true   ),
    doISRReweighting         ( true   ),
    doPileupReweighting      ( true   ),
    doAK8JetPtRescaling      ( false  ),
    applySmearing            ( true   ),
    applyScaleFactors        ( true   ),
    nSigmaScaleFactors       ( 17   ), // Count the number of sigmas you use in ScaleFactors.h
    varySystematics          ( SYST   ),
    systematicsFileName      ( "systematics/2024_12_18_1SigmaUpDown_NoPdf.txt" ),
    useJSON                  ( false  ), // by default: no need to apply, but can be useful if some lumisections need to be excluded additionally
    jsonFileName             ( "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/"
			       "Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt" ),
    useXSecFileForBkg        ( true   ) {}; // true: use file below, false: use value in the ntuple (evt_XSec)
  ~settings(){};

  const bool runOnSkim;
  const bool saveSkimmedNtuple;
  const bool doTopPtReweighting;
  const bool doISRReweighting;
  const bool doPileupReweighting;
  const bool doAK8JetPtRescaling;
  const bool applySmearing;
  const bool applyScaleFactors;
  const int  nSigmaScaleFactors;
  const bool varySystematics;
  const std::string systematicsFileName;
  const bool useJSON;
  const std::string jsonFileName;
  const bool useXSecFileForBkg;

} settings;
