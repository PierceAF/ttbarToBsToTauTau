#ifndef SCALEFACTORS_H
#define SCALEFACTORS_H

// Private headers
#include "eventBuffer.h" // make sure to set to same as in tnm.h
#include "Variables.h"
#include "EventSelections.h"

// 3rd party headers
#include "tnm.h" // for getplot
#include "BTagCalibrationStandalone.cpp" // From BTAG POG
#include "TFile.h"
#include "TH2.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TProfile.h"
#include "TStopwatch.h"
#include "TString.h"

// common libraries
#include <iostream>
#include <vector>

class ScaleFactors
{

public:
  typedef EventSelections::Regions Region;
  ScaleFactors(Variables& var) :
    v(var) {

    debug = 0;
    t_s0 = 0;
    t_s1 = 0;
    t_s2 = 0;
    t_s3 = 0;
    t_s4 = 0;
    t_s5 = 0;
    t_s6 = 0;
    t_s7 = 0;
    t_s8 = 0;
    t_s9 = 0;
    if (debug) {
      sw_s0 = new TStopwatch;
      sw_s1 = new TStopwatch;
      sw_s2 = new TStopwatch;
      sw_s3 = new TStopwatch;
      sw_s4 = new TStopwatch;
      sw_s5 = new TStopwatch;
      sw_s6 = new TStopwatch;
      sw_s7 = new TStopwatch;
      sw_s8 = new TStopwatch;
      sw_s9 = new TStopwatch;
      sw_s0->Reset();
      sw_s1->Reset();
      sw_s2->Reset();
      sw_s3->Reset();
      sw_s4->Reset();
      sw_s5->Reset();
      sw_s6->Reset();
      sw_s7->Reset();
      sw_s8->Reset();
      sw_s9->Reset();
    }
    
    // Select which scale/correction factors to use for each region
    scale_factors.resize(magic_enum::enum_count<Region>());
    
    scale_factors[Region::Pre_Lep].push_back(&sf_ele_medium);
    scale_factors[Region::Pre_Lep].push_back(&sf_muon_medium);
    scale_factors[Region::Pre_Lep].push_back(&sf_btag_medium);

    scale_factors[Region::Pre_e].push_back(&sf_ele_medium);
    scale_factors[Region::Pre_e].push_back(&sf_btag_medium);

    scale_factors[Region::Pre_u].push_back(&sf_muon_medium);
    scale_factors[Region::Pre_u].push_back(&sf_btag_medium);

    scale_factors[Region::Lep].push_back(&sf_ele_medium);
    scale_factors[Region::Lep].push_back(&sf_muon_medium);
    scale_factors[Region::Lep].push_back(&sf_btag_medium);

    scale_factors[Region::Lep_e].push_back(&sf_ele_medium);
    scale_factors[Region::Lep_e].push_back(&sf_btag_medium);

    scale_factors[Region::Lep_u].push_back(&sf_muon_medium);
    scale_factors[Region::Lep_u].push_back(&sf_btag_medium);

    scale_factors[Region::Pre_uu].push_back(&sf_muon_medium);
    scale_factors[Region::Pre_uu].push_back(&sf_btag_medium);

    scale_factors[Region::Pre_ee].push_back(&sf_ele_medium);
    scale_factors[Region::Pre_ee].push_back(&sf_btag_medium);

    scale_factors[Region::Pre_eu].push_back(&sf_ele_medium);
    scale_factors[Region::Pre_eu].push_back(&sf_muon_medium);
    scale_factors[Region::Pre_eu].push_back(&sf_btag_medium);

    scale_factors[Region::DiLep_uu].push_back(&sf_muon_medium);
    scale_factors[Region::DiLep_uu].push_back(&sf_btag_medium);

    scale_factors[Region::DiLep_ee].push_back(&sf_ele_medium);
    scale_factors[Region::DiLep_ee].push_back(&sf_btag_medium);

    scale_factors[Region::DiLep_eu].push_back(&sf_ele_medium);
    scale_factors[Region::DiLep_eu].push_back(&sf_muon_medium);
    scale_factors[Region::DiLep_eu].push_back(&sf_btag_medium);

  }
  ~ScaleFactors() {
    if (debug) {
      double t_s_sum = t_s0 + t_s1 + t_s2 + t_s3 + t_s4
        + t_s5 + t_s6 + t_s7 + t_s8 + t_s9;
      std::cout<<"ScaleFactors Time benchmarks:"<<std::endl;
      std::cout<<"- sf0: "<<t_s0<<"   ("<<(100*t_s0/t_s_sum)<<"%)"<<std::endl;
      std::cout<<"- sf1: "<<t_s1<<"   ("<<(100*t_s1/t_s_sum)<<"%)"<<std::endl;
      std::cout<<"- sf2: "<<t_s2<<"   ("<<(100*t_s2/t_s_sum)<<"%)"<<std::endl;
      std::cout<<"- sf3: "<<t_s3<<"   ("<<(100*t_s3/t_s_sum)<<"%)"<<std::endl;
      std::cout<<"- sf4: "<<t_s4<<"   ("<<(100*t_s4/t_s_sum)<<"%)"<<std::endl;
      std::cout<<"- sf5: "<<t_s5<<"   ("<<(100*t_s5/t_s_sum)<<"%)"<<std::endl;
      std::cout<<"- sf6: "<<t_s6<<"   ("<<(100*t_s6/t_s_sum)<<"%)"<<std::endl;
      std::cout<<"- sf7: "<<t_s7<<"   ("<<(100*t_s7/t_s_sum)<<"%)"<<std::endl;
      std::cout<<"- sf8: "<<t_s8<<"   ("<<(100*t_s8/t_s_sum)<<"%)"<<std::endl;
      std::cout<<"- sf9: "<<t_s9<<"   ("<<(100*t_s9/t_s_sum)<<"%)"<<std::endl;
      std::cout<<"- SUM:     "<<t_s_sum<<std::endl;
    }
  }

  void init_sf_histos();

  void init_input();

  void fill_sf_histos(const bool&, const bool&, const unsigned int&, const double&);

  double calc_pu_jet_sf(const double&);

  std::pair<double, double> calc_b_tagging_sf(const double&, const double&, const double&, const double&, const double&);

  std::tuple<double, double, double> calc_ele_sf(const double&, const double&);

  double calc_pho_sf();

  std::tuple<double, double, double> calc_muon_sf(const double&, const double&);
  
  void apply_scale_factors(const unsigned int&, std::vector<double>&, std::vector<std::vector<double> >&, 
                           const unsigned int&, const std::vector<std::vector<double> >&);

  std::vector<std::vector<double*> > scale_factors;

private:
  
  Variables& v;

  double get_syst_weight_(const double&, const double&, const double&, const double&);

  double get_syst_weight_(const double&, const double&, const double&);

  BTagCalibration* btag_calib_full_;
  BTagCalibration* btag_calib_fast_;
  BTagCalibrationReader* btag_sf_full_loose_;
  BTagCalibrationReader* btag_sf_fast_loose_;
  BTagCalibrationReader* btag_sf_full_medium_;
  BTagCalibrationReader* btag_sf_fast_medium_;

  //_______________________________________________________
  //                Input histograms

  TProfile* eff_btag_b_loose;
  TProfile* eff_btag_c_loose;
  TProfile* eff_btag_l_loose;
  TProfile* eff_btag_b_medium;
  TProfile* eff_btag_c_medium;
  TProfile* eff_btag_l_medium;
  
  TH2F* eff_pu;

	TH2F* eff_full_ele_veto;
	TH2F* eff_full_ele_razor;
	TH2F* eff_full_ele_nonIso;
	TH2F* eff_full_ele_nonIso_B2G ;
	TH2F* unc_full_ele_nonIso_B2G ;
	TH2F* eff_fast_ele_razor;
	TH2F* eff_fast_ele_nonIso;
	TH2F* eff_fast_ele_nonIso_B2G ;
	TH2F* unc_fast_ele_nonIso_B2G ;
	TH2F* eff_full_muon_veto;
	TH2F* eff_full_muon_razor;
	TH2F* eff_full_muon_nonIso;
	TH2F* eff_full_muon_nonIso_B2G;
	TH2F* unc_full_muon_nonIso_B2G;
	TH2F* eff_fast_muon_razor;
	TH2F* eff_fast_muon_nonIso;
	TH2F* eff_fast_muon_nonIso_B2G;
	TH2F* unc_fast_muon_nonIso_B2G;
  TH2F* eff_full_pho_mediumid;
  
  //_______________________________________________________
  //         List of scale and correction factors
  
  double sf_ele_veto, sf_ele_medium, sf_ele_nonIso;
  double sf_pho_medium;
  double sf_muon_veto, sf_muon_medium, sf_muon_nonIso;
  double sf_btag_loose, sf_btag_medium;

  //_______________________________________________________
  //           List of scale factor Histograms

  TH2D* h_btag_eff_b_loose;
  TH2D* h_btag_eff_c_loose;
  TH2D* h_btag_eff_l_loose;
  TH2D* h_btag_eff_b_medium;
  TH2D* h_btag_eff_c_medium;
  TH2D* h_btag_eff_l_medium;

  // Benchmarking
  int debug;
  double t_s0;
  double t_s1;
  double t_s2;
  double t_s3;
  double t_s4;
  double t_s5;
  double t_s6;
  double t_s7;
  double t_s8;
  double t_s9;
  TStopwatch *sw_s0;
  TStopwatch *sw_s1;
  TStopwatch *sw_s2;
  TStopwatch *sw_s3;
  TStopwatch *sw_s4;
  TStopwatch *sw_s5;
  TStopwatch *sw_s6;
  TStopwatch *sw_s7;
  TStopwatch *sw_s8;
  TStopwatch *sw_s9;

  void sw_(TStopwatch*, double&, bool start);
  
};


void ScaleFactors::init_input() {
  // B-tagging
  // Efficiencies (Oct31 - test)
  TFile* f;
  if (v.isWJets)
    f = TFile::Open("btag_eff/WJetsToLNu.root");
  else if (v.isTop)
    f = TFile::Open("btag_eff/TT.root");
  else
    f = TFile::Open("btag_eff/QCD.root");
  if (v.year==2018) {
    eff_btag_b_loose  = ((TH2D*)f->Get("btag_eff_b_loose_2018"))->ProfileX();
    eff_btag_c_loose  = ((TH2D*)f->Get("btag_eff_c_loose_2018"))->ProfileX();
    eff_btag_l_loose  = ((TH2D*)f->Get("btag_eff_l_loose_2018"))->ProfileX();
    eff_btag_b_medium = ((TH2D*)f->Get("btag_eff_b_medium_2018"))->ProfileX();
    eff_btag_c_medium = ((TH2D*)f->Get("btag_eff_c_medium_2018"))->ProfileX();
    eff_btag_l_medium = ((TH2D*)f->Get("btag_eff_l_medium_2018"))->ProfileX();
  } else if (v.year==2017) {
    eff_btag_b_loose  = ((TH2D*)f->Get("btag_eff_b_loose_2017"))->ProfileX();
    eff_btag_c_loose  = ((TH2D*)f->Get("btag_eff_c_loose_2017"))->ProfileX();
    eff_btag_l_loose  = ((TH2D*)f->Get("btag_eff_l_loose_2017"))->ProfileX();
    eff_btag_b_medium = ((TH2D*)f->Get("btag_eff_b_medium_2017"))->ProfileX();
    eff_btag_c_medium = ((TH2D*)f->Get("btag_eff_c_medium_2017"))->ProfileX();
    eff_btag_l_medium = ((TH2D*)f->Get("btag_eff_l_medium_2017"))->ProfileX();
  } else if (v.isAPV) {
    eff_btag_b_loose  = ((TH2D*)f->Get("btag_eff_b_loose_2016APV"))->ProfileX();
    eff_btag_c_loose  = ((TH2D*)f->Get("btag_eff_c_loose_2016APV"))->ProfileX();
    eff_btag_l_loose  = ((TH2D*)f->Get("btag_eff_l_loose_2016APV"))->ProfileX();
    eff_btag_b_medium = ((TH2D*)f->Get("btag_eff_b_medium_2016APV"))->ProfileX();
    eff_btag_c_medium = ((TH2D*)f->Get("btag_eff_c_medium_2016APV"))->ProfileX();
    eff_btag_l_medium = ((TH2D*)f->Get("btag_eff_l_medium_2016APV"))->ProfileX();
  } else {
    eff_btag_b_loose  = ((TH2D*)f->Get("btag_eff_b_loose_2016"))->ProfileX();
    eff_btag_c_loose  = ((TH2D*)f->Get("btag_eff_c_loose_2016"))->ProfileX();
    eff_btag_l_loose  = ((TH2D*)f->Get("btag_eff_l_loose_2016"))->ProfileX();
    eff_btag_b_medium = ((TH2D*)f->Get("btag_eff_b_medium_2016"))->ProfileX();
    eff_btag_c_medium = ((TH2D*)f->Get("btag_eff_c_medium_2016"))->ProfileX();
    eff_btag_l_medium = ((TH2D*)f->Get("btag_eff_l_medium_2016"))->ProfileX();
  }
  eff_btag_b_loose  ->SetDirectory(0);
  eff_btag_c_loose  ->SetDirectory(0);
  eff_btag_l_loose  ->SetDirectory(0);
  eff_btag_b_medium ->SetDirectory(0);
  eff_btag_c_medium ->SetDirectory(0);
  eff_btag_l_medium ->SetDirectory(0);
  f->Close();
  // Moriond17 SFs
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
  if (v.year==2018) {
    btag_calib_full_ = new BTagCalibration("DeepCSV", "scale_factors/btag/wp_deepJet_106XUL18_v2.csv");
  } else if (v.year==2017) {
    btag_calib_full_ = new BTagCalibration("DeepCSV", "scale_factors/btag/wp_deepJet_106XUL17_v3.csv");
  } else {
    if(v.isAPV) btag_calib_full_ = new BTagCalibration("DeepCSV", "scale_factors/btag/wp_deepJet_106XUL16preVFP_v2.csv");
    else        btag_calib_full_ = new BTagCalibration("DeepCSV", "scale_factors/btag/wp_deepJet_106XUL16postVFP_v3.csv");
  }
  // Loose WP
  btag_sf_full_loose_  = new BTagCalibrationReader(BTagEntry::OP_LOOSE, "central", {"up", "down"});
  btag_sf_full_loose_->load(*btag_calib_full_, BTagEntry::FLAV_B,    "comb");
  btag_sf_full_loose_->load(*btag_calib_full_, BTagEntry::FLAV_C,    "comb");
  btag_sf_full_loose_->load(*btag_calib_full_, BTagEntry::FLAV_UDSG, "incl");
  // Medium WP
  btag_sf_full_medium_ = new BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central", {"up", "down"});
  btag_sf_full_medium_->load(*btag_calib_full_, BTagEntry::FLAV_B,    "comb");
  btag_sf_full_medium_->load(*btag_calib_full_, BTagEntry::FLAV_C,    "comb");
  btag_sf_full_medium_->load(*btag_calib_full_, BTagEntry::FLAV_UDSG, "incl");
  // Spring16 FastSim
  if (v.year==2018) {
    btag_calib_fast_ =  new BTagCalibration("deepcsv", "scale_factors/btag/deepcsv_13TEV_18SL_7_5_2019.csv");
  } else if (v.year==2017) {
    // This file needed minor formatting to be readable
    // sed 's;^";;;s; "\;;;;s;"";";g;' scale_factors/btag/fastsim_csvv2_ttbar_26_1_2017.csv
    btag_calib_fast_ =  new BTagCalibration("deepcsv", "scale_factors/btag/deepcsv_13TEV_17SL_18_3_2019.csv");
  } else {
    btag_calib_fast_ =  new BTagCalibration("deepcsv", "scale_factors/btag/deepcsv_13TEV_16SL_18_3_2019.csv");
  }
  // Loose WP
  btag_sf_fast_loose_  = new BTagCalibrationReader(BTagEntry::OP_LOOSE, "central", {"up", "down"});
  btag_sf_fast_loose_->load(*btag_calib_fast_, BTagEntry::FLAV_B,    "fastsim");
  btag_sf_fast_loose_->load(*btag_calib_fast_, BTagEntry::FLAV_C,    "fastsim");
  btag_sf_fast_loose_->load(*btag_calib_fast_, BTagEntry::FLAV_UDSG, "fastsim");
  // Medium WP
  btag_sf_fast_medium_ = new BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central", {"up", "down"});
  btag_sf_fast_medium_->load(*btag_calib_fast_, BTagEntry::FLAV_B,    "fastsim");
  btag_sf_fast_medium_->load(*btag_calib_fast_, BTagEntry::FLAV_C,    "fastsim");
  btag_sf_fast_medium_->load(*btag_calib_fast_, BTagEntry::FLAV_UDSG, "fastsim");

  if (v.year==2018) eff_pu = getplot_TH2F("scale_factors/pu/PUID_106XTraining_ULRun2_EffSFandUncties_v1.root", "h2_eff_sfUL2018_L", "pu");
  else if (v.year==2017) eff_pu = getplot_TH2F("scale_factors/pu/PUID_106XTraining_ULRun2_EffSFandUncties_v1.root", "h2_eff_sfUL2017_L", "pu");
  else if (v.isAPV) eff_pu = getplot_TH2F("scale_factors/pu/PUID_106XTraining_ULRun2_EffSFandUncties_v1.root", "h2_eff_sfUL2016APV_L", "pu");
  else eff_pu = getplot_TH2F("scale_factors/pu/PUID_106XTraining_ULRun2_EffSFandUncties_v1.root", "h2_eff_sfUL2016_L", "pu");

  // Lepton scale factors
  if (v.year==2018) {
    eff_full_ele_veto                 = getplot_TH2F("scale_factors/electron/fullsim_electron_veto_UL2018.root",     						"EGamma_SF2D", 				"ele1");
    eff_full_ele_razor                = getplot_TH2F("scale_factors/electron/fullsim_electron_isolated_UL2018.root",     				"EGamma_SF2D", 				"ele2");
    eff_full_ele_nonIso               = getplot_TH2F("scale_factors/electron/fullsim_electron_nonIsolated_UL2018.root",     		"EGamma_SF2D", 				"ele3");
    eff_full_ele_nonIso_B2G           = getplot_TH2F("scale_factors/electron/fullsim_electron_nonIsolated_B2GCut_UL2018.root",  "h_eta_pT_electron",  "ele4");
    unc_full_ele_nonIso_B2G           = getplot_TH2F("scale_factors/electron/fullsim_electron_nonIsolated_B2GCut_UL2018.root",  "h_eta_pT_electron_stat",  "ele5");
    eff_fast_ele_razor                = getplot_TH2F("scale_factors/electron/fastsim_electron_isolated_UL2018.root",     				"EGamma_SF2D", 				"ele6");
    eff_fast_ele_nonIso               = getplot_TH2F("scale_factors/electron/fastsim_electron_nonIsolated_UL2018.root",     		"EGamma_SF2D", 				"ele7");
    eff_fast_ele_nonIso_B2G           = getplot_TH2F("scale_factors/electron/fastsim_electron_nonIsolated_B2GCut_UL2018.root",  "h_eta_pT_electron",  "ele8");
    unc_fast_ele_nonIso_B2G           = getplot_TH2F("scale_factors/electron/fastsim_electron_nonIsolated_B2GCut_UL2018.root",  "h_eta_pT_electron_stat",  "ele9");
    eff_full_muon_veto                = getplot_TH2F("scale_factors/muon/fullsim_muon_veto_UL2018.root",     						"NUM_RazorVeto_DEN_genTracks_abseta_pt", 				"mu1");
    eff_full_muon_razor               = getplot_TH2F("scale_factors/muon/fullsim_muon_isolated_UL2018.root",     				"NUM_RazorPass_DEN_genTracks_abseta_pt", 				"mu2");
    eff_full_muon_nonIso              = getplot_TH2F("scale_factors/muon/fullsim_muon_nonIsolated_UL2018.root",     		"NUM_RazorNoIso_DEN_genTracks_abseta_pt", 			"mu3");
    eff_full_muon_nonIso_B2G          = getplot_TH2F("scale_factors/muon/fullsim_muon_nonIsolated_B2GCut_UL2018.root",  "h_eta_pT_muon",  															"mu4");
    unc_full_muon_nonIso_B2G          = getplot_TH2F("scale_factors/muon/fullsim_muon_nonIsolated_B2GCut_UL2018.root",  "h_eta_pT_muon_stat",														"mu5");
    eff_fast_muon_razor               = getplot_TH2F("scale_factors/muon/fastsim_muon_isolated_UL2018.root",     				"NUM_RazorPass_DEN_genTracks_abseta_pt", 				"mu6");
    eff_fast_muon_nonIso              = getplot_TH2F("scale_factors/muon/fastsim_muon_nonIsolated_UL2018.root",     		"NUM_RazorNoIso_DEN_genTracks_abseta_pt", 			"mu7");
    eff_fast_muon_nonIso_B2G          = getplot_TH2F("scale_factors/muon/fastsim_muon_nonIsolated_B2GCut_UL2018.root",  "h_eta_pT_muon",  															"mu8");
    unc_fast_muon_nonIso_B2G          = getplot_TH2F("scale_factors/muon/fastsim_muon_nonIsolated_B2GCut_UL2018.root",  "h_eta_pT_muon_stat",														"mu9");
    eff_full_pho_mediumid             = getplot_TH2F("scale_factors/photon/fullsim/egammaEffi.txt_EGM2D_Pho_Tight_UL18.root","EGamma_SF2D","pho1");
  } else if (v.year==2017) {
    eff_full_ele_veto                 = getplot_TH2F("scale_factors/electron/fullsim_electron_veto_UL2017.root",     						"EGamma_SF2D", 				"ele1");
    eff_full_ele_razor                = getplot_TH2F("scale_factors/electron/fullsim_electron_isolated_UL2017.root",     				"EGamma_SF2D", 				"ele2");
    eff_full_ele_nonIso               = getplot_TH2F("scale_factors/electron/fullsim_electron_nonIsolated_UL2017.root",     		"EGamma_SF2D", 				"ele3");
    eff_full_ele_nonIso_B2G           = getplot_TH2F("scale_factors/electron/fullsim_electron_nonIsolated_B2GCut_UL2017.root",  "h_eta_pT_electron",  "ele4");
    unc_full_ele_nonIso_B2G           = getplot_TH2F("scale_factors/electron/fullsim_electron_nonIsolated_B2GCut_UL2017.root",  "h_eta_pT_electron_stat",  "ele5");
    eff_fast_ele_razor                = getplot_TH2F("scale_factors/electron/fastsim_electron_isolated_UL2017.root",     				"EGamma_SF2D", 				"ele6");
    eff_fast_ele_nonIso               = getplot_TH2F("scale_factors/electron/fastsim_electron_nonIsolated_UL2017.root",     		"EGamma_SF2D", 				"ele7");
    eff_fast_ele_nonIso_B2G           = getplot_TH2F("scale_factors/electron/fastsim_electron_nonIsolated_B2GCut_UL2017.root",  "h_eta_pT_electron",  "ele8");
    unc_fast_ele_nonIso_B2G           = getplot_TH2F("scale_factors/electron/fastsim_electron_nonIsolated_B2GCut_UL2017.root",  "h_eta_pT_electron_stat",  "ele9");
    eff_full_muon_veto                = getplot_TH2F("scale_factors/muon/fullsim_muon_veto_UL2017.root",     						"NUM_RazorVeto_DEN_genTracks_abseta_pt", 				"mu1");
    eff_full_muon_razor               = getplot_TH2F("scale_factors/muon/fullsim_muon_isolated_UL2017.root",     				"NUM_RazorPass_DEN_genTracks_abseta_pt", 				"mu2");
    eff_full_muon_nonIso              = getplot_TH2F("scale_factors/muon/fullsim_muon_nonIsolated_UL2017.root",     		"NUM_RazorNoIso_DEN_genTracks_abseta_pt", 			"mu3");
    eff_full_muon_nonIso_B2G          = getplot_TH2F("scale_factors/muon/fullsim_muon_nonIsolated_B2GCut_UL2017.root",  "h_eta_pT_muon",  															"mu4");
    unc_full_muon_nonIso_B2G          = getplot_TH2F("scale_factors/muon/fullsim_muon_nonIsolated_B2GCut_UL2017.root",  "h_eta_pT_muon_stat",														"mu5");
    eff_fast_muon_razor               = getplot_TH2F("scale_factors/muon/fastsim_muon_isolated_UL2017.root",     				"NUM_RazorPass_DEN_genTracks_abseta_pt", 				"mu6");
    eff_fast_muon_nonIso              = getplot_TH2F("scale_factors/muon/fastsim_muon_nonIsolated_UL2017.root",     		"NUM_RazorNoIso_DEN_genTracks_abseta_pt", 			"mu7");
    eff_fast_muon_nonIso_B2G          = getplot_TH2F("scale_factors/muon/fastsim_muon_nonIsolated_B2GCut_UL2017.root",  "h_eta_pT_muon",  															"mu8");
    unc_fast_muon_nonIso_B2G          = getplot_TH2F("scale_factors/muon/fastsim_muon_nonIsolated_B2GCut_UL2017.root",  "h_eta_pT_muon_stat",														"mu9");
    eff_full_pho_mediumid             = getplot_TH2F("scale_factors/photon/fullsim/egammaEffi.txt_EGM2D_PHO_Tight_UL17.root","EGamma_SF2D","pho1");
  } else if(v.isAPV) {
    eff_full_ele_veto                 = getplot_TH2F("scale_factors/electron/fullsim_electron_veto_UL2016_preVFP.root",     						"EGamma_SF2D", 				"ele1");
    eff_full_ele_razor                = getplot_TH2F("scale_factors/electron/fullsim_electron_isolated_UL2016_preVFP.root",     				"EGamma_SF2D", 				"ele2");
    eff_full_ele_nonIso               = getplot_TH2F("scale_factors/electron/fullsim_electron_nonIsolated_UL2016_preVFP.root",     		"EGamma_SF2D", 				"ele3");
    eff_full_ele_nonIso_B2G           = getplot_TH2F("scale_factors/electron/fullsim_electron_nonIsolated_B2GCut_UL2016_preVFP.root",  "h_eta_pT_electron",  "ele4");
    unc_full_ele_nonIso_B2G           = getplot_TH2F("scale_factors/electron/fullsim_electron_nonIsolated_B2GCut_UL2016_preVFP.root",  "h_eta_pT_electron_stat",  "ele5");
    eff_fast_ele_razor                = getplot_TH2F("scale_factors/electron/fastsim_electron_isolated_UL2016.root",     				"EGamma_SF2D", 				"ele6");
    eff_fast_ele_nonIso               = getplot_TH2F("scale_factors/electron/fastsim_electron_nonIsolated_UL2016.root",     		"EGamma_SF2D", 				"ele7");
    eff_fast_ele_nonIso_B2G           = getplot_TH2F("scale_factors/electron/fastsim_electron_nonIsolated_B2GCut_UL2016.root",  "h_eta_pT_electron",  "ele8");
    unc_fast_ele_nonIso_B2G           = getplot_TH2F("scale_factors/electron/fastsim_electron_nonIsolated_B2GCut_UL2016.root",  "h_eta_pT_electron_stat",  "ele9");
    eff_full_muon_veto                = getplot_TH2F("scale_factors/muon/fullsim_muon_veto_UL2016_preVFP.root",     						"NUM_RazorVeto_DEN_genTracks_abseta_pt", 				"mu1");
    eff_full_muon_razor               = getplot_TH2F("scale_factors/muon/fullsim_muon_isolated_UL2016_preVFP.root",     				"NUM_RazorPass_DEN_genTracks_abseta_pt", 				"mu2");
    eff_full_muon_nonIso              = getplot_TH2F("scale_factors/muon/fullsim_muon_nonIsolated_UL2016_preVFP.root",     		"NUM_RazorNoIso_DEN_genTracks_abseta_pt", 			"mu3");
    eff_full_muon_nonIso_B2G          = getplot_TH2F("scale_factors/muon/fullsim_muon_nonIsolated_B2GCut_UL2016_preVFP.root",  "h_eta_pT_muon",  															"mu4");
    unc_full_muon_nonIso_B2G          = getplot_TH2F("scale_factors/muon/fullsim_muon_nonIsolated_B2GCut_UL2016_preVFP.root",  "h_eta_pT_muon_stat",														"mu5");
    eff_fast_muon_razor               = getplot_TH2F("scale_factors/muon/fastsim_muon_isolated_UL2016.root",     				"NUM_RazorPass_DEN_genTracks_abseta_pt", 				"mu6");
    eff_fast_muon_nonIso              = getplot_TH2F("scale_factors/muon/fastsim_muon_nonIsolated_UL2016.root",     		"NUM_RazorNoIso_DEN_genTracks_abseta_pt", 			"mu7");
    eff_fast_muon_nonIso_B2G          = getplot_TH2F("scale_factors/muon/fastsim_muon_nonIsolated_B2GCut_UL2016.root",  "h_eta_pT_muon",  															"mu8");
    unc_fast_muon_nonIso_B2G          = getplot_TH2F("scale_factors/muon/fastsim_muon_nonIsolated_B2GCut_UL2016.root",  "h_eta_pT_muon_stat",														"mu9");
    eff_full_pho_mediumid             = getplot_TH2F("scale_factors/photon/fullsim/egammaEffi.txt_EGM2D_Pho_Tight_UL16.root","EGamma_SF2D","pho1");
  } else {
    eff_full_ele_veto                 = getplot_TH2F("scale_factors/electron/fullsim_electron_veto_UL2016_postVFP.root",     						"EGamma_SF2D", 				"ele1");
    eff_full_ele_razor                = getplot_TH2F("scale_factors/electron/fullsim_electron_isolated_UL2016_postVFP.root",     				"EGamma_SF2D", 				"ele2");
    eff_full_ele_nonIso               = getplot_TH2F("scale_factors/electron/fullsim_electron_nonIsolated_UL2016_postVFP.root",     		"EGamma_SF2D", 				"ele3");
    eff_full_ele_nonIso_B2G           = getplot_TH2F("scale_factors/electron/fullsim_electron_nonIsolated_B2GCut_UL2016_postVFP.root",  "h_eta_pT_electron",  "ele4");
    unc_full_ele_nonIso_B2G           = getplot_TH2F("scale_factors/electron/fullsim_electron_nonIsolated_B2GCut_UL2016_postVFP.root",  "h_eta_pT_electron_stat",  "ele5");
    eff_fast_ele_razor                = getplot_TH2F("scale_factors/electron/fastsim_electron_isolated_UL2016.root",     				"EGamma_SF2D", 				"ele6");
    eff_fast_ele_nonIso               = getplot_TH2F("scale_factors/electron/fastsim_electron_nonIsolated_UL2016.root",     		"EGamma_SF2D", 				"ele7");
    eff_fast_ele_nonIso_B2G           = getplot_TH2F("scale_factors/electron/fastsim_electron_nonIsolated_B2GCut_UL2016.root",  "h_eta_pT_electron",  "ele8");
    unc_fast_ele_nonIso_B2G           = getplot_TH2F("scale_factors/electron/fastsim_electron_nonIsolated_B2GCut_UL2016.root",  "h_eta_pT_electron_stat",  "ele9");
    eff_full_muon_veto                = getplot_TH2F("scale_factors/muon/fullsim_muon_veto_UL2016_postVFP.root",     						"NUM_RazorVeto_DEN_genTracks_abseta_pt", 				"mu1");
    eff_full_muon_razor               = getplot_TH2F("scale_factors/muon/fullsim_muon_isolated_UL2016_postVFP.root",     				"NUM_RazorPass_DEN_genTracks_abseta_pt", 				"mu2");
    eff_full_muon_nonIso              = getplot_TH2F("scale_factors/muon/fullsim_muon_nonIsolated_UL2016_postVFP.root",     		"NUM_RazorNoIso_DEN_genTracks_abseta_pt", 			"mu3");
    eff_full_muon_nonIso_B2G          = getplot_TH2F("scale_factors/muon/fullsim_muon_nonIsolated_B2GCut_UL2016_postVFP.root",  "h_eta_pT_muon",  															"mu4");
    unc_full_muon_nonIso_B2G          = getplot_TH2F("scale_factors/muon/fullsim_muon_nonIsolated_B2GCut_UL2016_postVFP.root",  "h_eta_pT_muon_stat",														"mu5");
    eff_fast_muon_razor               = getplot_TH2F("scale_factors/muon/fastsim_muon_isolated_UL2016.root",     				"NUM_RazorPass_DEN_genTracks_abseta_pt", 				"mu6");
    eff_fast_muon_nonIso              = getplot_TH2F("scale_factors/muon/fastsim_muon_nonIsolated_UL2016.root",     		"NUM_RazorNoIso_DEN_genTracks_abseta_pt", 			"mu7");
    eff_fast_muon_nonIso_B2G          = getplot_TH2F("scale_factors/muon/fastsim_muon_nonIsolated_B2GCut_UL2016.root",  "h_eta_pT_muon",  															"mu8");
    unc_fast_muon_nonIso_B2G          = getplot_TH2F("scale_factors/muon/fastsim_muon_nonIsolated_B2GCut_UL2016.root",  "h_eta_pT_muon_stat",														"mu9");
    eff_full_pho_mediumid             = getplot_TH2F("scale_factors/photon/fullsim/egammaEffi.txt_EGM2D_Pho_Tight_UL16_postVFP.root","EGamma_SF2D","pho1");
  }
}


//_______________________________________________________
//              Define Histograms here
void
ScaleFactors::init_sf_histos()
{
  // btagging efficiency
  double ptbins[11]  = { 20,30,50,70,100,140,200,300,600,1000,4000 };
  double effbins[3] = { -0.5,0.5,1.5 };
  h_btag_eff_b_loose            = new TH2D("btag_eff_b_loose",  ";AK4 Jet p_{T} (GeV);Pass b-tag", 10,ptbins, 2,effbins);
  h_btag_eff_c_loose            = new TH2D("btag_eff_c_loose",  ";AK4 Jet p_{T} (GeV);Pass b-tag", 10,ptbins, 2,effbins);
  h_btag_eff_l_loose            = new TH2D("btag_eff_l_loose",  ";AK4 Jet p_{T} (GeV);Pass b-tag", 10,ptbins, 2,effbins);
  h_btag_eff_b_medium           = new TH2D("btag_eff_b_medium", ";AK4 Jet p_{T} (GeV);Pass b-tag", 10,ptbins, 2,effbins);
  h_btag_eff_c_medium           = new TH2D("btag_eff_c_medium", ";AK4 Jet p_{T} (GeV);Pass b-tag", 10,ptbins, 2,effbins);
  h_btag_eff_l_medium           = new TH2D("btag_eff_l_medium", ";AK4 Jet p_{T} (GeV);Pass b-tag", 10,ptbins, 2,effbins);

}


//_______________________________________________________
//               Fill Histograms here
void
ScaleFactors::fill_sf_histos(const bool& varySystematics, const bool& runOnSkim, const unsigned int& syst_index, const double& weight)
{
  if (syst_index == 0) {
    // btag efficiency - No event selection cuts to be applied
    // When making this plot, should remove all baseline cuts
    while (v.Jet.Loop()) {
      if (v.Jet.Jet.pass[v.Jet.i]) {
        if (v.Jet().hadronFlavour==5) {
          h_btag_eff_b_loose ->Fill(v.Jet().pt, v.Jet.LooseBTag.pass[v.Jet.i]);
          h_btag_eff_b_medium->Fill(v.Jet().pt, v.Jet.MediumBTag.pass[v.Jet.i]);
        } else if (v.Jet().hadronFlavour==4) {
          h_btag_eff_c_loose ->Fill(v.Jet().pt, v.Jet.LooseBTag.pass[v.Jet.i]);
          h_btag_eff_c_medium->Fill(v.Jet().pt, v.Jet.MediumBTag.pass[v.Jet.i]);
        } else {
          h_btag_eff_l_loose ->Fill(v.Jet().pt, v.Jet.LooseBTag.pass[v.Jet.i]);
          h_btag_eff_l_medium->Fill(v.Jet().pt, v.Jet.MediumBTag.pass[v.Jet.i]);
        }
      }
    }
  }
}

//_______________________________________________________
//              function to get scaled weight
double
ScaleFactors::get_syst_weight_(const double& weight_nominal, const double& weight_up, const double& weight_down, const double& nSigma)
{
  double w = weight_nominal;
  if (nSigma == 0) {
    return w;
  } else {
    // Compute the weight according to the systematic variation considered
    // Use difference between nominal and up/down as 1 sigma variation
    double dw_up = weight_up - weight_nominal;
    double dw_down = weight_nominal - weight_down;
    if (nSigma >= 0.) {
      w += nSigma*dw_up;
    } else {
      w += nSigma*dw_down;
    }
    return w;
  }
}

double
ScaleFactors::get_syst_weight_(const double& weight_nominal, const double& uncertainty, const double& nSigma)
{
  double w = weight_nominal;
  // Use symmetrical difference for up/down variation
  if (nSigma!=0.) w *= 1.0 + nSigma * uncertainty;
  return w;
}

//_______________________________________________________
//                Calculate scale factors


double ScaleFactors::calc_pu_jet_sf(const double& nSigmaPUTagSF){
  double w = 1;
  double sf, sf_err;
  while (v.Jet.Jet.Loop()) {
		if(!v.Jet.Jet().matchJet) continue;
    geteff2D(eff_pu, v.Jet.Jet().pt, v.Jet.Jet().eta, sf, sf_err);
		w *= get_syst_weight_(sf, sf_err, nSigmaPUTagSF);
	}
  return w;
}

std::pair<double, double> ScaleFactors::calc_b_tagging_sf(const double& nSigmaBTagSFbc, const double& nSigmaBTagSFlight, const double& nSigmaBTagSF, const double& nSigmaBTagSFbcAPV, const double& nSigmaBTagSFlightAPV) {

	double nSigmaBTag = nSigmaBTagSF;
  double pMC_loose = 1, pData_loose = 1;
  double pMC_medium = 1, pData_medium = 1;
  while (v.Jet.Loop()) {
    // Jet ID
    if (v.Jet.Jet.pass[v.Jet.i]) {
      float pt = v.Jet().pt, eta = v.Jet().eta;
      // Btag efficiencies (quark flavour dependent)
      BTagEntry::JetFlavor FLAV;
      double eff_medium = 1.0, eff_loose = 1.0;
      if (v.Jet().hadronFlavour==5) {
        FLAV = BTagEntry::FLAV_B;
        eff_loose  = geteff1D(eff_btag_b_loose,  pt, false);
        eff_medium = geteff1D(eff_btag_b_medium, pt, false);
				if(v.isAPV == 1 && nSigmaBTagSFbcAPV != 0) nSigmaBTag = nSigmaBTagSFbcAPV;
				if(v.isAPV != 1 && nSigmaBTagSFbc != 0)    nSigmaBTag = nSigmaBTagSFbc;
      } else if (v.Jet().hadronFlavour==4) {
        FLAV = BTagEntry::FLAV_C;
        eff_loose  = geteff1D(eff_btag_c_loose,  pt, false);
        eff_medium = geteff1D(eff_btag_c_medium, pt, false);
				if(v.isAPV == 1 && nSigmaBTagSFbcAPV != 0) nSigmaBTag = nSigmaBTagSFbcAPV;
				if(v.isAPV != 1 && nSigmaBTagSFbc != 0)    nSigmaBTag = nSigmaBTagSFbc;
      } else {
        FLAV = BTagEntry::FLAV_UDSG;
        eff_loose  = geteff1D(eff_btag_l_loose,  pt, false);
        eff_medium = geteff1D(eff_btag_l_medium, pt, false);
				if(v.isAPV == 1 && nSigmaBTagSFlightAPV != 0) nSigmaBTag = nSigmaBTagSFlightAPV;
				if(v.isAPV != 1 && nSigmaBTagSFlight != 0)    nSigmaBTag = nSigmaBTagSFlight;
      }

      // Scale factors - FullSim
      double sf_loose_cen   = btag_sf_full_loose_ ->eval_auto_bounds("central", FLAV, eta, pt);
      double sf_loose_up    = btag_sf_full_loose_ ->eval_auto_bounds("up",      FLAV, eta, pt);
      double sf_loose_down  = btag_sf_full_loose_ ->eval_auto_bounds("down",    FLAV, eta, pt);
      double sf_medium_cen  = btag_sf_full_medium_->eval_auto_bounds("central", FLAV, eta, pt);
      double sf_medium_up   = btag_sf_full_medium_->eval_auto_bounds("up",      FLAV, eta, pt);
      double sf_medium_down = btag_sf_full_medium_->eval_auto_bounds("down",    FLAV, eta, pt);

      double sf_loose       = get_syst_weight_(sf_loose_cen,  sf_loose_up,  sf_loose_down,  nSigmaBTag);
      double sf_medium      = get_syst_weight_(sf_medium_cen, sf_medium_up, sf_medium_down, nSigmaBTag);

      // Working points
      if (v.Jet.LooseBTag.pass[v.Jet.i]) {
        pMC_loose   *= eff_loose;
        pData_loose *= eff_loose * sf_loose;
      } else {
        pMC_loose   *= 1 - eff_loose;
        pData_loose *= 1 - eff_loose * sf_loose;
      }

      if (v.Jet.MediumBTag.pass[v.Jet.i]) {
        pMC_medium   *= eff_medium;
        pData_medium *= eff_medium * sf_medium;
      } else {
        pMC_medium   *= 1 - eff_medium;
        pData_medium *= 1 - eff_medium * sf_medium;
      }
    }
  }
  double weight_loose = 1, weight_medium = 1;
  if(pMC_loose!=0) weight_loose = pData_loose/pMC_loose;
  if(pMC_medium!=0) weight_medium = pData_medium/pMC_medium;
  return std::make_pair(weight_loose, weight_medium);
}

std::tuple<double, double, double> ScaleFactors::calc_ele_sf(const double& nSigmaEleFullSimSF, const double& nSigmaEleFastSimSF) {
  double sf, sf_err,tmp;
  double weight_veto  = 1.0, weight_select = 1.0, weight_nonIso = 1.0;
  while (v.Electron.Loop()) {
    double pt  = v.Electron().pt;
    double eta = v.Electron().eta;

    if (v.Electron.Veto.pass[v.Electron.i]) {
      	geteff2D(eff_full_ele_veto, pt, eta, sf, sf_err);
     		weight_veto *= get_syst_weight_(sf, sf_err, nSigmaEleFullSimSF);
    }
    if (v.Electron.Select.pass[v.Electron.i]) {
      	geteff2D(eff_full_ele_razor, pt, eta, sf, sf_err);
      	weight_select *= get_syst_weight_(sf, sf_err, nSigmaEleFullSimSF);
    }
    if (v.Electron.NonIso.pass[v.Electron.i]) {
        geteff2D(eff_full_ele_nonIso, pt, eta, sf, sf_err);
      	weight_nonIso *= get_syst_weight_(sf, sf_err, nSigmaEleFullSimSF);
        geteff2D(eff_full_ele_nonIso_B2G, pt, abs(eta), sf, tmp);
        geteff2D(unc_full_ele_nonIso_B2G, pt, abs(eta), sf_err, tmp);
      	weight_nonIso *= get_syst_weight_(sf, sf_err, nSigmaEleFullSimSF);
    }
  }
  return std::make_tuple(weight_veto, weight_select, weight_nonIso);
}

double ScaleFactors::calc_pho_sf() {
  double sf, sf_err;
  double weight_select = 1.0;
  while (v.Photon.Loop()) {
    double pt  = v.Photon().pt;
    double eta = v.Photon().eta;

    // Selected Photons
    if (v.Photon.Select.pass[v.Photon.i]) {
      // ID
      geteff2D(eff_full_pho_mediumid, pt, eta, sf, sf_err);
      weight_select *= sf;
    }
  }
  return weight_select;
}

std::tuple<double, double, double> ScaleFactors::calc_muon_sf(const double& nSigmaMuonFullSimSF, const double& nSigmaMuonFastSimSF) {
	double sf, sf_err, tmp;
  double weight_veto  = 1.0, weight_select = 1.0, weight_nonIso = 1.0;
  while (v.Muon.Loop()) {
    double pt  = v.Muon().pt;
    double eta = abs(v.Muon().eta);

    // Veto Muons
    if (v.Muon.Veto.pass[v.Muon.i]) {
      	geteff2D(eff_full_muon_veto, pt, eta, sf, sf_err);
     		weight_veto *= get_syst_weight_(sf, sf_err, nSigmaMuonFullSimSF);
    }
    // Selected Muons
    if (v.Muon.Select.pass[v.Muon.i]) {
      	geteff2D(eff_full_muon_razor, pt, eta, sf, sf_err);
      	weight_select *= sf;
      	weight_select *= get_syst_weight_(sf, sf_err, nSigmaMuonFullSimSF);
    }
    // nonIso Muons
    if (v.Muon.NonIso.pass[v.Muon.i]) {
      	geteff2D(eff_full_muon_nonIso, pt, eta, sf, sf_err);
      	weight_nonIso *= get_syst_weight_(sf, sf_err, nSigmaMuonFullSimSF);
      	geteff2D(eff_full_muon_nonIso_B2G, pt, eta, sf, tmp);
      	geteff2D(unc_full_muon_nonIso_B2G, pt, eta, sf_err, tmp);
      	weight_nonIso *= get_syst_weight_(sf, sf_err, nSigmaMuonFullSimSF);
    }
  }

  return std::make_tuple(weight_veto, weight_select, weight_nonIso);
}


//____________________________________________________
//          Analysis Specific Scale factors
//    (Defined for each search region separately)

void
ScaleFactors::apply_scale_factors(const unsigned int& syst_index, std::vector<double>& all_weights, std::vector<std::vector<double> >& w_nm1,
                                  const unsigned int& s, const std::vector<std::vector<double> >& nSigmaSFs)
{
  size_t i = 0;

  // Don't forget to specify the total number of sigmas you use in settings_*.h !
  // Electron SFs (5 sigmas - reco, fullsim id/iso, fastsim)
  if (debug) sw_(sw_s1, t_s1, 1);
  std::tie(sf_ele_veto, sf_ele_medium, sf_ele_nonIso) = calc_ele_sf(nSigmaSFs[i][s], nSigmaSFs[i+1][s]);
	//sf_ele_veto =1; sf_ele_medium=1; sf_ele_nonIso=1;
  i+=2;
  if (debug) sw_(sw_s1, t_s1, 0);

  // Photon SFs
  if (debug) sw_(sw_s1, t_s1, 1);
  sf_pho_medium = calc_pho_sf();
  if (debug) sw_(sw_s1, t_s1, 0);

  // Muon SFs (3 sigmas - tracking, fullsim, fastsim)
  if (debug) sw_(sw_s2, t_s2, 1);
  std::tie(sf_muon_veto, sf_muon_medium, sf_muon_nonIso) =  calc_muon_sf(nSigmaSFs[i][s], nSigmaSFs[i+1][s]);
	//sf_muon_veto = 1; sf_muon_medium=1; sf_muon_nonIso=1;
  i+=2;
  if (debug) sw_(sw_s2, t_s2, 0);

  // b tagging SFs (2 sigma - fullsim, fastsim)
  if (debug) sw_(sw_s3, t_s3, 1);
  std::pair<double, double> sf_btag = calc_b_tagging_sf(nSigmaSFs[i][s], nSigmaSFs[i+1][s], nSigmaSFs[i+2][s], nSigmaSFs[i+3][s], nSigmaSFs[i+4][s]);
  sf_btag_loose = sf_btag.first, sf_btag_medium = sf_btag.second;
	//sf_btag_loose=1;sf_btag_medium=1;
  i+=5;
  if (debug) sw_(sw_s3, t_s3, 0);

    // N-1 weights
  // Calculate weight for all search regions, but without a specific weight
  // Do not allow syst variation, because it's very slow
  // And we probably don't care about systematics in N-1 weight plots
  if (syst_index==0) {
    //for (const auto& region : magic_enum::enum_values<Region>()) {
    for (size_t region=0; region<Region::DiLep_eu; ++region) {
      if (debug) sw_(sw_s6, t_s6, 1);
      size_t n=all_weights.size()+scale_factors[region].size();
      if (w_nm1[region].empty()) w_nm1[region].resize(n,1);
      else w_nm1[region].assign(n,1);
      if (debug) sw_(sw_s6, t_s6, 0);
      if (debug) sw_(sw_s7, t_s7, 1);
      // [i] is the index of the weight to exclude
      for (size_t i=0; i<n; ++i) {
        w_nm1[region][i] = 1;
        if (!v.isData) {
            // [j] all the rest is applied
          for (size_t j=0; j<n; ++j) if (j!=i) {
            if (j<all_weights.size()) w_nm1[region][i] *= all_weights[j];
            else  w_nm1[region][i] *= (*scale_factors[region][j-all_weights.size()]);
          }
        }
      }
      if (debug) sw_(sw_s7, t_s7, 0);
    }
  }
}

void ScaleFactors::sw_(TStopwatch* sw, double& t, bool start=true) {
  if (start) {
    sw->Start(kTRUE); 
  } else {
    sw->Stop();
    t += sw-> RealTime();
  }
}

#endif // End header guard
