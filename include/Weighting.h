#ifndef WEIGHTING_H
#define WEIGHTING_H

// Private headers
#include "Variables.h"
#include "EventSelections.h"

// 3rd party headers
#include "tnm.h" // for getplot
#include "TH2.h"
#include "TH3.h"
#include "TGraphAsymmErrors.h"
#include "TString.h"

// common libraries
#include <algorithm>
#include <iostream>
#include <vector>

class Weighting
{
public:
  Weighting(Variables& var) : v(var) {
    w_nm1.resize(magic_enum::enum_count<EventSelections::Regions>());
    sf_weight.resize(magic_enum::enum_count<EventSelections::Regions>());
    all_weights.resize(10,1);
  }
  ~Weighting() {}

  void init_weight_histos();

  void init_input();

  void fill_weight_histos(const bool&, const bool&, const unsigned int&, const double&);

  double get_xsec_from_ntuple(const std::vector<std::string>&, const bool&);

  std::pair<double, double> get_xsec_totweight_from_txt_file(const std::string&);

  double get_totweight_from_ntuple(const std::vector<std::string>&, const bool&);

  void init_pileup_reweighting(const bool&, const std::vector<std::string>&);

  double get_toppt_weight(const double&, const unsigned int&, const bool&);

  double get_isr_weight(const double&, const unsigned int&, const bool&);

  double get_pileup_weight(const double, const double&, const unsigned int&, const bool&);

  double get_l1_prefiring_weight(const double&);

  double get_alphas_weight(const double&, const int&);

  double get_scale_weight(const double&, const unsigned int&);

  double calc_lostlep_weight(const double&);

  double calc_trigger_efficiency(const double&);

  double weight;
  std::vector<double> all_weights;
  std::vector<double> sf_weight;
  // N-1 weights
  std::vector<std::vector<double> > w_nm1;

  double trigger_had = 1.0;
  double trigger_lep = 1.0;

  std::map<uint32_t, std::string> signal_bins;

  double get_syst_weight(const double&, const double&, const double&, const double&);

  double get_syst_weight(const double&, const double&, const double&);

  Variables& v;

  // make sure binning is the same as in PlottingBase
  std::vector<double> HT_2D_bins = {200,  450,  600,  700, 800, 900, 1000, 1200, 10000}; // 2D Trigger Eff Run2017-18
  std::vector<double> MET_2D_bins = {60, 100, 130, 160, 180, 200, 250, 300, 400, 4000}; // 2D Trigger Eff Run2017-18
  std::vector<int> merged_trigger_bins = {1,3,5,  11,   20, 37,44,  46,52,53, 55,57,59,61,62, 64,65,66,67,68,69,70,71};

  // ISR weights
  // https://indico.cern.ch/event/592621/contributions/2398559/attachments/1383909/2105089/16-12-05_ana_manuelf_isr.pdf
  // https://indico.cern.ch/event/616816/contributions/2489809/attachments/1418579/2174166/17-02-22_ana_isr_ewk.pdf
  int isr_type = 0;
  std::vector<float> isr_weights_strong = {1, 0.92,  0.821, 0.715, 0.662, 0.561, 0.511};
  std::vector<float> isr_weights_weak = {1, 1.052, 1.179, 1.150, 1.057, 1.000, 0.912, 0.783};
  double isr_normfact = 1;

  //_______________________________________________________
  //                Input histograms

  TGraphAsymmErrors* trig_had_mu;
  TGraphAsymmErrors* trig_had_ele;
  TH2F* trig_ele;
  TH2F* trig_mu;
  TH2F* trig_mu_d;
  TH2F* trig_mu_m;

  TH2F* h_prefmap_photon;
  TH2F* h_prefmap_jet;

  //_______________________________________________________
  //             List of weighting Histograms

  TH1D* h_totweight;
  TH1D* h_totweight_toppt;
  TH1D* h_totweight_pileup;
  TH1D* h_nisrjets;
  TH1D* h_totweight_isr;
  TH1D* h_npvLowHigh_data;
  std::vector<TH3D*> vh_npvLowHigh_signal;
  TH1D* h_pileup_data;
  TH1D* h_pileup_data_down;
  TH1D* h_pileup_data_up;
  TH1D* h_pileup_mc;
  TH1D* h_pileup_weight;
  TH1D* h_pileup_weight_down;
  TH1D* h_pileup_weight_up;
  TH1D* h_nvtx;
  TH1D* h_nvtx_rw;

  TH1D* h_trigger_pass;
  TH1D* h_trigger_total;
  TH2D* h_trigger2d_pass;
  TH2D* h_trigger2d_total;
};

void Weighting::init_input() {
  if (v.isSignal) {
    if (v.sample.Contains("TChi")) {
      isr_type = 1; 
    } else {
      isr_type = 2;
    }

    std::ifstream isrFile("data/isr_normfact.txt");
    // Read all nSigmas, nums
    int year = 0;
    std::string sample = "";
    double normfact = 0;
    std::string line;
    while (std::getline(isrFile, line)) {
      std::stringstream nth_line;
      nth_line<<line;
      nth_line>>year;
      nth_line>>sample;
      nth_line>>normfact;
      if (v.sample == TString(sample) && year == v.year) {
        isr_normfact = normfact;
      }
    }
    std::cout<<"Signal ISR normalization factor: "<<isr_normfact<<std::endl;
  }

  
  // L1 prefiring maps
  if (v.year==2016) {
    h_prefmap_photon = getplot_TH2F("data/L1PreFiring/L1PrefiringMaps_WithUL17.root", "L1prefiring_photonptvseta_2016BtoH", "l1prepho");
    h_prefmap_jet    = getplot_TH2F("data/L1PreFiring/L1PrefiringMaps_WithUL17.root", "L1prefiring_jetptvseta_2016BtoH",    "l1prejet");
  } else if (v.year==2017) {
    h_prefmap_photon = getplot_TH2F("data/L1PreFiring/L1PrefiringMaps_WithUL17.root", "L1prefiring_photonptvseta_2017BtoF", "l1prepho");
    h_prefmap_jet    = getplot_TH2F("data/L1PreFiring/L1PrefiringMaps_WithUL17.root", "L1prefiring_jetptvseta_2017BtoF",    "l1prejet");
  }
  
  if (v.year==2016) {
	  if(v.isAPV) {
    	trig_had_mu       = getplot_TGraphAsymmErrors("trigger_eff/240602/TriggerEffRun2.root", "2016APV_had_mu",       "trig1");
    	trig_had_ele      = getplot_TGraphAsymmErrors("trigger_eff/240602/TriggerEffRun2.root", "2016APV_had_ele",      "trig2");
    	trig_ele          = getplot_TH2F("trigger_eff/trig_2016preVFP.root", "EGamma_SF2D",          "trig3");
    	trig_mu           = getplot_TH2F("trigger_eff/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_SingleMuonTriggers.root", "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdMedium_and_PFIsoMedium_abseta_pt",           "trig6");
    	trig_mu_d         = getplot_TH2F("trigger_eff/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_SingleMuonTriggers.root", "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdMedium_and_PFIsoMedium_abseta_pt_efficiencyData",           "trig7");
    	trig_mu_m         = getplot_TH2F("trigger_eff/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_SingleMuonTriggers.root", "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdMedium_and_PFIsoMedium_abseta_pt_efficiencyData",           "trig8");
		} else {
    	trig_had_mu       = getplot_TGraphAsymmErrors("trigger_eff/240602/TriggerEffRun2.root", "2016_had_mu",       "trig1");
    	trig_had_ele      = getplot_TGraphAsymmErrors("trigger_eff/240602/TriggerEffRun2.root", "2016_had_ele",      "trig2");
    	trig_ele          = getplot_TH2F("trigger_eff/trig_2016postVFP.root", "EGamma_SF2D",          "trig3");
    	trig_mu           = getplot_TH2F("trigger_eff/Efficiencies_muon_generalTracks_Z_Run2016_UL_SingleMuonTriggers.root", "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdMedium_and_PFIsoMedium_abseta_pt",           "trig6");
    	trig_mu_d         = getplot_TH2F("trigger_eff/Efficiencies_muon_generalTracks_Z_Run2016_UL_SingleMuonTriggers.root", "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdMedium_and_PFIsoMedium_abseta_pt_efficiencyData",           "trig7");
    	trig_mu_m         = getplot_TH2F("trigger_eff/Efficiencies_muon_generalTracks_Z_Run2016_UL_SingleMuonTriggers.root", "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdMedium_and_PFIsoMedium_abseta_pt_efficiencyData",           "trig8");
		}
  } else if (v.year==2017) {
    trig_had_mu       = getplot_TGraphAsymmErrors("trigger_eff/240602/TriggerEffRun2.root", "2017_had_mu",       "trig1");
    trig_had_ele      = getplot_TGraphAsymmErrors("trigger_eff/240602/TriggerEffRun2.root", "2017_had_ele",      "trig2");
    trig_ele          = getplot_TH2F("trigger_eff/trig_2017.root", "EGamma_SF2D",          "trig3");
    trig_mu           = getplot_TH2F("trigger_eff/Efficiencies_muon_generalTracks_Z_Run2017_UL_SingleMuonTriggers.root", "NUM_IsoMu27_DEN_CutBasedIdMedium_and_PFIsoMedium_abseta_pt",           "trig6");
    trig_mu_d         = getplot_TH2F("trigger_eff/Efficiencies_muon_generalTracks_Z_Run2017_UL_SingleMuonTriggers.root", "NUM_IsoMu27_DEN_CutBasedIdMedium_and_PFIsoMedium_abseta_pt_efficiencyData",           "trig7");
    trig_mu_m         = getplot_TH2F("trigger_eff/Efficiencies_muon_generalTracks_Z_Run2017_UL_SingleMuonTriggers.root", "NUM_IsoMu27_DEN_CutBasedIdMedium_and_PFIsoMedium_abseta_pt_efficiencyData",           "trig8");
  } else if (v.year==2018) {
    trig_had_mu       = getplot_TGraphAsymmErrors("trigger_eff/240602/TriggerEffRun2.root", "2018_had_mu",       "trig1");
    trig_had_ele      = getplot_TGraphAsymmErrors("trigger_eff/240602/TriggerEffRun2.root", "2018_had_ele",      "trig2");
    trig_ele          = getplot_TH2F("trigger_eff/trig_2018.root", "EGamma_SF2D",          "trig3");
    trig_mu           = getplot_TH2F("trigger_eff/Efficiencies_muon_generalTracks_Z_Run2018_UL_SingleMuonTriggers.root", "NUM_IsoMu24_DEN_CutBasedIdMedium_and_PFIsoMedium_abseta_pt",           "trig6");
    trig_mu_d         = getplot_TH2F("trigger_eff/Efficiencies_muon_generalTracks_Z_Run2018_UL_SingleMuonTriggers.root", "NUM_IsoMu24_DEN_CutBasedIdMedium_and_PFIsoMedium_abseta_pt_efficiencyData",           "trig7");
    trig_mu_m         = getplot_TH2F("trigger_eff/Efficiencies_muon_generalTracks_Z_Run2018_UL_SingleMuonTriggers.root", "NUM_IsoMu24_DEN_CutBasedIdMedium_and_PFIsoMedium_abseta_pt_efficiencyData",           "trig8");
  }
}


//_______________________________________________________
//              Define Histograms here
void
Weighting::init_weight_histos()
{
  // bins
  Double_t gluinoBins[202]; for (int i=0; i<202; ++i) gluinoBins[i] = (i-0.5)*25;
  Double_t stopBins[402];   for (int i=0; i<402; ++i) stopBins[i] = (i-0.5)*5;
  Double_t npvLowHighBins[3] = { 0,20,100 };
  // total weight
  h_totweight                     = new TH1D("totweight",           "MC;;Total (generator) event weight", 1,0,1);
  h_totweight_toppt               = new TH1D("totweight_toppt",     "MC;;Total toppt weight",             2,0,2);
  h_totweight_pileup              = new TH1D("totweight_pileup",    "MC;;Total pileup weight",            2,0,2);
  // ISR reweighting
  h_nisrjets                      = new TH1D("nisrjets",            ";N_{ISR jets}", 16,-0.5,15.5);
  h_totweight_isr                 = new TH1D("totweight_isr",       "MC;;Total (generator) event weight", 2,0,2);
  // npv for extrapolations
  h_npvLowHigh_data               = new TH1D("npvLowHigh_data",     "Number of vertices - Data;N_{Vertices}", 2,npvLowHighBins);
  vh_npvLowHigh_signal   .push_back(new TH3D("npvLowHigh_T1tttt",   "T1tttt or T5ttcc or T5tttt;m_{#tilde{g}} (GeV);m_{#tilde{#chi}^{0}_{1}} (GeV);N_{vertex}", 201,gluinoBins, 201,gluinoBins, 2,npvLowHighBins));
  vh_npvLowHigh_signal   .push_back(new TH3D("npvLowHigh_T2tt",     "T2tt;m_{#tilde{t}} (GeV);m_{#tilde{#chi}^{0}_{1}} (GeV);N_{vertex}",                       401,stopBins, 401,stopBins,     2,npvLowHighBins));
  vh_npvLowHigh_signal   .push_back(new TH3D("npvLowHigh_TChiWZ",   "TChiWZ;m_{#tilde{#chi}^{#pm}_{0}=#tilde{#chi}^{0}_{2}} (GeV);m_{#tilde{#chi}^{0}_{1}} (GeV);N_{vertex}", 401,stopBins, 401,stopBins,     2,npvLowHighBins));
  vh_npvLowHigh_signal   .push_back(new TH3D("npvLowHigh_TChiHH",   "TChiHH;m_{#tilde{#chi}^{0}_{3}=#tilde{#chi}^{0}_{2}} (GeV);m_{#tilde{#chi}^{0}_{1}} (GeV);N_{vertex}", 401,stopBins, 401,stopBins,     2,npvLowHighBins));
  vh_npvLowHigh_signal   .push_back(new TH3D("npvLowHigh_T5qqqqZH", "T5qqqqZH;m_{#tilde{g}} (GeV);m_{#tilde{#chi}^{0}_{1}} (GeV);N_{vertex}", 201,gluinoBins, 201,gluinoBins, 2,npvLowHighBins));
  //vh_npvLowHigh_signal   .push_back(new TH3D("npvLowHigh_T6bbZH",   "T6bbZH;m_{#tilde{t}} (GeV);#tilde{#chi}^{0}_{2}} (GeV);m_{#tilde{#chi}^{0}_{1}} (GeV);N_{vertex}", 201,gluinoBins,201,gluinoBins, 201,gluinoBins, 2,npvLowHighBins));
  // pileup
  h_pileup_data                = new TH1D("pileup_data",        "Pile-up distribution - Data (Nominal);Pile-up", 100,0,100);
  h_pileup_data_down           = new TH1D("pileup_data_down",   "Pile-up distribution - Data (down);Pile-up",    100,0,100);
  h_pileup_data_up             = new TH1D("pileup_data_up",     "Pile-up distribution - Data (up);Pile-up",      100,0,100);
  h_pileup_mc                  = new TH1D("pileup_mc",          "Pile-up distribution - MC;Pile-up",             100,0,100);
  h_pileup_weight              = new TH1D("pileup_weight",      "Pile-up weights - Nominal MB X-sec (69 mb);Pile-up;Weight",    100,0,100);
  h_pileup_weight_down         = new TH1D("pileup_weight_down", "Pile-up weights - MB X-sec up 5% (72.45 mb);Pile-up;Weight",   100,0,100);
  h_pileup_weight_up           = new TH1D("pileup_weight_up",   "Pile-up weights - MB X-sec down 5% (65.55 mb);Pile-up;Weight", 100,0,100);
  h_nvtx                       = new TH1D("nvtx",               "Number of vertices - Nominal;N_{Vertices}",                      101,-0.5,100.5);
  h_nvtx_rw                    = new TH1D("nvtx_rw",            "Number of vertices - Pile-up reweighted (MC only);N_{Vertices}", 101,-0.5,100.5);

  // trigger efficiency
  double htbins[19]  = { 0, 200, 300, 400, 500, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1200, 1500, 2000, 4000, 10000 };
  double HTB[12] = {400, 500, 600, 700, 750, 800, 850, 900, 950, 1000, 1500, 10000};
  double ptB[9]  = {200, 300, 400, 450, 500, 550, 600, 1000, 10000};
  h_trigger_pass                = new TH1D("trigger_pass",    "Pass trigger;H_{T} (GeV)", 18,htbins);
  h_trigger_total               = new TH1D("trigger_total",          "Total;H_{T} (GeV)", 18,htbins);
  h_trigger2d_pass              = new TH2D("trigger2d_pass",  "Pass trigger;H_{T} (GeV);Leading AK8 jet p_{T} (GeV)", 11,HTB, 8,ptB);
  h_trigger2d_total             = new TH2D("trigger2d_total",        "Total;H_{T} (GeV);Leading AK8 jet p_{T} (GeV)", 11,HTB, 8,ptB);

}

//_______________________________________________________
//               Fill Histograms here
void
Weighting::fill_weight_histos(const bool& varySystematics, const bool& runOnSkim, const unsigned int& syst_index, const double& weight)
{
  // ISR jets counting
  // Taken from:
  // https://github.com/manuelfs/babymaker/blob/0136340602ee28caab14e3f6b064d1db81544a0a/bmaker/plugins/bmaker_full.cc#L1268-L1295
  if (!v.isData&&syst_index==0) {
    h_nisrjets->Fill(v.nJetISR, weight);
  }
}


//_______________________________________________________
//           Read cross-section from ntuple
double
Weighting::get_xsec_from_ntuple(const std::vector<std::string>& filenames, const bool& runOnSkim)
{
  float evt_XSec=0, prev_XSec=0;
  for (const auto& filename : filenames) {
    TFile *f = TFile::Open(filename.c_str());
    TTree* tree = (TTree*)f->Get(runOnSkim ? "B2GTree" : "B2GTTreeMaker/B2GTree");
    tree->GetBranch("evt_XSec")->SetAddress(&evt_XSec);
    tree->GetEntry(0);
    f->Close();
    if (prev_XSec!=0&&prev_XSec!=evt_XSec) {
      error("Weighting - Files added with different cross-sections. Please, add them separately!");
      return 0;
    }
    prev_XSec = evt_XSec;
  }
  return evt_XSec;
}

//_______________________________________________________
//           Read cross-section from txt file
std::pair<double, double>
Weighting::get_xsec_totweight_from_txt_file(const std::string& txt_file)
{
  double XSec = 0, Totweight = 0;
  std::ifstream xsecFile(txt_file.c_str());
  if ( !xsecFile.good() ) {
    return std::make_pair(0,0);
    std::cout<<"Unable to open cross-section file: "<<txt_file<<std::endl;
    error("Please provide the correct txt file for Cross-sections in settings.h!");
  } else {

    std::string line;
    std::string shortname, primary_dataset;
    double xsec, totweight;
    while ( std::getline(xsecFile, line) ) {
      std::stringstream nth_line;
      nth_line<<line;
      nth_line>>shortname;
      nth_line>>primary_dataset;
      nth_line>>xsec;
      nth_line>>totweight;
      if (v.sample.Contains(primary_dataset)) {
        XSec = xsec;
        Totweight = totweight;
      }

    }
  }
  if (XSec == 0) {
    std::cout<<"No crossection found for "<<v.sample<<" in cross section file: "<<txt_file<<std::endl;
    error("Please fix the cross-section file in settings.h!");
  }

  return std::make_pair(XSec, Totweight);
}

//_______________________________________________________
//          Read total weight from ntuple histos
double
Weighting::get_totweight_from_ntuple(const std::vector<std::string>& filenames, const bool& runOnSkim)
{
  // Merging totweight histos
  for (const auto& filename : filenames) {
    TFile* f = TFile::Open(filename.c_str());
    h_totweight->Add((TH1D*)f->Get(runOnSkim ? "totweight" : "EventCounter/totweight"));
    f->Close();
  }
  return h_totweight->GetBinContent(1);
}

//_______________________________________________________
//             Load pile-up reweighting infos
void
Weighting::init_pileup_reweighting(const bool& runOnSkim, const std::vector<std::string>& filenames)
{
  std::string pileupDir = "pileup/UltraLegacy2017/";
  if (v.year==2016) pileupDir = "pileup/UltraLegacy2016/";
  else if (v.year==2018) pileupDir = "pileup/UltraLegacy2018/";
	if(v.year != 2016) {
  	// Get data histogram (generated by pileupCalc.py script)
  	TFile* f_pileup_data = TFile::Open((pileupDir+"data_pileup.root").c_str());
 		h_pileup_data->Add((TH1D*)f_pileup_data->Get("pileup"));
	  f_pileup_data->Close();
	  // Also get up/down variations
	  TFile* f_pileup_data_down = TFile::Open((pileupDir+"data_pileup_down.root").c_str());
	  h_pileup_data_down->Add((TH1D*)f_pileup_data_down->Get("pileup"));
	  f_pileup_data_down->Close();
	  TFile* f_pileup_data_up = TFile::Open((pileupDir+"data_pileup_up.root").c_str());
	  h_pileup_data_up->Add((TH1D*)f_pileup_data_up->Get("pileup"));
 		f_pileup_data_up->Close();
	} else if(v.isAPV) {
  	// Get data histogram (generated by pileupCalc.py script)
  	TFile* f_pileup_data = TFile::Open((pileupDir+"data_pileup_preVFP.root").c_str());
 		h_pileup_data->Add((TH1D*)f_pileup_data->Get("pileup"));
	  f_pileup_data->Close();
	  // Also get up/down variations
	  TFile* f_pileup_data_down = TFile::Open((pileupDir+"data_pileup_down_preVFP.root").c_str());
	  h_pileup_data_down->Add((TH1D*)f_pileup_data_down->Get("pileup"));
	  f_pileup_data_down->Close();
	  TFile* f_pileup_data_up = TFile::Open((pileupDir+"data_pileup_up_preVFP.root").c_str());
	  h_pileup_data_up->Add((TH1D*)f_pileup_data_up->Get("pileup"));
 		f_pileup_data_up->Close();
	} else {
  	// Get data histogram (generated by pileupCalc.py script)
  	TFile* f_pileup_data = TFile::Open((pileupDir+"data_pileup_postVFP.root").c_str());
 		h_pileup_data->Add((TH1D*)f_pileup_data->Get("pileup"));
	  f_pileup_data->Close();
	  // Also get up/down variations
	  TFile* f_pileup_data_down = TFile::Open((pileupDir+"data_pileup_down_postVFP.root").c_str());
	  h_pileup_data_down->Add((TH1D*)f_pileup_data_down->Get("pileup"));
	  f_pileup_data_down->Close();
	  TFile* f_pileup_data_up = TFile::Open((pileupDir+"data_pileup_up_postVFP.root").c_str());
	  h_pileup_data_up->Add((TH1D*)f_pileup_data_up->Get("pileup"));
 		f_pileup_data_up->Close();
	}
  // get mc histogram (used to generate mc pile-up)
  TFile* f_pileup_mc = TFile::Open((pileupDir+"mc_pileup.root").c_str());
  h_pileup_mc->Add((TH1D*)f_pileup_mc->Get("pileup"));
  f_pileup_mc->Close();
  // Divide normalized data histo by normalized mc histo to get pileup weights for each bin
  h_pileup_weight     ->Divide(h_pileup_data,      h_pileup_mc, 1/h_pileup_data->Integral(),      1/h_pileup_mc->Integral());
  h_pileup_weight_down->Divide(h_pileup_data_down, h_pileup_mc, 1/h_pileup_data_down->Integral(), 1/h_pileup_mc->Integral());
  h_pileup_weight_up  ->Divide(h_pileup_data_up,   h_pileup_mc, 1/h_pileup_data_up->Integral(),   1/h_pileup_mc->Integral());
}

//_______________________________________________________
//              function to get scaled weight
double
Weighting::get_syst_weight(const double& weight_nominal, const double& weight_up, const double& weight_down, const double& nSigma)
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
Weighting::get_syst_weight(const double& weight_nominal, const double& uncertainty, const double& nSigma)
{
  double w = weight_nominal;
  // Use symmetrical difference for up/down variation
  if (nSigma!=0.) w *= 1.0 + nSigma * uncertainty;
  return w;
}


//_______________________________________________________
//                  Top pt reweighting
double
Weighting::get_toppt_weight(const double& nSigmaToppt, const unsigned int& syst_index, const bool& runOnSkim)
{
  double w_nom = 1;//, n=0;
  if (!v.isData) while (v.GenPart.Top.Loop()) {
		w_nom *= 0.103*std::exp(-0.0118*v.GenPart.Top().pt) -0.000134 * v.GenPart.Top().pt+0.973;
  }
  double w_toppt_up = 1;
  double w_toppt = std::sqrt(w_nom);
  double w_toppt_down = std::sqrt(w_nom);
  double w = get_syst_weight(w_toppt, w_toppt_up, w_toppt_down, nSigmaToppt);
  if (syst_index==0&&!runOnSkim) {
    h_totweight_toppt->Fill(0);
    h_totweight_toppt->Fill(1, w_toppt);
  }
  return w;
}

//_______________________________________________________
//                    ISR reweighting
int nJetISR;
double
Weighting::get_isr_weight(const double& nSigmaISR, const unsigned int& syst_index, const bool& runOnSkim)
{
  // Implementing the reweighting in this presentation:
  // https://indico.cern.ch/event/592621/contributions/2398559/attachments/1383909/2105089/16-12-05_ana_manuelf_isr.pdf
  // Using the values found on slide 8 (T2tt and T1tttt)
  double w = 1;
  isr_weights_strong = {1, 0.92,  0.821, 0.715, 0.662, 0.561, 0.511};
  isr_weights_weak = {1, 1.052, 1.179, 1.150, 1.057, 1.000, 0.912, 0.783};
  // ttbar ISR reweighting not needed, we do top pt reweighting!
  if (v.isSignal) {
    w = 0;
    if (isr_type==1) {
      double EWkino_pt = v.susy_mass[0];
      size_t bin;
      if      (EWkino_pt< 50) bin = 0;
      else if (EWkino_pt<100) bin = 1;
      else if (EWkino_pt<150) bin = 2;
      else if (EWkino_pt<200) bin = 3;
      else if (EWkino_pt<300) bin = 4;
      else if (EWkino_pt<400) bin = 5;
      else if (EWkino_pt<600) bin = 6;
      else bin = 7;
      w = isr_normfact * isr_weights_weak[bin];
    } else {
			if(syst_index==0) nJetISR = v.nJetISR;
			else        v.nJetISR = nJetISR;
      size_t bin = std::min(size_t(6), v.nJetISR);
      w = isr_normfact * isr_weights_strong[bin];
    }
    //double err = (1-w)/2;
    //double w_isr_up   = w + err;
    double w_isr      = w;
    //double w_isr_down = w - err;
    //w = get_syst_weight(w_isr, w_isr_up, w_isr_down, nSigmaISR);
    w = get_syst_weight(1, 1.01, 0.99, nSigmaISR);
    h_totweight_isr->Fill(0);
    h_totweight_isr->Fill(1, w_isr);
  }
  return w;
}

//_______________________________________________________
//                  Get pile-up weight
double
Weighting::get_pileup_weight(const double weight, const double& nSigmaPU, const unsigned int& syst_index, const bool& runOnSkim)
{
  // Background
  int pu_bin = v.Pileup_nTrueInt+1; // eg. pileup 0, is filled in bin 1
  double w_pileup = h_pileup_weight->GetBinContent(pu_bin);
  double w_pileup_up = h_pileup_weight_up->GetBinContent(pu_bin);
  double w_pileup_down = h_pileup_weight_down->GetBinContent(pu_bin);
  double w = get_syst_weight(w_pileup, w_pileup_up, w_pileup_down, nSigmaPU);
  if (syst_index==0&&!runOnSkim) {
    h_totweight_pileup->Fill(0);
    h_totweight_pileup->Fill(1, w_pileup);
  }
  if (syst_index == 0) {
    h_nvtx   ->Fill(v.PV_npvsGood, weight);
    h_nvtx_rw->Fill(v.PV_npvsGood, weight*w);
  }
  return w;
}


//_______________________________________________________
//          Get L1 prefiring weight for 2016/2017
double
Weighting::get_l1_prefiring_weight(const double& nSigmaL1PreFiring)
{
  // Background and Signal
  if (v.year==2016||v.year==2017) {
    if (v.isSignal) {
      // Implement by hand: PhysicsTools/PatUtils/plugins/L1ECALPrefiringWeightProducer.cc
      //Probability for the event NOT to prefire, computed with the prefiring maps per object.
      //Up and down values correspond to the resulting value when shifting up/down all prefiring rates in prefiring maps.
      bool useEMpt_ = false;
      double prefiringRateSystUnc_ = 0.2;
      double nonPrefiringProba[3] = {1., 1., 1.};  //0: central, 1: up, 2: down
      
      // Photons
      while (v.Photon.Loop()) {
        double pt_gam  = v.Photon().pt;
        double eta_gam = v.Photon().eta;
        if (pt_gam>=20 && std::abs(eta_gam)>=2 && std::abs(eta_gam)<=3) {
          // getPrefiringRate
          int nbinsy = h_prefmap_photon->GetNbinsY();
          double maxy = h_prefmap_photon->GetYaxis()->GetBinLowEdge(nbinsy + 1);
          if (pt_gam >= maxy) pt_gam = maxy - 0.01;
          int thebin = h_prefmap_photon->FindBin(eta_gam, pt_gam);
          double prefrate = h_prefmap_photon->GetBinContent(thebin);
          double statuncty = h_prefmap_photon->GetBinError(thebin);
          double systuncty = prefiringRateSystUnc_ * prefrate;
          double prefrate_up = std::min(1., prefrate + sqrt(pow(statuncty, 2) + pow(systuncty, 2)));
          double prefrate_dn = std::max(0., prefrate - sqrt(pow(statuncty, 2) + pow(systuncty, 2)));
          nonPrefiringProba[0] *= (1. - prefrate);
          nonPrefiringProba[1] *= (1. - prefrate_up);
          nonPrefiringProba[2] *= (1. - prefrate_dn);
        }
      } // end photon loop
      
      //Now applying the prefiring maps to jets in the affected regions.
      while (v.Jet.Loop()) {
        double pt_jet  = v.Jet().pt_nom;
        double eta_jet = v.Jet().eta;
        if (pt_jet>=20 && std::abs(eta_jet)>=2 && std::abs(eta_jet)<=3) {
      
          //Loop over photons to remove overlap
          double nonprefiringprobfromoverlappingphotons[3] = {1., 1., 1.};
          while (v.Photon.Loop()) {
            double pt_gam  = v.Photon().pt;
            double eta_gam = v.Photon().eta;
            if (pt_gam>=20 && std::abs(eta_gam)>=2 && std::abs(eta_gam)<=3) {
              double dR = DeltaR(v.Photon.v4(), v.Jet.v4());
              if (dR<=0.4) {
                // getPrefiringRate
                int nbinsy = h_prefmap_photon->GetNbinsY();
                double maxy = h_prefmap_photon->GetYaxis()->GetBinLowEdge(nbinsy + 1);
                if (pt_gam >= maxy) pt_gam = maxy - 0.01;
                int thebin = h_prefmap_photon->FindBin(eta_gam, pt_gam);
                double prefrate = h_prefmap_photon->GetBinContent(thebin);
                double statuncty = h_prefmap_photon->GetBinError(thebin);
                double systuncty = prefiringRateSystUnc_ * prefrate;
                double prefrate_up = std::min(1., prefrate + sqrt(pow(statuncty, 2) + pow(systuncty, 2)));
                double prefrate_dn = std::max(0., prefrate - sqrt(pow(statuncty, 2) + pow(systuncty, 2)));
                nonprefiringprobfromoverlappingphotons[0] = (1. - prefrate);
                nonprefiringprobfromoverlappingphotons[1] = (1. - prefrate_up);
                nonprefiringprobfromoverlappingphotons[2] = (1. - prefrate_dn);
              }
            }
          } // end photon loop within jet loop
          
          //useEMpt =true if one wants to use maps parametrized vs Jet EM pt instead of pt.
          if (useEMpt_) pt_jet *= (v.Jet().neEmEF + v.Jet().chEmEF);
          // getPrefiringRate
          int nbinsy = h_prefmap_jet->GetNbinsY();
          double maxy = h_prefmap_jet->GetYaxis()->GetBinLowEdge(nbinsy + 1);
          if (pt_jet >= maxy) pt_jet = maxy - 0.01;
          int thebin = h_prefmap_jet->FindBin(eta_jet, pt_jet);
          double prefrate = h_prefmap_jet->GetBinContent(thebin);
          double statuncty = h_prefmap_jet->GetBinError(thebin);
          double systuncty = prefiringRateSystUnc_ * prefrate;
          double prefrate_up = std::min(1., prefrate + sqrt(pow(statuncty, 2) + pow(systuncty, 2)));
          double prefrate_dn = std::max(0., prefrate - sqrt(pow(statuncty, 2) + pow(systuncty, 2)));
          double nonprefiringprobfromoverlappingjet[3] = {1. - prefrate, 1. - prefrate_up, 1. - prefrate_dn};
          for (int i=0; i<3; ++i) {
            if (nonprefiringprobfromoverlappingphotons[i] == 1.)
              nonPrefiringProba[i] *= nonprefiringprobfromoverlappingjet[i];
            //If overlapping photons have a non prefiring rate larger than the jet, then replace these weights by the jet one
            else if (nonprefiringprobfromoverlappingphotons[i] > nonprefiringprobfromoverlappingjet[i]) {
              if (nonprefiringprobfromoverlappingphotons[i] != 0.)
                nonPrefiringProba[i] *= nonprefiringprobfromoverlappingjet[i] / nonprefiringprobfromoverlappingphotons[i];
              else 
                nonPrefiringProba[i] = 0.;
            }
          }
          //Last case: if overlapping photons have a non prefiring rate smaller than the jet, don't consider the jet in the event weight, and do nothing.
        }
      } // end jet loop
      return get_syst_weight(nonPrefiringProba[0], nonPrefiringProba[1], nonPrefiringProba[2], nSigmaL1PreFiring);
    } else {
      // Get it from ntuple
      return get_syst_weight(v.L1PreFiringWeight_Nom, v.L1PreFiringWeight_Up, v.L1PreFiringWeight_Dn, nSigmaL1PreFiring);
    }
  }
  return 1;
}

//_______________________________________________________
//                  Get alpha_s weight
double
Weighting::get_alphas_weight(const double& nSigmaAlphaS, const int& LHA_PDF_ID)
{
  /*
  std::vector<float> alphas_Weights; // not available in NanoAOD

  // A set of two weights corresponding to
  // Powheg:  alpha_s = 0.118 -+ 0.002
  // aMC@NLO: alpha_s = 0.118 -+ 0.001
  // Recommendation is to use +- 0.0015 --> rescale difference by 0.75 or 1.5
  // Treat weight as usual, gaussian, rescale to desired nSigma
  double w_alphas = 1;
  //double w_alphas_up   = alphas_Weights[1];
  //double w_alphas_down = alphas_Weights[0];
  double w_alphas_up   = 0;
  double w_alphas_down = 0;
  double nSigma_0_0015 = nSigmaAlphaS;
  if (LHA_PDF_ID==260000||LHA_PDF_ID==260400) {
    // Powheg samples have -+ 0.001
    nSigma_0_0015 *= 1.5;
  } else {
    // aMC@NLO samples have -+ 0.002
    nSigma_0_0015 *= 0.75;
  }
  w_alphas = get_syst_weight(w_alphas, w_alphas_up, w_alphas_down, nSigma_0_0015);
  return w_alphas;
  */
  return 1;
}


//_______________________________________________________
//                  Get scale weight
double
//Weighting::get_scale_weight(const std::vector<double>& scale_weight_norm, const double& nSigmaScale, const unsigned int& numScale)
Weighting::get_scale_weight(const double& nSigmaScale, const unsigned int& numScale)
{
  //std::vector<float> scale_Weights;
  /*
    New LHEScaleWeight in NanoAOD:
    LHE scale variation weights (w_var / w_nominal); 
    v.LHEScaleWeight[0] is rensc=0.5d0 facsc=0.5d0
    v.LHEScaleWeight[1] is rensc=0.5d0 facsc=1d0
    v.LHEScaleWeight[2] is rensc=0.5d0 facsc=2d0
    v.LHEScaleWeight[3] is rensc=  1d0 facsc=0.5d0
    v.LHEScaleWeight[4] is rensc=  1d0 facsc=1d0
    v.LHEScaleWeight[5] is rensc=  1d0 facsc=2d0
    v.LHEScaleWeight[6] is rensc=  2d0 facsc=0.5d0
    v.LHEScaleWeight[7] is rensc=  2d0 facsc=1d0
    v.LHEScaleWeight[8] is rensc=  2d0 facsc=2d0 
  */
  if (nSigmaScale==0) return 1; // No systematics
  if (v.nLHEScaleWeight==0) return 1; // ST samples are known to miss scale weights
  double w_scale = 1;
  double w_scale_up = 1;   // Corresponds to 0.5 (More signal events)
  double w_scale_down = 1; // Corresponds to 2.0 (Less signal events)
  w_scale_up   = v.LHEScaleWeight[0];
  w_scale_down = v.LHEScaleWeight[8];
  w_scale = get_syst_weight(w_scale, w_scale_up, w_scale_down, nSigmaScale);
  return w_scale;
}

double Weighting::calc_lostlep_weight(const double& nSigmaLostLep) {
  // Lost Lepton event weight
  // First run this script to obtain the lost lepton uncertainties
  // python scripts/calc_lostlepsyst.py
  // Use the averages combination of all leptons for W/top
  double w = 1;
  if (!v.isData) if (v.Electron.Veto.n+v.Muon.Veto.n+v.Tau.Select.n==0) {
    // Final state leptons from W mother
    while (v.GenPart.LeptonFromW.Loop()) {
      double unc = 0;
      double abseta = std::abs(v.GenPart.LeptonFromW().eta);
      if      (abseta<0.5) unc = 0.125;
      else if (abseta<1.0) unc = 0.126;
      else if (abseta<1.5) unc = 0.129;
      else if (abseta<2.0) unc = 0.143;
      else if (abseta<2.5) unc = 0.175;
      w *= get_syst_weight(1.0, unc, nSigmaLostLep);
      //std::cout<<"Lost-lepton found: pt/eta/id = "<<v.GenPart.LeptonFromW().pt[i]<<" "
      //         <<v.GenPart.LeptonFromW().eta[i]<<" "<<v.GenPart.LeptonFromW().pdgId
      //         <<" E/M/T="<<v.Electron.Veto.n<<"/"<<v.Muon.Veto.n<<"/"<<v.Tau.Select.n<<std::endl;
    }
  }
  return w;
}


double Weighting::calc_trigger_efficiency(const double& nSigmaTrigger) {
  double sf = 0, err = 0;
  double eff1 = 0, eff2 = 0;
  // leptonic trigger efficiencies
  if (v.Muon.Select.n==1 && v.Electron.Select.n==1) {
		if (v.HLT_IsoMu24==1 || v.HLT_IsoTkMu24==1 || v.HLT_IsoMu27==1)
    	geteff2D(trig_mu, v.Muon.Select(0).pt, abs(v.Muon.Select(0).eta), sf, err);
		else if(v.HLT_Ele32_WPTight_Gsf==1 || v.HLT_Ele27_WPTight_Gsf==1 || v.HLT_Ele32_WPTight_Gsf_L1DoubleEG==1)
    	geteff2D(trig_ele, v.Electron.Select(0).pt, v.Electron.Select(0).eta, sf, err);
    trigger_lep = get_syst_weight(sf, err, nSigmaTrigger);
  } else if (v.Muon.Select.n>0) {
    geteff2D(trig_mu, v.Muon.Select(0).pt, abs(v.Muon.Select(0).eta), sf, err);
    trigger_lep = get_syst_weight(sf, err, nSigmaTrigger);
		if (v.Muon.Select.n == 2){
			geteff2D(trig_mu_d, v.Muon.Select(0).pt, abs(v.Muon.Select(0).eta), eff1, err);
			geteff2D(trig_mu_d, v.Muon.Select(1).pt, abs(v.Muon.Select(1).eta), eff2, err);
			trigger_lep = 1 - (get_syst_weight(eff1, err, nSigmaTrigger)*get_syst_weight(eff1, err, nSigmaTrigger));
			geteff2D(trig_mu_m, v.Muon.Select(0).pt, abs(v.Muon.Select(0).eta), eff1, err);
			geteff2D(trig_mu_m, v.Muon.Select(1).pt, abs(v.Muon.Select(1).eta), eff2, err);
			trigger_lep *= 1/(1 - (get_syst_weight(eff1, err, nSigmaTrigger)*get_syst_weight(eff1, err, nSigmaTrigger)));
		}
  } else if (v.Electron.Select.n>0) {
    geteff2D(trig_ele, v.Electron.Select(0).pt, v.Electron.Select(0).eta, sf, err);
    trigger_lep = get_syst_weight(sf, err, nSigmaTrigger);
  } else trigger_lep = 0;
  
	double w = trigger_lep;
  return w;
}

#endif // End header guard
