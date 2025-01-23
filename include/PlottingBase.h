#ifndef PLOTTINGBASE_H
#define PLOTTINGBASE_H

// Private headers
#include "Variables.h"
#include "Weighting.h"
#include "EventSelections.h"
#include "SmartHistos.h"

// 3rd party headers

// common libraries
#include <iostream>
#include <vector>

#define BM_GLU      2000
#define BM_GLU_NEU   500
#define BM_SQU       800
#define BM_SQU_NEU   100
#define BM_SQU_NEU2 1000
#define BM_EWK       800
#define BM_EWK_NEU   100

class PlottingBase
{
	public:
		PlottingBase(Variables& var) : v(var) {}
		~PlottingBase() {}

		SmartHistos sh;

		void define_histo_settings(const Weighting&, EventSelections&, const unsigned int&, const unsigned int&);

		void define_analysis_histos(const Weighting&, const unsigned int&, const unsigned int&);

		void calculate_plotting_variables(eventBuffer&, const unsigned int&);

		void fill_analysis_histos(EventSelections&, const Weighting&, const unsigned int&, const double&);

		void load_analysis_histos(std::string);

		void save_analysis_histos(bool);

		std::map<std::string, std::function<int()> > top_benchmarks;
		std::map<std::string, std::function<int()> > WZH_benchmarks;
		std::map<std::string, std::function<int()> > all_benchmarks;

		Variables& v;

		std::vector<std::string> all_cuts;

		typedef struct Sample { std::string postfix; std::string legend; std::string color; std::vector<std::string> dirs; } Sample;
		typedef struct PostfixOptions { size_t index; std::string postfixes; std::string legends; std::string colors; } PostfixOptions;

		typedef EventSelections::Regions Region;

		PostfixOptions get_pf_opts_(std::vector<std::vector<Sample> >, const TString&);

		void add_unrolled_bins(std::string, std::string, std::string,
				std::function<double()>, std::function<double()>,
				std::vector<double>, std::vector<double>,
				std::vector<int>, int, int, bool);



		std::map<int, std::vector<TH2D*> > m_vh_signal;
		std::map<int, std::vector<TH2D*> > m_vh_signal_new;
		//std::map<int, std::vector<TH2D*> > m_vh_signal_v2;

		// --------------------------------------------------------------------
		//                       Plotting options
		// --------------------------------------------------------------------

		std::string d = "HISTE1";
		// Control region stack plots (only data)
		std::string o_stk_d = "LogStack1AddRatioTwoCol98AddIntApproval15";
		std::string o_stk_s = "LogStack1AddRatioTwoCol98AddIntApproval45";
		std::string O_stk_d = "LogStack1AddRatioTwoCol98Approval15";
		std::string O_stk_s = "LogStack1AddRatioTwoCol98Approval45";
		// Signal region stack plots (data + N signal)
		std::string o_stk_d_S = "LogStack2AddRatioTwoCol108AddIntApproval15";
		std::string o_stk_s_S = "LogStack7AddRatioTwoCol108AddIntApproval45";
		std::string O_stk_d_S = "LogStack7AddRatioTwoCol108Approval15";
		std::string O_stk_s_S = "LogStack7AddRatioTwoCol108Approval45";
		std::string o_stk_d_T = "LogStack4AddRatioTwoCol48AddIntApproval15";
		std::string o_stk_s_T = "LogStack4AddRatioTwoCol48AddIntApproval45";
		std::string o_stk_d_V = "LogStack7AddRatioTwoCol78AddIntApproval15";
		std::string o_stk_s_V = "LogStack7AddRatioTwoCol78AddIntApproval45";
		std::string o_stk_d_H = "LogStack5AddRatioTwoCol58AddIntApproval15";
		std::string o_stk_s_H = "LogStack5AddRatioTwoCol58AddIntApproval45";
		std::string o_1or2d_d = "Approval15";
		std::string o_1or2d_s = "Approval45";
		std::string o_norm_d = "NormApproval15";
		std::string o_norm_s = "NormApproval45";

		std::vector<double> r_Stk4 = {0,0, 1.01e-2,1e4,  0.3,0.86};
		std::vector<double> r_Stk5 = {0,0, 1.01e-2,1e5,  0.3,0.86};
		std::vector<double> r_Stk6 = {0,0, 1.01e-2,1e6,  0.3,0.86};
		std::vector<double> r_Stk7 = {0,0, 1.01e-2,1e7,  0.3,0.86};
		std::vector<double> r_Stk8 = {0,0, 1.01e-2,1e8,  0.3,0.86};
		std::vector<double> r_Stk9 = {0,0, 1.01e-1,1e9,  0.3,0.86};
		std::string Stack = "StackPlot";

		// --------------------------------------------------------------------
		//                            Colors
		// --------------------------------------------------------------------


		/*
			 Complete list of main ROOT colors:
https://root.cern.ch/doc/master/classTColor.html#C02

kWhite  = 0,   kBlack  = 1,   kGray    = 920,  kRed    = 632,  kGreen  = 416,
kBlue   = 600, kYellow = 400, kMagenta = 616,  kCyan   = 432,  kOrange = 800,
kSpring = 820, kTeal   = 840, kAzure   = 860,  kViolet = 880,  kPink   = 900
*/

		std::string Black   = "1,";
		std::string DGray   = "922,";
		std::string DRed    = "634,";
		std::string Red     = "633,";
		std::string DPink   = "902,";
		std::string Pink    = "901,";
		std::string Purple  = "617,";
		std::string DMagen  = "619,";
		std::string FMagen  = "607,"; // Fainter
		std::string Violet  = "881,";
		std::string DViolet = "882,";
		std::string DBlue   = "602,";
		std::string Blue    = "601,";
		std::string FBlue   = "593,"; // Fainter
		std::string DAzure  = "862,";
		std::string Azure   = "861,";
		std::string DCyan   = "435,";
		std::string Cyan    = "433,";
		std::string DTeal   = "844,";
		std::string Teal    = "841,";
		std::string DGreen  = "418,";
		std::string Green   = "417,";
		std::string DSpring = "822,";
		std::string Spring  = "825,";
		std::string DYellow = "402,";
		std::string Yellow  = "401,"; // Beware - might still be faint
		std::string kOrange = "800,"; // Yellowish
		std::string Orange  = "797,";
		std::string DBrown  = "803,";
		std::string Brown   = "793,";
		std::string Carrot  = "807,"; // Reddish orange

		std::string col3_red_to_blue = "633,618,601,"; // red, purple, blue
		std::string col4_red_to_cyan = "633,618,601,434,"; // Red, purple, blue, cyan
		std::string col4_cyan_to_red = "434,601,618,633,"; // Cyan, blue, purple, red
		std::string col5_green_to_red = "418,434,601,618,633,"; // green, cyan, blue, purple, red
		std::string col5_red_to_green = "633,618,601,434,418,"; // red, , purple, blue, cyan, green
		std::string col6_rainbow_dark = "601,434,418,402,633,618,"; // blue, cyan, green, yellow, red, purple
		std::string col8 = "1,601,434,418,402,807,633,618,"; // above plus black and orange
		std::string col10 = "4,6,2,800,402,417,433,9,618,633,";
		std::string col12 = "1,4,6,2,800,402,417,433,9,618,633,924,"; // black, blue, magenta, red, orange, darker yellow, darker green, darker cyan, blue-purple, dark purple, dark red
		std::string col12_rainbow = "402,416,433,600,617,632,802,813,833,863,883,892,"; // Go around first bright and then dark colors
		std::string col16 = DRed+DPink+DMagen+DViolet+DBlue+DAzure+DCyan+DTeal+DGreen+DSpring+DYellow+Orange+Carrot+Brown+Black+DGray;

		// --------------------------------------------------------------------
		//                          Binning
		// --------------------------------------------------------------------

		// Bins
		std::vector<double> E;
		std::vector<double> Pt;
		std::vector<double> PtG;
		std::vector<double> PtF;
		std::vector<double> PtO;
		std::vector<double> PtT;
		std::vector<double> PtFine;
		std::vector<double> PtPho;
		std::vector<double> PtPho2;
		std::vector<double> M;
		std::vector<double> MFine;
		std::vector<double> MW;
		std::vector<double> DeepB;
		std::vector<double> MDP;
		std::vector<double> DP;
		std::vector<double> M2;
		std::vector<double> M3;
		std::vector<double> NVTX;
		std::vector<double> MET;
		std::vector<double> HT;
		std::vector<double> HTB;
		std::vector<double> PtB;
		std::vector<double> LepPt;
		std::vector<double> LepEta;
		std::vector<double> HT_2D_bins_old;
		std::vector<double> Pt_2D_bins;
		std::vector<double> PtLow_2D_bins;
		std::vector<double> PtHigh_2D_bins;
		std::vector<double> PtPho_2D_bins;
		std::vector<double> ElePt_2D_bins;
		std::vector<double> MuPt_2D_bins;
		std::vector<double> HT_2D_bins_lep;
		std::vector<double> HT1_2D_bins_lep;
		std::vector<double> HT2_2D_bins_lep;
		std::vector<double> HT3_2D_bins_lep;
		// HT/MET Trigger turnons
		std::vector<double> HT_2D_bins;
		std::vector<double> HT1_2D_bins;
		std::vector<double> HT2_2D_bins;
		std::vector<double> HT3_2D_bins;
		std::vector<double> MET_2D_bins;
		// unrolled bin mergin
		std::vector<int> merged_razor_bins;
		std::vector<int> merged_razor_bins_lep;
		std::vector<int> merged_trigger_bins;
		std::vector<int> merged_trigger1_bins;
		std::vector<int> merged_trigger2_bins;
		std::vector<int> merged_trigger3_bins;

};

class Plotting : public PlottingBase
{
	public:
		Plotting(Variables& var) : 
			PlottingBase(var),
			v(var)
	{}
		~Plotting() {}

		Variables& v;

		void define_additional_histo_settings(const Weighting&, EventSelections&, const unsigned int&, const unsigned int&);

		void define_additional_histos(const Weighting&, const unsigned int&, const unsigned int&);

};


//_______________________________________________________
//      Define Histo options: Filling/Postfixes


PlottingBase::PostfixOptions
PlottingBase::get_pf_opts_(std::vector<std::vector<Sample> > lists, const TString& sample) {
	std::vector<Sample> samples;
	for (auto list : lists) samples.insert(samples.end(), list.begin(), list.end());
	PostfixOptions opt{ (size_t)-1, "", "", "" };
	for (size_t i=0; i<samples.size(); ++i) {
		// Find index of matching directory
		for (size_t j=0; j<samples[i].dirs.size(); ++j)
			if (samples[i].dirs[j] == sample.Data()) opt.index = i;
		opt.postfixes += samples[i].postfix;
		opt.legends += samples[i].legend;
		opt.colors += samples[i].color;
		if (i+1!=samples.size()) {
			opt.postfixes +=  ";";
			opt.legends += ";";
			opt.colors += ",";
		}
	}
	return opt;
}

void
PlottingBase::add_unrolled_bins(std::string param_name, std::string axis1_title, std::string axis2_title,
		std::function<double()> var1, std::function<double()> var2,
		std::vector<double> bins1, std::vector<double> bins2,
		std::vector<int> merged_bins = {}, int precision1=0, int precision2=0, bool shrink_merged=false) {
	//std::cout<<param_name<<std::endl;
	std::vector<double> xbins = { 0 };
	std::map<int, std::string> bin1_labels, bin2_labels;
	std::stringstream label1, label2;
	size_t bin = 0, actual_bin = 0, big_bin = 0;
	for (size_t i=0, n=bins1.size(); i<n; ++i) {
		if (i!=0) {
			if (precision1) label1<<std::setprecision(precision1);
			label1<<bins1[i]<<"]";
			bin1_labels[big_bin] = label1.str();
		}
		if (i+1!=n){
			label1.str("");
			if (precision1) label1<<std::setprecision(precision1);
			label1<<"["<<bins1[i]<<", ";
			big_bin = shrink_merged ? actual_bin : bin;
			for (size_t j=0, m=bins2.size(); j<m; ++j) {
				if (j!=0) bin++;
				// Check if this bin is meant to be merged with the next one
				bool merge = false;
				for (auto merged_bin : merged_bins) if (bin==merged_bin) merge = true;
				if (!merge) {
					// Create bin label
					if (j!=0) {
						if (precision2) label2<<std::setprecision(precision2);
						label2<<bins2[j]<<"]";
						bin2_labels[++actual_bin] = label2.str();
						if (!merged_bins.empty()) xbins.push_back(shrink_merged ? actual_bin : bin);
					}
					if (j+1!=m){
						label2.str("");
						if (precision2) label2<<std::setprecision(precision2);
						label2<<"["<<bins2[j]<<", ";
					}
				}
			}
		}
	}
	if (merged_bins.empty()) xbins.push_back(bin);
	//for (auto xbin : xbins)
	//  std::cout<<xbin<<", ";
	//std::cout<<std::endl<<std::endl;
	if (shrink_merged) {
		sh.AddNewUnrolledFillParams(param_name, { .nbin=actual_bin, .bins=xbins, .bin_labels=bin2_labels, .unrolled_bin_labels=bin1_labels, .fill=[bins1,bins2,var1,var2,merged_bins] {
				double bin = -1, v1 = var1(), v2 = var2();
				for (size_t i=0, n=bins1.size(); i+1<n; ++i) if (v1>=bins1[i]&&v1<bins1[i+1])
				for (size_t j=0, m=bins2.size(); j+1<m; ++j) if (v2>=bins2[j]&&v2<bins2[j+1])
				bin = i*(m-1)+j;
				for (auto merged_bin : merged_bins) if (bin>=merged_bin) --bin;
				return bin;
				}, .axis_title=axis2_title, .unrolled_axis_title=axis1_title});
	} else {
		sh.AddNewUnrolledFillParams(param_name, { .nbin=actual_bin, .bins=xbins, .bin_labels=bin2_labels, .unrolled_bin_labels=bin1_labels, .fill=[bins1,bins2,var1,var2] {
				double bin = -1, v1 = var1(), v2 = var2();
				for (size_t i=0, n=bins1.size(); i+1<n; ++i) if (v1>=bins1[i]&&v1<bins1[i+1])
				for (size_t j=0, m=bins2.size(); j+1<m; ++j) if (v2>=bins2[j]&&v2<bins2[j+1])
				bin = i*(m-1)+j;
				return bin;
				}, .axis_title=axis2_title, .unrolled_axis_title=axis1_title});
	}
}

	void
PlottingBase::define_histo_settings(const Weighting& w, EventSelections& evt_sel,
		const unsigned int& syst_nSyst, const unsigned int& syst_index)
{
	const int debug = 0;

	if (debug) std::cout<<"PlottingBase::define_histo_settings: start"<<std::endl;

	sh.SetHistoWeights({ [&w] { return w.weight; } });

	// Keep this to be able to use analysis cuts
	//define_preselections(d);
	//define_selections(d);
	//for (const auto& region : evt_sel.analysis_cuts) w_nm1[region.first].resize(20,1);

	if (debug) std::cout<<"PlottingBase::define_histo_settings: weight, selections ok"<<std::endl;


	// --------------------------------------------------------------------
	//                            Cuts
	// --------------------------------------------------------------------

	// Pass each cut
	typedef EventSelections::Regions Region;
	//sh.AddNewCut("P", [&evt_sel] { return evt_sel.pass_all_cuts[static_cast<size_t>(Region::Pre)]; });







	if (debug) std::cout<<"PlottingBase::define_histo_settings: cuts ok"<<std::endl;

	//__________________________________
	//            Postfixes
	//__________________________________

	// Postfixes are vector definitions for histograms
	// They attach _<string> after histogram names
	// where <string> is chosen from a vector of strings
	// you need to give a natural number (0,1,2,3 ...) as a vector index
	// for the histogram to choose which histo to fill
	// The index should be size_t
	// In case the value is (size_t)-1 - the histogram won't be filled
	// This can effectively be used to define named cuts
	// The advantage over regular "Cuts" is that the postfix is attached to the histogram name
	// and the text is added to the legend (For this reason I prefer to use this over Cuts)

	// Notation:
	// AddNewPostfix("Name of Postfix", lambda function returning non-negative integer, "postfix 1;2;3; ...", "legend text 1;2;3; ...", "ROOT color 1,2,3, ...")

	// Define postfix - which can split the data to multiple plots
	// Parameter should be the index of the postfix vector
	// return (size_t)-1; --> don't fill histo, can be used to define a cut
	// 1st param: name of the postfix
	// 2nd param: lambda capture of the variable we want to fill
	// 3rd param: list of postixes that are attached to the histogram name, separate by ";"
	// 4th param: list of legend entrys, separate also by ";"
	// 3-4th param: can also use ranges for repetitive numbers: Jet[1to4] = "Jet1;Jet2;Jet3;Jet4", [1to10++2] = "1;3;5;7;9"
	// 5th param: list of TColors for each histogram, separated by ","

	// Sample postfixes for stack plots
	// Determine them from the directory names in which the input files are
	// Map directory names to postfix name, legend entry and color

	// ---------------------- Backgrounds ----------------------

	std::vector<Sample> bkg_ttbars;
	bkg_ttbars.push_back({ .postfix= "TT_powheg_pythia8",        .legend="t#bar{t}+jets", .color=DGreen, .dirs={ 
			// 2016
			"TT_TuneCUETP8M2T4_13TeV-powheg-pythia8",
			"TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8",
			"TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8",
			"TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8",
			// 2017/2018
			"TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8",
			"TTToHadronic_TuneCP5_13TeV-powheg-pythia8",
			"TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8" } });

	if (debug) std::cout<<"PlottingBase::define_histo_settings: ok2"<<std::endl;
	std::vector<Sample> bkg_nonttbars;
	bkg_nonttbars.push_back({ .postfix="Multijet",   .legend="Multijet",                  .color=FMagen, .dirs={ 
			// 2016
			"QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8",
			"QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8",
			"QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8",
			"QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8",
			"QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8",
			"QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8",
			"QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8",
			"QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8",
			"DYJetsToQQ_HT180_13TeV-madgraphMLM-pythia8",
			"WJetsToQQ_HT180_13TeV-madgraphMLM-pythia8",
			"WWTo4Q_13TeV-powheg",
			"ZJetsToQQ_HT600toInf_13TeV-madgraph",
			"ZZTo4Q_13TeV_amcatnloFXFX_madspin_pythia8",
			// 2017/2018
			"QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8",
			"QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8",
			"QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8",
			"QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8",
			"QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8",
			"QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8",
			"QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8",
			"QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8",
			"DYJetsToQQ_HT180_13TeV_TuneCP5-madgraphMLM-pythia8",
			"WJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8",
			"WJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8",
			"WJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8",
			"WJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8",
			"ZJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8",
			"ZJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8",
			"ZJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8",
			"ZJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8",
			"WWTo4Q_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8",
			// 2018
			"QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8",
			"QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8",
			"QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8",
			"QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8",
			"QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8",
			"QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8",
			"QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8",
			"QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8",
			"WJetsToQQ_HT400to600_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8",
			"WJetsToQQ_HT600to800_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8",
			"WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8",
			"WWTo4Q_NNPDF31_TuneCP5_13TeV-powheg-pythia8"
	} });

	bkg_nonttbars.push_back({ .postfix="WToLNu",     .legend="W+jets",                    .color=Red, .dirs={
			// 2016
			"WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8",
			"WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8",
			"WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8",
			"WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8",
			"WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8",
			"WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8",
			"WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8",
			"WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8",
			// 2017/2018
			"WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8",
			"WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8",
			"WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8",
			"WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8",
			"WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8",
			"WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8",
			"WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8"
			} });
	bkg_nonttbars.push_back({ .postfix="ZToNuNu",    .legend="Z#rightarrow#nu#nu", .color=Cyan, .dirs={
			// 2016/2017/2018
			"ZJetsToNuNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8",
			"ZJetsToNuNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8",
			"ZJetsToNuNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8",
			"ZJetsToNuNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8",
			"ZJetsToNuNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8",
			"ZJetsToNuNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8",
			"ZJetsToNuNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8" } });
	bkg_nonttbars.push_back({ .postfix="Higgs",    .legend="H+jets", .color=Brown, .dirs={
			// 2016/2017/2018
			"WminusH_HToBB_WToLNu_M-125_TuneCP5_13TeV_powheg_pythia8",
			"WminusH_HToBB_WToQQ_M-125_TuneCP5_13TeV_powheg_pythia8",
			"WplusH_HToBB_WToLNu_M-125_TuneCP5_13TeV_powheg_pythia8",
			"WplusH_HToBB_WToQQ_M-125_TuneCP5_13TeV_powheg_pythia8",
			"ZH_HToBB_ZToBB_M-125_TuneCP5_13TeV-powheg-pythia8",
			"ZH_HToBB_ZToLL_M-125_TuneCP5_13TeV_powheg_pythia8",
			"ZH_HToBB_ZToNuNu_M-125_TuneCP5_13TeV_powheg_pythia8",
			"ZH_HToBB_ZToQQ_M-125_TuneCP5_13TeV_powheg_pythia8",
			"ggZH_HToBB_ZToBB_M-125_TuneCP5_13TeV-powheg-pythia8",
			"ggZH_HToBB_ZToLL_M-125_TuneCP5_13TeV_powheg_pythia8",
			"ggZH_HToBB_ZToNuNu_M-125_TuneCP5_13TeV_powheg_pythia8",
			"ggZH_HToBB_ZToQQ_M-125_TuneCP5_13TeV_powheg_pythia8",
			"ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8",
			"ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8",
			"ttHJetToNonbb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8",
			"ttHJetTobb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8",
			// 2016
			"TTHH_TuneCP5_13TeV-amcatnlo-pythia8",
			//"VBFHToBB_M125_13TeV_amcatnlo_pythia8",
			"VBFHToBB_M-125_13TeV_powheg_pythia8_weightfix",
			"VBFHToWWTo2L2Nu_M125_13TeV_powheg_JHUgenv628_pythia8",
			"VBF_HToInvisible_M125_13TeV_powheg_pythia8",
			"TTWH_TuneCP5_13TeV-madgraph-pythia8",
			"TTZH_TuneCP5_13TeV-madgraph-pythia8",
			// 2017/2018
			"TTHH_TuneCP5_13TeV-madgraph-pythia8",
			"VBFHToBB_M-125_TuneCP5_13TeV-powheg-pythia8",
			"VBFHToWWTo2L2Nu_M-125_TuneCP5_13TeV-powheg-jhugen727-pythia8",
			// 2018
			"ZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8" } });
	bkg_nonttbars.push_back({ .postfix="Multiboson", .legend="VV(V)+t#bar{t}X",                  .color=DGray, .dirs={
			// 2016
			"WWTo1L1Nu2Q_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8",
			"WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8",
			"WWTo4Q_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8",
			"WZTo1L3Nu_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8",
			"WZTo1L1Nu2Q_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8",
			"WZTo2Q2Nu_13TeV_amcatnloFXFX_madspin_pythia8",
			"ZZTo2Nu2Q_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8",
			"ZZTo4Q_13TeV_amcatnloFXFX_madspin_pythia8",
			"ZZTo4L_TuneCP5_13TeV_powheg_pythia8",
			"TTTW_TuneCP5_13TeV-madgraph-pythia8",
			"TTWW_TuneCP5_13TeV-madgraph-pythia8",
			"TTWZ_TuneCP5_13TeV-madgraph-pythia8",
			// 2017
			"WWTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8",
			"WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8",
			"WWToLNuQQ_NNPDF31_TuneCP5_13TeV-powheg-pythia8",
			"WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8",
			"WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8_v2",
			"WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8",
			"WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8",
			"ZZTo2L2Nu_13TeV_powheg_pythia8",
			"ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8",
			"ZZTo2Q2Nu_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8",
			"ZZTo4L_13TeV_powheg_pythia8",
			"WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8",
			"WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8",
			"WZZ_TuneCP5_13TeV-amcatnlo-pythia8",
			"ZZZ_TuneCP5_13TeV-amcatnlo-pythia8",
			"ZZTo4L_M-1toInf_TuneCP5_13TeV_powheg_pythia8",
			"TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8",
			"TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8",
			"TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8",
			"TTZToNuNu_TuneCP5_13TeV-amcatnlo-pythia8",
			"TTZToQQ_TuneCP5_13TeV-amcatnlo-pythia8",
			"TTZZ_TuneCP5_13TeV-madgraph-pythia8",
			"TTTT_TuneCP5_13TeV-amcatnlo-pythia8",
			// 2018
			"WWZ_TuneCP5_13TeV-amcatnlo-pythia8",
			"WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8",
			"ZZTo2L2Nu_TuneCP5_13TeV_powheg_pythia8",
			"ZZTo4L_13TeV_powheg_pythia8_TuneCP5"
	} });
	bkg_nonttbars.push_back({ .postfix="Top",        .legend="Single t",                               .color=Orange,   .dirs={
			// 2016
			//"ST_s-channel_4f_InclusiveDecays_13TeV-amcatnlo-pythia8",
			"ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8",
			"ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1",
			"ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1",
			"ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1",
			"ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1",
			// 2017
			"ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8",
			"ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8",
			"ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8",
			"ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8",
			"ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8",
			// 2018
			"ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-madgraph-pythia8",
			"ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8",
			"ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8"
			} });
	//bkg_nonttbars.push_back({ .postfix="DYToLL",     .legend="Drell-Yan",                              .color=DBrown, .dirs={ 
	bkg_nonttbars.push_back({ .postfix="DYToLL",     .legend="Z#rightarrowll",                         .color=FBlue,  .dirs={
			// 2016
			"DYJetsToLL_M-5to50_HT-70to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
			"DYJetsToLL_M-5to50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
			"DYJetsToLL_M-5to50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
			"DYJetsToLL_M-5to50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
			"DYJetsToLL_M-5to50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
			"DYJetsToLL_M-50_HT-70to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
			"DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
			"DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
			"DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
			"DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
			"DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
			"DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
			"DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
			// 2017
			"DYJetsToLL_M-4to50_HT-70to100_TuneCP5_13TeV-madgraphMLM-pythia8",
			"DYJetsToLL_M-4to50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8",
			"DYJetsToLL_M-4to50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8",
			"DYJetsToLL_M-4to50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8",
			"DYJetsToLL_M-4to50_HT-600toInf_TuneCP5_13TeV-madgraphMLM-pythia8",
			"DYJetsToLL_M-50_HT-70to100_TuneCP5_13TeV-madgraphMLM-pythia8",
			"DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8",
			"DYJetsToLL_M-50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8",
			"DYJetsToLL_M-50_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8",
			"DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8",
			"DYJetsToLL_M-50_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8",
			"DYJetsToLL_M-50_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8",
			"DYJetsToLL_M-50_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8",
			// 2018
			"DYJetsToLL_M-4to50_HT-70to100_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8",
			"DYJetsToLL_M-4to50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8",
			"DYJetsToLL_M-4to50_HT-200to400_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8",
			"DYJetsToLL_M-4to50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8",
			"DYJetsToLL_M-4to50_HT-600toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8",
			"DYJetsToLL_M-50_HT-70to100_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8",
			"DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8",
			"DYJetsToLL_M-50_HT-200to400_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8",
			"DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8",
			"DYJetsToLL_M-50_HT-600to800_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8",
			"DYJetsToLL_M-50_HT-1200to2500_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8",
			"DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8",
			"DYJetsToLL_M-50_HT-800to1200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8"
	} });
	//bkg_nonttbars.push_back({ .postfix="GJets",      .legend="#gamma+jets",                            .color=DBrown,.dirs={ 
	bkg_nonttbars.push_back({ .postfix="GJets",      .legend="#gamma+jets",                            .color=kOrange, .dirs={
			// 2016
			"GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
			"GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
			"GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
			"GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
			"GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
			// 2017/2018
			"GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8",
			"GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8",
			"GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8",
			"GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8",
			"GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8"
			} });

	std::vector<Sample> bkg_all, bkg_selected, ttbar_selected, multijet, znunu, gjets;
	bkg_all.insert(bkg_all.end(), bkg_ttbars.begin(), bkg_ttbars.end());
	bkg_all.insert(bkg_all.end(), bkg_nonttbars.begin(), bkg_nonttbars.end());
	ttbar_selected.push_back(bkg_ttbars[0]); // powheg
	bkg_selected  .push_back(bkg_ttbars[0]); // powheg
	bkg_selected.insert(bkg_selected.end(), bkg_nonttbars.begin(), bkg_nonttbars.end());
	multijet.push_back(bkg_nonttbars[0]);
	znunu   .push_back(bkg_nonttbars[2]);
	gjets   .push_back(bkg_nonttbars[6]);

	std::vector<std::string> background_dirs;
	for (auto bkg : bkg_selected) for (auto dir : bkg.dirs) background_dirs.push_back(dir);
	std::vector<Sample> background, mc;
	background.push_back({ .postfix="Background", .legend="Background", .color=Black, .dirs=background_dirs });
	mc        .push_back({ .postfix="MC",         .legend="MC",         .color=Black, .dirs=background_dirs });
	static const PostfixOptions background_opt = get_pf_opts_({background}, v.sample);

	// ------------------------- Data --------------------------

	if (debug) std::cout<<"PlottingBase::define_histo_settings: ok3"<<std::endl;
	std::vector<Sample> data_all, data_selected, single_ele, single_mu, single_pho, jetht, met;
	data_all.push_back({ .postfix="Data",      .legend="Data",             .color=Black, .dirs={
			// Signal regions
			"JetHT_Run2016B",
			"JetHT_Run2016C",
			"JetHT_Run2016D",
			"JetHT_Run2016E",
			"JetHT_Run2016F",
			"JetHT_Run2016G",
			"JetHT_Run2016H",
			"JetHT_Run2017B",
			"JetHT_Run2017C",
			"JetHT_Run2017D",
			"JetHT_Run2017E",
			"JetHT_Run2017F",
			"JetHT_Run2018A",
			"JetHT_Run2018B",
			"JetHT_Run2018C",
			"JetHT_Run2018D",
			"MET_Run2016B",
			"MET_Run2016C",
			"MET_Run2016D",
			"MET_Run2016E",
			"MET_Run2016F",
			"MET_Run2016G",
			"MET_Run2016H",
			"MET_Run2017B",
			"MET_Run2017C",
			"MET_Run2017D",
			"MET_Run2017E",
			"MET_Run2017F",
			"MET_Run2018A",
			"MET_Run2018B",
			"MET_Run2018C",
			"MET_Run2018D",
			"HTMHT_Run2016B",
			"HTMHT_Run2016C",
			"HTMHT_Run2016D",
			"HTMHT_Run2016E",
			"HTMHT_Run2016F",
			"HTMHT_Run2016G",
			"HTMHT_Run2016H",
			// Control regions (and maybe leptonic signal)
			"SingleElectron_Run2016B",
			"SingleElectron_Run2016C",
			"SingleElectron_Run2016D",
			"SingleElectron_Run2016E",
			"SingleElectron_Run2016F",
			"SingleElectron_Run2016G",
			"SingleElectron_Run2016H",
			"SingleElectron_Run2017B",
			"SingleElectron_Run2017C",
			"SingleElectron_Run2017D",
			"SingleElectron_Run2017E",
			"SingleElectron_Run2017F",
			"EGamma_Run2018A",
			"EGamma_Run2018B",
			"EGamma_Run2018C",
			"EGamma_Run2018D",
			"SingleMuon_Run2016B",
			"SingleMuon_Run2016C",
			"SingleMuon_Run2016D",
			"SingleMuon_Run2016E",
			"SingleMuon_Run2016F",
			"SingleMuon_Run2016G",
			"SingleMuon_Run2016H",
			"SingleMuon_Run2017B",
			"SingleMuon_Run2017C",
			"SingleMuon_Run2017D",
			"SingleMuon_Run2017E",
			"SingleMuon_Run2017F",
			"SingleMuon_Run2018A",
			"SingleMuon_Run2018B",
			"SingleMuon_Run2018C",
			"SingleMuon_Run2018D"
				"SinglePhoton_Run2016B",
			"SinglePhoton_Run2016C",
			"SinglePhoton_Run2016D",
			"SinglePhoton_Run2016E",
			"SinglePhoton_Run2016F",
			"SinglePhoton_Run2016G",
			"SinglePhoton_Run2016H",
			"SinglePhoton_Run2017B",
			"SinglePhoton_Run2017C",
			"SinglePhoton_Run2017D",
			"SinglePhoton_Run2017E",
			"SinglePhoton_Run2017F"
	} });
	data_all.push_back({ .postfix="SingleEle", .legend="Data", .color=Black, .dirs={
			"SingleElectron_Run2016B",
			"SingleElectron_Run2016C",
			"SingleElectron_Run2016D",
			"SingleElectron_Run2016E",
			"SingleElectron_Run2016F",
			"SingleElectron_Run2016G",
			"SingleElectron_Run2016H",
			"SingleElectron_Run2017B",
			"SingleElectron_Run2017C",
			"SingleElectron_Run2017D",
			"SingleElectron_Run2017E",
			"SingleElectron_Run2017F",
			"EGamma_Run2018A",
			"EGamma_Run2018B",
			"EGamma_Run2018C",
			"EGamma_Run2018D"
			} });
	data_all.push_back({ .postfix="SingleMu",  .legend="Data",  .color=Black, .dirs={
			"SingleMuon_Run2016B",
			"SingleMuon_Run2016C",
			"SingleMuon_Run2016D",
			"SingleMuon_Run2016E",
			"SingleMuon_Run2016F",
			"SingleMuon_Run2016G",
			"SingleMuon_Run2016H",
			"SingleMuon_Run2017B",
			"SingleMuon_Run2017C",
			"SingleMuon_Run2017D",
			"SingleMuon_Run2017E",
			"SingleMuon_Run2017F",
			"SingleMuon_Run2018A",
			"SingleMuon_Run2018B",
			"SingleMuon_Run2018C",
			"SingleMuon_Run2018D"
			} });
	data_all.push_back({ .postfix="SinglePho",  .legend="Data",  .color=Black, .dirs={
			"SinglePhoton_Run2016B",
			"SinglePhoton_Run2016C",
			"SinglePhoton_Run2016D",
			"SinglePhoton_Run2016E",
			"SinglePhoton_Run2016F",
			"SinglePhoton_Run2016G",
			"SinglePhoton_Run2016H",
			"SinglePhoton_Run2017B",
			"SinglePhoton_Run2017C",
			"SinglePhoton_Run2017D",
			"SinglePhoton_Run2017E",
			"SinglePhoton_Run2017F"
			} });
	data_all.push_back({ .postfix="JetHT",       .legend="Data",     .color=Black, .dirs={
			"JetHT_Run2016B",
			"JetHT_Run2016C",
			"JetHT_Run2016D",
			"JetHT_Run2016E",
			"JetHT_Run2016F",
			"JetHT_Run2016G",
			"JetHT_Run2016H",
			"JetHT_Run2017B",
			"JetHT_Run2017C",
			"JetHT_Run2017D",
			"JetHT_Run2017E",
			"JetHT_Run2017F",
			"JetHT_Run2018A",
			"JetHT_Run2018B",
			"JetHT_Run2018C",
			"JetHT_Run2018D"
			} });
	data_all.push_back({ .postfix="MET",       .legend="Data",       .color=Black, .dirs={
			"MET_Run2016B",
			"MET_Run2016C",
			"MET_Run2016D",
			"MET_Run2016E",
			"MET_Run2016F",
			"MET_Run2016G",
			"MET_Run2016H",
			"MET_Run2017B",
			"MET_Run2017C",
			"MET_Run2017D",
			"MET_Run2017E",
			"MET_Run2017F",
			"MET_Run2018A",
			"MET_Run2018B",
			"MET_Run2018C",
			"MET_Run2018D"
			} });
	data_all.push_back({ .postfix="HTMHT",     .legend="Data",     .color=Black, .dirs={
			"HTMHT_Run2016B",
			"HTMHT_Run2016C",
			"HTMHT_Run2016D",
			"HTMHT_Run2016E",
			"HTMHT_Run2016F",
			"HTMHT_Run2016G",
			"HTMHT_Run2016H",
			"HTMHT_Run2017B",
			"HTMHT_Run2017C",
			"HTMHT_Run2017D",
			"HTMHT_Run2017E",
			"HTMHT_Run2017F",
			"HTMHT_Run2018A",
			"HTMHT_Run2018B",
			"HTMHT_Run2018C",
			"HTMHT_Run2018D"
			} });
	data_selected.push_back(data_all[0]);
	single_ele.push_back(data_all[1]);
	single_mu.push_back(data_all[2]);
	single_pho.push_back(data_all[3]);
	jetht.push_back(data_all[4]);
	met.push_back(data_all[5]);

	// Individual Cuts implemented as Postfixes
	sh.AddNewPostfix("Blind",      [this] { 
			if (v.isData) return (size_t)-1;
			return (size_t)0;
			}, "BlindData", "", Black);

	// ------------------------ Signals ------------------------

	if (debug) std::cout<<"PlottingBase::define_histo_settings: ok4"<<std::endl;
	std::vector<Sample> signal_all;

	signal_all.push_back({ .postfix="ttbarToBsToTauTau",           .legend="ttbarToBsToTauTau",        .color=DCyan,     .dirs={
			"ttbarToBsToTauTau_BsFilter_TauTauFilter_TuneCP5_13TeV-pythia8-evtgen" } });


	static const PostfixOptions plot_samples_opt=get_pf_opts_({data_selected, signal_all, bkg_selected}, v.sample);
	sh.AddNewPostfix("StackPlotSignal", [] { 
			return plot_samples_opt.index; 
			}, plot_samples_opt.postfixes, plot_samples_opt.legends, plot_samples_opt.colors);


	// ---------------------- Combinations ---------------------

	// Sample postfixes
	if (debug) std::cout<<"PlottingBase::define_histo_settings: ok5"<<std::endl;
	static const PostfixOptions all_samples_opt=get_pf_opts_({data_all, bkg_all, signal_all}, v.sample);
	sh.AddNewPostfix("AllSamples", [] { return all_samples_opt.index; }, all_samples_opt.postfixes, all_samples_opt.legends, all_samples_opt.colors);

	static const PostfixOptions plot_samples2_opt=get_pf_opts_({data_selected, bkg_selected}, v.sample);
	sh.AddNewPostfix("StackPlot", [] { 
			return plot_samples2_opt.index; 
			}, plot_samples2_opt.postfixes, plot_samples2_opt.legends, plot_samples2_opt.colors);

	static const PostfixOptions data_opt=get_pf_opts_({data_selected}, v.sample);
	sh.AddNewPostfix("Data", [] { return data_opt.index;  }, data_opt.postfixes, data_opt.legends, data_opt.colors);

	sh.AddNewPostfix("Background",  [] { return background_opt.index; }, background_opt.postfixes, background_opt.legends, background_opt.colors);
	static const PostfixOptions mc_opt = get_pf_opts_({mc}, v.sample);
	sh.AddNewPostfix("MC",  [] { return mc_opt.index; }, mc_opt.postfixes, mc_opt.legends, mc_opt.colors);

	static const PostfixOptions bkgs_opt=get_pf_opts_({bkg_selected}, v.sample);
	sh.AddNewPostfix("Backgrounds", [] { return bkgs_opt.index;  }, bkgs_opt.postfixes, bkgs_opt.legends, bkgs_opt.colors);

	static const PostfixOptions multijet_opt=get_pf_opts_({multijet}, v.sample);
	sh.AddNewPostfix("Multijet", [] { return multijet_opt.index;  }, multijet_opt.postfixes, multijet_opt.legends, multijet_opt.colors);
	static const PostfixOptions znunu_opt=get_pf_opts_({znunu}, v.sample);
	sh.AddNewPostfix("ZToNuNu", [] { return znunu_opt.index;  }, znunu_opt.postfixes, znunu_opt.legends, znunu_opt.colors);
	static const PostfixOptions gjets_opt=get_pf_opts_({gjets}, v.sample);
	sh.AddNewPostfix("GJets", [] { return gjets_opt.index;  }, gjets_opt.postfixes, gjets_opt.legends, gjets_opt.colors);

	if (debug) std::cout<<"PlottingBase::define_histo_settings: ok7"<<std::endl;
	static const PostfixOptions data_mc_opt = get_pf_opts_({data_selected, mc}, v.sample);
	sh.AddNewPostfix("Data_MC",  [] { return data_mc_opt.index; }, data_mc_opt.postfixes, data_mc_opt.legends, "1,633");
	//static const PostfixOptions data_fastsim_opt = get_pf_opts_({data_selected, signal_all}, v.sample);


	static const PostfixOptions triggers_opt = get_pf_opts_({single_ele, single_pho, single_mu, mc}, v.sample);

static const PostfixOptions triggers3_opt = get_pf_opts_({jetht, met, mc}, v.sample);
sh.AddNewPostfix("LeptonicMeasurements", [this]()
		{
		if (!(v.nLepTight==1)) return (size_t)-1;
		if (triggers3_opt.index==0) { // JetHT
		bool OR_HLT_Had = 
		v.HLT_PFHT125 ||
		v.HLT_PFHT180 ||
		v.HLT_PFHT200 ||
		v.HLT_PFHT250 ||
		v.HLT_PFHT300 ||
		v.HLT_PFHT350 ||
		v.HLT_PFHT370 ||
		v.HLT_PFHT400 ||
		v.HLT_PFHT430 ||
		v.HLT_PFHT475 ||
		v.HLT_PFHT510 ||
		v.HLT_PFHT590 ||
		v.HLT_PFHT650 ||
		v.HLT_PFHT680 ||
		v.HLT_PFHT780 ||
		v.HLT_PFHT800 ||
		v.HLT_PFHT890 ||
		v.HLT_PFHT900 ||
		v.HLT_PFHT1050 ||
		v.HLT_AK8PFJet40 ||
		v.HLT_AK8PFJet60 ||
		v.HLT_AK8PFJet80 ||
		v.HLT_AK8PFJet140 ||
		v.HLT_AK8PFJet200 ||
		v.HLT_AK8PFJet260 ||
		v.HLT_AK8PFJet320 ||
		v.HLT_AK8PFJet400 ||
		v.HLT_AK8PFJet450 ||
		v.HLT_AK8PFJet500 ||
		v.HLT_AK8PFJet550;
		if (!OR_HLT_Had) return (size_t)-1;
		if (v.Jet.Jet.n<3) return (size_t)-1;
		} else if (triggers3_opt.index==1) { // MET
			if (!(v.HLT_PFMET110_PFMHT110_IDTight==1||v.HLT_PFMETNoMu110_PFMHTNoMu110_IDTight==1||
						v.HLT_PFMET120_PFMHT120_IDTight==1||v.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight==1||
						v.HLT_PFMET120_PFMHT120_IDTight_PFHT60==1||v.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60==1)) return (size_t)-1;
			if (v.Jet.Jet.n<3) return (size_t)-1;
		} else if (triggers3_opt.index==2) { // Simulation
			if (v.Jet.Jet.n<3) return (size_t)-1;
		} else {
			return (size_t)-1;
		}
		return (size_t)triggers3_opt.index; 
		}, "JetHT_1Lep;MET_1Lep;MC_1Lep", "HT;MET;Simulation", Carrot+Green+Black);

static const PostfixOptions trigger_opt = get_pf_opts_({single_ele, met}, v.sample);
sh.AddNewPostfix("EleMETComb", [this]()
		{
		if (trigger_opt.index==0) { // SingleElectron
		if ((v.HLT_Ele27_WPTight_Gsf==1||v.HLT_Ele30_WPTight_Gsf==1 ||
					v.HLT_Ele32_WPTight_Gsf==1)&&
				v.Electron.Select.n==1&&v.Muon.Veto.n==0&&v.Tau.Select.n==0) return (size_t)0;
		} else if (trigger_opt.index==1) { // MET
		if (v.HLT_PFMET120_PFMHT120_IDTight==1&&v.nLepVeto==0&&v.Tau.Select.n==0) return (size_t)0;
		}
		return (size_t)-1; 
		}, "SingleEle_MET", "SingleEle + MET", Black);

static const PostfixOptions data_photruth_opt = get_pf_opts_({data_selected, mc}, v.sample);
sh.AddNewPostfix("Data_PhotonTruth",  [this] { 
		if (data_photruth_opt.index==0) return (size_t)0;
		else if (data_photruth_opt.index==1) {
		if      (v.Photon().matchGenPromptDirect) return (size_t)1;
		else if (v.Photon().matchGenPromptFrag)   return (size_t)2;
		else if (v.Photon().matchGenFake)         return (size_t)3;
		}
		return (size_t)-1; }, "Data;PromptDirect;PromptFrag;Fake", "Data;Prompt;Fragm.;Fakes", "1,418,601,633");
sh.AddNewPostfix("Data_Photon1Truth",  [this] { 
		if (v.Photon.Select.n!=1) return (size_t)-1;
		if (data_photruth_opt.index==0) return (size_t)0;
		else if (data_photruth_opt.index==1) {
		if      (v.Photon.Select(0).matchGenPromptDirect) return (size_t)1;
		else if (v.Photon.Select(0).matchGenPromptFrag)   return (size_t)2;
		else if (v.Photon.Select(0).matchGenFake)         return (size_t)3;
		}
		return (size_t)-1; }, "Data;PromptDirect;PromptFrag;Fake", "Data;Prompt;Fragm.;Fakes", "1,418,601,633");

sh.AddNewPostfix("run2", [] { return  (size_t)0; }, "run2;", "run2;", "418");
sh.AddNewPostfix("Year", [this]() { return (size_t)(v.year-2016); }, "2016;2017;2018", "2016;2017;2018", "418,601,633");
sh.AddNewPostfix("YearSF", [this]() { return v.isAPV ? (size_t)0 : (size_t)(v.year-2016+1); }, "2016APV;2016;2017;2018", "2016APV;2016;2017;2018", "1,418,601,633");
sh.AddNewPostfix("YearFS", [this]() { return (size_t)(v.year-2016); }, "2016;2017;2018", "2016;2017;2018", "418,601,633");


if (debug) std::cout<<"PlottingBase::define_histo_settings: sample postfixes ok"<<std::endl;

// ------------------ Signal/Control regions ---------------

// Cut names
std::map<std::string, std::string> legname;
legname["3Jet"]        = "Njet#geq3";
legname["HLT"]         = "HLT";
legname["0Ele"]        = "ele veto";
legname["0Mu"]         = "muon veto";
legname["0Tau"]        = "tau veto";
legname["1b"]          = "Nb#geq1";
legname["1W"]          = "NW#geq1";
legname["dPhi"]        = "#Delta#phi";
legname["InvdPhi"]     = "inv. #Delta#phi";
legname["0b"]          = "b-tag veto";
legname["1aW"]         = "NW(anti-tag)#geq1";
legname["1aTop"]       = "NTop(anti-tag)#geq1";
//legname["InvdPhi0p3"]  = "#Delta#phi<0.3";
legname["InvdPhi"]     = "inv. #Delta#phi";
legname["1Lep"]        = "Nlep=1";
legname["MT"]          = "m_{T}";
legname["1mW"]         = "NW(mass-tag)#geq1";
legname["1mTop"]       = "NTop(mass-tag)#geq1";
legname["1MTop"]       = "NTop(mass-tag)#geq1";
legname["R2ll"]        = "R^{2}";
legname["2Lep"]        = "Nlep=2";
legname["OppCharge"]   = "#sumq_{lep}=0";
legname["dPhi"]        = "#Delta#phi";
legname["Mll"]         = "|m_{ll} - m_{Z}| < 10 GeV";
legname["1Top"]        = "Ntop#geq1";
legname["1Pho"]        = "N#gamma=1";
std::map<Region, std::string> regionname;
regionname[Region::Pre_1Lep]          = "Semi-Leptonic Baseline selection";
regionname[Region::Pre_2Lep]          = "Di-Leptonic Baseline selection";
regionname[Region::Pre_Had]          = "Hadronic Baseline selection";
regionname[Region::Pre_1Lep_MT]     = "Semi-Leptonic Baseline selection + MT";
if (debug) std::cout<<"PlottingBase::define_histo_settings: region names ok"<<std::endl;

// Cut Postfixes
sh.AddNewPostfix("BaselineCuts", [] { return 0; }, "BaselineCuts", "Baseline cuts", Black);
all_cuts.push_back("BaselineCuts");
for (const auto& region : magic_enum::enum_entries<Region>()) {
	std::string cutflow_str = "";
	std::string name(region.second);
	if (debug>1) std::cout<<"PlottingBase::define_histo_settings: "<<region.first<<" ("<<name<<") start"<<std::endl;
	sh.AddNewPostfix(name, [region,&evt_sel] { return evt_sel.apply_all_cuts(region.first) ? 0 : (size_t)-1; },
			name, regionname[region.first], Black);
	if (debug>1) std::cout<<"PlottingBase::define_histo_settings: "<<region.first<<" ("<<name<<") all cuts ok"<<std::endl;
	for (size_t i=0, n=evt_sel.analysis_cuts[region.first].size(); i<n; ++i) {
		// Cuts in order 1-N: "PassNCuts[search region]"
		if (debug>1) std::cout<<"PlottingBase::define_histo_settings: "<<i<<" start"<<std::endl;
		sh.AddNewPostfix(name+"_"+std::to_string(i+1)+"Cuts", [i,region,&evt_sel] { return evt_sel.apply_ncut(region.first, i) ? 0 : (size_t)-1; },
				name+"_"+std::to_string(i+1)+"Cuts", name+" region, first "+std::to_string(i+1)+" cuts", Black);
		all_cuts.push_back(name+"_"+std::to_string(i+1)+"Cuts");
		if (debug>1) std::cout<<"PlottingBase::define_histo_settings: "<<i<<" n cuts ok"<<std::endl;
		cutflow_str += evt_sel.analysis_cuts[region.first][i].name+name+";";
		// N-1 Cuts: "[search region]_Excl[cut]"
		sh.AddNewPostfix(name+"_Excl"+evt_sel.analysis_cuts[region.first][i].name, [i,region,&evt_sel] { 
				unsigned int mask = (1<<evt_sel.analysis_cuts[region.first].size())-1 - (1<<i); 
				return ((evt_sel.cutbits[region.first] & mask) == mask) ? 0 : (size_t)-1; }, 
				//name+"_Excl"+evt_sel.analysis_cuts[region.first][i].name, regionname[region.first]+", no "+
				//legname[evt_sel.analysis_cuts[region.first][i].name]+" cut", Black);
			name+"_Excl"+evt_sel.analysis_cuts[region.first][i].name, regionname[region.first], Black);
		if (debug>1) std::cout<<"PlottingBase::define_histo_settings: "<<i<<" N-1 cuts ok"<<std::endl;
		// N-2 Cuts: "[search region]_Excl[cut1][cut2]"
		for (size_t j=i+1, n=evt_sel.analysis_cuts[region.first].size(); j<n; ++j)
			sh.AddNewPostfix(name+"_Excl"+evt_sel.analysis_cuts[region.first][i].name+evt_sel.analysis_cuts[region.first][j].name, [i,j,region,&evt_sel] { 
					unsigned int mask = (1<<evt_sel.analysis_cuts[region.first].size())-1 - (1<<i) - (1<<j); 
					return ((evt_sel.cutbits[region.first] & mask) == mask) ? 0 : (size_t)-1; }, 
					name+"_Excl"+evt_sel.analysis_cuts[region.first][i].name+evt_sel.analysis_cuts[region.first][j].name, 
					regionname[region.first]+", no "+legname[evt_sel.analysis_cuts[region.first][i].name]+", "+
					legname[evt_sel.analysis_cuts[region.first][j].name]+" cut", Black);
		if (debug>1) std::cout<<"PlottingBase::define_histo_settings: "<<i<<" N-2 cuts ok"<<std::endl;
		// N-3 Cuts: "[search region]_Excl[cut1][cut2][cut3]"
		for (size_t j=i+1, n=evt_sel.analysis_cuts[region.first].size(); j<n; ++j)
			for (size_t k=j+1; k<n; ++k)
				sh.AddNewPostfix(name+"_Excl"+evt_sel.analysis_cuts[region.first][i].name+
						evt_sel.analysis_cuts[region.first][j].name+evt_sel.analysis_cuts[region.first][k].name, 
						[i,j,k,region,&evt_sel] { 
						unsigned int mask = (1<<evt_sel.analysis_cuts[region.first].size())-1 - (1<<i) - (1<<j) - (1<<k); 
						return ((evt_sel.cutbits[region.first] & mask) == mask) ? 0 : (size_t)-1; }, 
						name+"_Excl"+evt_sel.analysis_cuts[region.first][i].name+
						evt_sel.analysis_cuts[region.first][j].name+evt_sel.analysis_cuts[region.first][k].name, 
						regionname[region.first]+", no "+legname[evt_sel.analysis_cuts[region.first][i].name]+", "+
						legname[evt_sel.analysis_cuts[region.first][j].name]+", "+legname[evt_sel.analysis_cuts[region.first][k].name]+" cut", Black);
		if (debug>1) std::cout<<"PlottingBase::define_histo_settings: "<<i<<" N-3 cuts ok"<<std::endl;
	}
	// Stackable Cut Histos: "CutFlow"
	sh.AddNewPostfix("CutFlow"+name, [region,&evt_sel] { 
			for (size_t i=0, n=evt_sel.analysis_cuts[region.first].size(); i<n; ++i) 
			if (!evt_sel.analysis_cuts[region.first][i].func()) return i; 
			return evt_sel.analysis_cuts[region.first].size();
			}, cutflow_str+"PassAll"+name, cutflow_str+regionname[region.first], col10+col10);
	if (debug>1) std::cout<<"PlottingBase::define_histo_settings: "<<region.first<<" ("<<name<<") end"<<std::endl<<std::endl;
}
//sh.AddNewPostfix("Baseline",       [&evt_sel] { return evt_sel.apply_cut(Region::Pre, 0) ? 0 : (size_t)-1; }, "Baseline",      "Baseline selection",               Black);
if (debug) std::cout<<"PlottingBase::define_histo_settings: regions ok"<<std::endl<<std::endl;


// ------------------------- Objects -----------------------

// Lepton Postfixes
// Photon Postfixes
sh.AddNewPostfix("SIEIE",          [this] {  return v.Photon().sieie < (v.Photon().isScEtaEB ?  0.01015 : 0.0272); }, "FailSIEIE;PassSIEIE", "Fail #sigma_{i#eta i#eta};Pass #sigma_{i#eta i#eta}", Red+Green);
sh.AddNewPostfix("CHIso",          [this] {  return v.Photon().pfRelIso03_chg < (v.Photon().isScEtaEB ? 0.441 : 0.442); }, "FailCHIso;PassCHIso", "Fail ch. iso.;Pass ch. iso.", Red+Green);
sh.AddNewPostfix("EB",             [this] {  return (size_t)(v.Photon().isScEtaEB ? 0 : -1); }, "Barrel", "Barrel", Black);
sh.AddNewPostfix("EE",             [this] {  return (size_t)(v.Photon().isScEtaEE ? 0 : -1); }, "Endcap", "Endcap", Black);
sh.AddNewPostfix("EB_EE",          [this] {  return v.Photon().isScEtaEE; }, "Barrel;Endcap", "Barrel;Endcap", "601,418");
sh.AddNewPostfix("Fake",           [this] {  return (size_t)(v.Photon().matchGenFake ? 0 : -1); }, "Fake", "Fake", Black);
sh.AddNewPostfix("Prompt",         [this] {  return (size_t)(v.Photon().matchGenPrompt ? 0 : -1); }, "Prompt", "Prompt", Black);
if (debug) std::cout<<"PlottingBase::define_histo_settings: photon pfs ok"<<std::endl<<std::endl;

// AK4 Jet Postfixes
sh.AddNewPostfix("Jets",    [this] { size_t i=v.Jet.Jet.viSel[v.Jet.i];        return (i<4)?i:(size_t)-1; }, "Jet[1to5]",  "1st Jet;2nd Jet;3rd Jet;[4to5]th Jet", col5_red_to_green);
sh.AddNewPostfix("BTags",   [this] { size_t i=v.Jet.MediumBTag.viSel[v.Jet.i]; return (i<4)?i:(size_t)-1; }, "BTag[1to5]", "1st b;2nd b;3rd b;[4to5]th b",         col5_red_to_green);
sh.AddNewPostfix("JetPhotonDR0.05", [this] {  return v.Jet().phoDR<0.05 ? (size_t)0 : (size_t)-1; }, "PhotonDR0p05",  "DR_{#gamma}<0.05", Black);
sh.AddNewPostfix("JetPhotonDR0.4",  [this] {  return v.Jet().phoDR<0.4  ? (size_t)0 : (size_t)-1; }, "PhotonDR0p4",   "DR_{#gamma}<0.4",  Black);
sh.AddNewPostfix("JetEleDR0.05",    [this] {  return v.Jet().eleDR<0.05 ? (size_t)0 : (size_t)-1; }, "EleDR0p05",     "DR_{ele}<0.05", Black);
sh.AddNewPostfix("JetEleDR0.4",     [this] {  return v.Jet().eleDR<0.4  ? (size_t)0 : (size_t)-1; }, "EleDR0p4",      "DR_{ele}<0.4",  Black);
sh.AddNewPostfix("JetMuonDR0.05",   [this] {  return v.Jet().muDR<0.05 ? (size_t)0 : (size_t)-1; },  "MuonDR0p05",    "DR_{muon}<0.05", Black);
sh.AddNewPostfix("JetMuonDR0.4",    [this] {  return v.Jet().muDR<0.4  ? (size_t)0 : (size_t)-1; },  "MuonDR0p4",     "DR_{muon}<0.4",  Black);
if (debug) std::cout<<"PlottingBase::define_histo_settings: AK4 pfs ok"<<std::endl<<std::endl;

sh.AddNewPostfix("LostLeptonFlavour", [this] {
		int id = std::abs(v.GenPart().pdgId);
		if (id==11) return (size_t)0;
		if (id==13) return (size_t)1;
		if (id==15) return (size_t)2;
		return (size_t)-1;
		}, "GenEle;GenMu;GenTau", "Electron;Muon;Tau", "417,601,633");
sh.AddNewPostfix("GenLeptonFlavour", [this] {
		int id = std::abs(v.GenPart().pdgId);
		if (id==11) return (size_t)0;
		if (id==13) return (size_t)1;
		if (id==15) return (size_t)2;
		return (size_t)-1;
		}, "GenEle;GenMu;GenTau", "Electron;Muon;Tau", "417,601,633");
sh.AddNewPostfix("GenTopLeptonFlavour", [this] {
		size_t iLep = v.GenPart().iGenLepGrandDaughter;
		if (iLep==(size_t)-1) return (size_t)-1;
		int id = std::abs(v.GenPart(iLep).pdgId);
		if      (id==11) return (size_t)0;
		else if (id==13) return (size_t)1;
		else if (id==15) return (size_t)2;
		return (size_t)-1;
		}, "GenEle;GenMu;GenTau", "Electron;Muon;Tau", "417,601,633");
sh.AddNewPostfix("GenWLeptonFlavour", [this] {
		size_t iLep = v.GenPart().iGenLepDaughter;
		if (iLep==(size_t)-1) return (size_t)-1;
		int id = std::abs(v.GenPart(iLep).pdgId);
		if      (id==11) return (size_t)0;
		else if (id==13) return (size_t)1;
		else if (id==15) return (size_t)2;
		return (size_t)-1;
		}, "GenEle;GenMu;GenTau", "Electron;Muon;Tau", "417,601,633");
sh.AddNewPostfix("GenLepMother", [this] {
		if     (!v.GenPart.Lepton.pass[v.GenPart.i]) return (size_t)-1;
		if      (v.GenPart.LeptonFromTop.pass[v.GenPart.i]) return (size_t)2;
		else if (v.GenPart.LeptonFromW.pass[v.GenPart.i]) return (size_t)1;
		else return (size_t)0;
		}, "Other;FromW;FromTop", "Other;from W decay;from top decay", Red+Green);
if (debug) std::cout<<"PlottingBase::define_histo_settings: gen pfs ok"<<std::endl<<std::endl;

// --------------------- Event Variables -------------------

sh.AddNewPostfix("OtherNonisoLep", [this] { return std::min(v.nLepVetoNoIso-v.nLepSelect,(size_t)1); }, "NoOtherNonisoLep;OtherNonisoLep", "0 other unisol. lepton;#geq1 other unisol. lepton", Green+Red);
sh.AddNewPostfix("OtherLooseLep",  [this] { return std::min(v.nLepVeto     -v.nLepSelect,(size_t)1); }, "NoOtherLep;OtherLep",           "0 other loose lepton;#geq1 other loose lepton", Red+Green);
sh.AddNewPostfix("Ele_Muon",       [this] { return (size_t)(v.Electron.Veto.n==1 ? 0 : v.Muon.Veto.n==1 ? 1 : -1); }, "1VetoEle;1MuVeto", "1e (veto);1#mu (veto)", "1,2");
sh.AddNewPostfix("Ele_or_Muon",    [this] { return (size_t)(v.Electron.Select.n==1 ? 0 : v.Muon.Select.n==1 ? 1 : -1); }, "EleOnly;MuOnly", "1e;1#mu", "1,2");
sh.AddNewPostfix("0Lep",           [this] { return (size_t)(v.Electron.Select.n+v.Muon.Select.n==1 ? 0 : -1); }, "1Lep",  "1 lep",  Black);
sh.AddNewPostfix("Lep",            [this] { return (size_t)((v.nLepVeto+v.Tau.Select.n)==0 ? 0 : v.Electron.Select.n+v.Muon.Select.n==1 ? 1 : -1); }, "0Lep;1Lep",  "0 lep;1 lep",  Black+Red);
sh.AddNewPostfix("Tau",            [this] { return (size_t)(v.Tau.Select.n==0 ? 0 : v.Tau.Select.n==1 ? 1 : -1); }, "0Tau;1Tau",  "0 #tau;1 #tau",  Black+Red);
sh.AddNewPostfix("1Ele",           [this] { return (size_t)(v.Electron.Select.n==1 ? 0 : -1); }, "1Ele",  "1e",  Black);
sh.AddNewPostfix("1Muon",          [this] { return (size_t)(v.Muon.Select.n==1 ?  0 : -1); }, "1Muon", "1#mu", Black);
sh.AddNewPostfix("Pho",            [this] { return (size_t)(v.Photon.Select.n==0 ? 0 : v.Photon.Select.n==1 ? 1 : -1); }, "0Pho;1Pho",  "0#gamma;1#gamma",  Black+Red);
sh.AddNewPostfix("1Pho",           [this] { return (size_t)(v.Photon.Select.n==1 ? 0 : -1); }, "1Pho",  "1#gamma",  Black);
sh.AddNewPostfix("1Fake",          [this] { return (size_t)(v.Photon.Fake.n==1 ? 0 : -1); }, "1FakePho",  "1#gamma",  Black);
sh.AddNewPostfix("1PrePho",        [this] { return (size_t)(v.Photon.PreSelect.n==1 ? 0 : -1); }, "1PrePho",  "1 pre-#gamma",  Black);
sh.AddNewPostfix("FakePhoton",     [this] { return (size_t)(v.Photon.Fake.pass[v.Photon.i] ? 0 : -1); }, "FakePhoton",  "Fake #gamma",  Black);
sh.AddNewPostfix("PromptPhoton",   [this] { return (size_t)(v.Photon.SelectNoIso.pass[v.Photon.i]&&v.Photon().matchGenPrompt ? 0 : -1); }, "PromptPhoton",  "Prompt #gamma",  Black);
sh.AddNewPostfix("SelectNoIsoPhoton", [this] { return (size_t)(v.Photon.SelectNoIso.pass[v.Photon.i] ? 0 : -1); }, "PhotonNoIso",  "#gamma (no iso.)",  Black);
sh.AddNewPostfix("2Ele_2Muon",     [this] { return (size_t)(v.Electron.Select.n==2 ? 0 : v.Muon.Select.n==2 ? 1 : -1); }, "2Ele;2Mu", "2e;2#mu", "1,2");
sh.AddNewPostfix("DiLep",          [this] { return (size_t)((v.Electron.Veto.n+v.Muon.Veto.n)!=2 ? -1 : v.Electron.Veto.n==2 ? 0 : v.Muon.Veto.n==2 ? 1 :  2); }, "2Ele;2Mu;1Ele1Mu", "2e;2#mu;1e+1#mu", Red+Blue+Purple);
//sh.AddNewPostfix("PromptDirect",   [this] { return (size_t)(v.Photon.Select.n==1 ? (v.Photon.Select(0).matchGenPromptDirect ? 0 : -1) : -1); }, "PromptDirect", "Prompt, direct", Black);
sh.AddNewPostfix("PromptDirect",   [this] { return (size_t)(v.Photon.PreSelect.n==1 ? (v.Photon().matchGenPromptDirect||v.Photon().matchGenPromptFrag ? 0 : -1) : -1); }, "PromptDirect", "Prompt, direct", Black);
sh.AddNewPostfix("1Pho_EB_EE",     [this] { return (size_t)(v.Photon.Select.n==1 ? v.Photon.Select(0).isScEtaEE : -1); }, "Barrel;Endcap", "Barrel;Endcap", "601,418");
sh.AddNewPostfix("1PrePho_EB_EE",  [this] { return (size_t)(v.Photon.PreSelect.n==1 ? v.Photon.PreSelect(0).isScEtaEE : -1); }, "Barrel;Endcap", "Barrel;Endcap", "601,418");
sh.AddNewPostfix("FailedJet",      [this] { return (size_t)(v.Jet.FailID.n>0); }, "NoFailedJet;FailedJet",  "0 bad jet;#geq1 bad jet",  Black+Red);
sh.AddNewPostfix("BadMuonJet",     [this] { return (size_t)(v.dPhiMuonJetMET>=2.74159); }, "NoBadMuonJet;BadMuonJet",  "0 bad #mu jet;#geq1 bad #mu jet",  Black+Red);
sh.AddNewPostfix("BadPFMET2",      [this] { return (size_t)(v.CaloMET_pt>0 ? v.MET_pt/v.CaloMET_pt>=2 : -1); }, "GoodPFMET2;BadPFMET2",  "PFMET/CaloMET<2;PFMET/CaloMET#geq2",  Black+Red);
sh.AddNewPostfix("BadPFMET3",      [this] { return (size_t)(v.CaloMET_pt>0 ? v.MET_pt/v.CaloMET_pt>=3 : -1); }, "GoodPFMET3;BadPFMET3",  "PFMET/CaloMET<3;PFMET/CaloMET#geq3",  Black+Red);
sh.AddNewPostfix("BadPFMET5",      [this] { return (size_t)(v.CaloMET_pt>0 ? v.MET_pt/v.CaloMET_pt>=5 : -1); }, "GoodPFMET5;BadPFMET5",  "PFMET/CaloMET<5;PFMET/CaloMET#geq5",  Black+Red);
sh.AddNewPostfix("NJet35",         [this] { return (size_t)(v.Jet.Jet.n<3 ? -1 : v.Jet.Jet.n>5); }, "NJet3to5;NJet6", "3#leqN_{jet}#leq5;6#leqN_{jet}", "1,2");
sh.AddNewPostfix("NJet46",         [this] { return (size_t)(v.Jet.Jet.n<4 ? -1 : v.Jet.Jet.n>6); }, "NJet4to6;NJet7", "4#leqN_{jet}#leq6;7#leqN_{jet}", "1,2");
sh.AddNewPostfix("NJet45",         [this] { return (size_t)(v.Jet.Jet.n==4||v.Jet.Jet.n==5 ? 0 : 1 ); }, "NJet4to5;NJet6", "4#leqN_{jet}#leq5;6#leqN_{jet}", "1,2");
sh.AddNewPostfix("NJetNoPho35",    [this] { return (size_t)(v.Jet.JetNoPho.n<3 ? -1 : v.Jet.JetNoPho.n>5); }, "NJet3to5;NJet6", "3#leqN_{jet}#leq5;6#leqN_{jet}", "1,2");
sh.AddNewPostfix("NJetNoPho46",    [this] { return (size_t)(v.Jet.JetNoPho.n<4 ? -1 : v.Jet.JetNoPho.n>6); }, "NJet4to6;NJet7", "4#leqN_{jet}#leq6;7#leqN_{jet}", "1,2");
sh.AddNewPostfix("NJetNoPho45",    [this] { return (size_t)(v.Jet.JetNoPho.n==4||v.Jet.JetNoPho.n==5 ? 0 : 1 ); }, "NJet4to5;NJet6", "4#leqN_{jet}#leq5;6#leqN_{jet}", "1,2");
sh.AddNewPostfix("MTBoost800",     [this] {  return (size_t)(v.MT_boost>=800 ? 0 : -1); }, "MTBoost800", "m_{T,Boost+MET}#geq800", Black);

sh.AddNewPostfix("Nb",             [this] {  return (size_t)(v.Jet.LooseBTag.n==0 ? 0 : v.Jet.MediumBTag.n>=1 ? 1: -1); },      "bveto;bincl",       "0b;#geq1b",               col4_cyan_to_red);
sh.AddNewPostfix("Nisob",          [this] {  return (size_t)(v.Jet.LooseIsoBTag.n==0 ? 0 : v.Jet.MediumIsoBTag.n>=1 ? 1: -1); },"isobveto;isobincl", "0b (iso.);#geq1b (iso.)", col4_cyan_to_red);
sh.AddNewPostfix("2b",             [this] {  return (size_t)(v.Jet.MediumBTag.n == 2 ? 0 : -1);          }, "2b",               "2b",                      col4_cyan_to_red);


// -------------------------- Other ------------------------

// Systematics Postfixes
sh.AddNewPostfix("Syst", [&syst_index] { return syst_index; }, std::string(";syst[1to")+std::to_string(syst_nSyst)+"]", std::string(";systematics [1to")+std::to_string(syst_nSyst)+"]", "1-999");
if (syst_nSyst>998) error("Error: Too large number of systematics, define more colors!");
// Weights
//sh.AddNewPostfix("NoPUWeight",     [] { return 0; }, "NoPUWeight",   "No pile-up reweighting", Black);
//sh.AddNewPostfix("NoTrigWeight",   [] { return 0; }, "NoTrigWeight", "No trigger weighting",   Black);
//sh.AddNewPostfix("NoEleSF",        [] { return 0; }, "NoEleSF",      "No ele SF",              Black);
//sh.AddNewPostfix("NoMuonSF",       [] { return 0; }, "NoMuonSF",     "No muon SF",             Black);
//sh.AddNewPostfix("NoBTagSF",       [] { return 0; }, "NoBTagSF",     "No b-tag SF",            Black);
//sh.AddNewPostfix("NoWTagSF",       [] { return 0; }, "NoWTagSF",     "No W-tag SF",            Black);
//sh.AddNewPostfix("NoTopTagSF",     [] { return 0; }, "NoTopTagSF",   "No top-tag SF",          Black);
//sh.AddNewPostfix("NoWeight",       [] { return 0; }, "NoWeight",     "",                       Black);
if (debug) std::cout<<"PlottingBase::define_histo_settings: all postfixes ok"<<std::endl;


// --------------------------------------------------------------------
//                          Binning
// --------------------------------------------------------------------

// Bins
E     = {0, 100, 200, 400, 600, 800, 1000, 1500, 2000, 10000};
//Pt    = {0, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 800, 900,       1200, 2000, 10000};
Pt    = {200,      300,      400,      500,      600,           800,                  2000, 10000};
M3    = {0, 50, 75, 100, 150, 200, 250, 800};
//PtG   = {0, 100, 150, 200, 250, 300, 350, 400,      500,                700,           1000,       2000, 10000};
PtG   = {200,      300,      400,                600,           800,                  2000, 10000};
//PtF   = {0,           200,      300,      400, 450, 500, 550, 600,      700,           1000,       2000, 10000};
PtF   = {200,      300,      400,                600,                                 2000, 10000};
PtO   = {200,                                                                                              10000};
PtT   = {0,                                         500,                                                 10000};
PtPho = {0, 100, 120, 140, 160, 180, 200, 225, 250, 300, 500, 1000, 4000};
PtPho2 = {0, 100, 120, 140, 160, 180, 200, 225, 250, 300, 4000};
for (double x=   0; x< 100; x+= 10) PtFine.push_back(x);
for (double x= 100; x< 500; x+= 20) PtFine.push_back(x);  
for (double x= 500; x<1000; x+= 50) PtFine.push_back(x);  
for (double x=1000; x<2000; x+=100) PtFine.push_back(x);  
for (double x=2000; x<5000; x+=500) PtFine.push_back(x);
M   = {0, 10, 20, 30, 40, 50, 65, 75, 85, 95, 105, 120, 135, 150, 165, 180, 195, 210, 230, 260, 300, 500, 1000};
for (double x=   0; x< 120; x+=  5) MFine.push_back(x);
for (double x= 120; x< 220; x+= 10) MFine.push_back(x);
for (double x= 220; x< 300; x+= 20) MFine.push_back(x);  
for (double x= 300; x< 500; x+= 50) MFine.push_back(x);  
for (double x= 500; x<1000; x+=100) MFine.push_back(x);  
for (double x=1000; x<5000; x+=500) MFine.push_back(x);  
MW  = {65, 75, 85, 95, 105};
DeepB = {0, 0.05, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0 };
for (double x=0.0; x< 1.8; x+=0.1) MDP.push_back(x);
for (double x=1.8; x< 2.4; x+=0.2) MDP.push_back(x);
for (double x=2.4; x<=3.2; x+=0.4) MDP.push_back(x);
for (double x=0.0; x<1.6; x+=0.4) DP.push_back(x);
for (double x=1.6; x<2.4; x+=0.2) DP.push_back(x);
for (double x=2.4; x<3.2; x+=0.1) DP.push_back(x);
for (double x=  0; x< 250; x+=10) M2.push_back(x);
for (double x=250; x< 500; x+=25) M2.push_back(x);
for (double x=500; x<1000; x+=50) M2.push_back(x);
NVTX.push_back(0);
for (double x=6;  x<  40;  x+=2) NVTX.push_back(x);
for (double x=40; x<=100;  x+=5) NVTX.push_back(x);
//MET = {0,     100,                     200,      300, 400, 500, 600, 700, 800, 1000, 1500, 2000};
//MET = {0, 80, 100, 120, 140, 160, 180, 200, 250, 300, 400, 500, 600, 700, 800, 1000, 1500, 2000};
MET = {0, 10, 20, 30, 40, 50, 75, 100, 125, 150, 200, 300};
HT = {0, 200, 300, 400, 500, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1100, 1200, 1500, 2000, 2500, 3000, 4000, 10000};
HTB = {400, 500, 600, 700, 750, 800, 850, 900, 950, 1000, 1500, 10000}; // 2D Trigger Eff
PtB = {200, 300, 400, 450, 500, 550, 600, 1000, 10000}; // 2D Trigger Eff
LepPt  = { 0, 5, 10, 15, 20, 30, 40, 50, 60, 80, 100, 120, 150, 200, 250, 300, 400, 500, 4000};
//   LepEta = { 0.0, 0.5, 1.0, 1.5, 2.0, 2.5 };
LepEta = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5 };
// Razor inclusive binning
HT_2D_bins_old  = {0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.2, 1.5, 10.0};
Pt_2D_bins      = {200, 300, 400, 450, 500, 550, 600, 700, 1000, 5000};
PtLow_2D_bins   = {200, 300, 400, 450, 500};
PtHigh_2D_bins  = {550, 600, 700, 1000, 5000};
PtPho_2D_bins   = {0, 120, 200, 300, 400, 500, 800, 5000};
ElePt_2D_bins   = {10, 15, 20, 30, 35, 40, 50, 100, 125, 250, 2000};
MuPt_2D_bins    = {10, 15, 20, 25, 30, 35, 40,  50,  60, 100, 250, 2000};
HT_2D_bins_lep  = {200, 300, 400, 500, 600, 10000};
HT1_2D_bins_lep = {200, 300, 400, 500};
HT2_2D_bins_lep = {500, 600, 10000};
// HT/MET Trigger turnons
// HT:   0  500   700   800  1050
// MHT:  0   75    85   100   120
// MET:  0  100   115   136   158
// rat:    1.33  1.35  1.36  1.32
//HT_2D_bins  = {0, 400, 500, 600, 700, 800, 1000, 1200, 1500, 5000}; // 2D Trigger Eff Run2017-18
HT_2D_bins  = {200,  450,  600,  700, 800, 900, 1000, 1200, 10000}; // 2D Trigger Eff Run2017-18
HT1_2D_bins = {200,  450,  600}; // 2D Trigger Eff Run2017-18
HT2_2D_bins = {600,  700,  800,  900}; // 2D Trigger Eff Run2017-18
HT3_2D_bins = {900, 1000, 1200, 10000}; // 2D Trigger Eff Run2017-18
MET_2D_bins = {0, 100, 130, 160, 180, 200, 250, 300, 400, 4000}; // 2D Trigger Eff Run2017-18
// HT,MET Changgi {300, 500, 550, 600, 800, 1500, 10000}, {50, 150, 200, 400, 10000}, {}, 0, 0);
// unrolled bin mergin
merged_razor_bins = {19,23,24};
merged_razor_bins_lep = {24,28,29};
merged_trigger_bins  = {1,3,5,  11,   20, 37,44,  46,52,53, 55,57,59,61,62, 64,65,66,67,68,69,70,71};
merged_trigger1_bins = {1,3,5,  11                                                                 };
merged_trigger2_bins = {               2, 19,26,  28,34,35                                         };
merged_trigger3_bins = {                                    10,12,14,16,17, 19,20,21,22,23,24,25,26};

// ------------------- Binning with Postfix ----------------

std::stringstream HT_pf, HT_leg;
for (size_t i=0, n=HTB.size(); i<n-1; ++i) {
	HT_pf<<"HT"<<HTB[i]<<"to"<<HTB[i+1];
	HT_leg<<"H_{T} #subset ["<<HTB[i]<<","<<HTB[i+1]<<"[";
	if (i!=n-2) { HT_pf<<";"; HT_leg<<";"; }
}
sh.AddNewPostfix("HTBins", [this] { for (size_t i=0, n=HTB.size(); i<n-1; ++i) if (v.AK4_Ht>=HTB[i]&&v.AK4_Ht<HTB[i+1]) return i; return (size_t)-1; },
		HT_pf.str(), HT_leg.str(), col12+col12);
std::stringstream AK8Pt_pf, AK8Pt_leg;
for (size_t i=0, n=PtB.size(); i<n-1; ++i) {
	AK8Pt_pf<<"Jet1AK8Pt"<<PtB[i]<<"to"<<PtB[i+1];
	AK8Pt_leg<<"AK8 jet1 p_{T} #subset ["<<PtB[i]<<","<<PtB[i+1]<<"[";
	if (i!=n-2) { AK8Pt_pf<<";"; AK8Pt_leg<<";"; }
}

// --------------------------------------------------------------------
//                         Fill Parameters
// --------------------------------------------------------------------

// Define the filling parameters: bins, variable to fill, default plotting range etc.
// Usage:
// AddNewFillParams("Name of fill param", { .nbin= int(number of bins), .bins={vector of bin edges} OR {min,max}, .fill=lambda function which returns the value you want to fill, .axis_title="title and (unit)", .def_range={min,max range to draw}});

sh.AddNewFillParams("Bin", { .nbin=1,   .bins={0,1}, .fill=[] { return 0; }, .axis_title="Bin", .def_range={0,1}}); // Just one bin for averages/counts

// ------------------------- Leptons -----------------------

// Veto Leptons
sh.AddNewFillParams("VetoElePt",       { .nbin=LepPt.size()-1, .bins=LepPt,   .fill=[this] { return v.Electron().pt;                  }, .axis_title="Loose Electron p_{T} (GeV)", .def_range={0,500}});
sh.AddNewFillParams("VetoEleEta",      { .nbin=LepEta.size()-1, .bins=LepEta, .fill=[this] { return std::abs(v.Electron().eta);                 }, .axis_title="Loose Electron |#eta|",  .def_range={-2.5,2.5}});
sh.AddNewFillParams("VetoMuPt",        { .nbin=LepPt.size()-1, .bins=LepPt,   .fill=[this] { return v.Muon().pt;                    }, .axis_title="Loose Muon p_{T} (GeV)",     .def_range={0,500}});
sh.AddNewFillParams("VetoMuEta",       { .nbin=LepEta.size()-1, .bins=LepEta, .fill=[this] { return std::abs(v.Muon().eta);                   }, .axis_title="Loose Muon |#eta|",      .def_range={-2.4,2.4}});
sh.AddNewFillParams("VetoTauPt",       { .nbin=LepPt.size()-1, .bins=LepPt,   .fill=[this] { return v.Tau().pt;                  }, .axis_title="Tau p_{T} (GeV)",     .def_range={0,500}});
sh.AddNewFillParams("VetoTauEta",      { .nbin=LepEta.size()-1, .bins=LepEta, .fill=[this] { return std::abs(v.Tau().eta);                 }, .axis_title="Tau |#eta|",      .def_range={-2.6,2.6}});

// Selected Leptons
// Electrons
sh.AddNewFillParams("ElePt",           { .nbin=LepPt.size()-1, .bins=LepPt,  .fill=[this] { return v.Electron.Select.n<1 ? -9999 : v.Electron.Select(0).pt;  }, .axis_title="Electron p_{T} (GeV)", .def_range={ELE_SELECT_PT_CUT,500}});
sh.AddNewFillParams("EleEta",          { .nbin=LepEta.size()-1, .bins=LepEta,.fill=[this] { return v.Electron.Select.n<1 ? -9999 : std::abs(v.Electron.Select(0).eta); }, .axis_title="Electron |#eta|",  .def_range={-ELE_SELECT_ETA_CUT,ELE_SELECT_ETA_CUT}});
sh.AddNewFillParams("EleJetDR",        { .nbin=  60, .bins={     0,      6}, .fill=[this] { return v.Electron.Select.n<1 ? -9999 : v.Electron.Select(0).jetDR;         }, .axis_title="#DeltaR (ele, jet)",         .def_range={0,4}});
sh.AddNewFillParams("EleJetPt",        { .nbin= 100, .bins={     0,    500}, .fill=[this] { return v.Electron.Select.n<1 ? -9999 : v.Electron.Select(0).iMatchedAK4==(size_t)-1 ? -9999 : v.Jet(v.Electron.Select(0).iMatchedAK4).pt;         }, .axis_title="p_{T, nearest jet to ele}"});
sh.AddNewFillParams("EleJetDPhi",      { .nbin=MDP.size()-1, .bins=MDP,      .fill=[this] { return v.Electron.Select.n<1 ? -9999 : v.Electron.Select(0).jetDPhi;       }, .axis_title="#Delta#phi (ele, jet)"});
sh.AddNewFillParams("Ele1JetDPhi",     { .nbin=MDP.size()-1, .bins=MDP,      .fill=[this] { return v.Electron.Select.n<1 ? -9999 : v.Electron.Select(0).jetDPhi; }, .axis_title="#Delta#phi (1st ele, jet)"});
sh.AddNewFillParams("Ele2JetDPhi",     { .nbin=MDP.size()-1, .bins=MDP,      .fill=[this] { return v.Electron.Select.n<2 ? -9999 : v.Electron.Select(1).jetDPhi; }, .axis_title="#Delta#phi (2nd ele, jet)"});
sh.AddNewFillParams("EleCleanJetDRmin",{ .nbin=  60, .bins={     0,      6}, .fill=[this] { return v.Electron().jetDRmin; }, .axis_title="e, jet #DeltaR_{min}"});
sh.AddNewFillParams("EleCleanJetPtrel",{ .nbin=  40, .bins={     0,    100}, .fill=[this] { return v.Electron().cleanJetPtrel; }, .axis_title="e, cleaned jet p_{T,rel} (GeV)"});
sh.AddNewFillParams("EleTightPt",       { .nbin= 100, .bins={     0,    500}, .fill=[this] { return v.Electron.Tight.n<1 ? -9999 : v.Electron.Tight(0).pt;  }, .axis_title="Electron p_{T} (GeV)", .def_range={ELE_TIGHT_PT_CUT,250}});
sh.AddNewFillParams("EleTightEta",      { .nbin=  40, .bins={    -4,      4}, .fill=[this] { return v.Electron.Tight.n<1 ? -9999 : v.Electron.Tight(0).eta; }, .axis_title="Electron #eta",  .def_range={-ELE_TIGHT_ETA_CUT,ELE_TIGHT_ETA_CUT}});
sh.AddNewFillParams("EleNoIsoPt",       { .nbin= 100, .bins={     0,    500}, .fill=[this] { return v.Electron.NoIso.n<1 ? -9999 : v.Electron.NoIso(0).pt;  }, .axis_title="Electron p_{T} (GeV)", .def_range={ELE_TIGHT_PT_CUT,250}});
sh.AddNewFillParams("EleNoIsoEta",      { .nbin=  40, .bins={    -4,      4}, .fill=[this] { return v.Electron.NoIso.n<1 ? -9999 : v.Electron.NoIso(0).eta; }, .axis_title="Electron #eta",  .def_range={-ELE_TIGHT_ETA_CUT,ELE_TIGHT_ETA_CUT}});
sh.AddNewFillParams("EleNonIsoPt",       { .nbin= 100, .bins={     0,    500}, .fill=[this] { return v.Electron.NonIso.n<1 ? -9999 : v.Electron.NonIso(0).pt;  }, .axis_title="Electron p_{T} (GeV)", .def_range={ELE_TIGHT_PT_CUT,250}});
sh.AddNewFillParams("EleNonIsoEta",      { .nbin=  40, .bins={    -4,      4}, .fill=[this] { return v.Electron.NonIso.n<1 ? -9999 : v.Electron.NonIso(0).eta; }, .axis_title="Electron #eta",  .def_range={-ELE_TIGHT_ETA_CUT,ELE_TIGHT_ETA_CUT}});
// Muons
sh.AddNewFillParams("MuPt",            { .nbin=LepPt.size()-1, .bins=LepPt,   .fill=[this] { return v.Muon.Select.n<1 ? -9999 : v.Muon.Select(0).pt;  }, .axis_title="Muon p_{T} (GeV)", .def_range={MU_SELECT_PT_CUT,500}});
sh.AddNewFillParams("MuEta",           { .nbin=LepEta.size()-1, .bins=LepEta, .fill=[this] { return v.Muon.Select.n<1 ? -9999 : std::abs(v.Muon.Select(0).eta); }, .axis_title="Muon |#eta|",        .def_range={-MU_SELECT_ETA_CUT,MU_SELECT_ETA_CUT}});
sh.AddNewFillParams("MuJetDR",         { .nbin=  60, .bins={     0,      6}, .fill=[this] { return v.Muon.Select.n<1 ? -9999 : v.Muon.Select(0).jetDR;      }, .axis_title="#DeltaR (muon, jet)",        .def_range={0,4}});
sh.AddNewFillParams("MuJetPt",         { .nbin= 100, .bins={     0,    500}, .fill=[this] { return v.Muon.Select.n<1 ? -9999 : v.Muon.Select(0).iMatchedAK4==(size_t)-1 ? -9999 : v.Jet(v.Muon.Select(0).iMatchedAK4).pt;      }, .axis_title="p_{T, nearest jet to muon}"});
sh.AddNewFillParams("MuJetDPhi",       { .nbin=MDP.size()-1, .bins=MDP,      .fill=[this] { return v.Muon.Select.n<1 ? -9999 : v.Muon.Select(0).jetDPhi;    }, .axis_title="#Delta#phi (muon, jet)"});
sh.AddNewFillParams("Mu1JetDPhi",      { .nbin=MDP.size()-1, .bins=MDP,      .fill=[this] { return v.Muon.Select.n<1 ? -9999 : v.Muon.Select(0).jetDPhi; }, .axis_title="#Delta#phi (1st muon, jet)"});
sh.AddNewFillParams("Mu2JetDPhi",      { .nbin=MDP.size()-1, .bins=MDP,      .fill=[this] { return v.Muon.Select.n<2 ? -9999 : v.Muon.Select(1).jetDPhi; }, .axis_title="#Delta#phi (2nd muon, jet)"});
sh.AddNewFillParams("MuCleanJetDRmin",{ .nbin=  60, .bins={     0,      6}, .fill=[this] { return v.Muon().jetDRmin; }, .axis_title="#mu, jet #DeltaR_{min}"});
sh.AddNewFillParams("MuCleanJetPtrel",{ .nbin=  40, .bins={     0,    100}, .fill=[this] { return v.Muon().cleanJetPtrel; }, .axis_title="#mu, cleaned jet p_{T,rel} (GeV)"});
sh.AddNewFillParams("MuTightPt",            { .nbin= 100, .bins={     0,    500}, .fill=[this] { return v.Muon.Tight.n<1 ? -9999 : v.Muon.Tight(0).pt;  }, .axis_title="Muon p_{T} (GeV)", .def_range={MU_TIGHT_PT_CUT,250}});
sh.AddNewFillParams("MuTightEta",           { .nbin=  40, .bins={    -4,      4}, .fill=[this] { return v.Muon.Tight.n<1 ? -9999 : std::abs(v.Muon.Tight(0).eta); }, .axis_title="Muon |#eta|",        .def_range={-MU_TIGHT_ETA_CUT,MU_TIGHT_ETA_CUT}});
sh.AddNewFillParams("MuNoIsoPt",            { .nbin= 100, .bins={     0,    500}, .fill=[this] { return v.Muon.NoIso.n<1 ? -9999 : v.Muon.NoIso(0).pt;  }, .axis_title="Muon p_{T} (GeV)", .def_range={MU_TIGHT_PT_CUT,250}});
sh.AddNewFillParams("MuNoIsoEta",           { .nbin=  40, .bins={    -4,      4}, .fill=[this] { return v.Muon.NoIso.n<1 ? -9999 : std::abs(v.Muon.NoIso(0).eta); }, .axis_title="Muon |#eta|",        .def_range={-MU_TIGHT_ETA_CUT,MU_TIGHT_ETA_CUT}});
sh.AddNewFillParams("MuNonIsoPt",            { .nbin= 100, .bins={     0,    500}, .fill=[this] { return v.Muon.NonIso.n<1 ? -9999 : v.Muon.NonIso(0).pt;  }, .axis_title="Muon p_{T} (GeV)", .def_range={MU_TIGHT_PT_CUT,250}});
sh.AddNewFillParams("MuNonIsoEta",           { .nbin=  40, .bins={    -4,      4}, .fill=[this] { return v.Muon.NonIso.n<1 ? -9999 : std::abs(v.Muon.NonIso(0).eta); }, .axis_title="Muon |#eta|",        .def_range={-MU_TIGHT_ETA_CUT,MU_TIGHT_ETA_CUT}});

// ------------------------- Photons -----------------------

// Photons
sh.AddNewFillParams("PhotonPt",        { .nbin=  25, .bins={     0,   2000}, .fill=[this] { return v.Photon.Select.n==1 ? v.Photon.Select(0).pt  : -9999;  }, .axis_title="Photon p_{T} (GeV)", .def_range={PHOTON_SELECT_PT_CUT,1000}});
sh.AddNewFillParams("PrePhotonPt",     { .nbin=  25, .bins={     0,   2000}, .fill=[this] { return v.Photon.PreSelect.n==1 ? v.Photon.PreSelect(0).pt  : -9999;  }, .axis_title="Preselected #gamma p_{T} (GeV)", .def_range={PHOTON_SELECT_PT_CUT,1000}});
sh.AddNewFillParams("PhotonPtBins",    { .nbin=PtPho.size()-1, .bins=PtPho,  .fill=[this] { return v.Photon.Select.n==1 ? v.Photon.Select(0).pt  : -9999;  }, .axis_title="Photon p_{T} (GeV)", .def_range={0,1000}});
sh.AddNewFillParams("PhotonPtFewBins", { .nbin=PtPho2.size()-1, .bins=PtPho2,  .fill=[this] { return v.Photon.Select.n==1 ? v.Photon.Select(0).pt  : -9999;  }, .axis_title="Photon p_{T} (GeV)", .def_range={0,300}});
sh.AddNewFillParams("PrePhotonPtBins", { .nbin=PtPho.size()-1, .bins=PtPho,  .fill=[this] { return v.Photon.PreSelect.n==1 ? v.Photon.PreSelect(0).pt  : -9999;  }, .axis_title="Preselected #gamma p_{T} (GeV)", .def_range={0,1000}});
sh.AddNewFillParams("PhotonEta",       { .nbin=  40, .bins={    -4,      4}, .fill=[this] { return v.Photon.Select.n==1 ? v.Photon.Select(0).eta : -9999; }, .axis_title="Photon #eta",        .def_range={-PHOTON_SELECT_ETA_CUT,PHOTON_SELECT_PT_CUT}});
sh.AddNewFillParams("PhotonSIEIE",     { .nbin= 100, .bins={     0,    0.1}, .fill=[this] { return v.Photon().sieie;        }, .axis_title="Photon #sigma_{i#eta i#eta}", .def_range={0,0.1}});
sh.AddNewFillParams("PhotonCHIso",     { .nbin=  40, .bins={     0,     20}, .fill=[this] { return v.Photon().pfRelIso03_chg;     }, .axis_title="Photon Charged Isolation (GeV)", .def_range={0,10}});
sh.AddNewFillParams("Photon1Eta",      { .nbin=  25, .bins={  -2.5,    2.5}, .fill=[this] { return v.Photon.Select.n<1 ? -9999 : v.Photon.Select(0).eta; }, .axis_title="Photon #eta"});
sh.AddNewFillParams("PhotonCHIso_preslectphopt",       { .nbin=  40, .bins={     0,     20}, .fill=[this] { return v.Photon().pt*v.Photon().pfRelIso03_chg;     }, .axis_title="Photon Charged Isolation (GeV)", .def_range={0,20}});
sh.AddNewFillParams("PhotonCHIso_preslectphoptREBin",  { .nbin=  10, .bins={0, 0.00001, 1, 1.5, 2, 3, 4, 5, 7, 10, 1000}, .fill=[this] { return v.Photon().pt*v.Photon().pfRelIso03_chg;     }, .axis_title="PhoCHIso_pt (GeV)"});
sh.AddNewFillParams("PhotonCHIso_preslectphoptLog",    { .nbin=  16, .bins={-4, 4}, .fill=[this] { return std::log10(v.Photon().pt*v.Photon().pfRelIso03_chg) ;     }, .axis_title="PhoCHIso_pt (GeV)", .def_range={-4,2}});
sh.AddNewFillParams("PhotonCHIso_preslectphoptLN",     { .nbin=  16, .bins={-8, 8}, .fill=[this] { return std::log(v.Photon().pt*v.Photon().pfRelIso03_chg) ;     }, .axis_title="PhoCHIso_pt (GeV)"});

// --------------------------- Jets ------------------------

// AK4 Jets
sh.AddNewFillParams("JetPtBins",            { .nbin=PtFine.size()-1,.bins=PtFine, .fill=[this] { return v.Jet().pt;           }, .axis_title="Jet p_{T} (GeV)", .def_range={  0,2000} });
sh.AddNewFillParams("JetPtFewBins",         { .nbin=PtF.size()-1,  .bins=PtF,     .fill=[this] { return v.Jet().pt;           }, .axis_title="Jet p_{T} (GeV)", .def_range={200,2000} });
sh.AddNewFillParams("JetPtOneBin",          { .nbin=PtO.size()-1,  .bins=PtO,     .fill=[this] { return v.Jet().pt;           }, .axis_title="Jet p_{T} (GeV)", .def_range={400,2000} });
sh.AddNewFillParams("JetPt",                { .nbin=  80, .bins={     0,   4000}, .fill=[this] { return v.Jet().pt;           }, .axis_title="Jet p_{T} (GeV)", .def_range={0,2000} });
sh.AddNewFillParams("JetEta",               { .nbin=  80, .bins={    -4,      4}, .fill=[this] { return v.Jet().eta;          }, .axis_title="Jet #eta",        .def_range={-2.4,2.4}});
sh.AddNewFillParams("JetPhi",               { .nbin=  16, .bins={-3.142,  3.142}, .fill=[this] { return v.Jet().phi;          }, .axis_title="Jet #phi"});
sh.AddNewFillParams("JetDeepB",             { .nbin=  20, .bins={     0,   1.00}, .fill=[this] { return std::min(v.Jet().btagDeepB,float(0.999)); }, .axis_title="Jet DeepB"});
sh.AddNewFillParams("JetPhotonDR",          { .nbin= 120, .bins={     0,      6}, .fill=[this] { return v.Jet().phoDR;           }, .axis_title="#DeltaR (jet, photon)", .def_range={0,0.8}});
sh.AddNewFillParams("JetPhotonPtRatio",     { .nbin=  40, .bins={     0,      4}, .fill=[this] { return v.Jet().phoPtRatio;      }, .axis_title="p_{T}^{photon}/p_{T}^{jet}", .def_range={0,2}});
sh.AddNewFillParams("JetEleDR",             { .nbin= 120, .bins={     0,      6}, .fill=[this] { return v.Jet().eleDR;           }, .axis_title="#DeltaR (jet, electron)", .def_range={0,0.8}});
sh.AddNewFillParams("JetElePtRatio",        { .nbin=  40, .bins={     0,      4}, .fill=[this] { return v.Jet().elePtRatio;      }, .axis_title="p_{T}^{electron}/p_{T}^{jet}", .def_range={0,2}});
sh.AddNewFillParams("JetMuonDR",            { .nbin= 120, .bins={     0,      6}, .fill=[this] { return v.Jet().muDR;            }, .axis_title="#DeltaR (jet, muon)", .def_range={0,0.8}});
sh.AddNewFillParams("JetMuonPtRatio",       { .nbin=  40, .bins={     0,      4}, .fill=[this] { return v.Jet().muPtRatio;       }, .axis_title="p_{T}^{muon}/p_{T}^{jet}", .def_range={0,2}});
// BJets
sh.AddNewFillParams("BJetPtBins",           { .nbin=PtFine.size()-1,.bins=PtFine, .fill=[this] { return v.Jet().pt;           }, .axis_title="B-jet p_{T} (GeV)", .def_range={0,2000} });
sh.AddNewFillParams("BJetPt",               { .nbin=  80, .bins={     0,   4000}, .fill=[this] { return v.Jet().pt;           }, .axis_title="B-jet p_{T} (GeV)", .def_range={0,2000} });
sh.AddNewFillParams("BJetEta",              { .nbin=  80, .bins={    -4,      4}, .fill=[this] { return v.Jet().eta;          }, .axis_title="B-jet #eta",        .def_range={-2.4,2.4}});
sh.AddNewFillParams("BJetPhi",              { .nbin=  16, .bins={-3.142,  3.142}, .fill=[this] { return v.Jet().phi;          }, .axis_title="B-jet #phi"});
sh.AddNewFillParams("BJetDeepB",              { .nbin=  20, .bins={     0,   1.00}, .fill=[this] { return std::min(v.Jet().btagDeepB,float(0.999)); }, .axis_title="B-jet DeepB"});
// Megajets
sh.AddNewFillParams("MegaJetPt",            { .nbin=  80, .bins={     0,   4000}, .fill=[this] { return v.iMegaJet<v.megajets.size() ? v.megajets[v.iMegaJet].Pt()  : -9999;      }, .axis_title="Megajet p_{T} (GeV)", .def_range={0,2000} });
sh.AddNewFillParams("MegaJetEta",           { .nbin=  80, .bins={    -4,      4}, .fill=[this] { return v.iMegaJet<v.megajets.size() ? v.megajets[v.iMegaJet].Eta() : -9999;     }, .axis_title="Megajet #eta",        .def_range={-2.4,2.4}});
sh.AddNewFillParams("MegaJetPhi",           { .nbin=  16, .bins={-3.142,  3.142}, .fill=[this] { return v.iMegaJet<v.megajets.size() ? v.megajets[v.iMegaJet].Phi() : -9999;     }, .axis_title="Megajet #phi"});

// gen
sh.AddNewFillParams("GenLepPt",             { .nbin= 100, .bins={     0,    500},  .fill=[this] { return v.GenPart().pt;            }, .axis_title="Gen. Lepton p_{T} (GeV)", .def_range={5,200}});
sh.AddNewFillParams("GenLepPtBins",         { .nbin=LepPt.size()-1,  .bins=LepPt,  .fill=[this] { return v.GenPart().pt;            }, .axis_title="Gen. Lepton p_{T} (GeV)", .def_range={5,200}});

// Lepton reco+id+iso efficiencies
sh.AddNewFillParams("GenElePt",             { .nbin=LepPt.size()-1,  .bins=LepPt,  .fill=[this] { return abs(v.GenPart().pdgId) == 11 ? v.GenPart().pt : 9999; }, .axis_title="Gen. Ele p_{T} (GeV)", .def_range={5,200}});
sh.AddNewFillParams("GenEleEta",            { .nbin=LepEta.size()-1, .bins=LepEta, .fill=[this] { return abs(v.GenPart().pdgId) == 11 ? fabs(v.GenPart().eta) : 9999; }, .axis_title="Gen. Ele |\eta|", .def_range={0,2.5}});
sh.AddNewFillParams("GenMuPt",              { .nbin=LepPt.size()-1,  .bins=LepPt,  .fill=[this] { return abs(v.GenPart().pdgId) == 13 ? v.GenPart().pt : 9999; }, .axis_title="Gen. Mu p_{T} (GeV)", .def_range={5,200}});
sh.AddNewFillParams("GenMuEta",             { .nbin=LepEta.size()-1, .bins=LepEta, .fill=[this] { return abs(v.GenPart().pdgId) == 13 ? fabs(v.GenPart().eta) : 9999; }, .axis_title="Gen. Mu |\eta|", .def_range={0,2.5}});

sh.AddNewFillParams("GenLepWMatchedGenLepPt",       { .nbin= 100, .bins={     0,    500},  .fill=[this] { return v.GenPart().iGenLepDaughter  == (size_t)-1 ? -9999 : v.GenPart(v.GenPart().iGenLepDaughter).pt; }, .axis_title="Gen. Lepton (from W) p_{T} (GeV)", .def_range={5,200}});
sh.AddNewFillParams("GenLepWMatchedGenLepPtBins",   { .nbin=LepPt.size()-1,  .bins=LepPt,  .fill=[this] { return v.GenPart().iGenLepDaughter  == (size_t)-1 ? -9999 : v.GenPart(v.GenPart().iGenLepDaughter).pt; }, .axis_title="Gen. Lepton (from W) p_{T} (GeV)", .def_range={5,200}});

sh.AddNewFillParams("GenLepTopMatchedGenLepPt",     { .nbin= 100, .bins={     0,    500},  .fill=[this] { return v.GenPart().iGenLepGrandDaughter == (size_t)-1 ? -9999 : v.GenPart(v.GenPart().iGenLepGrandDaughter).pt; }, .axis_title="Gen. Lepton (from top) p_{T} (GeV)", .def_range={5,200}});
sh.AddNewFillParams("GenLepTopMatchedGenLepPtBins", { .nbin=LepPt.size()-1,  .bins=LepPt,  .fill=[this] { return v.GenPart().iGenLepGrandDaughter == (size_t)-1 ? -9999 : v.GenPart(v.GenPart().iGenLepGrandDaughter).pt; }, .axis_title="Gen. Lepton (from top) p_{T} (GeV)", .def_range={5,200}});

sh.AddNewFillParams("GenLepEta",                    { .nbin=LepEta.size()-1, .bins=LepEta, .fill=[this] { return std::abs(v.GenPart().eta); }, .axis_title="Gen. Lepton |#eta|"});
sh.AddNewFillParams("GenLepTopMatchedGenLepEta",    { .nbin=LepEta.size()-1, .bins=LepEta, .fill=[this] { return v.GenPart().iGenLepGrandDaughter == (size_t)-1 ? -9999 : std::abs(v.GenPart(v.GenPart().iGenLepGrandDaughter).eta); }, .axis_title="Gen. Lepton (from top) |#eta|"});
//   sh.AddNewFillParams("NGenLepTopMatchedGenLep",    { .nbin=3, .bins={  -0.5,    2.5}, .fill=[this] { return v.GenPart().iGenLepGrandDaughter == (size_t)-1 ? -9999 : v.GenPart(v.GenPart().iGenLepGrandDaughter).n; }, .axis_title="Number of Gen. Lepton (from top)"});
sh.AddNewFillParams("GenHadWPt",                    { .nbin=   80, .bins={    0,    4000}, .fill=[this] { return v.GenPart().pt;  }, .axis_title="Gen-W (had.) p_{T} (GeV)",   .def_range={0, 2000}});
sh.AddNewFillParams("GenHadWPtBins",                { .nbin=Pt.size()-1, .bins=Pt,         .fill=[this] { return v.GenPart().pt;  }, .axis_title="Gen-W (had.) p_{T} (GeV)",   .def_range={0, 2000}});
sh.AddNewFillParams("GenHadWEta",                   { .nbin=LepEta.size()-1, .bins=LepEta, .fill=[this] { return std::abs(v.GenPart().eta); }, .axis_title="Gen. Had W |#eta|"});

sh.AddNewFillParams("GenTopPt",                     { .nbin=   80, .bins={    0,    4000}, .fill=[this] { return v.GenPart().pt;  }, .axis_title="Gen-top p_{T} (GeV)", .def_range={0, 2000}});
sh.AddNewFillParams("GenTopPtBins",                 { .nbin=Pt.size()-1, .bins=Pt,         .fill=[this] { return v.GenPart().pt;  }, .axis_title="Gen-top p_{T} (GeV)", .def_range={0, 2000}});
sh.AddNewFillParams("GenTopEta",                    { .nbin=LepEta.size()-1, .bins=LepEta, .fill=[this] { return std::abs(v.GenPart().eta); }, .axis_title="Gen. top |#eta|"});

sh.AddNewFillParams("GenHadTopPt",                  { .nbin=   80, .bins={    0,    4000}, .fill=[this] { return v.GenPart().pt;  }, .axis_title="Gen-top (had.) p_{T} (GeV)", .def_range={0, 2000}});
sh.AddNewFillParams("GenHadTopPtBins",              { .nbin=Pt.size()-1, .bins=Pt,         .fill=[this] { return v.GenPart().pt;  }, .axis_title="Gen-top (had.) p_{T} (GeV)", .def_range={0, 2000}});

sh.AddNewFillParams("GenLepTopPt",                  { .nbin=   80, .bins={    0,    4000}, .fill=[this] { return v.GenPart().pt;  }, .axis_title="Gen-top (lep.) p_{T} (GeV)", .def_range={0, 2000}});
sh.AddNewFillParams("GenLepTopPtBins",              { .nbin=Pt.size()-1, .bins=Pt,         .fill=[this] { return v.GenPart().pt;  }, .axis_title="Gen-top (lep.) p_{T} (GeV)", .def_range={0, 2000}});


// --------------------- Event Variables -------------------

// Object counts
sh.AddNewFillParams("NVtx",                 { .nbin=NVTX.size()-1, .bins=NVTX,    .fill=[this] { return v.PV_npvsGood;                  }, .axis_title="N_{Vertices}",         .def_range={0,50}});
sh.AddNewFillParams("NJet",                 { .nbin=  20, .bins={    0,      20}, .fill=[this] { return v.Jet.Jet.n;                    }, .axis_title="N_{Jet}",              .def_range={2,20}});
sh.AddNewFillParams("NJetBins",             { .nbin=   5, .bins={2,3,4,5,6,20},   .fill=[this] { return v.Jet.Jet.n;                    }, .axis_title="N_{Jet}",              .def_range={2,20}});
sh.AddNewFillParams("NJetNoPho",            { .nbin=  20, .bins={    0,      20}, .fill=[this] { return v.Jet.JetNoPho.n;               }, .axis_title="N_{Jet}",              .def_range={2,20}});
sh.AddNewFillParams("NJetNoPhoBins",        { .nbin=  3,  .bins={2,4,6,20},   .fill=[this] { return v.Jet.JetNoPho.n;               }, .axis_title="N_{Jet}",              .def_range={2,20}});  
sh.AddNewFillParams("NBTag",                { .nbin=  10, .bins={    0,      10}, .fill=[this] { return v.Jet.MediumBTag.n;             }, .axis_title="N_{b}",                .def_range={0,10}});
sh.AddNewFillParams("NLooseBTag",           { .nbin=   5, .bins={    0,       5}, .fill=[this] { return v.Jet.LooseBTag.n;              }, .axis_title="N_{b, loose tag}",     .def_range={0,4}});
sh.AddNewFillParams("NTightBTag",           { .nbin=   5, .bins={    0,       5}, .fill=[this] { return v.Jet.TightBTag.n;              }, .axis_title="N_{b, tight tag}",     .def_range={0,4}});
sh.AddNewFillParams("NBTagNoPho",           { .nbin=   5, .bins={    0,       5}, .fill=[this] { return v.Jet.MediumBTagNoPho.n;        }, .axis_title="N_{b, no #gamma}",     .def_range={0,4}});
sh.AddNewFillParams("NLepVeto",             { .nbin=   5, .bins={    0,       5}, .fill=[this] { return v.nLepVeto;                     }, .axis_title="N_{lepton, Veto}",     .def_range={0,4}});
sh.AddNewFillParams("NEleVeto",             { .nbin=   5, .bins={    0,       5}, .fill=[this] { return v.Electron.Veto.n;              }, .axis_title="N_{ele, Veto}",        .def_range={0,4}});
sh.AddNewFillParams("NMuVeto",              { .nbin=   5, .bins={    0,       5}, .fill=[this] { return v.Muon.Veto.n;                  }, .axis_title="N_{muon, Veto}",       .def_range={0,4}});
sh.AddNewFillParams("NTauSelect",             { .nbin=  10, .bins={    0,      10}, .fill=[this] { return v.Tau.Select.n;                   }, .axis_title="N_{tau, Veto}",        .def_range={0,10}});
sh.AddNewFillParams("NIsoTrk",              { .nbin=   5, .bins={    0,       5}, .fill=[this] { return v.nIsoTrack;                    }, .axis_title="N_{iso trk}",          .def_range={0,4}});
sh.AddNewFillParams("NLep",                 { .nbin=   5, .bins={    0,       5}, .fill=[this] { return v.nLepSelect;                   }, .axis_title="N_{lepton}",           .def_range={0,4}});
sh.AddNewFillParams("NLepTight",            { .nbin=   5, .bins={    0,       5}, .fill=[this] { return v.nLepTight;                    }, .axis_title="N_{e/#mu, tight}",         .def_range={0,5}});
sh.AddNewFillParams("NEleNoIso",            { .nbin=   5, .bins={    0,       5}, .fill=[this] { return v.Electron.NoIso.n;             }, .axis_title="N_{e, tight, no iso}", .def_range={0,5}});
sh.AddNewFillParams("NMuNoIso",             { .nbin=   5, .bins={    0,       5}, .fill=[this] { return v.Muon.NoIso.n;                 }, .axis_title="N_{#mu, tight, no iso}", .def_range={0,5}});
sh.AddNewFillParams("NLepNoIso",            { .nbin=   5, .bins={    0,       5}, .fill=[this] { return v.nLepNoIso;                    }, .axis_title="N_{e/#mu, tight, no iso}", .def_range={0,5}});
sh.AddNewFillParams("NEleNonIso",           { .nbin=   5, .bins={    0,       5}, .fill=[this] { return v.Electron.NonIso.n;            }, .axis_title="N_{e, non-iso.}", .def_range={0,5}});
sh.AddNewFillParams("NMuNonIso",            { .nbin=   5, .bins={    0,       5}, .fill=[this] { return v.Muon.NonIso.n;                }, .axis_title="N_{#mu, non-iso.}", .def_range={0,5}});
sh.AddNewFillParams("NLepNonIso",           { .nbin=   5, .bins={    0,       5}, .fill=[this] { return v.nLepNonIso;                   }, .axis_title="N_{e/#mu, non-iso.}", .def_range={0,5}});
sh.AddNewFillParams("NEle",                 { .nbin=   5, .bins={    0,       5}, .fill=[this] { return v.Electron.Select.n;            }, .axis_title="N_{electron}",         .def_range={0,4}});
sh.AddNewFillParams("NMu",                  { .nbin=   5, .bins={    0,       5}, .fill=[this] { return v.Muon.Select.n;                }, .axis_title="N_{muon}",             .def_range={0,4}});
sh.AddNewFillParams("NPhoton",              { .nbin=   5, .bins={    0,       5}, .fill=[this] { return v.Photon.Select.n;              }, .axis_title="N_{photon}",           .def_range={0,4}});
// Gen truth
sh.AddNewFillParams("NGenEleFromW",         { .nbin=   5, .bins={    0,       5}, .fill=[this] { return v.GenPart.EleFromW.n;      }, .axis_title="N_{e, gen}",   .def_range={0,5}});
sh.AddNewFillParams("NGenMuFromW",          { .nbin=   5, .bins={    0,       5}, .fill=[this] { return v.GenPart.MuFromW.n;       }, .axis_title="N_{#mu, gen}", .def_range={0,5}});
sh.AddNewFillParams("NGenLepFromW",         { .nbin=   5, .bins={    0,       5}, .fill=[this] { return v.GenPart.LeptonFromW.n;   }, .axis_title="N_{lep, gen}", .def_range={0,5}});
sh.AddNewFillParams("NGenEleFromZ",         { .nbin=   5, .bins={    0,       5}, .fill=[this] { return v.GenPart.EleFromZ.n;      }, .axis_title="N_{e, gen}",   .def_range={0,5}});
sh.AddNewFillParams("NGenMuFromZ",          { .nbin=   5, .bins={    0,       5}, .fill=[this] { return v.GenPart.MuFromZ.n;       }, .axis_title="N_{#mu, gen}", .def_range={0,5}});
sh.AddNewFillParams("NGenLepFromZ",         { .nbin=   5, .bins={    0,       5}, .fill=[this] { return v.GenPart.LeptonFromZ.n;   }, .axis_title="N_{lep, gen}", .def_range={0,5}});
sh.AddNewFillParams("NGenEleFromTop",       { .nbin=   5, .bins={    0,       5}, .fill=[this] { return v.GenPart.EleFromTop.n;    }, .axis_title="N_{e, gen}",   .def_range={0,5}});
sh.AddNewFillParams("NGenMuFromTop",        { .nbin=   5, .bins={    0,       5}, .fill=[this] { return v.GenPart.MuFromTop.n;     }, .axis_title="N_{#mu, gen}", .def_range={0,5}});
sh.AddNewFillParams("NGenLepFromTop",       { .nbin=   5, .bins={    0,       5}, .fill=[this] { return v.GenPart.LeptonFromTop.n; }, .axis_title="N_{lep, gen}", .def_range={0,5}});

sh.AddNewFillParams("NGenEleAndMuFrom2Tops",        { .nbin=   3, .bins={    0,       3}, .fill=[this] { return (v.GenPart.EleFromTop.n+v.GenPart.MuFromTop.n); }, .axis_title="N_{ele, gen top}+N_{mu, gen top}", .def_range={}});           
sh.AddNewFillParams("NGenEleAndMuAndTauFrom2Tops",  { .nbin=   3, .bins={    0,       3}, .fill=[this] { return (v.GenPart.EleFromTop.n+v.GenPart.MuFromTop.n+v.GenPart.TauFromTop.n); }, .axis_title="N_{ele, gen top}+N_{mu, gen top}+N_{tau, gen top}", .def_range={}});
sh.AddNewFillParams("NGenEleAndMuFrom1Top",         { .nbin=   2, .bins={    0,       2}, .fill=[this] { return (v.GenPart.EleFromTop.n+v.GenPart.MuFromTop.n); }, .axis_title="N_{ele, gen top}+N_{mu, gen top}", .def_range={}});
sh.AddNewFillParams("NGenEleAndMuAndTauFrom1Top",   { .nbin=   2, .bins={    0,       2}, .fill=[this] { return (v.GenPart.EleFromTop.n+v.GenPart.MuFromTop.n+v.GenPart.TauFromTop.n); }, .axis_title="N_{ele, gen top}+N_{mu, gen top}+N_{tau, gen top}", .def_range={}});
sh.AddNewFillParams("NGenEleAndMuFrom1W",           { .nbin=   2, .bins={    0,       2}, .fill=[this] { return (v.GenPart.EleFromW.n+v.GenPart.MuFromW.n); }, .axis_title="N_{ele, gen W}+N_{mu, gen W}", .def_range={}});
sh.AddNewFillParams("NGenEleAndMuAndTauFrom1W",     { .nbin=   2, .bins={    0,       2}, .fill=[this] { return (v.GenPart.LeptonFromW.n); }, .axis_title="N_{ele, gen W}+N_{mu, gen W}+N_{tau, gen W}", .def_range={}});
sh.AddNewFillParams("NGenHadW",             { .nbin=   5, .bins={    0,       5}, .fill=[this] { return v.GenPart.HadW.n;          }, .axis_title="N_{W (had.), gen}", .def_range={0,5}});
sh.AddNewFillParams("NGenHadZ",             { .nbin=   5, .bins={    0,       5}, .fill=[this] { return v.GenPart.HadZ.n;          }, .axis_title="N_{Z (had.), gen}", .def_range={0,5}});
sh.AddNewFillParams("NGenHadTop",           { .nbin=   5, .bins={    0,       5}, .fill=[this] { return v.GenPart.HadTop.n;        }, .axis_title="N_{top (had.), gen}", .def_range={0,5}});
sh.AddNewFillParams("NGenLepTop",           { .nbin=   5, .bins={    0,       5}, .fill=[this] { return v.GenPart.LepTop.n;        }, .axis_title="N_{top (lep.), gen}", .def_range={0,5}});

// HT
//sh.AddNewFillParams("HT",                   { .nbin=HT.size()-1, .bins=HT,        .fill=[this] { return v.AK4_Ht;             }, .axis_title="H_{T} (GeV)",                 .def_range={200, 2000}});
sh.AddNewFillParams("HTFine",               { .nbin=  50, .bins={    0,    5000}, .fill=[this] { return v.AK4_Ht;             }, .axis_title="H_{T} (GeV)",                   .def_range={200, 2000}});
sh.AddNewFillParams("HT",                   { .nbin=HT.size()-1, .bins=HT,        .fill=[this] { return v.AK4_Ht;             }, .axis_title="H_{T} (GeV)",                   .def_range={200, 2000}});
sh.AddNewFillParams("HT2DBins",             { .nbin=HT_2D_bins.size()-1, .bins=HT_2D_bins, .fill=[this] { return v.AK4_Ht;             }, .axis_title="H_{T} (GeV)",          .def_range={200, 1500}});
sh.AddNewFillParams("OnlineHT",             { .nbin= 100, .bins={    0,    5000}, .fill=[this] { return v.AK4_HtOnline;       }, .axis_title="H_{T}^{HLT} (GeV)",             .def_range={200, 2000}});
sh.AddNewFillParams("AK8HT",                { .nbin=HT.size()-1, .bins=HT,        .fill=[this] { return v.AK8_Ht;             }, .axis_title="H_{T}^{AK8} (GeV)",             .def_range={0, 2000}});
// MET
sh.AddNewFillParams("METPhi",               { .nbin=  24, .bins={-3.142,  3.142}, .fill=[this] { return v.MET_phi;             }, .axis_title="MET #phi"});
sh.AddNewFillParams("MET",                  { .nbin=MET.size()-1, .bins=MET,      .fill=[this] { return v.MET_pt;              }, .axis_title="#slash{E}_{T} (GeV)",              .def_range={0,300}});
sh.AddNewFillParams("MET2DBins",            { .nbin=MET_2D_bins.size()-1, .bins=MET_2D_bins, .fill=[this] { return v.MET_pt;   }, .axis_title="#slash{E}_{T} (GeV)",              .def_range={0,600}});
sh.AddNewFillParams("METNoPho",             { .nbin=MET.size()-1, .bins=MET,      .fill=[this] { return std::sqrt(v.MET_pho.Perp2());   }, .axis_title="#slash{E}_{T, no #gamma} (GeV)", .def_range={0,2000}});
sh.AddNewFillParams("METNo1Lep",            { .nbin=MET.size()-1, .bins=MET,      .fill=[this] { return std::sqrt(v.MET_1l.Perp2());    }, .axis_title="#slash{E}_{T,no lep} (GeV)",     .def_range={0,2000}});
sh.AddNewFillParams("METNo1VLep",           { .nbin=MET.size()-1, .bins=MET,      .fill=[this] { return std::sqrt(v.MET_1vl.Perp2());   }, .axis_title="#slash{E}_{T,no lep} (GeV)",     .def_range={0,2000}});
sh.AddNewFillParams("METNo2Lep",            { .nbin=MET.size()-1, .bins=MET,      .fill=[this] { return std::sqrt(v.MET_2l.Perp2());    }, .axis_title="#slash{E}_{T,ll} (GeV)",         .def_range={0,2000}});
sh.AddNewFillParams("METNoDiLep",           { .nbin=MET.size()-1, .bins=MET,      .fill=[this] { return std::sqrt(v.MET_dilep.Perp2()); }, .axis_title="#slash{E}_{T,1l missing} (GeV)",         .def_range={0,2000}});
sh.AddNewFillParams("METFine",              { .nbin=  40, .bins={    0,    2000}, .fill=[this] { return v.MET_pt;                       }, .axis_title="#slash{E}_{T} (GeV)",            .def_range={0,2000}});
sh.AddNewFillParams("PuppiMET",             { .nbin=MET.size()-1, .bins=MET,      .fill=[this] { return v.PuppiMET_pt;                  }, .axis_title="#slash{E}_{T} (GeV)",            .def_range={0,2000}});
sh.AddNewFillParams("RawMET",               { .nbin=MET.size()-1, .bins=MET,      .fill=[this] { return v.RawMET_pt;                    }, .axis_title="#slash{E}_{T} (GeV)",            .def_range={0,2000}});
//sh.AddNewFillParams("TkMET",                { .nbin=MET.size()-1, .bins=MET,      .fill=[this] { return v.TkMET_pt;                     }, .axis_title="#slash{E}_{T} (GeV)",            .def_range={0,2000}});
sh.AddNewFillParams("PFMETOverCaloMET",     { .nbin=100, .bins={-10,10}, .fill=[this] { return v.CaloMET_pt>0 ? std::log(v.MET_pt/v.CaloMET_pt) : 9999; }, .axis_title="ln(PFMET / CaloMET)"});

// Unrolled HT vs AK8Pt1 - For trigger efficiency
add_unrolled_bins("HTMET",           "H_{T} (GeV)", "#slash{E}_{T} (GeV)", [this] { return v.AK4_Ht; }, [this] { return v.MET_pt; }, HT_2D_bins,  MET_2D_bins, merged_trigger_bins);
add_unrolled_bins("HT1MET",          "H_{T} (GeV)", "#slash{E}_{T} (GeV)", [this] { return v.AK4_Ht; }, [this] { return v.MET_pt; }, HT1_2D_bins, MET_2D_bins, merged_trigger1_bins);
add_unrolled_bins("HT2MET",          "H_{T} (GeV)", "#slash{E}_{T} (GeV)", [this] { return v.AK4_Ht; }, [this] { return v.MET_pt; }, HT2_2D_bins, MET_2D_bins, merged_trigger2_bins);
add_unrolled_bins("HT3MET",          "H_{T} (GeV)", "#slash{E}_{T} (GeV)", [this] { return v.AK4_Ht; }, [this] { return v.MET_pt; }, HT3_2D_bins, MET_2D_bins, merged_trigger3_bins);
add_unrolled_bins("HTMETPrev",       "H_{T} (GeV)", "#slash{E}_{T} (GeV)", [this] { return v.AK4_Ht; }, [this] { return v.MET_pt; }, {300, 500, 550, 600, 800, 1500, 10000}, {50, 150, 200, 400, 10000});
// HT vs Photon Pt
add_unrolled_bins("HTPhotonPt",      "H_{T} (GeV)", "Photon p_{T} (GeV)", [this] { return v.AK4_Ht; }, [this] { return v.Photon.Select.n<1? -9999 : v.Photon.Select(0).pt; }, HT_2D_bins,  PtPho_2D_bins, {}, 0, 0);
add_unrolled_bins("HT1PhotonPt",     "H_{T} (GeV)", "Photon p_{T} (GeV)", [this] { return v.AK4_Ht; }, [this] { return v.Photon.Select.n<1? -9999 : v.Photon.Select(0).pt; }, HT1_2D_bins,  PtPho_2D_bins, {}, 0, 0);
add_unrolled_bins("HT2PhotonPt",     "H_{T} (GeV)", "Photon p_{T} (GeV)", [this] { return v.AK4_Ht; }, [this] { return v.Photon.Select.n<1? -9999 : v.Photon.Select(0).pt; }, HT2_2D_bins,  PtPho_2D_bins, {}, 0, 0);
add_unrolled_bins("HT3PhotonPt",     "H_{T} (GeV)", "Photon p_{T} (GeV)", [this] { return v.AK4_Ht; }, [this] { return v.Photon.Select.n<1? -9999 : v.Photon.Select(0).pt; }, HT3_2D_bins,  PtPho_2D_bins, {}, 0, 0);
// Unrolled HT vs LeptonPt - For 2017-18 leptonic trigger efficiencies
add_unrolled_bins("HTElePt",         "H_{T} (GeV)", "Ele p_{T} (GeV)", [this] { return v.AK4_Ht; }, [this] { return v.Electron.Select.n<1 ? -9999 : v.Electron.Select(0).pt; }, HT_2D_bins_lep,  ElePt_2D_bins);
add_unrolled_bins("HT1ElePt",        "H_{T} (GeV)", "Ele p_{T} (GeV)", [this] { return v.AK4_Ht; }, [this] { return v.Electron.Select.n<1 ? -9999 : v.Electron.Select(0).pt; }, HT1_2D_bins_lep, ElePt_2D_bins);
add_unrolled_bins("HT2ElePt",        "H_{T} (GeV)", "Ele p_{T} (GeV)", [this] { return v.AK4_Ht; }, [this] { return v.Electron.Select.n<1 ? -9999 : v.Electron.Select(0).pt; }, HT2_2D_bins_lep, ElePt_2D_bins);
add_unrolled_bins("HTMuPt",          "H_{T} (GeV)", "Mu p_{T} (GeV)",  [this] { return v.AK4_Ht; }, [this] { return v.Muon.Select.n<1 ? -9999 : v.Muon.Select(0).pt;       }, HT_2D_bins_lep,  MuPt_2D_bins);
add_unrolled_bins("HT1MuPt",         "H_{T} (GeV)", "Mu p_{T} (GeV)",  [this] { return v.AK4_Ht; }, [this] { return v.Muon.Select.n<1 ? -9999 : v.Muon.Select(0).pt;       }, HT1_2D_bins_lep, MuPt_2D_bins);
add_unrolled_bins("HT2MuPt",         "H_{T} (GeV)", "Mu p_{T} (GeV)",  [this] { return v.AK4_Ht; }, [this] { return v.Muon.Select.n<1 ? -9999 : v.Muon.Select(0).pt;       }, HT2_2D_bins_lep, MuPt_2D_bins);

// DPhi
sh.AddNewFillParams("DeltaPhi",                 { .nbin=DP.size()-1,  .bins=DP,       .fill=[this] { return v.dPhiRazor;               }, .axis_title="#Delta#phi_{megajets}"});
sh.AddNewFillParams("DeltaPhiNoPho",            { .nbin=DP.size()-1,  .bins=DP,       .fill=[this] { return v.dPhiRazorNoPho;          }, .axis_title="#Delta#phi_{megajets, no #gamma}"});
sh.AddNewFillParams("MinDeltaPhi",              { .nbin=  64, .bins={    0,     3.2}, .fill=[this] { return v.minDeltaPhi;             }, .axis_title="#Delta#phi_{min}"});
//sh.AddNewFillParams("MinDeltaPhi",              { .nbin=MDP.size()-1, .bins=MDP,      .fill=[this] { return v.minDeltaPhi;             }, .axis_title="#Delta#phi_{min}"});
sh.AddNewFillParams("MinDeltaPhiNo1Lep",        { .nbin=MDP.size()-1, .bins=MDP,      .fill=[this] { return v.minDeltaPhi_1l;          }, .axis_title="#Delta#phi_{min,no lep}"});
sh.AddNewFillParams("MinDeltaPhiNo1VLep",       { .nbin=MDP.size()-1, .bins=MDP,      .fill=[this] { return v.minDeltaPhi_1vl;         }, .axis_title="#Delta#phi_{min,no lep}"});
sh.AddNewFillParams("MinDeltaPhiNo2Lep",        { .nbin=MDP.size()-1, .bins=MDP,      .fill=[this] { return v.minDeltaPhi_2l;          }, .axis_title="#Delta#phi_{min,ll}"});
sh.AddNewFillParams("MinDeltaPhiNoPho",         { .nbin=MDP.size()-1, .bins=MDP,      .fill=[this] { return v.minDeltaPhi_pho;         }, .axis_title="#Delta#phi_{min, no #gamma}"});
sh.AddNewFillParams("DeltaPhiLLMET",            { .nbin=MDP.size()-1, .bins=MDP,      .fill=[this] { return v.dPhi_2l_met;             }, .axis_title="#Delta#phi (ll, MET)"});
sh.AddNewFillParams("DeltaPhiLLJet",            { .nbin=MDP.size()-1, .bins=MDP,      .fill=[this] { return v.dPhi_2l_jet;             }, .axis_title="#Delta#phi_{min} (ll, jet)"});
sh.AddNewFillParams("DeltaPhiBoostedJetMET",    { .nbin=DP.size()-1,  .bins=DP,       .fill=[this] { return v.dPhiBoostedJetMET;       }, .axis_title="#Delta#phi (AK8 jet, MET)"});
sh.AddNewFillParams("DeltaPhiBoostedJetLep",    { .nbin=DP.size()-1,  .bins=DP,       .fill=[this] { return v.dPhiBoostedJetLep;       }, .axis_title="#Delta#phi (AK8 jet, lep)"});
sh.AddNewFillParams("DeltaPhiBoostedJetLepMET", { .nbin=DP.size()-1,  .bins=DP,       .fill=[this] { return v.dPhiBoostedJetLepMET;    }, .axis_title="#Delta#phi (AK8 jet, lep+MET)"});
sh.AddNewFillParams("DeltaPhiMuonJetMET",       { .nbin=8, .bins={0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.74159, 3.1416}, .fill=[this] { return v.dPhiMuonJetMET; }, .axis_title="#Delta#phi (#mu jet, MET)"});
// MT/Mll
sh.AddNewFillParams("MT",                   { .nbin=  50, .bins={    0,    1000}, .fill=[this] { return v.MT_lepVeto;              }, .axis_title="m_{T} (GeV)",  .def_range={0,500}});
sh.AddNewFillParams("MTSelect",             { .nbin=  50, .bins={    0,    1000}, .fill=[this] { return v.MT;                      }, .axis_title="m_{T} (GeV)",  .def_range={0,500}});
sh.AddNewFillParams("MTTight",              { .nbin=  50, .bins={    0,    1000}, .fill=[this] { return v.MT_lepTight;             }, .axis_title="m_{T} (GeV)",  .def_range={0,500}});
sh.AddNewFillParams("MTNonIso",             { .nbin=  50, .bins={    0,    1000}, .fill=[this] { return v.MT_lepNonIso;            }, .axis_title="m_{T} (GeV)",  .def_range={0,500}});
sh.AddNewFillParams("MTBoost",              { .nbin=  20, .bins={    0,    4000}, .fill=[this] { return v.MT_boost;                }, .axis_title="m_{T,Boost+MET} (GeV)",  .def_range={0,2000}});
sh.AddNewFillParams("Mll",                  { .nbin=  50, .bins={    0,     500}, .fill=[this] { return v.M_2l;                    }, .axis_title="m_{ll} (GeV)", .def_range={0,200}});
// SUSY
sh.AddNewFillParams("MGluino",              { .nbin= 121, .bins={-12.5, 3012.5 }, .fill=[this] { return v.susy_mass[0];      }, .axis_title="m_{#tilde{g}} (GeV)",        .def_range={550,2350}});
sh.AddNewFillParams("MSquark",              { .nbin=  81, .bins={-12.5, 2012.5 }, .fill=[this] { return v.susy_mass[0];      }, .axis_title="m_{#tilde{t}} (GeV)",        .def_range={  0,1650}});
sh.AddNewFillParams("MEWK",                 { .nbin=  81, .bins={-12.5, 2012.5 }, .fill=[this] { return v.susy_mass[0];      }, .axis_title="m_{#tilde{#chi}^{#pm}_{1}} (GeV)",        .def_range={  0,1650}});
sh.AddNewFillParams("MLSP",                 { .nbin=  81, .bins={-12.5, 2012.5 }, .fill=[this] { return v.susy_mass[1];     }, .axis_title="m_{#tilde{#chi}^{0}_{1}} (GeV)", .def_range={  0,1650}});
sh.AddNewFillParams("SquarkLSPMassDiff",    { .nbin= 400, .bins={0, 2000 },       .fill=[this] { return v.susy_mass[0]-v.susy_mass[1]; }, .axis_title="m_{#tilde{t}}-m_{#tilde{#chi}^{0}_{1}} (GeV)"});
if (debug) std::cout<<"PlottingBase::define_histo_settings: non-special fillparams ok"<<std::endl;


// --------------------------------------------------------------------
//              Special Fill Parameters (not for binning)
// --------------------------------------------------------------------

// ------------------------- Special -----------------------
//     e.g. Efficiencies, averages, fake rates etc.

// Special Y/Z axis parameters:
// Systematics
sh.AddSpecial({ .name="Counts", .name_plus_1d="Syst", .axis="Events / bin", .axis_plus_1d="Systematics variation index"});
//sh.AddSpecial({ .name="SignalSelectionEfficiency",    .name_plus_1d="PassSignalSelection",    .axis="Signal Selection Efficiency", .axis_plus_1d="Pass Signal Selection"});
sh.AddSpecial({ .name="SignalSignificance_T5ttcc",    .name_plus_1d="Bkg_T5ttcc",             .axis="S/#sqrt{S+B} - T5ttcc",        .axis_plus_1d="Background, Signal - T5ttcc"});
sh.AddSpecial({ .name="SignalSignificance_TChiZH",    .name_plus_1d="Bkg_TChiZH",                 .axis="S/#sqrt{S+B} - TChiZH",        .axis_plus_1d="Background, Signal - TChiZH"});







sh.AddNewFillParams("Counts",                         { .nbin= 1+syst_nSyst, .bins={-0.5, syst_nSyst+0.5}, .fill=[&syst_index] { return syst_index; }, .axis_title="Events / bin"});

sh.AddNewSpecialFillParams("HLTEff_HT_MET_or_HT",           { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[this] {
		return 
		v.HLT_PFHT1050==1 ||
		//v.HLT_PFMET120_PFMHT120_IDTight==1 ||
		//v.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight==1 ||
		//v.HLT_PFMETTypeOne120_PFMHT120_IDTight==1 ||
		v.HLT_PFHT500_PFMET100_PFMHT100_IDTight==1 || 
		v.HLT_PFHT700_PFMET85_PFMHT85_IDTight==1 ||
		v.HLT_PFHT800_PFMET75_PFMHT75_IDTight==1;
		}, .axis_title="#epsilon_{HT*_MET* OR HT1050}",               .def_range={0,1} });
sh.AddNewSpecialFillParams("HLTEff_AllMET_or_HT_MET_or_HT", { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[this] { 
		return 
		v.HLT_PFHT1050==1 ||
		v.HLT_PFMET120_PFMHT120_IDTight==1 ||
		v.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight==1 ||
		v.HLT_PFMETTypeOne120_PFMHT120_IDTight==1 ||
		v.HLT_PFHT500_PFMET100_PFMHT100_IDTight==1 || 
		v.HLT_PFHT700_PFMET85_PFMHT85_IDTight==1 ||
		v.HLT_PFHT800_PFMET75_PFMHT75_IDTight==1;    
		}, .axis_title="#epsilon_{AllMET120 OR HT*_MET* OR HT1050}",  .def_range={0,1} });
// Leptonic triggers
sh.AddNewSpecialFillParams("HLTEff_Mu15_HT450",      { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[this] { return v.HLT_Mu15_IsoVVVL_PFHT450==1;      }, .axis_title="#epsilon_{Mu15_HT450}",  .def_range={0,1} });
sh.AddNewSpecialFillParams("HLTEff_IsoMu27",         { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[this] { return v.HLT_IsoMu27==1;                   }, .axis_title="#epsilon_{IsoMu27}",     .def_range={0,1} });
sh.AddNewSpecialFillParams("HLTEff_Mu50",            { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[this] { return v.HLT_Mu50==1;                      }, .axis_title="#epsilon_{Mu50}",        .def_range={0,1} });
sh.AddNewSpecialFillParams("HLTEff_Ele15_HT450",     { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[this] { return v.HLT_Ele15_IsoVVVL_PFHT450==1;     }, .axis_title="#epsilon_{Ele15_HT450}", .def_range={0,1} });
sh.AddNewSpecialFillParams("HLTEff_Ele35",           { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[this] { return v.HLT_Ele35_WPTight_Gsf==1;         }, .axis_title="#epsilon_{Ele35}",       .def_range={0,1} });
sh.AddNewSpecialFillParams("HLTEff_Ele115",          { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[this] { return v.HLT_Ele115_CaloIdVT_GsfTrkIdT==1; }, .axis_title="#epsilon_{Ele115}",      .def_range={0,1} });
sh.AddNewSpecialFillParams("HLTEff_IsoMu27_or_Mu50", { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[this] { return 
		v.HLT_IsoMu27==1 ||
		v.HLT_Mu50==1;
		}, .axis_title="#epsilon_{IsoMu27 OR Mu50}",                .def_range={0,1} });
sh.AddNewSpecialFillParams("HLTEff_IsoMu27_or_Mu50_or_Mu15_HT450",  { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[this] { return
		v.HLT_IsoMu27==1 ||
		v.HLT_Mu50==1 ||
		v.HLT_Mu15_IsoVVVL_PFHT450==1;
		}, .axis_title="#epsilon_{IsoMu27 OR Mu50 OR Mu15_HT450}",  .def_range={0,1} });
sh.AddNewSpecialFillParams("HLTEff_Ele35_or_Ele115",                { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[this] { return 
		v.HLT_Ele35_WPTight_Gsf==1 ||
		v.HLT_Ele115_CaloIdVT_GsfTrkIdT==1;
		}, .axis_title="#epsilon_{Ele35 OR Ele115}",                .def_range={0,1} });
sh.AddNewSpecialFillParams("HLTEff_Ele35_or_Ele115_or_Ele15_HT450", { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[this] { return 
		v.HLT_Ele35_WPTight_Gsf==1 ||
		v.HLT_Ele115_CaloIdVT_GsfTrkIdT==1 ||
		v.HLT_Ele15_IsoVVVL_PFHT450==1;
		}, .axis_title="#epsilon_{Ele35 OR Ele115 OR Ele15_HT450}", .def_range={0,1} });
sh.AddNewSpecialFillParams("HLTEff_Hadronic", { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[this] {
		if (v.year==2016) {
		return
		// JetHT
		v.HLT_PFHT800==1 ||
		v.HLT_PFHT900==1 ||
		v.HLT_AK8PFJet450==1 ||
		// MET
		v.HLT_PFMET110_PFMHT110_IDTight==1 ||
		v.HLT_PFMETNoMu110_PFMHTNoMu110_IDTight==1 ||
		// HTMHT
		v.HLT_PFHT300_PFMET110==1;
		} else {
		return
		// JetHT
		v.HLT_PFHT1050==1 ||
		v.HLT_AK8PFJet500==1 ||
		// MET
		v.HLT_PFMET120_PFMHT120_IDTight==1 ||
		v.HLT_PFMET120_PFMHT120_IDTight_PFHT60==1 ||
		v.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight==1 ||
		v.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60==1 ||
		v.HLT_PFHT500_PFMET100_PFMHT100_IDTight==1 || 
		v.HLT_PFHT700_PFMET85_PFMHT85_IDTight==1 ||
		v.HLT_PFHT800_PFMET75_PFMHT75_IDTight==1;
		}
}, .axis_title="#epsilon_{MET120 OR HT*_MET* OR HT*}",     .def_range={0,1} });
sh.AddNewSpecialFillParams("HLTEff_SingleMu", { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[this] { 
		if (v.year==2016) return v.HLT_IsoMu24==1 || v.HLT_IsoTkMu24==1;
		else return
		v.HLT_IsoMu27==1 || v.HLT_IsoTkMu27==1;
}, .axis_title="#epsilon_{Single muon}", .def_range={0,1} });
sh.AddNewSpecialFillParams("HLTEff_SingleEle", { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[this] { 
		if (v.year==2016) return v.HLT_Ele27_WPTight_Gsf==1;
		else return v.HLT_Ele32_WPTight_Gsf==1;
}, .axis_title="#epsilon_{Single electron}", .def_range={0,1} });
sh.AddNewSpecialFillParams("HLTEff_SinglePho", { .nbin=    2, .bins={ -0.5,     1.5}, .fill=[this] { 
		if (v.year==2016) return v.HLT_Photon175==1;
		else return v.HLT_Photon200==1;
		}, .axis_title="#epsilon_{Single photon}", .def_range={0,1} });

if (debug) std::cout<<"PlottingBase::define_histo_settings: special fillparams ok"<<std::endl;

// ----------------------- ROC Curves ----------------------

// Leptons
// GenParticles for GenTruth matching
std::map<std::string, std::function<int()> > gen_leptons;
gen_leptons["EleFromW"]           = [this] { return (int)v.Electron().matchGenEleFromW; };
gen_leptons["MuFromW"]            = [this] { return (int)v.Muon().matchGenMuFromW; };
gen_leptons["EleFromTop"]         = [this] { return (int)v.Electron().matchGenEleFromTop; };
gen_leptons["EleNuFromTop"]       = [this] { return (int)v.Electron().nuMatchGenEleNuFromTop; };
gen_leptons["MuFromTop"]          = [this] { return (int)v.Muon().matchGenMuFromTop; };
//   gen_leptons["TauFromTop"]         = [this] { return (int)v.Tau().matchGenTauFromTop; };
gen_leptons["MuNuFromTop"]        = [this] { return (int)v.Muon().nuMatchGenMuNuFromTop; };
gen_leptons["EleFromHardProcess"] = [this] { return (int)v.Electron().matchGenEleFromHardProcess; };
gen_leptons["MuFromHardProcess"]  = [this] { return (int)v.Muon().matchGenMuFromHardProcess; };
// GenTruth for Leptons
for (const auto& genp : gen_leptons) for (const auto& bm : all_benchmarks)
	sh.AddNewFillParams(genp.first+"_"+bm.first, { .nbin= 2, .bins={-0.5,1.5}, .fill=[genp,bm] { int signal = bm.second(); return signal==1 ? genp.second() : signal; }});

	// cuts not for comparison
	sh.AddNewFillParams("elept5_cut",           { .nbin=1,  .bins={5, 5000}, .fill=[this] { return v.Electron().pt; }, .axis_title="p_{T}^{e}", .def_range={5, 5000}});
	sh.AddNewFillParams("elept10_cut",          { .nbin=1,  .bins={10,5000}, .fill=[this] { return v.Electron().pt; }, .axis_title="p_{T}^{e}", .def_range={10,5000}});
	sh.AddNewFillParams("elept30_cut",          { .nbin=1,  .bins={30,5000}, .fill=[this] { return v.Electron().pt; }, .axis_title="p_{T}^{e}", .def_range={30,5000}});
	sh.AddNewFillParams("mupt5_cut",            { .nbin=1,  .bins={5 ,5000}, .fill=[this] { return v.Muon().pt; }, .axis_title="p_{T}^{#mu}", .def_range={5, 5000}});
	sh.AddNewFillParams("mupt10_cut",           { .nbin=1,  .bins={10,5000}, .fill=[this] { return v.Muon().pt; }, .axis_title="p_{T}^{#mu}", .def_range={10,5000}});
	sh.AddNewFillParams("mupt30_cut",           { .nbin=1,  .bins={30,5000}, .fill=[this] { return v.Muon().pt; }, .axis_title="p_{T}^{#mu}", .def_range={30,5000}});
	sh.AddNewFillParams("eleeta_veto_cut",      { .nbin=1,  .bins={1,2}, .fill=[this]     { return std::abs(v.Electron().eta)<2.5; }, .axis_title="#eta_{e}", .def_range={1,2}});
	sh.AddNewFillParams("eleeta_cut",           { .nbin=1,  .bins={2,3}, .fill=[this] { 
			double abseta = std::abs(v.Electron().eta);
			if      (abseta>=2.5)                 return 0;
			else if (abseta>=1.442&&abseta<1.556) return 1;
			else return 2; }, .axis_title="#eta_{e}", .def_range={2,3}});
	sh.AddNewFillParams("mueta_cut",            { .nbin=1,  .bins={1,2}, .fill=[this] { return std::abs(v.Muon().eta)<2.4; }, .axis_title="#eta_{#mu}", .def_range={1,2}});
	sh.AddNewFillParams("eleip_loose_cut",      { .nbin=1,  .bins={1,2}, .fill=[this] { return std::abs(v.Electron().dz)<0.50&&std::abs(v.Electron().dxy)<0.2; }, .axis_title="vtx. d_{xy,e}/d_{z,e}", .def_range={1,2}});
	sh.AddNewFillParams("eleip_medium_cut",     { .nbin=1,  .bins={1,2}, .fill=[this] { return std::abs(v.Electron().dz)<0.10&&std::abs(v.Electron().dxy)<0.05; }, .axis_title="vtx. d_{xy,e}/d_{z,e}", .def_range={1,2}});
	sh.AddNewFillParams("eleip_tight_cut",      { .nbin=1,  .bins={1,2}, .fill=[this] { return std::abs(v.Electron().dz)<0.10&&std::abs(v.Electron().dxy)<0.05&&v.Electron().sip3d<4; }, .axis_title="vtx. d_{xy,e}/d_{z,e}", .def_range={1,2}});
	sh.AddNewFillParams("muip_loose_cut",       { .nbin=1,  .bins={1,2}, .fill=[this] { return std::abs(v.Muon().dz)<0.50&&std::abs(v.Muon().dxy)<0.2; }, .axis_title="vtx. d_{xy,#mu}/d_{z,#mu}", .def_range={1,2}});
	sh.AddNewFillParams("muip_medium_cut",      { .nbin=1,  .bins={1,2}, .fill=[this] { return std::abs(v.Muon().dz)<0.10&&std::abs(v.Muon().dxy)<0.05; }, .axis_title="vtx. d_{xy,#mu}/d_{z,#mu}", .def_range={1,2}});
	sh.AddNewFillParams("muip_tight_cut",       { .nbin=1,  .bins={1,2}, .fill=[this] { return std::abs(v.Muon().dz)<0.10&&std::abs(v.Muon().dxy)<0.05&&v.Muon().sip3d<4; }, .axis_title="vtx. d_{xy,#mu}/d_{z,#mu}", .def_range={1,2}});
	sh.AddNewFillParams("muminiiso_loose_cut",  { .nbin=1,  .bins={1,2}, .fill=[this] { return v.Muon().miniPFRelIso_all<0.4; }, .axis_title="#mu Mini iso.", .def_range={1,2}});
	sh.AddNewFillParams("muminiiso_medium_cut", { .nbin=1,  .bins={1,2}, .fill=[this] { return v.Muon().miniPFRelIso_all<0.2; }, .axis_title="#mu Mini iso.", .def_range={1,2}});
	sh.AddNewFillParams("mu2diso15_cut",        { .nbin=1,  .bins={1,2}, .fill=[this] { return !(v.Muon().jetDRmin<0.4&&v.Muon().cleanJetPtrel<15); }, .axis_title="#mu 2D iso.", .def_range={1,2}});
	sh.AddNewFillParams("muid_soft_cut",        { .nbin=1,  .bins={1,2}, .fill=[this] { return v.Muon().softId; }, .axis_title="Soft (Cut-based) ID", .def_range={1,2}});
	sh.AddNewFillParams("muid_medium_cut",      { .nbin=1,  .bins={1,2}, .fill=[this] { return v.Muon().mediumId; }, .axis_title="Medium (Cut-based) ID", .def_range={1,2}});
	sh.AddNewFillParams("muid_mediumprompt_cut",{ .nbin=1,  .bins={1,2}, .fill=[this] { return v.Muon().mediumPromptId; }, .axis_title="Medium (Cut-based) ID", .def_range={1,2}});
	sh.AddNewFillParams("muid_tight_cut",       { .nbin=1,  .bins={1,2}, .fill=[this] { return v.Muon().tightId; }, .axis_title="Tight (Cut-based) ID", .def_range={1,2}});
	sh.AddNewFillParams("muid_mva_loose_cut",   { .nbin=1,  .bins={1,2}, .fill=[this] { return v.Muon().mvaId>0; }, .axis_title="#mu MVA Loose ID", .def_range={1,2}});
	sh.AddNewFillParams("muid_mva_medium_cut",  { .nbin=1,  .bins={1,2}, .fill=[this] { return v.Muon().mvaId>1; }, .axis_title="#mu MVA Medium ID", .def_range={1,2}});
	sh.AddNewFillParams("eleid_mva_loose_cut",  { .nbin=1,  .bins={1,2}, .fill=[this] { return v.Electron().mvaFall17V2noIso_WPL; }, .axis_title="e MVA Loose ID (no iso.)", .def_range={1,2}});
	sh.AddNewFillParams("eleid_mva_medium_cut", { .nbin=1,  .bins={1,2}, .fill=[this] { return v.Electron().mvaFall17V2noIso_WP90; }, .axis_title="e MVA Medium ID (no iso.)", .def_range={1,2}});
	// pt/eta - can easily be tuned
	// pt      - black 1
	// eta     - green 418
	sh.AddNewFillParams("elept5",                  { .nbin=23,  .bins={0,5,10,15,20,25,30,35,40,45,50,60,70,80,100,125,150,200,250,300,400,500,1000,5000}, .fill=[this] { return v.Electron().pt; }, .axis_title="p_{T}^{e}", .def_range={5, 5000, 1}});
	sh.AddNewFillParams("elept10",                 { .nbin=23,  .bins={0,5,10,15,20,25,30,35,40,45,50,60,70,80,100,125,150,200,250,300,400,500,1000,5000}, .fill=[this] { return v.Electron().pt; }, .axis_title="p_{T}^{e}", .def_range={10,5000, 1}});
	sh.AddNewFillParams("elept30",                 { .nbin=23,  .bins={0,5,10,15,20,25,30,35,40,45,50,60,70,80,100,125,150,200,250,300,400,500,1000,5000}, .fill=[this] { return v.Electron().pt; }, .axis_title="p_{T}^{e}", .def_range={30,5000, 1}});
	sh.AddNewFillParams("mupt5",                   { .nbin=23,  .bins={0,5,10,15,20,25,30,35,40,45,50,60,70,80,100,125,150,200,250,300,400,500,1000,5000}, .fill=[this] { return v.Muon().pt; }, .axis_title="p_{T}^{#mu}", .def_range={5, 5000, 1}});
	sh.AddNewFillParams("mupt10",                  { .nbin=23,  .bins={0,5,10,15,20,25,30,35,40,45,50,60,70,80,100,125,150,200,250,300,400,500,1000,5000}, .fill=[this] { return v.Muon().pt; }, .axis_title="p_{T}^{#mu}", .def_range={10,5000, 1}});
	sh.AddNewFillParams("mupt30",                  { .nbin=23,  .bins={0,5,10,15,20,25,30,35,40,45,50,60,70,80,100,125,150,200,250,300,400,500,1000,5000}, .fill=[this] { return v.Muon().pt; }, .axis_title="p_{T}^{#mu}", .def_range={30,5000, 1}});
	sh.AddNewFillParams("eleeta_veto",             { .nbin=4, .bins={0,4}, .fill=[this] {
			double abseta = std::abs(v.Electron().eta);
			if      (abseta>=2.5)                 return 0;
			else if (abseta>=1.442&&abseta<1.556) return 1;
			else if (abseta>=2.1)                 return 2;
			else                                  return 3;
			}, .axis_title="#eta_{e}", .def_range={1,4, 1}});
sh.AddNewFillParams("eleeta",                  { .nbin=4,  .bins={0,4}, .fill=[this] { 
		double abseta = std::abs(v.Electron().eta);
		if      (abseta>=2.5)                 return 0;
		else if (abseta>=1.442&&abseta<1.556) return 1;
		else if (abseta>=2.1)                 return 2;
		else                                  return 3;
		}, .axis_title="#eta_{#mu}", .def_range={2,4, 1}});
sh.AddNewFillParams("mueta",                   { .nbin=3,   .bins={0,3}, .fill=[this] { 
		double abseta = std::abs(v.Muon().eta);
		if      (abseta>=2.4) return 0;
		else if (abseta>=2.1) return 1;
		else                  return 2; 
		}, .axis_title="#eta_{#mu}", .def_range={1,3, 1}});
// impact parameters - simple pick
// ip      - green 418;
sh.AddNewFillParams("eleip_loose",             { .nbin=3, .bins={1,4}, .fill=[this] {
		if      (std::abs(v.Electron().dz)>=0.50||std::abs(v.Electron().dxy)>=0.20) return 0;
		else if (std::abs(v.Electron().dz)>=0.10||std::abs(v.Electron().dxy)>=0.05) return 1;
		else if (v.Electron().sip3d>=4) return 2;
		else                                      return 3; 
		}, .axis_title="Impact param.", .def_range={0.9,4, 418}});
sh.AddNewFillParams("eleip_medium",            { .nbin=3, .bins={1,4}, .fill=[this] {
		if      (std::abs(v.Electron().dz)>=0.50||std::abs(v.Electron().dxy)>=0.20) return 0;
		else if (std::abs(v.Electron().dz)>=0.10||std::abs(v.Electron().dxy)>=0.05) return 1;
		else if (v.Electron().sip3d>=4) return 2;
		else                                      return 3; 
		}, .axis_title="Impact param.", .def_range={2,4, 418}});
sh.AddNewFillParams("eleip_tight",             { .nbin=3, .bins={1,4}, .fill=[this] {
		if      (std::abs(v.Electron().dz)>=0.50||std::abs(v.Electron().dxy)>=0.20) return 0;
		else if (std::abs(v.Electron().dz)>=0.10||std::abs(v.Electron().dxy)>=0.05) return 1;
		else if (v.Electron().sip3d>=4) return 2;
		else                                      return 3; 
		}, .axis_title="Impact param.", .def_range={3,4, 418}});
sh.AddNewFillParams("muip_loose",              { .nbin=3, .bins={1,4}, .fill=[this] {
		if      (std::abs(v.Muon().dz)>=0.50||std::abs(v.Muon().dxy)>=0.20) return 0;
		else if (std::abs(v.Muon().dz)>=0.10||std::abs(v.Muon().dxy)>=0.05) return 1;
		else if (v.Muon().sip3d>=4) return 2;
		else                              return 3; 
		}, .axis_title="Impact param.", .def_range={0.9,4, 418}});
sh.AddNewFillParams("muip_medium",             { .nbin=3, .bins={1,4}, .fill=[this] {
		if      (std::abs(v.Muon().dz)>=0.50||std::abs(v.Muon().dxy)>=0.20) return 0;
		else if (std::abs(v.Muon().dz)>=0.10||std::abs(v.Muon().dxy)>=0.05) return 1;
		else if (v.Muon().sip3d>=4) return 2;
		else                              return 3; 
		}, .axis_title="Impact param.", .def_range={2,4, 418}});
sh.AddNewFillParams("muip_tight",              { .nbin=3, .bins={1,4}, .fill=[this] {
		if      (std::abs(v.Muon().dz)>=0.50||std::abs(v.Muon().dxy)>=0.20) return 0;
		else if (std::abs(v.Muon().dz)>=0.10||std::abs(v.Muon().dxy)>=0.05) return 1;
		else if (v.Muon().sip3d>=4) return 2;
		else                              return 3; 
		}, .axis_title="Impact param.", .def_range={3,4, 418}});
// isolations - mini iso recommended, but we can cross-check others
// pfiso   - yellow 402
// miniiso - orange 801
// 2diso   - gray   921
sh.AddNewFillParams("elepfiso",                { .nbin=4, .bins={1,5}, .fill=[this] { 
		// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Offline_selection_criteria_for_V
		float pt = v.Electron().pt, pfiso = v.Electron().pfRelIso03_all;
		if (std::abs(v.Electron().eta+v.Electron().deltaEtaSC)<=1.479) {
		if      (pfiso>=0.1980+0.506/pt) return 0;
		else if (pfiso>=0.1120+0.506/pt) return 1;
		else if (pfiso>=0.0478+0.506/pt) return 2;
		else if (pfiso>=0.0287+0.506/pt) return 3;
		else                             return 4;
		} else {
		if      (pfiso>=0.2030+0.963/pt) return 0;
		else if (pfiso>=0.1080+0.963/pt) return 1;
		else if (pfiso>=0.0658+0.963/pt) return 2;
		else if (pfiso>=0.0445+0.963/pt) return 3;
		else                             return 4;
		} }, .axis_title="PF iso.", .def_range={0,5, 402}});
sh.AddNewFillParams("eleminiiso",              { .nbin=4, .bins={1,5}, .fill=[this] {
		float miniIso = v.Electron().miniPFRelIso_all;
		if      (miniIso>=0.40) return 0;
		else if (miniIso>=0.20) return 1;
		else if (miniIso>=0.10) return 2;
		else if (miniIso>=0.05) return 3;
		else                    return 4;
		}, .axis_title="Mini iso.", .def_range={0,5, 801}});
sh.AddNewFillParams("eleminiiso_loose",        { .nbin=4, .bins={1,5}, .fill=[this] {
		float miniIso = v.Electron().miniPFRelIso_all;
		if      (miniIso>=0.40) return 0;
		else if (miniIso>=0.20) return 1;
		else if (miniIso>=0.10) return 2;
		else if (miniIso>=0.05) return 3;
		else                    return 4;
		}, .axis_title="Mini iso.", .def_range={0.9,5, 801}});
sh.AddNewFillParams("eleminiiso_tight",        { .nbin=4, .bins={1,5}, .fill=[this] {
		float miniIso = v.Electron().miniPFRelIso_all;
		if      (miniIso>=0.40) return 0;
		else if (miniIso>=0.20) return 1;
		else if (miniIso>=0.10) return 2;
		else if (miniIso>=0.05) return 3;
		else                    return 4;
		}, .axis_title="Mini iso.", .def_range={3,5, 402}});
sh.AddNewFillParams("ele2diso15",              { .nbin=10, .bins={1,11}, .fill=[this] {
		if (v.Electron().jetDRmin<0.4) for (int i=0; i<10; ++i) 
		if (v.Electron().cleanJetPtrel<(i+1)*5) return i;
		return 10;
		}, .axis_title="2D iso.", .def_range={3,11, 921}});
sh.AddNewFillParams("ele2diso",                { .nbin=10, .bins={1,11}, .fill=[this] {
		if (v.Electron().jetDRmin<0.4) for (int i=0; i<10; ++i) 
		if (v.Electron().cleanJetPtrel<(i+1)*5) return i;
		return 10;
		}, .axis_title="2D iso.", .def_range={0,11, 921}});
sh.AddNewFillParams("mupfiso",                 { .nbin=8, .bins={-1,7}, .fill=[this] { return v.Muon().pfIsoId; }, .axis_title="#mu PF iso.", .def_range={0,7}}); 
//sh.AddNewFillParams("mumultiiso",              { .nbin=4, .bins={-1,3}, .fill=[this] { return v.Muon().multiIsoId; }, .axis_title="#mu Multi iso.", .def_range={0,7}});// currently unavailable
sh.AddNewFillParams("muminiiso",               { .nbin=4, .bins={1,5}, .fill=[this] {
		float miniIso = v.Muon().miniPFRelIso_all;
		if      (miniIso>=0.40) return 0;
		else if (miniIso>=0.20) return 1;
		else if (miniIso>=0.10) return 2;
		else if (miniIso>=0.05) return 3;
		else                    return 4;
		}, .axis_title="Mini iso.", .def_range={0,5, 801}});
sh.AddNewFillParams("muminiiso_loose",         { .nbin=4, .bins={1,5}, .fill=[this] {
		float miniIso = v.Muon().miniPFRelIso_all;
		if      (miniIso>=0.40) return 0;
		else if (miniIso>=0.20) return 1;
		else if (miniIso>=0.10) return 2;
		else if (miniIso>=0.05) return 3;
		else                    return 4;
		}, .axis_title="Mini iso.", .def_range={0.9,5, 801}});
sh.AddNewFillParams("muminiiso_medium",        { .nbin=4, .bins={1,5}, .fill=[this] {
		float miniIso = v.Muon().miniPFRelIso_all;
		if      (miniIso>=0.40) return 0;
		else if (miniIso>=0.20) return 1;
		else if (miniIso>=0.10) return 2;
		else if (miniIso>=0.05) return 3;
		else                    return 4;
		}, .axis_title="Mini iso.", .def_range={2,5, 801}});
sh.AddNewFillParams("mu2diso",                 { .nbin=10, .bins={1,11}, .fill=[this] {
		if (v.Muon().jetDRmin<0.4) for (int i=0; i<10; ++i) 
		if (v.Muon().cleanJetPtrel<(i+1)*5) return i;
		return 10;
		}, .axis_title="2D iso.", .def_range={0,11, 921}});
sh.AddNewFillParams("mu2diso15",               { .nbin=10, .bins={1,11}, .fill=[this] {
		if (v.Muon().jetDRmin<0.4) for (int i=0; i<10; ++i) 
		if (v.Muon().cleanJetPtrel<(i+1)*5) return i;
		return 10;
		}, .axis_title="2D iso.", .def_range={3,11, 921}});
// id with/without isolation
// id      - cyan to red 434,601,618,633
sh.AddNewFillParams("eleid_cut_iso",           { .nbin=4, .bins={1,5}, .fill=[this] { return v.Electron().cutBased; }, .axis_title="Cut-based ID", .def_range={0,  5, 601}});
sh.AddNewFillParams("eleid_cut_iso_veto",      { .nbin=4, .bins={1,5}, .fill=[this] { return v.Electron().cutBased; }, .axis_title="Cut-based ID", .def_range={0.9,5, 601}});
sh.AddNewFillParams("eleid_cut_iso_loose",     { .nbin=4, .bins={1,5}, .fill=[this] { return v.Electron().cutBased; }, .axis_title="Cut-based ID", .def_range={2,  5, 601}});
sh.AddNewFillParams("eleid_cut_iso_medium",    { .nbin=4, .bins={1,5}, .fill=[this] { return v.Electron().cutBased; }, .axis_title="Cut-based ID", .def_range={3,  5, 601}});
sh.AddNewFillParams("eleid_cut_iso_tight",     { .nbin=4, .bins={1,5}, .fill=[this] { return v.Electron().cutBased; }, .axis_title="Cut-based ID", .def_range={4,  5, 601}});
sh.AddNewFillParams("eleid_mva",               { .nbin=3, .bins={1,4}, .fill=[this] {
		if (!v.Electron().mvaFall17V2noIso_WPL)  return 0;
		if (!v.Electron().mvaFall17V2noIso_WP90) return 1;
		if (!v.Electron().mvaFall17V2noIso_WP80) return 2;
		else                                               return 3;
		}, .axis_title="MVA ID (no iso.)", .def_range={0,4, 633}});
sh.AddNewFillParams("eleid_mva_loose",         { .nbin=3, .bins={1,4}, .fill=[this] {
		if (!v.Electron().mvaFall17V2noIso_WPL)  return 0;
		if (!v.Electron().mvaFall17V2noIso_WP90) return 1;
		if (!v.Electron().mvaFall17V2noIso_WP80) return 2;
		else                                               return 3;
		}, .axis_title="MVA ID (no iso.)", .def_range={0.9,4, 633}});
sh.AddNewFillParams("eleid_mva_medium",        { .nbin=3, .bins={1,4}, .fill=[this] {
		if (!v.Electron().mvaFall17V2noIso_WPL)  return 0;
		if (!v.Electron().mvaFall17V2noIso_WP90) return 1;
		if (!v.Electron().mvaFall17V2noIso_WP80) return 2;
		else                                               return 3;
		}, .axis_title="MVA ID (no iso.)", .def_range={2,4, 633}});
sh.AddNewFillParams("eleid_mva_tight",         { .nbin=3, .bins={1,4}, .fill=[this] {
		if (!v.Electron().mvaFall17V2noIso_WPL)  return 0;
		if (!v.Electron().mvaFall17V2noIso_WP90) return 1;
		if (!v.Electron().mvaFall17V2noIso_WP80) return 2;
		else                                               return 3;
		}, .axis_title="MVA ID (no iso.)", .def_range={3,4, 633}});
sh.AddNewFillParams("eleid_mva_iso",           { .nbin=3, .bins={1,4}, .fill=[this] {
		if (!v.Electron().mvaFall17V2Iso_WPL)  return 0;
		if (!v.Electron().mvaFall17V2Iso_WP90) return 1;
		if (!v.Electron().mvaFall17V2Iso_WP80) return 2;
		else                                             return 3;
		}, .axis_title="MVA ID (w/ iso.)", .def_range={0,4, 618}});
sh.AddNewFillParams("eleid_mva_iso_loose",     { .nbin=3, .bins={1,4}, .fill=[this] {
		if (!v.Electron().mvaFall17V2Iso_WPL)  return 0;
		if (!v.Electron().mvaFall17V2Iso_WP90) return 1;
		if (!v.Electron().mvaFall17V2Iso_WP80) return 2;
		else                                             return 3;
		}, .axis_title="MVA ID (w/ iso.)", .def_range={0.9,4, 618}});
sh.AddNewFillParams("eleid_mva_iso_medium",    { .nbin=3, .bins={1,4}, .fill=[this] {
		if (!v.Electron().mvaFall17V2Iso_WPL)  return 0;
		if (!v.Electron().mvaFall17V2Iso_WP90) return 1;
		if (!v.Electron().mvaFall17V2Iso_WP80) return 2;
		else                                             return 3;
		}, .axis_title="MVA ID (w/ iso.)", .def_range={2,4, 618}});
sh.AddNewFillParams("eleid_mva_iso_tight",     { .nbin=3, .bins={1,4}, .fill=[this] {
		if (!v.Electron().mvaFall17V2Iso_WPL)  return 0;
		if (!v.Electron().mvaFall17V2Iso_WP90) return 1;
		if (!v.Electron().mvaFall17V2Iso_WP80) return 2;
		else                                             return 3;
		}, .axis_title="MVA ID (w/ iso.)", .def_range={3,4, 618}});
sh.AddNewFillParams("eleid_mva_miniiso",       { .nbin=6, .bins={2,8}, .fill=[this] {
		float miniIso = v.Electron().miniPFRelIso_all;
		if      (!v.Electron().mvaFall17V2noIso_WPL)  return 0;
		else if (miniIso>=0.40)                                 return 1;
		else if (!v.Electron().mvaFall17V2noIso_WP90) return 2;
		else if (miniIso>=0.20)                                 return 3;
		else if (miniIso>=0.10)                                 return 4;
		else if (!v.Electron().mvaFall17V2noIso_WP80) return 5;
		else if (miniIso>=0.05)                                 return 6;
		else                                                    return 7;
		}, .axis_title="MVA ID (w/ Mini iso.)", .def_range={0,8, 633}});
sh.AddNewFillParams("eleid_mva_2diso",         { .nbin=11, .bins={1,12}, .fill=[this] {
		if      (!v.Electron().mvaFall17V2noIso_WPL)                      return 0;
		else if (v.Electron().jetDRmin<0.4&&v.Electron().cleanJetPtrel<5)  return 1;
		else if (v.Electron().jetDRmin<0.4&&v.Electron().cleanJetPtrel<10) return 2;
		else if (!v.Electron().mvaFall17V2noIso_WP90)                     return 3;
		else if (v.Electron().jetDRmin<0.4&&v.Electron().cleanJetPtrel<15) return 4;
		else if (v.Electron().jetDRmin<0.4&&v.Electron().cleanJetPtrel<20) return 5;
		else if (v.Electron().jetDRmin<0.4&&v.Electron().cleanJetPtrel<25) return 6;
		else if (!v.Electron().mvaFall17V2noIso_WP80)                     return 7;
		else if (v.Electron().jetDRmin<0.4&&v.Electron().cleanJetPtrel<30) return 8;
		else if (v.Electron().jetDRmin<0.4&&v.Electron().cleanJetPtrel<40) return 9;
		else if (v.Electron().jetDRmin<0.4&&v.Electron().cleanJetPtrel<50) return 10;
		else                                                                        return 11;
		}, .axis_title="MVA ID (w/ 2D iso.)", .def_range={0,12, 921}}); // Gray
sh.AddNewFillParams("eleid_mva_miniiso_loose", { .nbin=6, .bins={2,8}, .fill=[this] {
		float miniIso = v.Electron().miniPFRelIso_all;
		if (!v.Electron().mvaFall17V2noIso_WPL)       return 0;
		else if (miniIso>=0.40)                                 return 1;
		else if (!v.Electron().mvaFall17V2noIso_WP90) return 2;
		else if (miniIso>=0.20)                                 return 3;
		else if (miniIso>=0.10)                                 return 4;
		else if (!v.Electron().mvaFall17V2noIso_WP80) return 5;
		else if (miniIso>=0.05)                                 return 6;
		else                                                    return 7;
		}, .axis_title="MVA ID (w/ Mini iso.)", .def_range={1.9,8, 633}});
sh.AddNewFillParams("muid_soft_comp",          { .nbin=1, .bins={1,2}, .fill=[this] { return v.Muon().softId; }, .axis_title="Soft (Cut-based) ID", .def_range={0,  2, 434}});
sh.AddNewFillParams("muid_soft",               { .nbin=1, .bins={1,2}, .fill=[this] { return v.Muon().softId; }, .axis_title="Soft (Cut-based) ID", .def_range={0.9,2, 434}});
sh.AddNewFillParams("muid_medium_comp",        { .nbin=2, .bins={1,3}, .fill=[this] { 
		if      (!v.Muon().mediumId)       return 0;
		else if (!v.Muon().mediumPromptId) return 1;
		else                                     return 2;
		}, .axis_title="Medium (Cut-based) ID", .def_range={0,3, 601}});
sh.AddNewFillParams("muid_medium",             { .nbin=2, .bins={1,3}, .fill=[this] { 
		if      (!v.Muon().mediumId)       return 0;
		else if (!v.Muon().mediumPromptId) return 1;
		else                                     return 2;
		}, .axis_title="Medium (Cut-based) ID", .def_range={0.9,3, 601}});
sh.AddNewFillParams("muid_mediumprompt",       { .nbin=2, .bins={1,3}, .fill=[this] { 
		if      (!v.Muon().mediumId)       return 0;
		else if (!v.Muon().mediumPromptId) return 1;
		else                                     return 2;
		}, .axis_title="Medium (Cut-based) ID", .def_range={2,3, 601}});
sh.AddNewFillParams("muid_tight_comp",         { .nbin=1, .bins={1,2}, .fill=[this] { return v.Muon().tightId; }, .axis_title="Tight (Cut-based) ID", .def_range={0,  2, 618}});
sh.AddNewFillParams("muid_tight",              { .nbin=1, .bins={1,2}, .fill=[this] { return v.Muon().tightId; }, .axis_title="Tight (Cut-based) ID", .def_range={0.9,2, 618}});
sh.AddNewFillParams("muid_mva",                { .nbin=3, .bins={1,4}, .fill=[this] { return v.Muon().mvaId; }, .axis_title="MVA ID", .def_range={0,  4, 633}});
sh.AddNewFillParams("muid_mva_loose",          { .nbin=3, .bins={1,4}, .fill=[this] { return v.Muon().mvaId; }, .axis_title="MVA ID", .def_range={0.9,4, 633}});
sh.AddNewFillParams("muid_mva_medium",         { .nbin=3, .bins={1,4}, .fill=[this] { return v.Muon().mvaId; }, .axis_title="MVA ID", .def_range={2,  4, 633}});
sh.AddNewFillParams("muid_mva_tight",          { .nbin=3, .bins={1,4}, .fill=[this] { return v.Muon().mvaId; }, .axis_title="MVA ID", .def_range={3,  4, 633}});


if (debug) std::cout<<"PlottingBase::define_histo_settings: ROC curve fillparams ok"<<std::endl;
if (debug) std::cout<<"PlottingBase::define_histo_settings: all fillparams ok"<<std::endl;

// Define histo types (for different object to loop on, and different cuts to apply)
// A new type should be defined for all different objects, depending on what we loop on
sh.AddHistoType("ele",                 "Electrons");
sh.AddHistoType("ele select",          "Electrons");
sh.AddHistoType("ele veto",            "Veto electrons");
sh.AddHistoType("ele tight noiso",     "Electrons");
sh.AddHistoType("syst ele",            "Electrons");
sh.AddHistoType("syst ele select",     "Electrons");
sh.AddHistoType("syst ele veto",       "Veto electrons");
sh.AddHistoType("mu",                  "Muons");
sh.AddHistoType("mu select",           "Muons");
sh.AddHistoType("mu veto",             "Veto muons");
sh.AddHistoType("mu tight noiso",      "Muons");
sh.AddHistoType("syst mu",             "Muons");
sh.AddHistoType("syst mu select",      "Muons");
sh.AddHistoType("syst mu veto",        "Veto muons");
sh.AddHistoType("syst mu tight noiso", "Muons");
sh.AddHistoType("tau",            		 "taus");
sh.AddHistoType("syst tau",       		 "taus");
sh.AddHistoType("pho",                 "Photons");
sh.AddHistoType("prepho",              "Pre-selected photons");
sh.AddHistoType("syst pho",            "Photons");
sh.AddHistoType("syst prepho",         "Pre-selected photons");
sh.AddHistoType("AK4",                 "Jets");
sh.AddHistoType("syst AK4",            "Jets");
sh.AddHistoType("b",                   "b-tagged jets");
sh.AddHistoType("b loose",             "Loose b-tagged jets");
sh.AddHistoType("syst b",              "b-tagged jets");
sh.AddHistoType("syst b loose",        "Loose b-tagged jets");
sh.AddHistoType("megajet",             "Megajets");
sh.AddHistoType("syst megajet",        "Megajets");
sh.AddHistoType("AK8",                 "AK8 jets");
sh.AddHistoType("AK8Mass",                 "AK8Mass jets");
sh.AddHistoType("syst AK8",            "AK8 jets");
sh.AddHistoType("hadtop",              "Tops (had.)");
sh.AddHistoType("leptop",              "Tops (lep.)");
sh.AddHistoType("leptop cand",         "Top candidates (lep.)");
sh.AddHistoType("hadw",                "Ws (had.)");
sh.AddHistoType("hadz",                "Zs (had.)");
sh.AddHistoType("hadh",                "Hs (H#rightarrowb#bar{b})");
sh.AddHistoType("syst hadtop",         "Tops (had.)");
sh.AddHistoType("syst leptop",         "Tops (lep.)");
sh.AddHistoType("syst leptop cand",    "Top candidates (lep.)");
sh.AddHistoType("syst hadw",           "Ws (had.)");
sh.AddHistoType("syst hadz",           "Zs (had.)");
sh.AddHistoType("syst hadh",           "Hs (H#rightarrowb#bar{b})");
sh.AddHistoType("gen lep",             "Gen-leptons");
sh.AddHistoType("gen ele",             "Gen-electrons");
sh.AddHistoType("gen mu",              "Gen-muons");
sh.AddHistoType("gen hadW",            "Gen-Ws (had.)");
sh.AddHistoType("gen hadZ",            "Gen-Zs (had.)");
sh.AddHistoType("gen hadV",            "Gen-Vs (had.)");
sh.AddHistoType("gen hadH",            "Gen-Hs (had.)");
sh.AddHistoType("gen top",             "Gen-tops");
sh.AddHistoType("gen hadtop",          "Gen-tops (had.)");
sh.AddHistoType("gen leptop",          "Gen-tops (lep.)");
sh.AddHistoType("evt",                 "Events / bin");
sh.AddHistoType("syst evt",            "Events / bin");

}

//_______________________________________________________
//              Define Histograms here
void
PlottingBase::define_analysis_histos(const Weighting& w, const unsigned int& syst_nSyst, const unsigned int& syst_index) {

	//__________________________________
	//        Define Smarthistos
	//__________________________________

	// Usage of AddHistos:
	// sh.AddHistos(type, { .fill=fill parameter(s), .pfs={list of posfix names}, .cuts={list of cuts},.draw="ROOT Draw() options", .opt="SmartHistos options", .ranges={plot ranges and overridden legend location}});
	// 1st param: string, the object/event type defined by the first parameter in AddHistoType()
	// 2nd param (.fill): string, either the name of a FillParam using the first parameter of AddNeFillParams() for a 1D histo, or two fillparams separated by "_vs_" for a 2D one
	// 3rd param (.pfs):  vector of strings for each Postfix name defined by the first parameter in AddNewPostfix(), supports up to 5 Postfix
	//                    N.B: for 1D plots, the first postfix is used to draw different plots on the same canvas, while creating multiple plots with the rest
	// 4th param (.cuts): vector of cuts before filling (I do not encourage the usage, because usually postfixes work nicer by attaching a string to the plot name)
	// 5th param (.draw): The option to use when using the TH1D/TH2D.Draw() method
	// 6th param (.opt):  SmartHistos specific options, some common examples:
	//   "Stack"+"N"   : Make a stack histo for all postfix, except the first N (can be data, signals, etc.)
	//   "Log"         : Log Y or Z scale
	//   "LogX"        : Log X scale
	//   "TwoCol"+"XY" : Split legend to two columns, and specify how many rows should each hav maximum
	//   "AddInt"      : Add integral of the histogram to the legend
	//   "AddRatio"    : Create a ratio of either 1st pf / 2nd pf, or 1st and pf/stacked histo
	//   There are many more, please peek into source code for complete list or see other examples
	// 7th .ranges: List of axis ranges in order: xlow,xhigh, ylow,yhigh, then  zlow,zhigh for 2D plot 
	//              OR the location of the top left edge of the legend (x1,y2), the legend size is auto adjusted 
	//              (by setting x2 and y1 depending on the number of entries)

	/*

Examples:

	// Single stack plot
	sh.AddHistos("evt", { .fill="MRR2", .pfs={"StackPlotTopSignal","Blind","SR_Had_1htop"}, .cuts={},.draw="",.opt="LogStack4AddRatioTwoCol48AddIntApproval45",.ranges={0,0, 1.01e-2,1e6, 0.3,0.86}});

	// Multiple stack plots in multiple regions (hadronic Higgs)
	for (const auto& enum_and_name : magic_enum::enum_entries<Region>()) {
	Region region = region.first;
	std::string cut(region.second);
	if (TString(cut).BeginsWith("SR_Had_")&&TString(cut).Constains("H")) {
	sh.SetHistoWeights({ [&w,region] { return w.sf_weight[region]; } });
	for (auto std_plot : {"HT", "MET", "MR", "R2", "MRR2"})
	sh.AddHistos(s+"evt", { .fill=c+std_plot, .pfs={"StackPlotHSignal","Blind",cut},  .cuts={},.draw=d,.opt=opt,.ranges=r_Stk6});
	}
	}

*/

	// Systematics:

	// The systematics of the plot are hidden in an additional dimension added to the histogram
	// eg.  a TH1D and it's systematic variations are stored in a TH2D
	// Each extra dimensional bin contain a variation of the underlying distribution
	// To enable this, a special dummy axis needs to be added which is called "Counts"
	// The string "Counts_vs_" is added to the plot name
	// and the "syst " is attached to the histogram type so the code knows to fill this
	// additional axis

	bool doSyst = (syst_nSyst>0);
	bool systematics = 1;

	std::string s = "";
	std::string c = "";
	if (systematics) {
		s = "syst ";
		c = "Counts_vs_";
	}

	std::vector<std::string> standard_plots;
	if (doSyst) standard_plots = {"HT", "MET", "HTMET"};
	else        standard_plots = {"HT", "HTFine", "METPhi", "METFine", "MET", "HTMET"};


	// ----------------------------------------------------------------------------------------------
	//                                    New preselection segions
	//-----------------------------------------------------------------------------------------------

	if (doSyst) {
		standard_plots = {"HT", "METPhi", "MET"};
	} else {
		standard_plots = {"HT", "METPhi", "METFine", "MET"};
	}
	//standard_plots.push_back("HTFine");
	//standard_plots.push_back("HTMET");

	for (auto region : {Region::Pre_1Lep, Region::Pre_1Lep_MT, Region::Pre_2Lep}) {
		sh.SetHistoWeights({ [&w,region] { return w.sf_weight[region]; } });
		std::string cut(magic_enum::enum_name(region));
		// Stack plots
		std::string opt  = o_stk_d_S;
		sh.AddHistos(  "evt",     { .fill="NJet",                    .pfs={"StackPlotSignal","Year",cut}, .cuts={},.draw=d,.opt=opt,.ranges=r_Stk9});
		sh.AddHistos(  "evt",     { .fill="NBTag",                   .pfs={"StackPlotSignal","Year",cut}, .cuts={},.draw=d,.opt=opt,.ranges=r_Stk9});
		sh.AddHistos(  "evt",     { .fill="MET",                     .pfs={"StackPlotSignal","Year",cut}, .cuts={},.draw=d,.opt=opt,.ranges=r_Stk9});
		sh.AddHistos(  "evt",     { .fill="NTauSelect",                .pfs={"StackPlotSignal","Year",cut}, .cuts={},.draw=d,.opt=opt,.ranges=r_Stk9});
		sh.AddHistos(  "evt",     { .fill="MTSelect",                .pfs={"StackPlotSignal","Year",cut}, .cuts={},.draw=d,.opt=opt,.ranges=r_Stk9});
		sh.AddHistos(  "evt",     { .fill="MinDeltaPhi",             .pfs={"StackPlotSignal","Year",cut}, .cuts={},.draw=d,.opt=opt,.ranges=r_Stk9});
		sh.AddHistos(  "evt",     { .fill="ElePt",                   .pfs={"StackPlotSignal","Year",cut}, .cuts={},.draw=d,.opt=opt,.ranges=r_Stk9});
		sh.AddHistos(  "evt",     { .fill="MuPt",                    .pfs={"StackPlotSignal","Year",cut}, .cuts={},.draw=d,.opt=opt,.ranges=r_Stk9});
		sh.AddHistos(  "AK4",     { .fill="JetPtBins",               .pfs={"StackPlotSignal","Year",cut}, .cuts={},.draw=d,.opt=opt,.ranges=r_Stk9});
		sh.AddHistos(  "AK4",     { .fill="JetPt",                   .pfs={"StackPlotSignal","Year",cut}, .cuts={},.draw=d,.opt=opt,.ranges=r_Stk9});
		sh.AddHistos(  "AK4",     { .fill="JetEta",                  .pfs={"StackPlotSignal","Year",cut}, .cuts={},.draw=d,.opt=opt,.ranges=r_Stk9});
		sh.AddHistos(  "AK4",     { .fill="JetPhi",                  .pfs={"StackPlotSignal","Year",cut}, .cuts={},.draw=d,.opt=opt,.ranges=r_Stk9});
	}
}

	void
PlottingBase::fill_analysis_histos(EventSelections& evt_sel, const Weighting& w, const unsigned int& syst_index, const double& weight)
{
	// SmartHistos
	if (syst_index == 0) { // Default (no systematic variation)
		//__________________________________ 
		//         Fill Smarthistos
		//__________________________________
		while(v.Electron.Loop())                                           sh.Fill("ele");
		while(v.Electron.Loop()) if (v.Electron.Select.pass[v.Electron.i]) sh.Fill("ele select");
		while(v.Electron.Loop()) if (v.Electron.Veto.  pass[v.Electron.i]) sh.Fill("ele veto");
		while(v.Electron.Loop()) if (v.Electron.NoIso. pass[v.Electron.i]) sh.Fill("ele tight noiso");
		while(v.Muon.Loop())                                       sh.Fill("mu");
		while(v.Muon.Loop())     if (v.Muon.Select.pass[v.Muon.i]) sh.Fill("mu select");
		while(v.Muon.Loop())     if (v.Muon.Veto.  pass[v.Muon.i]) sh.Fill("mu veto");
		while(v.Muon.Loop())     if (v.Muon.NoIso. pass[v.Muon.i]) sh.Fill("mu tight noiso");
		while(v.Tau.Loop())      if (v.Tau.Select.   pass[v.Tau.i])  sh.Fill("tau");
		while(v.Photon.Loop())   if (v.Photon.Select.   pass[v.Photon.i]) sh.Fill("pho");
		while(v.Jet.Loop())      if (v.Jet.Jet.       pass[v.Jet.i]) sh.Fill("AK4");
		while(v.Jet.Loop())      if (v.Jet.MediumBTag.pass[v.Jet.i]) sh.Fill("b");
		while(v.Jet.Loop())      if (v.Jet.LooseBTag. pass[v.Jet.i]) sh.Fill("b loose");
		if (!v.isData) {
			while(v.GenPart.Loop())  if (v.GenPart.Lepton.pass[v.GenPart.i])  sh.Fill("gen lep");
			while(v.GenPart.Loop())  if (v.GenPart.Ele.   pass[v.GenPart.i])  sh.Fill("gen ele");
			while(v.GenPart.Loop())  if (v.GenPart.Mu.    pass[v.GenPart.i])  sh.Fill("gen mu");
		}
		sh.Fill("evt");

	}

	// Do the same for systematics plots:
	while(v.Electron.Loop())                                           sh.Fill("syst ele");
	while(v.Electron.Loop()) if (v.Electron.Select.pass[v.Electron.i]) sh.Fill("syst ele select");
	while(v.Electron.Loop()) if (v.Electron.Veto.  pass[v.Electron.i]) sh.Fill("syst ele veto");
	while(v.Electron.Loop()) if (v.Electron.NoIso. pass[v.Electron.i]) sh.Fill("syst ele tight noiso");
	while(v.Muon.Loop())                                       sh.Fill("syst mu");
	while(v.Muon.Loop())     if (v.Muon.Select.pass[v.Muon.i]) sh.Fill("syst mu select");
	while(v.Muon.Loop())     if (v.Muon.Veto.  pass[v.Muon.i]) sh.Fill("syst mu veto");
	while(v.Muon.Loop())     if (v.Muon.NoIso. pass[v.Muon.i]) sh.Fill("syst mu tight noiso");
	while(v.Tau.Loop())      if (v.Tau.Select.   pass[v.Tau.i])  sh.Fill("syst tau");
	while(v.Photon.Loop())   if (v.Photon.Select.   pass[v.Photon.i]) sh.Fill("syst pho");
	while(v.Jet.Loop())      if (v.Jet.Jet.       pass[v.Jet.i]) sh.Fill("syst AK4");
	while(v.Jet.Loop())      if (v.Jet.MediumBTag.pass[v.Jet.i]) sh.Fill("syst b");
	while(v.Jet.Loop())      if (v.Jet.LooseBTag. pass[v.Jet.i]) sh.Fill("syst b loose");
	sh.Fill("syst evt");
}


	void
PlottingBase::load_analysis_histos(std::string inputfile)
{
	sh.Add(inputfile.c_str());
}

	void
PlottingBase::save_analysis_histos(bool draw=0)
{
	if (draw) sh.DrawPlots();
	else sh.Write();
}

#endif // End header guard
