#ifndef EVENTSELECTIONS_H
#define EVENTSELECTIONS_H

// Private headers
#include "Variables.h"
#include "Weighting.h"

// 3rd party headers
#include "tnm.h"
#include "magic_enum.h"
#include "TString.h"
#include "TRandom3.h" // For HEM failure in MC

// common libraries
#include <iostream>
#include <functional>
#include <map>
#include <vector>

class EventSelections
{
public:
  EventSelections(Variables& var) :
    v(var) {}
  ~EventSelections() {}

  typedef struct Cut { std::string name; std::function<bool()> func; } Cut;
  std::vector<Cut> baseline_cuts;

  void define_preselections();
  
  void define_event_selections();

  enum Regions : size_t {
    // Preection regions
    Pre_1Lep,
    Pre_Had,
    Pre_2Lep,
		Pre_1Lep_MT,
      };
  typedef Regions Region;

  void define_region(Region, Region, std::vector<Cut>);

  std::vector<std::vector<Cut> > analysis_cuts;
  std::vector<bool> pass_all_cuts;
  std::vector<unsigned int> cutbits;

  bool apply_cut(Region, std::string);
  bool apply_cut(Region, unsigned int);
  bool apply_ncut(Region, std::string);
  bool apply_ncut(Region, unsigned int);
  bool apply_cuts(Region, std::vector<std::string>);
  bool apply_cuts(Region, std::vector<unsigned int>);
  bool apply_all_cuts(Region);
  bool apply_all_cuts_except(Region, std::string);
  bool apply_all_cuts_except(Region, unsigned int);
  bool apply_all_cuts_except(Region, std::vector<std::string>);
  bool apply_all_cuts_except(Region, std::vector<unsigned int>);

  bool pass_skimming();

  void apply_event_selections();

private:

  Variables& v;  
  TRandom3 rnd_;

};




//_______________________________________________________
//  Apply analysis cuts in the specified search region

bool
EventSelections::apply_all_cuts(Region region) {
  return apply_ncut(region, analysis_cuts[region].size());
}

bool
EventSelections::apply_ncut(Region region, unsigned int ncut) {
  if (ncut>analysis_cuts[region].size()) return 0;
  for (unsigned int i=0; i<ncut; ++i) if ( ! analysis_cuts[region][i].func() ) return 0;
  return 1;
}

// Cuts to apply/exclude by cut name
bool
EventSelections::apply_cut(Region region, std::string cut_name) {
  for (const auto& cut : analysis_cuts[region]) if (cut_name == cut.name) return cut.func();
  return 0;
}

bool
EventSelections::apply_cuts(Region region, std::vector<std::string> cuts) {
  for (const auto& cut_in_region : analysis_cuts[region]) for (const auto& cut : cuts)
    if (cut == cut_in_region.name) if (!cut_in_region.func()) return 0;
  return 1;
}

bool
EventSelections::apply_all_cuts_except(Region region, std::string cut_to_skip) {
  bool result = true, found = false;
  for (const auto& cut : analysis_cuts[region]) {
    if (cut.name == cut_to_skip) {
      found = true;
      continue;
    }
    if (!cut.func()) result = false;
  }
  // If a certain cut meant to be skipped (N-1) is not found for some reason
  // eg. mistyped, then end the job with ar error
  // This is for safety: We do not want to fill histograms wrongly by mistake
  if (!found) {
    std::cout<<"No cut to be skipped exsists in search region \""<<magic_enum::enum_name(region)<<"\" with name: \""<<cut_to_skip<<"\""<<std::endl;
    error("EventSelections - the second argument for apply_all_cuts_except() is a non-sensical cut");
  }
  return result;
}

bool
EventSelections::apply_all_cuts_except(Region region, std::vector<std::string> cuts_to_skip) {
  bool result = true;
  unsigned int found = 0;
  for (const auto& cut : analysis_cuts[region]) {
    for (const auto& cut_to_skip : cuts_to_skip) if (cut.name==cut_to_skip) {
      ++found;
      continue;
    }
    if (!cut.func()) result = false;
  }
  // If a certain cut meant to be skipped is not found for some reason
  // eg. mistyped, then end the job with ar error
  // This is for safety: We do not want to fill histograms wrongly by mistake
  if (found!=cuts_to_skip.size()) {
    std::cout<<"A cut to be skipped does not exsist in seaerch region \""<<magic_enum::enum_name(region)<<"\" with names: ";
    for (const auto& cut : cuts_to_skip) std::cout<<cut<<", "; std::cout<<std::endl;
    error("EventSelections - the second argument for apply_all_cuts_except() contains at least one non-sensical cut");
  }
  return result;
}


// Same functions but with cut index which is faster (can use an enum, to make it nicer)
bool
EventSelections::apply_cut(Region region, unsigned int cut_index) { return analysis_cuts[region][cut_index].func(); }

bool
EventSelections::apply_cuts(Region region, std::vector<unsigned int> cuts) {
  for (const unsigned int& cut : cuts) if ( ! analysis_cuts[region][cut].func() ) return 0;
  return 1;
}

bool
EventSelections::apply_all_cuts_except(Region region, unsigned int cut_to_skip) {
  if (cut_to_skip>=analysis_cuts[region].size()) {
    std::cout<<"Index ("<<cut_to_skip<<") is too high for the cut to be skipped in search region '"<<magic_enum::enum_name(region)<<"'"<<std::endl;
    error("EventSelections::apply_all_cuts_except(Region region, unsigned int cut_to_skip)");
  }
  for (unsigned int i=0, n=analysis_cuts[region].size(); i<n; ++i) {
    if (i==cut_to_skip) continue;
    if ( ! analysis_cuts[region][i].func() ) return 0;
  }
  return 1;
}

bool
EventSelections::apply_all_cuts_except(Region region, std::vector<unsigned int> cuts_to_skip) {
  for (unsigned int i=0, n=analysis_cuts[region].size(); i<n; ++i) {
    for (const unsigned int& cut_to_skip : cuts_to_skip) if (i!=cut_to_skip)
      if ( ! analysis_cuts[region][i].func() ) return 0;
  }
  return 1;
}

void
EventSelections::apply_event_selections()
{
  // Calculate decision of each individual cut
  for (size_t i=0, n=analysis_cuts.size(); i<n; ++i) {
    cutbits[i] = 0;
    for (size_t j=0, m=analysis_cuts[i].size(); j<m; ++j)
      if (analysis_cuts[i][j].func()) cutbits[i] += 1<<j;
    pass_all_cuts[i] = (cutbits[i]==(unsigned int)((1<<analysis_cuts[i].size())-1));
  }
}
//_______________________________________________________
//                Define Skimming cuts
//      Needed, if you want to skim the ntuple

bool
EventSelections::pass_skimming()
{
  return 1;
}



//_______________________________________________________
//              Define event cleaning cuts
void
EventSelections::define_preselections()
{
  baseline_cuts.clear();

  // Apply the same cuts as it is in the ntuple - Only for check
  // cut is an std::function, which we can define easily with a lambda function

  // Recommended event filters by MET group
  // https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2?rev=137
  baseline_cuts.push_back({ .name="Clean_goodVertices",        .func = [this] { return v.Flag_goodVertices; } });
  baseline_cuts.push_back({ .name="Clean_CSC_Halo_SuperTight", .func = [this] { return v.isSignal ? 1 : v.Flag_globalSuperTightHalo2016Filter; } });
  baseline_cuts.push_back({ .name="Clean_HBHE_Noise",          .func = [this] { return v.Flag_HBHENoiseFilter; } });
  baseline_cuts.push_back({ .name="Clean_HBHE_IsoNoise",       .func = [this] { return v.Flag_HBHENoiseIsoFilter; } });
  baseline_cuts.push_back({ .name="Clean_Ecal_Dead_Cell_TP",   .func = [this] { return v.Flag_EcalDeadCellTriggerPrimitiveFilter; } });
  baseline_cuts.push_back({ .name="Clean_Bad_PF_Muon",         .func = [this] { return v.Flag_BadPFMuonFilter; } });
  baseline_cuts.push_back({ .name="Clean_Bad_PF_Muon_Dz",      .func = [this] { return v.Flag_BadPFMuonDzFilter; } });
  baseline_cuts.push_back({ .name="Clean_hf_NoisyHits",        .func = [this] { return v.Flag_hfNoisyHitsFilter; } });
  //baseline_cuts.push_back({ .name="Clean_Bad_Charged",         .func = [this] { return v.Flag_BadChargedCandidateFilter; } });
  baseline_cuts.push_back({ .name="Clean_EE_Bad_Sc",           .func = [this] { return v.isData ? v.Flag_eeBadScFilter : 1; } });
  if (v.year>2016) {
    baseline_cuts.push_back({ .name="Clean_Ecal_Bad_Calib",    .func = [this] { return v.Flag_ecalBadCalibFilter; } });
  }
  // HEM 15/16 failure
  if (v.year==2018) {
    baseline_cuts.push_back({ .name="Clean_HEM_failure", .func = [this] {
                                if (v.isData ? v.run>=319077 : rnd_.Rndm()<0.645) {
																	while (v.Jet.Loop()){
                                  	if (v.Jet().pt  >= 30 &&
                                      //v.Jet().eta > -4.7 && v.Jet().eta < -1.4 &&
                                      //v.Jet().phi > -1.6 && v.Jet().phi < -0.8)
                                      v.Jet().eta > -3.2 && v.Jet().eta < -1.2 &&
                                      v.Jet().phi > -1.77 && v.Jet().phi < -0.67)
                                    return 0;
																	}
																	while (v.Electron.Loop()){
                                  	if (v.Electron().pt  >= 30 &&
                                      v.Electron().eta > -3.0 && v.Electron().eta < -1.4 &&
                                      v.Electron().phi > -1.57 && v.Electron().phi < -0.87)
                                    return 0;
																	}
																}
                                return 1;
                              } });
  }
  //https://indico.cern.ch/event/884106/contributions/3725350/attachments/1977822/3292861/MET_news_slides.pdf#page=3&zoom=auto,-58,405
  baseline_cuts.push_back({ .name="Clean_MET_Extra", .func = [this] { 
                              if (v.isData) {
                                if (v.run==321149) {
                                  if (v.event==91061433   ||
                                      v.event==2202873820 ||
                                      v.event==2202827292 ||
                                      v.event==2354264119 ||
                                      v.event==2354264306) return 0;
                                } else if (v.run==321730) {
                                  if (v.event==457120628 ||
                                      v.event==457027731 ||
                                      v.event==457120560 ||
                                      v.event==51878652 ||
                                      v.event==51878651 ||
                                      v.event==51878615) return 0;
                                }
                              }
                              return 1;
                            } });
  baseline_cuts.push_back({ .name="Clean_Bad_Muon_Jet", .func = [this] { return v.dPhiMuonJetMET<2.74159; } });
  baseline_cuts.push_back({ .name="Clean_Bad_PFMET",    .func = [this] { 
                              double ratio = v.CaloMET_pt>0 ? v.MET_pt/v.CaloMET_pt : 9999;
                              if (v.year==2016) {
                                return (ratio>= 0.5 && ratio<5.0);
                              } else if (v.year==2017) {
                                return (ratio>= 0.5 && ratio<3.0);
                              } else if (v.year==2018) {
                                return (ratio>= 0.5 && ratio<3.0);
                              } else return bool(1);
                            } });
}

//_______________________________________________________
//          Define Analysis event selection cuts
//     Can define all sorts of Signal/Control regions

void
EventSelections::define_region(Region region, Region presel, std::vector<Cut> cuts) {
  for (const auto& presel_cut : analysis_cuts[presel])
    analysis_cuts[region].push_back(presel_cut);
  for (const auto& cut : cuts)
    analysis_cuts[region].push_back(cut);
}

void
EventSelections::define_event_selections()
{
  analysis_cuts.clear();
  analysis_cuts.resize(magic_enum::enum_count<Region>());
  pass_all_cuts.resize(magic_enum::enum_count<Region>());
  cutbits      .resize(magic_enum::enum_count<Region>());

  // Define here cuts that are include in all Signal/Control regions
  // MET Filters, etc. are already applied with baseline_cuts

  // If you change the trigger selection
  // Make sure to change also to match trigger efficiency selection in Plottingbase
  // then when new efficiency is calculated, load the new trigger efficiency in Weighting
  std::function<bool()> hadronic_triggers;
  if (v.sample.Contains("HTMHT")) {
    hadronic_triggers = [this] { 
      if (v.year==2016) return v.HLT_PFHT300_PFMET110==1;
      else return (bool)0;
    };
  } else if (v.sample.Contains("MET")) {
    hadronic_triggers = [this] { 
      if (v.year==2016) return !(v.HLT_PFHT300_PFMET110==1) && // veto HTMHT
        (v.HLT_PFMET110_PFMHT110_IDTight==1||v.HLT_PFMETNoMu110_PFMHTNoMu110_IDTight==1);
      else return 
        (v.HLT_PFMET120_PFMHT120_IDTight==1 ||
         v.HLT_PFMET120_PFMHT120_IDTight_PFHT60==1 ||
         v.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight==1 ||
         v.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60==1 ||
         v.HLT_PFHT500_PFMET100_PFMHT100_IDTight==1 || 
         v.HLT_PFHT700_PFMET85_PFMHT85_IDTight==1 ||
         v.HLT_PFHT800_PFMET75_PFMHT75_IDTight==1); };
  } else if (v.sample.Contains("JetHT")) {
    hadronic_triggers = [this] { 
      if (v.year==2016) return 
        !(v.HLT_PFHT300_PFMET110==1|| // veto HTMHT
          v.HLT_PFMET110_PFMHT110_IDTight==1||v.HLT_PFMETNoMu110_PFMHTNoMu110_IDTight==1) && // veto MET
        (v.HLT_AK8PFJet450==1||v.HLT_PFHT800==1||v.HLT_PFHT900==1);
      else return !(v.HLT_PFMET120_PFMHT120_IDTight==1 || // veto MET
                    v.HLT_PFMET120_PFMHT120_IDTight_PFHT60==1 ||
                    v.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight==1 ||
                    v.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60==1 ||
                    v.HLT_PFHT500_PFMET100_PFMHT100_IDTight==1 || 
                    v.HLT_PFHT700_PFMET85_PFMHT85_IDTight==1 ||
                    v.HLT_PFHT800_PFMET75_PFMHT75_IDTight==1) &&
        v.HLT_PFHT1050==1;
        //(v.HLT_PFHT1050==1||v.HLT_AK8PFJet500==1); // Add this when new efficiencies ready
      
    };
  } else {
    // Data histos should not contain events from other datasets
    hadronic_triggers = [this] { return !v.isData; };
		//Use trigger selection to Data
/*
    hadronic_triggers = [this] { 
      if (v.year==2016) return 
        (v.HLT_PFHT300_PFMET110==1|| // veto HTMHT
          v.HLT_PFMET110_PFMHT110_IDTight==1||v.HLT_PFMETNoMu110_PFMHTNoMu110_IDTight==1) ||
        (v.HLT_AK8PFJet450==1||v.HLT_PFHT800==1||v.HLT_PFHT900==1);
      else return (v.HLT_PFMET120_PFMHT120_IDTight==1 || // veto MET
                    v.HLT_PFMET120_PFMHT120_IDTight_PFHT60==1 ||
                    v.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight==1 ||
                    v.HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60==1 ||
                    v.HLT_PFHT500_PFMET100_PFMHT100_IDTight==1 || 
                    v.HLT_PFHT700_PFMET85_PFMHT85_IDTight==1 ||
                    v.HLT_PFHT800_PFMET75_PFMHT75_IDTight==1) ||
        v.HLT_PFHT1050==1;
    };
*/
  }

  std::function<bool()> leptonic_triggers;
  if (v.sample.Contains("SingleElectron")||v.sample.Contains("EGamma")) {
    leptonic_triggers = [this] {
      if (v.year==2016) return
        v.HLT_Ele27_WPTight_Gsf==1 ||
      else return
        v.HLT_Ele32_WPTight_Gsf==1 ||
    };
  } else if (v.sample.Contains("SingleMuon")) {
    leptonic_triggers = [this] {
      // Veto events already collected by Single Electron trigger
      if (v.year==2016) {
        return
        v.HLT_IsoMu24==1 ||
        v.HLT_IsoTkMu24==1 ||
      } else {
        return
        v.HLT_IsoMu27==1 ||
        v.HLT_IsoTkMu27==1 ||
      }
    };
  } else {
    // Data histos should not contain events from other datasets
    leptonic_triggers = [this] { return !v.isData; };
  }

  std::function<bool()> photonic_triggers;
  if (v.sample.Contains("SinglePhoton")||v.sample.Contains("EGamma")) {
    photonic_triggers = [this] {
      if (v.year==2016) return v.HLT_Photon175==1;
      else return v.HLT_Photon200==1;
    };
  } else {
    // Data histos should not contain events from other datasets
    photonic_triggers = [this] { return !v.isData; };
  }
  
  v.get_signal_mass();

  // Preselection leptonic regions
  analysis_cuts[Region::Pre_1Lep] = {
    { .name="NJetPre",    .func = [this] { return v.Jet.Jet.n>=4;              }},
    { .name="1Lep",       .func = [this] { return v.nLepSelect==1;             }},
    { .name="HLT",        .func =                leptonic_triggers              },
    { .name="1b",         .func = [this] { return v.Jet.MediumBTag.n>=1;            }},
  };
  analysis_cuts[Region::Pre_2Lep] = {
    { .name="NJetPre",    .func = [this] { return v.Jet.Jet.n>=4;              }},
    { .name="1Lep",       .func = [this] { return v.nLepSelect==2;             }},
    { .name="HLT",        .func =                leptonic_triggers              },
    { .name="1b",         .func = [this] { return v.Jet.MediumBTag.n>=1;            }},
  };
  // Preselection hadronic regions
  analysis_cuts[Region::Pre_Had] = {
    { .name="NJetPre",    .func = [this] { return v.Jet.Jet.n>=6;              }},
    { .name="1Lep",       .func = [this] { return v.nLepVeto==0;             }},
    { .name="HLT",        .func =                hadronic_triggers              },
    { .name="1b",         .func = [this] { return v.Jet.MediumBTag.n>=1;            }},
  };
  // Preselection leptonic regions
  analysis_cuts[Region::Pre_1Lep_MT] = {
    { .name="NJetPre",    .func = [this] { return v.Jet.Jet.n>=4;              }},
    { .name="1Lep",       .func = [this] { return v.nLepSelect==1;             }},
    { .name="HLT",        .func =                leptonic_triggers              },
    { .name="1b",         .func = [this] { return v.Jet.MediumBTag.n>=1;            }},
    { .name="MT",         .func = [this] { return v.MT>=80;                   }},
  };

}


#endif // End header guard
