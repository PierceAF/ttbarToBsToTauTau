#ifndef VARIABLES_H
#define VARIABLES_H

// Private headers
#include "eventBuffer.h" // make sure to set to same as in tnm.h
#include "Razor.h"
#include "MT2.h"
#include "XYMETCorrection.h"

// 3rd party headers
#include "tnm.h"
//#include "TLorentzVector.h" // should be replaced to ROOT::Math::LorentzVector
#include "Math/LorentzVector.h"
#include "TRandom3.h" // For selecting which lepton to add to MET in DiLep CR
#include "TString.h"

// common libraries
#include <iostream>
#include <functional>
#include <algorithm>

#define NOVAL_F 9999.0
#define NOVAL_I 9999

//_______________________________________________________
//                 Object selections

/*
  Latest Electron IDs:
  [1] Cut Based  - https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2?rev=41#Working_points_for_2016_data_for
  [2] MVA        - https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentificationRun2?rev=30#Recommended_MVA_recipes_for_2016
  [3] SUSY (Use) - https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF?rev=172#Electrons

  Latest Isolation WPs:
  [4] SUSY MiniIso Loose/Tight - https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF?rev=172#ID_IP_ISO_AN1
  [5] RazorInclusive - https://github.com/RazorCMS/RazorAnalyzer/blob/master/src/RazorAnalyzer.cc#L1210

  Latest Impact Point Cut:
  [6] SUSY Loose/Tight IP2D (Use) - https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF?rev=172#ID_IP_ISO_AN1
  [7] POG  Tight - https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2?rev=41#Offline_selection_criteria

  For Veto (Regions except Z) Choose:
  - Spring15 Cut based Veto ID without relIso (EA) cut
  - for pt<20  absolute iso03 (EA) < 5
  - Mini-Isolation (EA)/pt < 0.1 (Medium WP [4])
  - pt >= 5
  - |eta| < 2.5
  - |d0| < 0.2, |dz| < 0.5 (Loose IP2D [6])
  OR when using MVA:
  - Spring16 MVA Loose ID
  - for pt>=20 Mini-Isolation (EA)/pt < 0.2 (RazorInclusive cut [5])
  - pt >= 5
  - |eta| < 2.5
  - 3D IP sig < 4

  For Selection (Z) Choose:
  - Spring15 Cut based Medium ID without relIso (EA) cut
  - Mini-Isolation (EA)/pt < 0.1 (Tight WP [4])
  - pt >= 10
  - |eta| < 2.5, also exclude barrel-endcap gap [1.442,1556]
  - |d0| < 0.05, |dz| < 0.1 (Tight IP2D [6])
  OR when using MVA:
  - Spring16 MVA Loose ID
  - Mini-Isolation (EA)/pt < 0.1 (Tight WP [4])
  - pt >= 10
  - |eta| < 2.5, also exclude barrel-endcap gap [1.442,1556]
  - |d0| < 0.05, |dz| < 0.1 (Tight IP2D [6])


  For Tight Selection (TriggerEff Only) Choose:
  - Spring15 Cut based Tight ID (including relIso (EA) cut)
  - pt >= 30
  - |eta| < 2.5, also exclude barrel-endcap gap [1.442,1556]
  - |d0| < 0.05, |dz| < 0.1 (Tight IP2D and IP3D [6])

  For Loose selection (deprecated) Choose:
  - Spring15 Cut based Loose ID without relIso (EA) cut
  - Mini-Isolation (EA)/pt < 0.1 (Tight WP [4])
  - pt >= 10
  - |eta| < 2.5, also exclude barrel-endcap gap [1.442,1556]
  - |d0| < 0.2, |dz| < 0.5 (Loose IP2D [6])

*/

#define ELE_VETO_PT_CUT        5
#define ELE_VETO_ETA_CUT       2.5
#define ELE_VETO_MINIISO_CUT   0.4
#define ELE_VETO_IP_D0_CUT     0.05
#define ELE_VETO_IP_DZ_CUT     0.1
#define ELE_VETO_IP_3D_CUT     4   // For skim only
#define ELE_VETO_ABSISO_CUT    5   // For skim only

#define ELE_SELECT_PT_CUT      10
#define ELE_SELECT_ETA_CUT     2.5
#define ELE_SELECT_MINIISO_CUT 0.1
#define ELE_SELECT_IP_D0_CUT   0.05
#define ELE_SELECT_IP_DZ_CUT   0.1

#define ELE_TIGHT_PT_CUT       30
#define ELE_TIGHT_ETA_CUT      2.5
#define ELE_TIGHT_IP_D0_CUT    0.05
#define ELE_TIGHT_IP_DZ_CUT    0.1
//#define ELE_TIGHT_IP_SIG_CUT   4

/*
  Latest Muon IDs (Loose/Medium):
  [1] POG Loose/Medium - https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2?rev=28#Short_Term_Instructions_for_Mori

  Latest Isolation WPs:
  [2] SUSY MiniISo Loose/Tight - https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF?rev=172#ID_IP_ISO
  [3] POG Tight RelIso - https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2?rev=28#Muon_Isolation
  [4] RazorInclusive - https://github.com/RazorCMS/RazorAnalyzer/blob/master/src/RazorAnalyzer.cc#L1602

  Latest Impact Point Cut (Loose/Tight):
  [5] SUSY Loose/Tight IP2D - https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF?rev=172#ID_IP_ISO

  For Veto (All except Z) Choose:
  - POG recommended Loose ID (No Iso/IP)
  - for pt>=20 Mini-Isolation (EA)/pt < 0.2 (RazorInclusive cut [4])
  - for pt<20  absolute iso04 (EA) < 10
  - pt >= 5
  - |eta| < 2.4
  - 3D IP sig < 4

  For Loose Choose:
  - POG recommended Loose ID (No Iso/IP)
  - Mini-Isolation (EA)/pt < 0.2 (Tight WP [2])
  - pt >= 10
  - |eta| < 2.4
  - |d0| < 0.2, |dz| < 0.5 (Loose IP2D [5])

  For Selection (Z) Choose:
  - POG recommended Medium ID (No Iso/IP)
  - Mini-Isolation (EA)/pt < 0.2 (tight WP [2])
  - pt >= 5
  - |eta| < 2.4
  - |d0| < 0.05, |dz| < 0.1 (Tight IP2D [5])

  For Tight Selection (TriggerEff Only) Choose:
  - POG recommended Tight ID (No Iso/IP)
  - comb. rel. Isolation (R=0.4) < 0.15 (tight WP [3])
  - pt >= 30
  - |eta| < 2.4
  - |d0| < 0.05, |dz| < 0.1 (Tight IP2D and IP3D [5])

*/

#define MU_VETO_PT_CUT         5
#define MU_VETO_ETA_CUT        2.4
#define MU_VETO_MINIISO_CUT    0.4
#define MU_VETO_IP_D0_CUT      0.2
#define MU_VETO_IP_DZ_CUT      0.5
#define MU_VETO_ABSISO_CUT     10  // For skim only
#define MU_VETO_IP_3D_CUT      4   // For skim only

#define MU_SELECT_PT_CUT       10
#define MU_SELECT_ETA_CUT      2.4
#define MU_SELECT_MINIISO_CUT  0.2
#define MU_SELECT_IP_D0_CUT    0.05
#define MU_SELECT_IP_DZ_CUT    0.1

#define MU_TIGHT_PT_CUT        30
#define MU_TIGHT_ETA_CUT       2.4
#define MU_TIGHT_RELISO_CUT    0.15
#define MU_TIGHT_IP_D0_CUT     0.05
#define MU_TIGHT_IP_DZ_CUT     0.1
#define MU_TIGHT_IP_SIG_CUT    4

/*
  Latest Tau IDs:
  https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV?rev=35

  For Veto Choose:
  - POG Loose ID:
  + byLooseCombinedIsolationDeltaBetaCorr3Hits"
  - pt >= 18 (Same as MINIAOD)

  OR use:
  - Isolated charged trk = 0
  + pt>=5 (ele/mu), 10 (hadrons)
  + isolation (dR=0.3) / pt < 0.2 (ele/mu), 0.1 (hadrons)
  + |dz| < 0.1

*/

#define  TAU_VETO_PT_CUT  20
#define  TAU_VETO_ETA_CUT 2.3


/*
  Latest Photon IDs:
  https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2?rev=36

  For Selection Choose:
  - Spring16 Cut based Medium photon ID
  - Pass electron veto
  - pt >= 80
  - |eta| < 2.5

*/

#define PHOTON_SELECT_PT_CUT        80
#define PHOTON_SELECT_ETA_CUT       2.5

/*
  Jet ID (Oct31/Jan12 ntuple):
  https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2017

  For AK4 Jet Selection Choose:
  - Tight jet ID
  - pt > 30
  - |eta| < 2.4 //Maybe we should change the higher

  For AK8 Jet Selection Choose:
  - Tight jet ID
  - pt > 200 (this cut was lowered to 170 for skimming)
  - |eta| < 2.4 //Maybe we should change the higher

*/

#define JET_AK4_PT_CUT  30 //30 is original
#define JET_AK4_ETA_CUT 2.4
#define JET_AK8_PT_CUT  200
#define JET_AK8_ETA_CUT 2.4


/*
  Latest b-tagging WPs/SFs:
  https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X

  Choose:
  - CombinedSecondaryVertex v2
  - DeepCSV >= 0.1522 (Loose - for Veto)
  - DeepCSV >= 0.4941 (Medium - for Tag)

*/

#define B_SUBJET_DEPP_LOOSE_CUT 0.1522
#define B_DEEP_LOOSE_CUT        0.1522
#define B_DEEP_MEDIUM_CUT       0.4941
#define B_DEEP_TIGHT_CUT        0.8001

//_______________________________________________________
//                    Event Variables

class EventData : public eventBuffer::Event_s {
public:
  EventData(eventBuffer& d) : d_(d) {}
  
private:
  eventBuffer& d_;
  
public:
  void initEventData() {
    eventBuffer::Event_s& base = *this;
    base = std::move(d_);
    resetEventData();
  }

  // Selected object counts
  // ----------- All Leptons ----------
  size_t nLepVetoNoIso;
  size_t nLepVeto;
  size_t nLepSelect;
  size_t nLepTight;
  size_t nLepNoIso;
  size_t nLepNonIso;

  double lepNeutrinoDR;

  // -------------- Jets --------------
  size_t nJetISR;

  // --------------- MET --------------
  Vector3 MET;
  Vector3 MET_1vl;
  Vector3 MET_1l;
  Vector3 MET_2l;
  Vector3 MET_dilep;
  Vector3 MET_fakepho;
  Vector3 MET_pho;
  
  // ---------- Discriminators --------
  // Signal discriminators
  double AK4_Ht;
  double AK4_HtOnline;
  double AK8_Ht;
  double MT_lepVeto;
  double MT;
  double MT_lepTight;
  double MT_lepNonIso;
  double MT_boost;
  double M_2l;
  double dPhi_2l_met;
  double dPhi_2l_jet;
  double minDeltaPhi; // Min(DeltaPhi(Jet_i, MET)), i=1,2,3,4
  double minDeltaPhi_1vl;
  double minDeltaPhi_1l;
  double minDeltaPhi_2l;
  double minDeltaPhi_pho;

  // Other
  std::vector<LorentzVector> megajets;
  std::vector<LorentzVector> megajets_nophoton;
  size_t iMegaJet;
  std::vector<LorentzVector> hemJets;

  void resetEventData() {
    nLepVetoNoIso = 0;
    nLepVeto      = 0;
    nLepSelect    = 0;
    nLepTight     = 0;
    nLepNoIso     = 0;
    nLepNonIso    = 0;
    lepNeutrinoDR = NOVAL_F;
    
    nJetISR       = 0;

    MET           .SetXYZ(0,0,0);
    MET_1vl       .SetXYZ(0,0,0);
    MET_1l        .SetXYZ(0,0,0);
    MET_2l        .SetXYZ(0,0,0);
    MET_fakepho   .SetXYZ(0,0,0);
    MET_pho       .SetXYZ(0,0,0);
    MET_dilep     .SetXYZ(0,0,0);
    
    AK4_Ht       = 0;
    AK4_HtOnline = 0;
    AK8_Ht       = 0;
    MT_lepVeto    = NOVAL_F;
    MT            = NOVAL_F;
    MT_lepTight   = NOVAL_F;
    MT_lepNonIso  = NOVAL_F;
    MT_boost      = NOVAL_F;
    M_2l          = -NOVAL_F;
    dPhi_2l_met   = NOVAL_F;
    dPhi_2l_jet   = NOVAL_F;
    minDeltaPhi     = NOVAL_F;
    minDeltaPhi_1l  = NOVAL_F;
    minDeltaPhi_1vl = NOVAL_F;
    minDeltaPhi_2l  = NOVAL_F;
    minDeltaPhi_pho = NOVAL_F;

    megajets.clear();
    megajets_nophoton.clear();
    iMegaJet = -1;
    hemJets.clear();
  }

};


class Variables : public EventData {

public:

  Variables(eventBuffer& d, const int& y, const bool& data, const bool& signal, const std::string& dirname, const bool& isapv) :
    EventData(d),
    year(y),
    isData(data),
    isSignal(signal),
    isBackground(!data&&!signal),
    sample(dirname),
    Electron(d),
    Muon(d),
    Tau(d),
    Photon(d),
    Jet(d),
    GenJet(d),
    SubJet(d),
    FatJet(d),
    GenPart(d)
  { 
    resetEventData(); 
    isFastSim = sample.Contains("SMS");
    isQCD = sample.Contains("QCD_HT") || sample.Contains("QQ_HT-");
    isGJets = sample.BeginsWith("GJets");
    isWJets = sample.Contains("WJetsToLNu");
    isTop = sample.Contains("TTTo")||sample.Contains("ST");
    isZInv = sample.Contains("ZJetsToNuNu");
    //isAPV = sample.Contains("HIPM") || sample.Contains("APV");
    isAPV = isapv;
  }
  ~Variables() {}

  const int  year;
  const bool isData;
  const bool isSignal;
  const bool isBackground;
  const TString sample;
  bool isFastSim;
  bool isQCD;
  bool isGJets;
  bool isWJets;
  bool isTop;
  bool isZInv;
	bool isAPV;

  std::vector<double> susy_mass;
  int signal_index = -1;

  TRandom3 rnd_;

  double get_syst_weight(const double& weight_nominal, const double& weight_up, const double& weight_down, const double& nSigma) {
    if (nSigma == 0) {
      return weight_nominal;
    } else {
      double w = weight_nominal;
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

  double get_syst_weight(const double& weight_nominal, const double& uncertainty, const double& nSigma) {
    double w = weight_nominal;
    // Use symmetrical difference for up/down variation
    if (nSigma!=0.) w *= 1.0 + nSigma * uncertainty;
    return w;
  }
  

  //_______________________________________________________
  //                   Analysis Variables

  //_______________________________________________________
  //                    Object containers

  // Selected object container class
  // It also stores variables for a derived collection
  // which inherits its data structure from a base type
  template<typename BaseType, typename DerivedType>
  class ObjectVectorContainer {
    
  public:
    ObjectVectorContainer<BaseType, DerivedType>
    (std::vector<BaseType>& baseColl, 
     size_t default_buffer_size) : 
      baseColl_(baseColl),
      default_buffer_size_(default_buffer_size) {
      init(); 
    };
    ~ObjectVectorContainer<BaseType, DerivedType>() {};
    
  private:
    // reference to the base collection
    std::vector<BaseType>& baseColl_;
    size_t default_buffer_size_; // default buffer size from eventBuffer class

    // new data storage containers
    std::vector<DerivedType> derivedColl_;
    std::vector<LorentzVector> v4_;

  public:
    size_t i;
    size_t n;
    
    void init() {
      i = -1;
      n = 0;
    }

    // Reset only the newly declared variables
    // Lazy implementation --> use the move constructor (FAST) for the base class
    // AND the class member initialization for the rest (new feature of C++11)
    // This way, we don't have to implement and use an initializer at all :)
    // variables can be declared and initialized in one line in the class declaration
    // this way of resetting is benefitting from this automatic initialization
    // while not compromising speed much excessive copying
    // Please, keep the base collection intact (do not pass it's reference anywhere)
    // and make sure the object survives so this method can (re-)access it's contents
    void moveData() {
      init();
      derivedColl_.clear();
      v4_.clear();
      for (auto& obj : baseColl_) {
        derivedColl_.emplace_back(std::move(obj));
        LorentzVector v4;
        float pt = std::abs(obj.pt);
        float m = obj.mass;
        float px = pt*std::cos(obj.phi);
        float py = pt*std::sin(obj.phi);
        float pz = pt*sinh(obj.eta);
        float e = m>=0 ? std::sqrt(px*px+py*py+pz*pz+m*m) : std::sqrt(std::max((px*px+py*py+pz*pz-m*m), (float)0.));
        v4.SetPxPyPzE(px, py, pz, e);
        v4_.emplace_back(v4);
      }
      n=derivedColl_.size();
    }

    bool Loop() {
      if (++i<n) { 
        return 1; 
      } else {
        i=-1;
        return 0;
      }
    }

    DerivedType& operator() () {
      if (i!=(size_t)-1) return derivedColl_[i]; 
      else {
        std::cout<<"ERROR - Variables::ObjectVectorContainer::operator(): "
                 <<"method called outside while(ObjectVectorContainer.Loop())"<<std::endl;
        std::exit(1);
      }
    };
    LorentzVector&    v4   () { 
      if (i!=(size_t)-1) return v4_[i];
      else {
        std::cout<<"ERROR - Variables::ObjectVectorContainer::v4(): "
                 <<"method called outside while(ObjectVectorContainer.Loop())"<<std::endl;
        std::exit(1);
      }
    };
    DerivedType& operator() (size_t& iRef) { 
      if (iRef<derivedColl_.size()) return derivedColl_[iRef];
      else {
        std::cout<<"ERROR - Variables::ObjectVectorContainer::operator(size_t): Object index out of range "
                 <<"iRef="<<iRef<<" >= "<<derivedColl_.size()<<"(collection size)"<<std::endl;
        std::exit(1);
      }
    };
    LorentzVector&    v4   (size_t& iRef) { 
      if (iRef<v4_.size()) return v4_[iRef];
      else {
        std::cout<<"ERROR - Variables::ObjectVectorContainer::v4(size_t): Object index out of range "
                 <<"iRef="<<iRef<<" >= "<<v4_.size()<<"(collection size)"<<std::endl;
        std::exit(1);
      }
    };
    DerivedType& operator() (int     iRef) { 
      if (iRef<derivedColl_.size()&&iRef>=0) return derivedColl_[iRef];
      else {
        if (iRef<0) std::cout<<"ERROR - Variables::ObjectVectorContainer::operator(int): Object index negative "
                             <<"iRef="<<iRef<<". Please use a non-negative number!"<<std::endl;

        else std::cout<<"ERROR - Variables::ObjectVectorContainer::operator(int): Object index out of range "
                      <<"iRef="<<iRef<<" >= "<<derivedColl_.size()<<"(collection size)"<<std::endl;
        std::exit(1);
      }
    };
    LorentzVector&    v4   (int     iRef) { 
      if (iRef<v4_.size()&&iRef>=0) return v4_[iRef];
      else {
        if (iRef<0) std::cout<<"ERROR - Variables::ObjectVectorContainer::v4(int): Object index negative "
                             <<"iRef="<<iRef<<". Please use a non-negative number!"<<std::endl;
        
        else std::cout<<"ERROR - Variables::ObjectVectorContainer::v4(int): Object index out of range "
                      <<"iRef="<<iRef<<" >= "<<v4_.size()<<"(collection size)"<<std::endl;
        std::exit(1);
      }
    };
    size_t       get_default_buffer_size() { return default_buffer_size_; }

  };
  
  // Class to store multiple selected objects
  template<typename BaseType, typename DerivedType> class Object {
  public:
    Object(ObjectVectorContainer<BaseType, DerivedType>* x) : 
      parent_(x) {
      init(x->get_default_buffer_size());
    };
    Object(ObjectVectorContainer<BaseType, DerivedType>* x, std::function<bool()> s) : 
      parent_(x), selector_(s) {
      init(x->get_default_buffer_size());
    };
    ~Object() {};

    std::vector<bool> pass;            // tells if the object passes the selection
    std::vector<size_t> viSel;         // indices of the selected objects
    std::vector<size_t> viRef;         // indices of the reference objects (in the whole collection)

    size_t  i;         // index of the selected object
    size_t  iRef;      // index of the reference object
    size_t  n;         // number of selected objects

    // Initialize buffers
    void init(size_t nreserve) {
      pass .reserve(nreserve);
      viSel.reserve(nreserve);
      viRef.reserve(nreserve);
      reset();
    }
    void reset() {
      pass.assign(parent_->n, 0);
      viSel.assign(parent_->n, (size_t)-1);
      viRef.clear();
      i = -1;
      n = 0;
    }

    bool auto_define() { return define(selector_()); }
    
    bool define(bool passObj) {
      if ((pass[parent_->i] = passObj)) {
        viSel[parent_->i] = n++;
        viRef.emplace_back(parent_->i);
      }
      return passObj;
    }

    bool Loop() {
      if (++i<n) { 
        iRef    = viRef[i];
        return 1;
      } else {
        i=-1;
        return 0;
      }
    }

    const DerivedType& operator() ()             { return parent_->operator()(viRef[i]);    };
    const LorentzVector&    v4   ()             { return parent_->v4(viRef[i]);            };
    const DerivedType& operator() (size_t& iRef) { return parent_->operator()(viRef[iRef]); };
    const LorentzVector&    v4   (size_t& iRef) { return parent_->v4(viRef[iRef]);         };
    const DerivedType& operator() (int     iRef) { return parent_->operator()(viRef[iRef]); };
    const LorentzVector&    v4   (int     iRef) { return parent_->v4(viRef[iRef]);         };

  private:
    ObjectVectorContainer<BaseType, DerivedType>* parent_;
    std::function<bool()> selector_;
  };

  //_______________________________________________________
  //                      Electrons

  struct ElectronData : eventBuffer::Electron_s {
    ElectronData(eventBuffer::Electron_s&& arg) : 
      eventBuffer::Electron_s(std::move(arg)) {}
  
    // Matched AK4 jet
    size_t iMatchedAK4    = -1;
    double jetDR          =  NOVAL_F;
    double jetDPhi        =  NOVAL_F;
    double jetDRmin       =  NOVAL_F;
    double cleanJetPtrel  =  NOVAL_F;
    bool   isPartOfJet    =  0;
    // Associated neutrino
    LorentzVector nu     {0,0,0,0};
    double nuDR           =  NOVAL_F;
    // Matched AK8 jet
    bool   matchAK8       =  0;
    size_t iMatchedAK8    = -1;
    double ak8JetDR       =  NOVAL_F;
    double ak8JetDPhi     =  NOVAL_F;
    double nuMatchedAK8DR =  NOVAL_F;
    // GenTruth
    bool matchGenEle                = 0;
    bool matchGenEleFromHardProcess = 0;
    bool matchGenEleFromW           = 0;
    bool matchGenEleFromZ           = 0;
    bool matchGenEleFromH           = 0;
    bool matchGenEleFromTop         = 0;
    bool nuMatchGenEleNu            = 0;
    bool nuMatchGenEleNuFromW       = 0;
    bool nuMatchGenEleNuFromTop     = 0;
  };

  class ElectronSelection : 
    public ObjectVectorContainer<eventBuffer::Electron_s, ElectronData> {
    
  public:
    ElectronSelection(eventBuffer& d) :
      ObjectVectorContainer<eventBuffer::Electron_s, ElectronData>
    (d.Electron, d.Electron_pt.size()) { init(); };
    ~ElectronSelection() {};
    typedef Object<eventBuffer::Electron_s, ElectronData> Electron_c;
    
    Electron_c CBVeto      {this};
    Electron_c CBVetoNoIso {this};
    Electron_c Veto      {this};
    Electron_c VetoNoIso {this};
    Electron_c Select    {this};
    Electron_c Tight     {this};
    Electron_c NoIso     {this};
    Electron_c NonIso    {this};
    
    void initObjects() {
      moveData();
      CBVeto     .reset();
      CBVetoNoIso.reset();
      Veto     .reset();
      VetoNoIso.reset();
      Select   .reset();
      Tight    .reset();
      NoIso    .reset();
      NonIso   .reset();
    }

  } Electron;

  //_______________________________________________________
  //                        Muons

  struct MuonData : eventBuffer::Muon_s {
    MuonData(eventBuffer::Muon_s&& arg) : 
      eventBuffer::Muon_s(std::move(arg)) {}
  
    // Matched AK4 jet
    size_t iMatchedAK4    = -1;
    double jetDR          =  NOVAL_F;
    double jetDPhi        =  NOVAL_F;
    double jetDRmin       =  NOVAL_F;
    double cleanJetPtrel  =  NOVAL_F;
    bool   isPartOfJet    =  0;
    // Associated neutrino
    LorentzVector nu     {0,0,0,0};
    double nuDR           =  NOVAL_F;
    // Matched AK8 jet
    bool   matchAK8       =  0;
    size_t iMatchedAK8    = -1;
    double ak8JetDR       =  NOVAL_F;
    double ak8JetDPhi     =  NOVAL_F;
    double nuMatchedAK8DR =  NOVAL_F;
    // GenTruth
    bool matchGenMu                = 0;
    bool matchGenMuFromHardProcess = 0;
    bool matchGenMuFromW           = 0;
    bool matchGenMuFromZ           = 0;
    bool matchGenMuFromH           = 0;
    bool matchGenMuFromTop         = 0;
    bool nuMatchGenMuNu            = 0;
    bool nuMatchGenMuNuFromW       = 0;
    bool nuMatchGenMuNuFromTop     = 0;
  };

  class MuonSelection : public ObjectVectorContainer
  <eventBuffer::Muon_s, MuonData> {
    
  public:
    MuonSelection(eventBuffer& d) :
      ObjectVectorContainer<eventBuffer::Muon_s, MuonData>
    (d.Muon, d.Muon_pt.size()) { init(); };
    ~MuonSelection() {};
    typedef Object<eventBuffer::Muon_s, MuonData> Muon_c;
    
    Muon_c CBLoose   {this};
    Muon_c CBLooseNoIso   {this};
    Muon_c CBMedium  {this};
    Muon_c CBMediumNoIso  {this};
    Muon_c Veto      {this};
    Muon_c VetoNoIso {this};
    Muon_c Select    {this};
    Muon_c Tight     {this};
    Muon_c NoIso     {this};
    Muon_c NonIso    {this};

    void initObjects() {
      moveData();
      CBLoose  .reset();
      CBMedium .reset();
      Veto     .reset();
      VetoNoIso.reset();
      Select   .reset();
      Tight    .reset();
      NoIso    .reset();
      NonIso   .reset();
    }
  } Muon;

  //_______________________________________________________
  //                        Taus

  struct TauData : eventBuffer::Tau_s {
    TauData(eventBuffer::Tau_s&& arg) : 
      eventBuffer::Tau_s(std::move(arg)) {}
    // GenTruth
    bool matchGenTau = 0;
  };

  class TauSelection : public ObjectVectorContainer
  <eventBuffer::Tau_s, TauData> {
    
  public:
    TauSelection(eventBuffer& d) :
      ObjectVectorContainer<eventBuffer::Tau_s, TauData>
    (d.Tau, d.Tau_pt.size()) { init(); };
    ~TauSelection() {};
    typedef Object<eventBuffer::Tau_s, TauData> Tau_c;
    
    Tau_c Veto{this};
    
    void initObjects() {
      moveData();
      Veto     .reset();
    }
  } Tau;

  //_______________________________________________________
  //                       Photons

  struct PhotonData : eventBuffer::Photon_s {
    PhotonData(eventBuffer::Photon_s&& arg) : 
      eventBuffer::Photon_s(std::move(arg)) {}
    // GenTruth
    bool fromFrag             = 0;
    bool matchGenFake         = 1;
    bool matchGenPrompt       = 0;
    bool matchGenPromptDirect = 0;
    bool matchGenPromptFrag   = 0;
  };

  class PhotonSelection : public ObjectVectorContainer
  <eventBuffer::Photon_s, PhotonData> {
    
  public:
    PhotonSelection(eventBuffer& d) :
      ObjectVectorContainer<eventBuffer::Photon_s, PhotonData>
    (d.Photon, d.Photon_pt.size()) { init(); };
    ~PhotonSelection() {};
    typedef Object<eventBuffer::Photon_s, PhotonData> Photon_c;
    
    Photon_c PreSelect  {this};
    Photon_c Select     {this};
    Photon_c Fake       {this};
    Photon_c SelectNoIso{this};

    void initObjects() {
      moveData();
      PreSelect  .reset();
      Select     .reset();
      Fake       .reset();
      SelectNoIso.reset();
    }
  } Photon;


  //_______________________________________________________
  //                        Jets 

  struct JetData : eventBuffer::Jet_s {
    JetData(eventBuffer::Jet_s&& arg) : 
      eventBuffer::Jet_s(std::move(arg)) {}
  
    double eleDR       =  NOVAL_F;
    double elePtRatio  =  0;
    double muDR        =  NOVAL_F;
    double muPtRatio   =  0;
    double phoDR       =  NOVAL_F;
    double phoPtRatio  =  0;
    bool matchJet =  0;
    // GenTruth
  };

  class JetSelection : public ObjectVectorContainer
  <eventBuffer::Jet_s, JetData> {
    
  public:
    JetSelection(eventBuffer& d) :
      ObjectVectorContainer<eventBuffer::Jet_s, JetData>
    (d.Jet, d.Jet_pt.size()) { init(); };
    ~JetSelection() {};
    typedef Object<eventBuffer::Jet_s, JetData> Jet_c;
    
    Jet_c Jet            {this};
    Jet_c FailID         {this};
    Jet_c MuonJet        {this};
    Jet_c JetLepNoIso    {this};
    Jet_c JetLepNonIso   {this};
    Jet_c JetNoLep       {this};
    Jet_c JetNoPho       {this};
    Jet_c LooseBTag      {this};
    Jet_c MediumBTag     {this};
    Jet_c TightBTag      {this};
    Jet_c MediumBTagNoPho{this};
    Jet_c LooseIsoBTag   {this};
    Jet_c MediumIsoBTag  {this};

    void initObjects() {
      moveData();
      Jet          .reset();
      FailID       .reset();
      MuonJet      .reset();
      JetLepNoIso  .reset();
      JetLepNonIso .reset();
      JetNoLep     .reset();
      JetNoPho     .reset();
      LooseBTag    .reset();
      MediumBTag   .reset();
      TightBTag    .reset();
      LooseIsoBTag .reset();
      MediumIsoBTag.reset();
    }
  } Jet;

  //_______________________________________________________
  //                     GenJets 

  struct GenJetData : eventBuffer::GenJet_s {
    GenJetData(eventBuffer::GenJet_s&& arg) : 
      eventBuffer::GenJet_s(std::move(arg)) {}
  
  };

  class GenJetSelection : public ObjectVectorContainer
  <eventBuffer::GenJet_s, GenJetData> {
    
  public:
    GenJetSelection(eventBuffer& d) :
      ObjectVectorContainer<eventBuffer::GenJet_s, GenJetData>
    (d.GenJet, d.GenJet_pt.size()) { init(); };
    ~GenJetSelection() {};
    typedef Object<eventBuffer::GenJet_s, GenJetData> GenJet_c;
    
    void initObjects() {
      moveData();
    }
  } GenJet;

  //_______________________________________________________
  //                  FatJets and SubJets

  struct SubJetData : eventBuffer::SubJet_s {
    SubJetData(eventBuffer::SubJet_s&& arg) : 
      eventBuffer::SubJet_s(std::move(arg)) {}
    // GenTruth
  };

  class SubJetSelection : public ObjectVectorContainer
  <eventBuffer::SubJet_s, SubJetData> {
    
  public:
    SubJetSelection(eventBuffer& d) :
      ObjectVectorContainer<eventBuffer::SubJet_s, SubJetData>
    (d.SubJet, d.SubJet_pt.size()) { init(); };
    ~SubJetSelection() {};
    typedef Object<eventBuffer::SubJet_s, SubJetData> SubJet_c;
    
    void initObjects() {
      moveData();
    }
  } SubJet;

  struct FatJetData : eventBuffer::FatJet_s {
    FatJetData(eventBuffer::FatJet_s&& arg) : 
      eventBuffer::FatJet_s(std::move(arg)) {}
  
    bool   passSubJetBTag =  0;
    size_t nSubJet        =  0;
    size_t nSubJetBTag    =  0;
    double maxSubJetDeepB = -NOVAL_F;
    double tau21          =  NOVAL_F;
    double tau31          =  NOVAL_F;
    double tau32          =  NOVAL_F;
    double eleDR          =  NOVAL_F;
    double elePtRatio     =  0;
    double muDR           =  NOVAL_F;
    double muPtRatio      =  0;
    double phoDR          =  NOVAL_F;
    double phoPtRatio     =  0;
    bool   matchLepNoIso  =  0;
    bool   matchLepNonIso =  0;
    double lepNonIsoNuDR  =  NOVAL_F;
    double LSF            = -NOVAL_F;
    double LSF_NoIso      = -NOVAL_F;
    double matchedNoIsoLepJetDRmin      = -NOVAL_F;
    double matchedNoIsoLepCleanJetPtrel = -NOVAL_F;
    // GenTruth
    bool matchGenHadW     =  0;
    bool matchGenHadZ     =  0;
    bool matchGenHadH     =  0;
    bool matchGenHadTop   =  0;
    bool matchGenLepTop   =  0;
    bool matchGenLepton   =  0;
    int  matchedGenLeptonMotherID = -NOVAL_F;
  };

  class FatJetSelection : public ObjectVectorContainer
  <eventBuffer::FatJet_s, FatJetData> {
    
  public:
    FatJetSelection(eventBuffer& d) :
      ObjectVectorContainer<eventBuffer::FatJet_s, FatJetData>
    (d.FatJet, d.FatJet_pt.size()) { init(); };
    ~FatJetSelection() {};
    typedef Object<eventBuffer::FatJet_s, FatJetData> FatJet_c;
    
    FatJet_c JetAK8               {this};
    FatJet_c JetAK8Mass           {this};
    // All old and new WPs
    FatJet_c W1                   {this};
    FatJet_c W2                   {this};
    FatJet_c W3                   {this};
    FatJet_c Top1                 {this};
    FatJet_c Top2                 {this};
    FatJet_c Top3                 {this};
    FatJet_c Top4                 {this};
    FatJet_c Top5                 {this};
    FatJet_c WDeepMD1             {this};
    FatJet_c WDeepMD2             {this};
    FatJet_c WDeepMD3             {this};
    FatJet_c WDeepMD4             {this};
    FatJet_c WDeep1               {this};
    FatJet_c WDeep2               {this};
    FatJet_c WDeep3               {this};
    FatJet_c WDeep4               {this};
    FatJet_c ZDeepMD1             {this};
    FatJet_c ZDeepMD2             {this};
    FatJet_c ZDeepMD3             {this};
    FatJet_c ZDeep1               {this};
    FatJet_c ZDeep2               {this};
    FatJet_c ZDeep3               {this};
    FatJet_c VDeep1               {this};
    FatJet_c VDeep2               {this};
    FatJet_c VDeep3               {this};
    FatJet_c HDeepMD1             {this};
    FatJet_c HDeepMD2             {this};
    FatJet_c HDeepMD3             {this};
    FatJet_c HDeep1               {this};
    FatJet_c HDeep2               {this};
    FatJet_c HDeep3               {this};
    FatJet_c HDeep4               {this};
    FatJet_c TopDeepMD1           {this};
    FatJet_c TopDeepMD2           {this};
    FatJet_c TopDeepMD3           {this};
    FatJet_c TopDeepMD4           {this};
    FatJet_c TopDeep1             {this};
    FatJet_c TopDeep2             {this};
    FatJet_c TopDeep3             {this};
    FatJet_c TopDeep4             {this};
    FatJet_c WparticleNet1        {this};
    FatJet_c WparticleNet2        {this};
    FatJet_c WparticleNet3        {this};
    FatJet_c ZparticleNet1        {this};
    FatJet_c ZparticleNet2        {this};
    FatJet_c ZparticleNet3        {this};
    FatJet_c VparticleNet1        {this};
    FatJet_c VparticleNet2        {this};
    FatJet_c VparticleNet3        {this};
    FatJet_c HparticleNet1        {this};
    FatJet_c HparticleNet2        {this};
    FatJet_c HparticleNet3        {this};
    FatJet_c HparticleNet4        {this};
    FatJet_c TopparticleNet1      {this};
    FatJet_c TopparticleNet2      {this};
    FatJet_c TopparticleNet3      {this};
    // New definitions
    FatJet_c HadW                 {this};
    FatJet_c HadZ                 {this};
    FatJet_c HadV                 {this};
    FatJet_c HadH                 {this};
    FatJet_c HadTop               {this};
    FatJet_c LepJetCand           {this};
    FatJet_c HadW_MD              {this};
    FatJet_c HadZ_MD              {this};
    FatJet_c HadTop_MD            {this};
    FatJet_c LepJet               {this};
    FatJet_c LepJetNoIso          {this};
    FatJet_c LepJetNoPt           {this};
    FatJet_c LepTop               {this};
    FatJet_c LepJetHighMass       {this};
    FatJet_c LepTopHighMass       {this};
    FatJet_c LepTopNoIso          {this};
    FatJet_c LepTopNoPt           {this};
    FatJet_c LepTopNoMass         {this};
    FatJet_c LepTopNoSubJetB      {this};

    void initObjects() {
      moveData();
      JetAK8              .reset();
      JetAK8Mass          .reset();
      W1                  .reset();
      W2                  .reset();
      W3                  .reset();
      Top1                .reset();
      Top2                .reset();
      Top3                .reset();
      Top4                .reset();
      Top5                .reset();
      WDeepMD1            .reset();
      WDeepMD2            .reset();
      WDeepMD3            .reset();
      WDeepMD4            .reset();
      WDeep1              .reset();
      WDeep2              .reset();
      WDeep3              .reset();
      WDeep4              .reset();
      ZDeepMD1            .reset();
      ZDeepMD2            .reset();
      ZDeepMD3            .reset();
      ZDeep1              .reset();
      ZDeep2              .reset();
      ZDeep3              .reset();
      VDeep1              .reset();
      VDeep2              .reset();
      VDeep3              .reset();
      HDeepMD1            .reset();
      HDeepMD2            .reset();
      HDeepMD3            .reset();
      HDeep1              .reset();
      HDeep2              .reset();
      HDeep3              .reset();
      HDeep4              .reset();
      TopDeepMD1          .reset();
      TopDeepMD2          .reset();
      TopDeepMD3          .reset();
      TopDeepMD4          .reset();
      TopDeep1            .reset();
      TopDeep2            .reset();
      TopDeep3            .reset();
      TopDeep4            .reset();
      WparticleNet1       .reset();
      WparticleNet2       .reset();
      WparticleNet3       .reset();
      ZparticleNet1       .reset();
      ZparticleNet2       .reset();
      ZparticleNet3       .reset();
      VparticleNet1       .reset();
      VparticleNet2       .reset();
      VparticleNet3       .reset();
      HparticleNet1       .reset();
      HparticleNet2       .reset();
      HparticleNet3       .reset();
      HparticleNet4       .reset();
      TopparticleNet1     .reset();
      TopparticleNet2     .reset();
      TopparticleNet3     .reset();
      HadW                .reset();
      HadZ                .reset();
      HadV                .reset();
      HadH                .reset();
      HadTop              .reset();
      LepJetCand          .reset();
      HadW_MD             .reset();
      HadZ_MD             .reset();
      HadTop_MD           .reset();
      LepJet              .reset();
      LepJetNoIso         .reset();
      LepJetNoPt          .reset();
      LepTop              .reset();
      LepJetHighMass      .reset();
      LepTopHighMass      .reset();
      LepTopNoIso         .reset();
      LepTopNoPt          .reset();
      LepTopNoMass        .reset();
      LepTopNoSubJetB     .reset();
    }
  } FatJet;
  
  //_______________________________________________________
  //                     Gen Particles

  struct GenPartData : eventBuffer::GenPart_s {
    GenPartData(eventBuffer::GenPart_s&& arg) : 
      eventBuffer::GenPart_s(std::move(arg)) {}
  
    bool LastCopyCand   = 0;
    bool NoSameDaughter = 1;
    bool FinalState     = 0;
    // For Object efficiency - N.B: Not for Mistag!
    bool passEleCBVetoNoIso = 0;
    bool passEleCBVeto     = 0;
    bool passMuoCBLooseNoIso = 0;
    bool passMuoCBLoose    = 0;
    bool passMuoCBMediumNoIso = 0;
    bool passMuoCBMedium   = 0;
    bool passLepVeto       = 0;
    bool passLepNoIso      = 0;
    bool passLepNonIso     = 0;
    bool passHadWTag       = 0;
    bool passHadZTag       = 0;
    bool passHadHTag       = 0;
    bool passHadTopTag     = 0;
    bool passLepTopTag     = 0;
    // Pointers to matched object indices
    bool   matchAK8                  =  0;
    size_t iMatchedAK8               = -1;
    size_t iMother                   = -1; // eg. (W->nu e)
    size_t iGrandMother              = -1; // eg. (t->b (W->nu mu))
    size_t iGreatGrandMother         = -1; // eg. (t->b (W->nu (tau->nu e)))
    size_t iGenLepDaughter           = -1;
    size_t iGenLepGrandDaughter      = -1;
    size_t iGenLepGreatGrandDaughter = -1;
    std::vector<size_t> iDaughters;
    //size_t iGenHadWMatchedAK8; --> replace with above
    //size_t iGenLepWMatchedGenLep;
    //size_t iGenLepTopMatchedGenLep;
    // Ancestor ID
    int motherId           = -NOVAL_I;
    int grandMotherId      = -NOVAL_I;
    int greatGrandMotherId = -NOVAL_I;
  };

  class GenPartSelection : public ObjectVectorContainer
  <eventBuffer::GenPart_s, GenPartData> {
    
  public:
    GenPartSelection(eventBuffer& d) :
      ObjectVectorContainer<eventBuffer::GenPart_s, GenPartData>
    (d.GenPart, d.GenPart_pt.size()) { init(); };
    ~GenPartSelection() {};
    typedef Object<eventBuffer::GenPart_s, GenPartData> GenPart_c;
    
    GenPart_c Lepton                {this};
    GenPart_c Ele                   {this};
    GenPart_c Mu                    {this};
    GenPart_c Tau                   {this};
    GenPart_c LeptonFromHardProcess {this};
    GenPart_c EleFromHardProcess    {this};
    GenPart_c MuFromHardProcess     {this};
    GenPart_c TauFromHardProcess    {this};
    GenPart_c LeptonFromW           {this};
    GenPart_c EleFromW              {this};
    GenPart_c MuFromW               {this};
    GenPart_c TauFromW              {this};
    GenPart_c LeptonFromZ           {this};
    GenPart_c EleFromZ              {this};
    GenPart_c MuFromZ               {this};
    GenPart_c TauFromZ              {this};
    GenPart_c LeptonFromH           {this};
    GenPart_c EleFromH              {this};
    GenPart_c MuFromH               {this};
    GenPart_c TauFromH              {this};
    GenPart_c LeptonFromTop         {this};
    GenPart_c EleFromTop            {this};
    GenPart_c MuFromTop             {this};
    GenPart_c TauFromTop            {this};
    GenPart_c Nu                    {this};
    GenPart_c EleNu                 {this};
    GenPart_c MuNu                  {this};
    GenPart_c TauNu                 {this};
    GenPart_c NuFromW               {this};
    GenPart_c EleNuFromW            {this};
    GenPart_c MuNuFromW             {this};
    GenPart_c TauNuFromW            {this};
    GenPart_c NuFromTop             {this};
    GenPart_c EleNuFromTop          {this};
    GenPart_c MuNuFromTop           {this};
    GenPart_c TauNuFromTop          {this};
    GenPart_c PromptPhoton          {this};
    GenPart_c PromptDirectPhoton    {this};
    GenPart_c PromptFragPhoton      {this};
    GenPart_c b                     {this};
    GenPart_c W                     {this};
    GenPart_c LepW                  {this};
    GenPart_c HadW                  {this};
    GenPart_c Z                     {this};
    GenPart_c LepZ                  {this};
    GenPart_c HadZ                  {this};
    GenPart_c H                     {this};
    GenPart_c LepH                  {this};
    GenPart_c HadH                  {this};
    GenPart_c Top                   {this};
    GenPart_c LepTop                {this};
    GenPart_c HadTop                {this};
    
    void initObjects() {
      moveData();
      Lepton               .reset();
      Ele                  .reset();
      Mu                   .reset();
      Tau                  .reset();
      LeptonFromHardProcess.reset();
      EleFromHardProcess   .reset();
      MuFromHardProcess    .reset();
      TauFromHardProcess   .reset();
      LeptonFromW          .reset();
      EleFromW             .reset();
      MuFromW              .reset();
      TauFromW             .reset();
      LeptonFromZ          .reset();
      EleFromZ             .reset();
      MuFromZ              .reset();
      TauFromZ             .reset();
      LeptonFromH          .reset();
      EleFromH             .reset();
      MuFromH              .reset();
      TauFromH             .reset();
      LeptonFromTop        .reset();
      EleFromTop           .reset();
      MuFromTop            .reset();
      TauFromTop           .reset();
      Nu                   .reset();
      EleNu                .reset();
      MuNu                 .reset();
      TauNu                .reset();
      NuFromW              .reset();
      EleNuFromW           .reset();
      MuNuFromW            .reset();
      TauNuFromW           .reset();
      NuFromTop            .reset();
      EleNuFromTop         .reset();
      MuNuFromTop          .reset();
      TauNuFromTop         .reset();
      PromptPhoton         .reset();
      PromptDirectPhoton   .reset();
      PromptFragPhoton     .reset();
      b                    .reset();
      W                    .reset();
      LepW                 .reset();
      HadW                 .reset();
      Z                    .reset();
      LepZ                 .reset();
      HadZ                 .reset();
      H                    .reset();
      LepH                 .reset();
      HadH                 .reset();
      Top                  .reset();
      LepTop               .reset();
      HadTop               .reset();
    }
  } GenPart;

  int recalc_jets     = 0;
  int recalc_megajets = 0;
  int recalc_met      = 0;

  // New methods instead of calculate_common_variables
  void get_signal_mass() {
    susy_mass.assign(3, -NOVAL_F);
    if (isSignal) while (GenPart.Loop()) {
      if (signal_index==0) { if(std::abs(GenPart().pdgId) == 1000021) susy_mass[0] = GenPart().mass; } // gluino
      else if (signal_index==1||signal_index==4){ if(std::abs(GenPart().pdgId) == 1000006 || std::abs(GenPart().pdgId) == 1000005) susy_mass[0] = GenPart().mass; } // stop
      else if (signal_index==2||signal_index==3){ if(std::abs(GenPart().pdgId) == 1000024 || std::abs(GenPart().pdgId) == 1000023 || std::abs(GenPart().pdgId) == 1000025) susy_mass[0] = GenPart().mass; } // chargino
      else if (signal_index==6) { if(std::abs(GenPart().pdgId) == 1000021) susy_mass[0] = GenPart().mass; } // gluino

      if (signal_index==6)               { if(std::abs(GenPart().pdgId) == 1000006) susy_mass[1] = GenPart().mass; } // stop
      else                                    { if(std::abs(GenPart().pdgId) == 1000022) susy_mass[1] = GenPart().mass; } // LSP neutralino (chi^0_1)
    }
   	susy_mass[0] = int(std::round(susy_mass[0]/5)*5);
   	susy_mass[1] = int(std::round(susy_mass[1]/25)*25);
  }

  void define_lepton_and_photon_variables() {
    // Lepton selections do not depend on any systematics
    // --> event weights are assigned instead (for SF, lost lepton etc)
    define_leptons_and_photons_();
  }

  void define_jet_variables(const unsigned int& syst_index) {
    // Vairables need to be recalculated each time the jet energy is changed
    // eg. Jet selection, W/top tags etc. that depends on jet pt
    define_jets_(syst_index);
  }
  
  void define_genparticle_variables() {
    define_genparticles_();
  }
  
  void define_event_variables(const unsigned int& syst_index) {
    define_event_(syst_index);
  }

  //_______________________________________________________
  //              Rescale jet 4-momenta

  std::vector<LorentzVector> saved_Jet_v4;
  std::vector<float> saved_Jet_pt;
  std::vector<LorentzVector> saved_FatJet_v4;
  std::vector<float> saved_FatJet_pt;
  std::vector<float> saved_FatJet_msoftdrop;
  float saved_MET_pt;
  float saved_MET_phi;
  void rescale_smear_jet_met(const bool& applySmearing, const unsigned int& syst_index,
                             const double& nSigmaJES, const double& nSigmaJER)
  {
		if (syst_index!=0) {
      // Load
      while (Jet.Loop()) {
        Jet.v4() = saved_Jet_v4[Jet.i];
        Jet().pt = saved_Jet_pt[Jet.i];
      }
      while (FatJet.Loop()) {
        FatJet.v4() = saved_FatJet_v4[Jet.i];
        FatJet().pt = saved_FatJet_pt[FatJet.i];
        FatJet().msoftdrop = saved_FatJet_msoftdrop[FatJet.i];
      }
      MET_pt  = saved_MET_pt;
      MET_phi = saved_MET_phi;
    }
    // Aplly JES/JER
    // Replace pt with the nominal value after smearing
		//cout << syst_index << ": ";
		double temp=0;
    while (Jet.Loop()) {
			temp = Jet().pt_nom;
			Jet().pt = Jet().pt_nom;
      //double scaleJES = get_syst_weight(Jet().pt_nom, Jet().pt_jesTotalUp, Jet().pt_jesTotalDown, nSigmaJES)/Jet().pt;
      double scaleJES = get_syst_weight(Jet().pt_nom, Jet().pt_jesTotalUp, Jet().pt_jesTotalDown, nSigmaJES)/temp;
      Jet().pt *= scaleJES;
      Jet.v4() *= scaleJES;
      if (applySmearing) {
        //double scaleJER = get_syst_weight(Jet().pt_nom, Jet().pt_jerUp, Jet().pt_jerDown, nSigmaJER)/Jet().pt;
        double scaleJER = get_syst_weight(Jet().pt_nom, Jet().pt_jerUp, Jet().pt_jerDown, nSigmaJER)/temp;
        Jet().pt *= scaleJER;
        Jet.v4() *= scaleJER;
      }
			//cout << Jet().pt << ", ";
    }
    while (FatJet.Loop()) {
			temp = FatJet().pt_nom;
			FatJet().pt = FatJet().pt_nom;
      //double scaleJES = get_syst_weight(FatJet().pt_nom, FatJet().pt_jesTotalUp, FatJet().pt_jesTotalDown, nSigmaJES)/FatJet().pt;
      double scaleJES = get_syst_weight(FatJet().pt_nom, FatJet().pt_jesTotalUp, FatJet().pt_jesTotalDown, nSigmaJES)/temp;
      FatJet().pt *= scaleJES;
      FatJet.v4() *= scaleJES;
      if (applySmearing) {
        //double scaleJER = get_syst_weight(FatJet().pt_nom, FatJet().pt_jerUp, FatJet().pt_jerDown, nSigmaJER)/FatJet().pt;
        double scaleJER = get_syst_weight(FatJet().pt_nom, FatJet().pt_jerUp, FatJet().pt_jerDown, nSigmaJER)/temp;
        FatJet().pt *= scaleJER;
        FatJet.v4() *= scaleJER;
      }
			//cout << FatJet().pt << ", ";
      // SoftDrop mass uncertainty
      // JMS/JMR added to JES/JER in quadrature
			temp = FatJet().msoftdrop_nom;
			FatJet().msoftdrop = FatJet().msoftdrop_nom;
      double jesUNC = get_syst_weight(FatJet().msoftdrop_nom, FatJet().msoftdrop_jesTotalUp, FatJet().msoftdrop_jesTotalDown, nSigmaJES)/FatJet().msoftdrop_nom - 1;
      double jmsUNC = get_syst_weight(FatJet().msoftdrop_nom, FatJet().msoftdrop_jmsUp, FatJet().msoftdrop_jmsDown, nSigmaJES)/FatJet().msoftdrop_nom - 1;
      //scaleJES = get_syst_weight(FatJet().msoftdrop_nom, sqrt(jesUNC*jesUNC + jmsUNC*jmsUNC), nSigmaJES)/FatJet().msoftdrop;
      scaleJES = get_syst_weight(FatJet().msoftdrop_nom, sqrt(jesUNC*jesUNC + jmsUNC*jmsUNC), nSigmaJES)/temp;
      FatJet().msoftdrop *= scaleJES;
      if (applySmearing) {
        double jerUNC = get_syst_weight(FatJet().msoftdrop_nom, FatJet().msoftdrop_jerUp, FatJet().msoftdrop_jerDown, nSigmaJER)/FatJet().msoftdrop_nom - 1;
        double jmrUNC = get_syst_weight(FatJet().msoftdrop_nom, FatJet().msoftdrop_jmrUp, FatJet().msoftdrop_jmrDown, nSigmaJER)/FatJet().msoftdrop_nom - 1;
        //double scaleJER = get_syst_weight(FatJet().msoftdrop_nom, sqrt(jerUNC*jerUNC + jmrUNC*jmrUNC), nSigmaJES)/FatJet().msoftdrop;
        double scaleJER = get_syst_weight(FatJet().msoftdrop_nom, sqrt(jerUNC*jerUNC + jmrUNC*jmrUNC), nSigmaJES)/temp;
        FatJet().msoftdrop *= scaleJER;
      }
    }

		//cout << "MET : ";
    // MET uncertainties
    // JES/JER and unclustered energy variations
		if (!isSignal) {
      if (applySmearing) {
        MET_pt   = get_syst_weight(MET_T1Smear_pt,  MET_T1Smear_pt_jesTotalUp,  MET_T1Smear_pt_jesTotalDown,  nSigmaJES);
        MET_phi  = get_syst_weight(MET_T1Smear_phi, MET_T1Smear_phi_jesTotalUp, MET_T1Smear_phi_jesTotalDown, nSigmaJES);
        MET_pt  *= get_syst_weight(MET_T1Smear_pt,  MET_T1Smear_pt_jerUp,  MET_T1Smear_pt_jerDown,  nSigmaJER)/MET_T1Smear_pt;
        MET_phi *= get_syst_weight(MET_T1Smear_phi, MET_T1Smear_phi_jerUp, MET_T1Smear_phi_jerDown, nSigmaJER)/MET_T1Smear_phi;
      } else {
        MET_pt   = get_syst_weight(MET_T1_pt,  MET_T1_pt_jesTotalUp,  MET_T1_pt_jesTotalDown,  nSigmaJES);
        MET_phi  = get_syst_weight(MET_T1_phi, MET_T1_phi_jesTotalUp, MET_T1_phi_jesTotalDown, nSigmaJES);
        MET_pt  *= get_syst_weight(MET_T1_pt,  MET_T1_pt_jerUp,  MET_T1_pt_jerDown,  nSigmaJER)/MET_T1_pt;
        MET_phi *= get_syst_weight(MET_T1_phi, MET_T1_phi_jerUp, MET_T1_phi_jerDown, nSigmaJER)/MET_T1_phi;
      }
		} else {
			MET_pt = MET_T1Smear_pt;
			MET_phi = MET_T1Smear_phi;
		}

    // MET correction
    // https://twiki.cern.ch/twiki/bin/view/CMS/MissingETRun2Corrections#xy_Shift_Correction_MET_phi_modu
    std::pair<double,double> corr_MET = METXYCorr_Met_MetPhi(MET_pt, MET_phi, run, year, !isData, PV_npvs);
    MET_pt  = corr_MET.first;
    MET_phi = corr_MET.second;

		//cout << MET_pt << endl;
    // Save and load nominal values
    if (syst_index==0) {
      saved_MET_pt  = MET_pt;
      saved_MET_phi = MET_phi;
      // Save
      saved_Jet_v4.clear();
      saved_Jet_pt.clear();
      saved_FatJet_v4.clear();
      saved_FatJet_pt.clear();
      saved_FatJet_msoftdrop.clear();
      while (Jet.Loop()) {
        saved_Jet_v4.push_back(Jet.v4());
        saved_Jet_pt.push_back(Jet().pt);
      }
      while (FatJet.Loop()) {
        saved_FatJet_v4.push_back(FatJet.v4());
        saved_FatJet_pt.push_back(FatJet().pt);
        saved_FatJet_msoftdrop.push_back(FatJet().msoftdrop);
      }
    }
  }

  void initObjects() {
    Electron.initObjects();
    Muon.initObjects();
    Tau.initObjects();
    Photon  .initObjects();
    Jet.initObjects();
    FatJet.initObjects();
    GenJet.initObjects();
    SubJet.initObjects();
    GenPart.initObjects();
    initEventData();
  }

  
private:
  
  // Definitions of the internal methods defining variables and object selections

  void define_leptons_and_photons_(int debug = 0) {
    if (debug) std::cout<<"Variables::define_leptons_and_photons_: start (new event)"<<std::endl;

    // Initial loop on AK4 jets 
    // only needed for ele/muon 2D cut (which relies on the associated jet)
    while (Jet.Loop()) {
      double abseta = std::abs(Jet().eta);
      double NHF  = Jet().neHEF;
      double CHF  = Jet().chHEF;
      double NEMF = Jet().neEmEF;
      double CEMF = Jet().chEmEF;
      double CF   = CHF+CEMF;
      double MUF  = Jet().muEF;
      int NumConst = Jet().nConstituents;
      bool tightLepVetoJetID = 0;
      if (year==2016) {
        tightLepVetoJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((abseta<=2.4 && CHF>0 && CF>0 && CEMF<0.90) || abseta>2.4) && abseta<=2.7;
      } else if (year==2017) {
        tightLepVetoJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((abseta<=2.4 && CHF>0 && CF>0 && CEMF<0.90) || abseta>2.4) && abseta<=2.7;
      } else if (year==2018) {
        tightLepVetoJetID = (abseta<=2.6 && CEMF<0.8 && CF>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 ); 
      }
			bool pujet = true;
      if (year==2016) {
				if(Jet().pt < 40 && Jet().pt > 30 && Jet().puIdDisc < -0.71) pujet = false;
				else if(Jet().pt < 50 && Jet().pt > 40 && Jet().puIdDisc < -0.42) pujet = false;
			} else if (year==2017) {
				if(Jet().pt < 40 && Jet().pt > 30 && Jet().puIdDisc < -0.63) pujet = false;
				else if(Jet().pt < 50 && Jet().pt > 40 && Jet().puIdDisc < -0.19) pujet = false;
			} else {
				if(Jet().pt < 40 && Jet().pt > 30 && Jet().puIdDisc < -0.63) pujet = false;
				else if(Jet().pt < 50 && Jet().pt > 40 && Jet().puIdDisc < -0.19) pujet = false;
			}
      
      Jet.Jet.define( tightLepVetoJetID &&
                      //pujet &&
                      Jet().pt            >= JET_AK4_PT_CUT &&
                      std::abs(Jet().eta)  < JET_AK4_ETA_CUT);
    }
    if (debug) std::cout<<"Variables::define_leptons_and_photons_: end AK4 jet definition"<<std::endl;
  
    if (debug) std::cout<<"Variables::define_leptons_and_photons_: init AK8 jet collection"<<std::endl;

    // Electrons - full definitions
    while (Electron.Loop()) {
      float pt      = Electron().pt;
      float abseta  = std::abs(Electron().eta);
      float absd0   = std::abs(Electron().dxy);
      float absdz   = std::abs(Electron().dz);
      float miniIso = Electron().miniPFRelIso_all;
      //float ipsig   = std::abs(Electron().sip3d);
      bool id_noIso_WPL = 0;
      if (year==2016||year==2018) id_noIso_WPL = Electron().mvaFall17V2noIso_WPL;
      else if (year==2017)        id_noIso_WPL = Electron().mvaFall17V2noIso_WPL;
      bool id_noIso_WP90 = 0;
      if (year==2016||year==2018) id_noIso_WP90 = Electron().mvaFall17V2noIso_WP90;
      else if (year==2017)        id_noIso_WP90 = Electron().mvaFall17V2noIso_WP90;
      bool id_Iso_WP90 = 0;
      if (year==2016||year==2018) id_Iso_WP90 = Electron().mvaFall17V2Iso_WP90;
      else if (year==2017)        id_Iso_WP90 = Electron().mvaFall17V2Iso_WP90;

      // Calculate the B2G 2D cut variables
      // [Lepton]_jetPtRelv2 is available in v5 NanoAOD, but it is bugged :(
      // Try to do what (should be) done in NanoAOD (revert JEC, subtract ele, reapply JEC) based on the original NanoAOD recipe:
      // https://github.com/cms-nanoAOD/cmssw/blob/master-cmsswmaster/PhysicsTools/NanoAOD/plugins/LeptonJetVarProducer.cc#L152-L162
      // 1) if no jetIdx, or no jet id --> match nearest jet passing the ID for DR, leave Ptrel = NOVAL_F
      // 2) if    jetIdx -> calculate DR
      // 3)              -> subtract lep pt from jet, check if passs pt cut (if not Ptrel=NOVAL_F), calculate Ptrel
      // Calculate dRmin for the associated/nearest jet
      int jetIdx = Electron().jetIdx;
      Electron().isPartOfJet = true;
      if (jetIdx==-1||jetIdx>=Jet.n) Electron().isPartOfJet = false;
      else if (!Jet.Jet.pass[jetIdx]) Electron().isPartOfJet = false;
      if (debug>1) std::cout<<"Start calculating 2D cut for - iEle="<<Electron.i<<" isPartOfJet="<<Electron().isPartOfJet<<" matched jetIdx="<<jetIdx<<std::endl;
      if (Electron().isPartOfJet) {
        Electron().jetDRmin = DeltaR(Electron.v4(), Jet.v4(jetIdx));
      } else while (Jet.Jet.Loop()) {
        double dR = DeltaR(Electron.v4(),Jet.Jet.v4());
        if (dR<Electron().jetDRmin) {
          Electron().jetDRmin = dR;
          if (debug>1) std::cout<<" - matched to jetdIdx="<<Jet.Jet.iRef<<" dR="<<dR<<std::endl; 
        }
      }
      if (debug>1) std::cout<<"jetIdx="<<jetIdx<<std::endl;
      // Calculate pTrel for the associated/nearest jet
      // First subtract the ele from the jet
      // The JEC reverted ele energy subtraction should work in the mean time
      if (jetIdx!=-1&&jetIdx<Jet.n) {
        LorentzVector cleanjet_v4 = Jet.v4(jetIdx);
        if (debug>1) std::cout<<"-       jet pt="<<cleanjet_v4.Pt()<<std::endl;
        cleanjet_v4 -= Electron.v4()*(1/(1-Jet(jetIdx).rawFactor));
        // Then, for version 2 check also that the cleaned jet passes the minimum pt threshold cut that is considered in the analysis
        // If not, revert back to the previous closest jet
        if (debug>1) std::cout<<"- clean jet pt="<<cleanjet_v4.Pt()<<std::endl;
        if (cleanjet_v4.Pt() >= JET_AK4_PT_CUT) {
          Electron().cleanJetPtrel = Perp(Electron.v4().Vect(), cleanjet_v4.Vect());
          if (debug>1) std::cout<<"  + pass pt cut, Ptrel="<<Electron().cleanJetPtrel<<std::endl;
        }
      }
      if (debug>1) std::cout<<"-------------------------------------------------------------"<<std::endl;
      if (debug>1) std::cout<<"ele jet DR="<<Electron().jetDRmin<<" Ptrel="<<Electron().cleanJetPtrel<<std::endl<<std::endl;

      // Reconstruct the neutrino 4 momentum using W mass constraint
      // https://twiki.cern.ch/twiki/bin/view/Main/TopPairProduction
      // https://github.com/BoostedScalefactors/WTopScalefactorProducer/blob/master/Skimmer/python/variables.py
      // Do the calculation for all ellectron and muon candidates
      double MET_px = MET_pt*std::cos(MET_phi);
      double MET_py = MET_pt*std::sin(MET_phi);
      double a = 80.379*80.379 - Electron.v4().M()*Electron.v4().M() + 2*Electron.v4().Px()*MET_px + 2*Electron.v4().Py()*MET_py;
      double A =  4 * (Electron.v4().E()*Electron.v4().E() - Electron.v4().Pz()*Electron.v4().Pz());
      double B = -4 * a * Electron.v4().Pz();
      double C =  4 * (Electron.v4().E()*Electron.v4().E()) * (MET_px*MET_px + MET_py*MET_py) - a*a;
      double D = B*B - 4*A*C;
      // If there are real solutions, use the one with lowest Pz                                            
      double MET_pz = D>=0 ? std::min((-B+std::sqrt(D))/(2*A), (-B-std::sqrt(D))/(2*A)) : -B/(2*A);
      Electron().nu.SetXYZT(MET_px, MET_py, MET_pz, std::sqrt(MET_px*MET_px+MET_py*MET_py+MET_pz*MET_pz));
      Electron().nuDR = DeltaR(Electron().nu, Electron.v4());

      // Assoicated AK8 jet
      double mindR_AK8 = NOVAL_F;
      size_t iMatchAK8 = -1;
      while(FatJet.Loop()) {
        double dR = DeltaR(Electron.v4(), FatJet.v4());
        if (dR<0.8 && dR<mindR_AK8) {
          mindR_AK8 = dR;
          Electron().matchAK8 = true;
          iMatchAK8 = FatJet.i;
        }
      }
      if (Electron().matchAK8) {
        Electron().iMatchedAK8 = iMatchAK8;
        Electron().nuMatchedAK8DR = DeltaR(Electron().nu, FatJet.v4(iMatchAK8));
      }

      // Selected objects
      // Cut-based veto
      if (Electron.CBVetoNoIso.define(Electron().cutBased &&
                                    pt      >= ELE_VETO_PT_CUT &&
                                    abseta  <  ELE_VETO_ETA_CUT && //!(abseta>=1.442 && abseta< 1.556) &&
                                    absd0   <  0.2 &&
                                    absdz   <  0.5))
        Electron.CBVeto       .define((miniIso < 0.1));

      // Veto
      if (Electron.VetoNoIso.define(id_noIso_WPL &&
                                    pt      >= ELE_VETO_PT_CUT &&
                                    abseta  <  ELE_VETO_ETA_CUT && //!(abseta>=1.442 && abseta< 1.556) &&
                                    absd0   <  ELE_VETO_IP_D0_CUT &&
                                    absdz   <  ELE_VETO_IP_DZ_CUT))
        Electron.Veto       .define((miniIso <  ELE_VETO_MINIISO_CUT));

      // Select
      Electron.Select.define( id_noIso_WP90 &&
                              pt        >= ELE_SELECT_PT_CUT &&
                              abseta    <  ELE_SELECT_ETA_CUT && !(abseta>=1.442 && abseta< 1.556) &&
                              miniIso   <  ELE_SELECT_MINIISO_CUT &&
                              absd0     <  ELE_SELECT_IP_D0_CUT &&
                              absdz     <  ELE_SELECT_IP_DZ_CUT);
      
      if ( pt        >= ELE_TIGHT_PT_CUT &&
           abseta    <  ELE_TIGHT_ETA_CUT && !(abseta>=1.442 && abseta< 1.556) &&
           //ipsig     <  ELE_TIGHT_IP_SIG_CUT &&
           absd0     <  ELE_TIGHT_IP_D0_CUT &&
           absdz     <  ELE_TIGHT_IP_DZ_CUT) {
        // Tight - ID with isolation - Medium MVA ID
        Electron.Tight.define(id_Iso_WP90);
        // Non-isolated lepton - Loose MVA ID
        if (Electron.NoIso.define(id_noIso_WPL)) {
          // + Loose 2D isolation
          Electron.NonIso.define(!(Electron().jetDRmin<0.4 && Electron().jetPtRelv2<15));
        }
      }      
    }
    if (debug) std::cout<<"Variables::define_leptons_and_photons_: end electron definitions"<<std::endl;


    // Muons - full definitions
    while (Muon.Loop()) {
      float pt      = Muon().pt;
      float abseta  = std::abs(Muon().eta);
      float absd0   = std::abs(Muon().dxy);
      float absdz   = std::abs(Muon().dz);
      float miniIso = Muon().miniPFRelIso_all;
      float ipsig   = std::abs(Muon().sip3d);
      float relIso = Muon().pfRelIso04_all;

      // Calculate the B2G 2D cut variables
      // [Lepton]_jetPtRelv2 is available in v5 NanoAOD, but it is bugged :(
      // Try to do what (should be) done in NanoAOD (revert JEC, subtract ele, reapply JEC) based on the original NanoAOD recipe:
      // https://github.com/cms-nanoAOD/cmssw/blob/master-cmsswmaster/PhysicsTools/NanoAOD/plugins/LeptonJetVarProducer.cc#L152-L162
      // 1) if no jetIdx, or no jet id --> match nearest jet passing the ID for DR, leave Ptrel = NOVAL_F
      // 2) if    jetIdx -> calculate DR
      // 3)              -> subtract lep pt from jet, check if passs pt cut (if not Ptrel=NOVAL_F), calculate Ptrel
      // Calculate dRmin for the associated/nearest jet
      int jetIdx = Muon().jetIdx;
      Muon().isPartOfJet = true;
      if (jetIdx==-1||jetIdx>=Jet.n) Muon().isPartOfJet = false;
      else if (!Jet.Jet.pass[jetIdx]) Muon().isPartOfJet = false;
      if (debug>1) std::cout<<"Start calculating 2D cut for - iMu="<<Muon.i<<" isPartOfJet="<<Muon().isPartOfJet<<" matched jetIdx="<<jetIdx<<std::endl;
      if (Muon().isPartOfJet) {
        Muon().jetDRmin = DeltaR(Muon.v4(), Jet.v4(jetIdx));
      } else while (Jet.Jet.Loop()) {
        double dR = DeltaR(Muon.v4(), Jet.Jet.v4());
        if (dR<Muon().jetDRmin) {
          Muon().jetDRmin = dR;
          if (debug>1) std::cout<<" - matched to jetdIdx="<<Jet.Jet.iRef<<" dR="<<dR<<std::endl; 
        }
      }
      if (debug>1) std::cout<<"jetIdx="<<jetIdx<<std::endl;
      // Calculate pTrel for the associated/nearest jet
      // First subtract the muon from the jet
      // Similar to [Lepton]_jetPtRelv2, a useful variable Jet_muonSubtrFactor is in NanoAOD v5
      // Problem is that this variable is also bugged :(
      // The JEC reverted muon energy subtraction should work in the mean time (similar to electrons)
      if (jetIdx!=-1&&jetIdx<Jet.n) {
        LorentzVector cleanjet_v4 = Jet.v4(jetIdx);
        if (debug>1) std::cout<<"-       jet pt="<<cleanjet_v4.Pt()<<std::endl;
        //cleanjet_v4 *= Jet(jetIdx).muonSubtrFactor; // Variable bugged
        cleanjet_v4 -= Muon.v4()*(1/(1-Jet(jetIdx).rawFactor));
        // Then, for version 2 check also that the cleaned jet passes the minimum pt threshold cut that is considered in the analysis
        // If not, revert back to the previous closest jet
        if (debug>1) std::cout<<"- clean jet pt="<<cleanjet_v4.Pt()<<std::endl;
        if (cleanjet_v4.Pt() >= JET_AK4_PT_CUT) {
          Muon().cleanJetPtrel = Perp(Muon.v4().Vect(), cleanjet_v4.Vect());
          if (debug>1) std::cout<<"  + pass pt cut, Ptrel="<<Muon().cleanJetPtrel<<std::endl;
        }
      }
      if (debug>1) std::cout<<"-------------------------------------------------------------"<<std::endl;
      if (debug>1) std::cout<<"mu jet DRmin="<<Muon().jetDRmin<<" Ptrel="<<Muon().cleanJetPtrel<<std::endl<<std::endl;

      // Reconstruct the neutrino 4 momentum using W mass constraint
      // https://twiki.cern.ch/twiki/bin/view/Main/TopPairProduction
      // https://github.com/BoostedScalefactors/WTopScalefactorProducer/blob/master/Skimmer/python/variables.py
      // Do the calculation for all ellectron and muon candidates
      double MET_px = MET_pt*std::cos(MET_phi);
      double MET_py = MET_pt*std::sin(MET_phi);
      double a = 80.379*80.379 - Muon.v4().M()*Muon.v4().M() + 2*Muon.v4().Px()*MET_px + 2*Muon.v4().Py()*MET_py;
      double A =  4 * (Muon.v4().E()*Muon.v4().E() - Muon.v4().Pz()*Muon.v4().Pz());
      double B = -4 * a * Muon.v4().Pz();
      double C =  4 * (Muon.v4().E()*Muon.v4().E()) * (MET_px*MET_px + MET_py*MET_py) - a*a;
      double D = B*B - 4*A*C;
      // If there are real solutions, use the one with lowest Pz                                            
      double MET_pz = D>=0 ? std::min((-B+std::sqrt(D))/(2*A), (-B-std::sqrt(D))/(2*A)) : -B/(2*A);
      Muon().nu.SetXYZT(MET_px, MET_py, MET_pz, std::sqrt(MET_px*MET_px+MET_py*MET_py+MET_pz*MET_pz));
      Muon().nuDR = DeltaR(Muon().nu, Muon.v4());
    
      // Assoicated AK8 jet
      double mindR_AK8 = NOVAL_F;
      size_t iMatchAK8 = -1;
      while(FatJet.Loop()) {
        double dR = DeltaR(Muon.v4(), FatJet.v4());
        if (dR<0.8 && dR<mindR_AK8) {
          mindR_AK8 = dR;
          Muon().matchAK8 = true;
          iMatchAK8 = FatJet.i;
        }
      }
      if (Muon().matchAK8) {
        Muon().iMatchedAK8 = iMatchAK8;
        Muon().nuMatchedAK8DR = DeltaR(Muon().nu, FatJet.v4(iMatchAK8));
      }
    
      // Selected objects
      // Veto
      if (Muon.VetoNoIso.define(Muon().softId &&
                                pt      >= MU_VETO_PT_CUT &&
                                abseta  <  MU_VETO_ETA_CUT &&
                                absd0   <  MU_VETO_IP_D0_CUT &&
                                absdz   <  MU_VETO_IP_DZ_CUT))
        Muon.Veto       .define(miniIso <  MU_VETO_MINIISO_CUT);
      // Cut-based loose
      if (Muon.CBLooseNoIso.define(Muon().looseId &&
                                pt      >= MU_VETO_PT_CUT &&
                                abseta  <  MU_VETO_ETA_CUT &&
                                absd0   <  MU_VETO_IP_D0_CUT &&
                                absdz   <  MU_VETO_IP_DZ_CUT))
        //Muon.CBLoose       .define(miniIso < 0.2);
        Muon.CBLoose       .define(miniIso < 1.2);
      // Cut-based medium
      if (Muon.CBMediumNoIso.define(Muon().mediumId &&
                                pt      >= MU_VETO_PT_CUT &&
                                abseta  <  MU_VETO_ETA_CUT &&
                                absd0   <  MU_VETO_IP_D0_CUT &&
                                absdz   <  MU_VETO_IP_DZ_CUT))
        //Muon.CBMedium      .define(miniIso < 0.2);
        Muon.CBMedium      .define(miniIso < 1.2);
      // Select
      Muon.Select.define( Muon().mediumPromptId &&
                          pt        >= MU_SELECT_PT_CUT &&
                          abseta    <  MU_SELECT_ETA_CUT &&
                          miniIso   <  MU_SELECT_MINIISO_CUT &&
                          absd0     <  MU_SELECT_IP_D0_CUT &&
                          absdz     <  MU_SELECT_IP_DZ_CUT);

      if (pt        >= MU_TIGHT_PT_CUT &&
          abseta    <  MU_TIGHT_ETA_CUT &&
          ipsig     <  MU_TIGHT_IP_SIG_CUT &&
          absd0     <  MU_TIGHT_IP_D0_CUT &&
          absdz     <  MU_TIGHT_IP_DZ_CUT) {
        // Tight - ID without isolation + relIso
        Muon.Tight.define(Muon().tightId &&
                          relIso  <  MU_TIGHT_RELISO_CUT);
        // Non-isolated lepton - Loose MVA ID
        if (Muon.NoIso.define(Muon().softId)) {
        //if (Muon.NoIso.define(Muon().mvaId>0)) {
          // + Loose 2D isolation
          Muon.NonIso.define(!(Muon().jetDRmin<0.4 && Muon().jetPtRelv2<15));
        }
      }
    }
    if (debug) std::cout<<"Variables::define_leptons_and_photons_: end muon definitions"<<std::endl;

    // Taus - full definitions
    while (Tau.Loop()) {
      // Use Loose working point
      int WP = 3; // 0: VLoose, 1: Loose, 2: Medium, 3: Tight
      int bm_jet = 4<<WP;
      WP = 1; // 0: VLoose, 1: Loose, 2: Medium, 3: Tight
      int bm_e   = 4<<WP;
      int bm_mu  = 1<<WP;
      Tau.Veto.define(( Tau().pt            >= TAU_VETO_PT_CUT &&
                        std::abs(Tau().eta) <  TAU_VETO_ETA_CUT &&
                        ( ((Tau().idDeepTau2017v2p1VSjet & bm_jet) == bm_jet) ||
                          ((Tau().idDeepTau2017v2p1VSe   & bm_e  ) == bm_e  ) ||
                          ((Tau().idDeepTau2017v2p1VSmu  & bm_mu ) == bm_mu ) ) ));
    }
    // No pt/eta cut
    if (debug) std::cout<<"Variables::define_leptons_and_photons_: end tau definitions"<<std::endl;
    
    // Photons - full definitions
    while (Photon.Loop()) {
      float pt = Photon().pt;
      float abseta = std::abs(Photon().eta);
      bool ele_veto  = Photon().electronVeto;
      bool pixel_veto = Photon().pixelSeed;
      //bool id_select = (Photon().Photon_cutBasedBitmap&2)==2;
  
      // Implementing cuts for Fall17 V2 ID
      // https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2?rev=49#Working_points_for_94X_and_later
      // https://github.com/cms-sw/cmssw/blob/master/RecoEgamma/PhotonIdentification/python/Identification/cutBasedPhotonID_Fall17_94X_V2_cff.py#L53-L81
      // Revert Sigma_ietaieta cut in order to be able to measure QCD background
      // Fit shape of the charge distribution
      // Do not apply both cuts for preselection
      // select fakes with reverted sieie cut
      
      bool id_preselect = false, id_select = false;
      // Cut based ID, without
      if (Photon().isScEtaEB){
        // Barrel cuts (EB)
        id_preselect =
          Photon().hoe                                      < 0.02197 &&
          Photon().pfRelIso03_all - Photon().pfRelIso03_chg < 1.189+0.01512*pt+2.259e-05*pt*pt &&
          Photon().pfRelIso03_all                           < 2.08+0.004017*pt;
        id_select = id_preselect &&
          Photon().sieie                                    < 0.01015 &&
          Photon().pfRelIso03_chg                           < 1.141;
      } else {
        // Encap cuts (EE)
        id_preselect =
          Photon().hoe                                      < 0.0326 &&
          Photon().pfRelIso03_all - Photon().pfRelIso03_chg < 2.718+0.0117*pt+2.3e-05*pt*pt &&
          Photon().pfRelIso03_all                           < 3.867+0.0037*pt;
        id_select = id_preselect &&
          Photon().sieie                                    < 0.0272 &&
          Photon().pfRelIso03_chg                           < 1.051;
      }
      // Medium ID without Sigma_ietaieta cut
      if (Photon.PreSelect.define(id_preselect &&
                                  ele_veto &&
                                  !pixel_veto &&
                                  pt        >= PHOTON_SELECT_PT_CUT &&
                                  abseta    <  PHOTON_SELECT_ETA_CUT )) {
        // Fake photons (those that fail the SigmaIeteIeta cut)
        Photon.Fake.define(Photon().sieie >= (Photon().isScEtaEB ? 0.01015 : 0.0272));
        Photon.SelectNoIso.define(Photon().sieie < (Photon().isScEtaEB ? 0.01015 : 0.0272));
      }
  
      // Photons passing full ID
      Photon.Select.define(id_select &&
                           ele_veto &&
                           !pixel_veto &&
                           pt        >= PHOTON_SELECT_PT_CUT &&
                           abseta    <  PHOTON_SELECT_ETA_CUT );
    }
    if (debug) std::cout<<"Variables::define_leptons_and_photons_: end photon definitions"<<std::endl;
  }
 
  void define_jets_(const unsigned int& syst_index, int debug = 0) {
    Jet.initObjects();
    FatJet.initObjects();
    SubJet.initObjects();
    GenJet.initObjects();
		if (debug < 0) std::cout << int((syst_index-1)/2) << ", " ;
    // AK4 jets - Definitions except those depending on AK8 jets
    while (Jet.Loop()) {
      if (debug>1) std::cout<<"Variables::define_jets_: AK4 "<<Jet.i<<" start"<<std::endl;
		  if (debug < 0) std::cout << Jet().pt << ", " ;
			//if(syst_index == 0) Jet().pt = Jet().pt_nom;
      // Jet ID
      double abseta = std::abs(Jet().eta);
      double NHF  = Jet().neHEF;
      double CHF  = Jet().chHEF;
      double NEMF = Jet().neEmEF;
      double CEMF = Jet().chEmEF;
      double CF   = CHF+CEMF;
      double MUF  = Jet().muEF;
      int NumConst = Jet().nConstituents;
      bool tightLepVetoJetID = 0;
      if (year==2016) {
        tightLepVetoJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((abseta<=2.4 && CHF>0 && CF>0 && CEMF<0.80) || abseta>2.4) && abseta<=2.7;
      } else if (year==2017 || year==2018) {
        tightLepVetoJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((abseta<=2.4 && CHF>0 && CF>0 && CEMF<0.80) || abseta>2.4) && abseta<=2.7;
      }
			bool pujet = true;
      if (year==2016) {
				if(Jet().pt < 40 && Jet().pt > 30 && Jet().puIdDisc < -0.71) pujet = false;
				else if(Jet().pt < 50 && Jet().pt > 40 && Jet().puIdDisc < -0.42) pujet = false;
			} else if (year==2017) {
				if(Jet().pt < 40 && Jet().pt > 30 && Jet().puIdDisc < -0.63) pujet = false;
				else if(Jet().pt < 50 && Jet().pt > 40 && Jet().puIdDisc < -0.19) pujet = false;
			} else {
				if(Jet().pt < 40 && Jet().pt > 30 && Jet().puIdDisc < -0.63) pujet = false;
				else if(Jet().pt < 50 && Jet().pt > 40 && Jet().puIdDisc < -0.19) pujet = false;
			}
      Jet.FailID.define( !tightLepVetoJetID &&
                         Jet().pt            >= JET_AK4_PT_CUT &&
                         std::abs(Jet().eta)  < JET_AK4_ETA_CUT);
      if (Jet.Jet.define( tightLepVetoJetID &&
                          //pujet &&
                          Jet().pt            >= JET_AK4_PT_CUT &&
                          std::abs(Jet().eta)  < JET_AK4_ETA_CUT)) {
        
  	while (GenJet.Loop()) {
			if(sqrt((Jet().phi-GenJet().phi)*(Jet().phi-GenJet().phi) + (Jet().eta-GenJet().eta)*(Jet().eta-GenJet().eta)) < 0.4) Jet().matchJet=true;
		}
        if (debug>1) std::cout<<"Variables::define_jets_: AK4 "<<Jet.i<<" id ok"<<std::endl;
        
        // b tagging
        if (year==2018) {
          Jet.LooseBTag .define(Jet().btagDeepFlavB >= 0.0532);
          Jet.MediumBTag.define(Jet().btagDeepFlavB >= 0.3040);
          Jet.TightBTag .define(Jet().btagDeepFlavB >= 0.7476);
        } else if (year==2017) {
          Jet.LooseBTag .define(Jet().btagDeepFlavB >= 0.0490);
          Jet.MediumBTag.define(Jet().btagDeepFlavB >= 0.2783);
          Jet.TightBTag .define(Jet().btagDeepFlavB >= 0.7100);
        } else {
        	if (isAPV) {
          	Jet.LooseBTag .define(Jet().btagDeepFlavB >= 0.0508);
          	Jet.MediumBTag.define(Jet().btagDeepFlavB >= 0.2598);
          	Jet.TightBTag .define(Jet().btagDeepFlavB >= 0.6502);
        	} else {
          	Jet.LooseBTag .define(Jet().btagDeepFlavB >= 0.0480);
          	Jet.MediumBTag.define(Jet().btagDeepFlavB >= 0.2489);
          	Jet.TightBTag .define(Jet().btagDeepFlavB >= 0.6377);
					}
        }
        if (debug>1) std::cout<<"Variables::define_jets_: AK4 "<<Jet.i<<" b-tagging ok"<<std::endl;
        
        // Lepton-jet overlap
        while (Electron.Loop()) {
          double dR = DeltaR(Electron.v4(), Jet.v4());
          if (dR<Electron().jetDR) {
            Electron().iMatchedAK4 = Jet.i;
            Electron().jetDR   = dR;
            Electron().jetDPhi = std::abs(DeltaPhi(Electron().phi, Jet().phi));
          }
          if (Electron.Veto.pass[Electron.i]&&dR<Jet().eleDR) {
            Jet().eleDR = dR;
            Jet().elePtRatio = Electron().pt/Jet().pt;
          }
        }
        while (Muon.Loop()) {
          double dR = DeltaR(Muon.v4(), Jet.v4());
          if (dR<Muon().jetDR) {
            Muon().iMatchedAK4 = Jet.i;
            Muon().jetDR   = dR;
            Muon().jetDPhi = std::abs(DeltaPhi(Muon().phi, Jet().phi));
          }
          if (Muon.Veto.pass[Muon.i]&&dR<Jet().muDR) {
            Jet().muDR = dR;
            Jet().muPtRatio = Muon().pt/Jet().pt;
          }
        }
        if (debug>1) std::cout<<"Variables::define_jets_: AK4 "<<Jet.i<<" lepton-overlap ok"<<std::endl;
        // Photon-jet overlap
        while (Photon.Select.Loop()) {
          double dR = DeltaR(Photon.Select.v4(), Jet.v4());
          if (dR<Jet().phoDR) {
            Jet().phoDR = dR;
            Jet().phoPtRatio = Photon.Select().pt/Jet().pt;
          }
        }
        // Exclude jets that have overlapping photon
        if (Photon.Select.n==1) {
          if (Jet.JetNoPho.define(Photon.Select(0).jetIdx!=Jet.i))
            Jet.MediumBTagNoPho.define(Jet().btagDeepB >= B_DEEP_MEDIUM_CUT);
        } else if (Photon.PreSelect.n==1) {
          if (Jet.JetNoPho.define(Photon.PreSelect(0).jetIdx!=Jet.i))
            Jet.MediumBTagNoPho.define(Jet().btagDeepB >= B_DEEP_MEDIUM_CUT);
        }
        if (debug>1) std::cout<<"Variables::define_jets_: AK4 "<<Jet.i<<" photon-overlap ok"<<std::endl;
      }
      
    } // End AK4 Jet Loop
    if (debug) std::cout<<"Variables::define_leptons_and_photons_: AK4 ordinary jets end"<<std::endl;


    /*
    // Match non isolated tight lepton to nearest AK4 jet to be able to tell which Razor megajet it belongs
    for (size_t i=0; i<noiso_leptons.size(); ++i) {
      bool match = false;
      if (iJetLepNoIso[i]==-1) match = true;
      else if (!Jet.Jet.pass[iJetLepNoIso[i]]) match = true;
      if (match) {
        double mindR = NOVAL_F;
        for(size_t j=0; j<selected_jets.size(); ++j) {
          double dR = DeltaR(noiso_leptons[i], selected_jets[j]);
          if (dR<mindR) {
            mindR = dR;
            iJetLepNoIso[i] = iJet[j];
          }
        }
      }
    }
    for (size_t i=0; i<noniso_leptons.size(); ++i) {
      bool match = false;
      if (iJetLepNonIso[i]==-1) match = true;
      else if (!Jet.Jet.pass[iJetLepNonIso[i]]) match = true;
      if (match) {
        double mindR = NOVAL_F;
        for(size_t j=0; j<selected_jets.size(); ++j) {
          double dR = DeltaR(noniso_leptons[i], selected_jets[j]);
          if (dR<mindR) {
            mindR = dR;
            iJetLepNonIso[i] = iJet[j];
          }
        }
      }
    }
    */
  
  void define_genparticles_(int debug = 0) {
    if (debug) std::cout<<"Variables::define_genparticles_: start"<<std::endl;

    // Loop on generator particles
    // 1st Loop is a recursive search for (grand/)mothers
    // set flags based on (grand/)daughters (eg. leptonic decays for W/top)
    // TODO: debug Data, which says the collection size is not 0
    while (GenPart.Loop()) {
      if (debug>1) std::cout<<"Variables::define_genparticles_: gen 1st loop "<<GenPart.i<<" start"<<std::endl;
			if(GenPart().pdgId == 0) continue;
      int motherId = -NOVAL_I, grandMotherId = -NOVAL_I, greatGrandMotherId = -NOVAL_I;
      int iMother = GenPart().genPartIdxMother, iGrandMother = -1, iGreatGrandMother = -1;
      // Set the mother to the first ancestor that has different pdg id
      if (iMother!=-1&&iMother<GenPart.n) {
        while (GenPart(iMother).pdgId == GenPart().pdgId) {
          GenPart(iMother).NoSameDaughter = 0;
          iMother = GenPart(iMother).genPartIdxMother;
          if (iMother==-1) break;
        }
				if(debug>1) std::cout << "Test " << std::endl;
        if (iMother!=-1) {
          motherId = GenPart(iMother).pdgId;
          // The mother is a candidate for final non-decayed state
          // But have to check that it does not have a daughter with the same pdgId
          // found in the previous iteration
          GenPart(iMother).LastCopyCand = 1;
          iGrandMother = GenPart(iMother).genPartIdxMother;
          if (iGrandMother!=-1) {
            // Remove intermediate states from history (where id matches that of the mother)
            while (GenPart(iGrandMother).pdgId == GenPart(iMother).pdgId) {
              iGrandMother = GenPart(iGrandMother).genPartIdxMother;
              if (iGrandMother==-1) break;
            }
            if (iGrandMother!=-1) {
              grandMotherId = GenPart(iGrandMother).pdgId;
              // This seems a bit too much, but there are cases like: top->W->tau->ele/mu
              iGreatGrandMother = GenPart(iGrandMother).genPartIdxMother;
              if (iGreatGrandMother!=-1) {
                while (GenPart(iGreatGrandMother).pdgId == GenPart(iGrandMother).pdgId) {
                  iGreatGrandMother = GenPart(iGreatGrandMother).genPartIdxMother;
                  if (iGreatGrandMother==-1) break;
                }
                if (iGreatGrandMother!=-1) greatGrandMotherId = GenPart(iGreatGrandMother).pdgId;
              }
            }
          }
        }
      }
      // Save info
      GenPart().iMother = iMother;
      GenPart().iGrandMother = iGrandMother;
      GenPart().iGreatGrandMother = iGreatGrandMother;
      GenPart().motherId = motherId;
      GenPart().grandMotherId = grandMotherId;
      GenPart().greatGrandMotherId = greatGrandMotherId;
      if (iMother!=-1) {
        GenPart(iMother).iDaughters.push_back(GenPart.i);
      }
      if (debug>1) std::cout<<"Variables::define_genparticles_: gen 1st loop "<<GenPart.i<<" (grand/)mother ok"<<std::endl;
      
      // Check if particle is the last copy
      // gen status flags stored bitwise, bits are: 
      // 0 : isPrompt, 1 : isDecayedLeptonHadron, 2 : isTauDecayProduct, 3 : isPromptTauDecayProduct, 
      // 4 : isDirectTauDecayProduct, 5 : isDirectPromptTauDecayProduct, 6 : isDirectHadronDecayProduct, 
      // 7 : isHardProcess, 8 : fromHardProcess, 9 : isHardProcessTauDecayProduct, 
      // 10 : isDirectHardProcessTauDecayProduct, 11 : fromHardProcessBeforeFSR, 
      // 12 : isFirstCopy, 13 : isLastCopy, 14 : isLastCopyBeforeFSR
      // Good presentation from Josh Bendavid (slide 13-17)
      // https://indico.cern.ch/event/402279/contributions/960421/attachments/805964/1104514/mcaod-Jun17-2015.pdf#search=isLastCopy%20GenPart
      // isLastCopy: last copy of a given particle.  More likely, but notguaranteed, to have "final" physical kinematics
      // --> Correctly give status=1 for e/mu, and status=2 for tau
      //     But is unreliable for Ws, there we need to recursively
      //     eg.: lep -> W (final) -> W (interm) -> top
      //     isLastCopy fails O(%) for the final state W, especially for electrons (as high as 3%)
      bool isLastCopy = (GenPart().statusFlags>>13)&1;
      if (isLastCopy) {
        // Generator final state leptons
        if (debug>1) std::cout<<"Variables::define_genparticles_: lastcopy ok"<<std::endl;
        if (GenPart.Lepton.define( std::abs(GenPart().pdgId)==11 ||
                                     std::abs(GenPart().pdgId)==13 ||
                                     std::abs(GenPart().pdgId)==15 )) {
          if (debug>1) std::cout<<"Variables::define_genparticles_: lepton start"<<std::endl;
          // Leptons from the hard process
          GenPart.LeptonFromHardProcess.define((GenPart().statusFlags>>8)&1);
          if (debug>1) std::cout<<"Variables::define_genparticles_: lepton from hard process ok"<<std::endl;
          // Leptons from W
          if (GenPart.LeptonFromW.define(std::abs(motherId)==24)) {
            if (debug>1) std::cout<<"Variables::define_genparticles_: lepton from W start"<<std::endl;
            GenPart.LepW.pass[iMother] = 1;
            GenPart(iMother).iGenLepDaughter = GenPart.i;
            // also from leptonic top
            if (GenPart.LeptonFromTop.define(std::abs(grandMotherId)==6)) {
              GenPart.LepTop.pass[iGrandMother] = 1;
              GenPart(iGrandMother).iGenLepGrandDaughter = GenPart.i;
            }
            // also from Higgs
            if (GenPart.LeptonFromH.define(std::abs(grandMotherId)==25)) {
              GenPart.LepH.pass[iGrandMother] = 1;
              GenPart(iGrandMother).iGenLepGrandDaughter = GenPart.i;
            }
            if (debug>1) std::cout<<"Variables::define_genparticles_: lepton from W ok"<<std::endl;
          }
          // Leptons from Z
          else if (GenPart.LeptonFromZ.define(std::abs(motherId)==23)) {
            if (debug>1) std::cout<<"Variables::define_genparticles_: lepton from Z start"<<std::endl;
            GenPart.LepZ.pass[iMother] = 1;
            GenPart(iMother).iGenLepDaughter = GenPart.i;
            // also from Higgs
            if (GenPart.LeptonFromH.define(std::abs(grandMotherId)==25)) {
              GenPart.LepH.pass[iGrandMother] = 1;
              GenPart(iGrandMother).iGenLepGrandDaughter = GenPart.i;
            }
            if (debug>1) std::cout<<"Variables::define_genparticles_: lepton from Z ok"<<std::endl;
          }
          // Leptons from h
          else if (GenPart.LeptonFromH.define(std::abs(motherId)==25)) {
            if (debug>1) std::cout<<"Variables::define_genparticles_: lepton from h start"<<std::endl;
            GenPart.LepH.pass[iMother] = 1;
            GenPart(iMother).iGenLepDaughter = GenPart.i;
            if (debug>1) std::cout<<"Variables::define_genparticles_: lepton from h ok"<<std::endl;
          }
          // Lepton from taus, here we go ... :)
          else if ( ( std::abs(GenPart().pdgId)==11 ||
                      std::abs(GenPart().pdgId)==13 ) &&
                    std::abs(motherId)==15) {
            if (debug>1) std::cout<<"Variables::define_genparticles_: lepton from taus start"<<std::endl;
            // also lepton from W
            if (GenPart.LeptonFromW.define(std::abs(grandMotherId)==24)) {
              GenPart.LepW.pass[iGrandMother] = 1;
              GenPart(iGrandMother).iGenLepGrandDaughter = GenPart.i;
              // also lepton from top
              if (GenPart.LeptonFromTop.define(std::abs(greatGrandMotherId)==6)) {
                GenPart.LepTop.pass[iGreatGrandMother] = 1;
                GenPart(iGreatGrandMother).iGenLepGreatGrandDaughter = GenPart.i;
              }
            }
          }
          if (debug>1) std::cout<<"Variables::define_genparticles_: lepton ok"<<std::endl;
        }
        // Generator neutrinos
        else if (GenPart.Nu.define 
              ( std::abs(GenPart().pdgId)==12 ||
                std::abs(GenPart().pdgId)==14 ||
                std::abs(GenPart().pdgId)==16 )) {
          if (debug>1) std::cout<<"Variables::define_genparticles_: nu start"<<std::endl;
          // Nus from W
          if (GenPart.NuFromW.define(std::abs(motherId)==24)) {
            // also from leptonic top
            GenPart.NuFromTop.define(std::abs(grandMotherId)==6);
          }
          // Nu from taus
          else if ( ( std::abs(GenPart().pdgId)==12 ||
                      std::abs(GenPart().pdgId)==14 ) &&
                    std::abs(motherId)==15) {
            // also lepton from W
            if (GenPart.NuFromW.define(std::abs(grandMotherId)==24)) {
              // also lepton from top
              GenPart.NuFromTop.define(std::abs(greatGrandMotherId)==6);
            }
          }
          if (debug>1) std::cout<<"Variables::define_genparticles_: nu ok"<<std::endl;
        }
        // Checking distance of q/g to photons
        // This is needed for flagging fragmentation photons
        if( GenPart().status==23 && ( std::abs(GenPart().pdgId)==1 ||
                                      std::abs(GenPart().pdgId)==2 ||
                                      std::abs(GenPart().pdgId)==3 ||
                                      std::abs(GenPart().pdgId)==4 ||
                                      std::abs(GenPart().pdgId)==5 ||
                                      std::abs(GenPart().pdgId)==6 ||
                                      std::abs(GenPart().pdgId)==21) ) {
          if (debug>1) std::cout<<"Variables::define_genparticles_: q/g start"<<std::endl;
          while (Photon.Loop()) if (Photon().genPartIdx!=-1 && Photon().genPartIdx<GenPart.n) {
            double dR_genqg = DeltaR(GenPart.v4(Photon().genPartIdx), GenPart.v4());
            if (dR_genqg<0.4) Photon().fromFrag = true;
          }
          if (debug>1) std::cout<<"Variables::define_genparticles_: q/g ok"<<std::endl;
        }
        if (debug>1) std::cout<<"Variables::define_genparticles_: end direct/frag photon matching"<<std::endl;

      }
      if (debug>1) std::cout<<"Variables::define_genparticles_: gen 1st loop "<<GenPart.i<<" end"<<std::endl;
    } // End GenPart 1st loop
    if (debug) std::cout<<"Variables::define_genparticles_: end gen 1st loop"<<std::endl;
    
    // 2nd Loop, to calculate gen matched objects
    while (GenPart.Loop()) {
      if (debug>1) std::cout<<"Variables::define_genparticles_: gen 2nd loop "<<GenPart.i<<" start"<<std::endl;
      float m = GenPart().mass;
      // mass is not stored in NanoAOD for particle masses below 10 GeV, need to look it up in PDG
      if      (std::abs(GenPart().pdgId)==1)  m = 0.005;
      else if (std::abs(GenPart().pdgId)==2)  m = 0.002;
      else if (std::abs(GenPart().pdgId)==3)  m = 0.095;
      else if (std::abs(GenPart().pdgId)==4)  m = 1.275;
      else if (std::abs(GenPart().pdgId)==5)  m = 4.180;
      else if (std::abs(GenPart().pdgId)==11) m = 0.0005;
      else if (std::abs(GenPart().pdgId)==13) m = 0.106;
      else if (std::abs(GenPart().pdgId)==15) m = 1.777;
      float pt = std::abs(GenPart().pt);
      float px = pt*std::cos(GenPart().phi);
      float py = pt*std::sin(GenPart().phi);
      float pz = pt*sinh(GenPart().eta);
      float e = m>=0 ? std::sqrt(px*px+py*py+pz*pz+m*m) : std::sqrt(std::max((px*px+py*py+pz*pz-m*m), (float)0.));
      GenPart.v4().SetPxPyPzE(px, py, pz, e);
      if (debug>1) std::cout<<"Variables::define_genparticles_: gen 2nd loop "<<GenPart.i<<" masses ok"<<std::endl;

      // Match to AK8
      double mindR = NOVAL_F;
      while (FatJet.JetAK8.Loop()) {
        double dR = DeltaR(GenPart.v4(), FatJet.JetAK8.v4());
        if (dR<0.8 && dR<mindR) {
          mindR = dR;
          GenPart().matchAK8 = true;
          GenPart().iMatchedAK8 = FatJet.JetAK8.iRef;
        }
      }
  
      // Generator leptons
      if (GenPart.Lepton.pass[GenPart.i]) {
        if (debug>1) std::cout<<"Variables::define_genparticles_: gen 2nd loop "<<GenPart.i<<" gen lep"<<std::endl;
        // Propagate lepton info to AK8
        while (FatJet.Loop()) if (DeltaR(GenPart.v4(), FatJet.v4())<0.8) {
          FatJet().matchGenLepton = true;
          if      (GenPart.LeptonFromTop.pass[GenPart.i]) FatJet().matchedGenLeptonMotherID = 6;
          else if (GenPart.LeptonFromW.  pass[GenPart.i]) FatJet().matchedGenLeptonMotherID = 24;
          else if (GenPart.LeptonFromZ  .pass[GenPart.i]) FatJet().matchedGenLeptonMotherID = 23;
          else if (GenPart.LeptonFromH  .pass[GenPart.i]) FatJet().matchedGenLeptonMotherID = 25;
          else                                            FatJet().matchedGenLeptonMotherID = GenPart().motherId;
        }

        // RECO object matching
        if (GenPart.Ele.define(std::abs(GenPart().pdgId)==11)) {
          while (Electron.Loop()) if (Electron().genPartIdx==GenPart.i) {
            Electron().matchGenEle = true;
            if (Electron.CBVeto  .pass[Electron.i]) GenPart().passEleCBVeto   = true;
            if (Electron.CBVetoNoIso  .pass[Electron.i]) GenPart().passEleCBVetoNoIso  = true;
            if (Electron.Veto  .pass[Electron.i]) GenPart().passLepVeto   = true;
            if (Electron.NoIso .pass[Electron.i]) GenPart().passLepNoIso  = true;
            if (Electron.NonIso.pass[Electron.i]) GenPart().passLepNonIso = true;
          }
        } else if (GenPart.Mu.define(std::abs(GenPart().pdgId)==13)) {
          while (Muon.Loop()) if (Muon().genPartIdx==GenPart.i) {
            Muon().matchGenMu = true;
            if (Muon.CBLooseNoIso  .pass[Muon.i]) GenPart().passMuoCBLooseNoIso = true;
            if (Muon.CBLoose  .pass[Muon.i]) GenPart().passMuoCBLoose = true;
            if (Muon.CBMediumNoIso  .pass[Muon.i]) GenPart().passMuoCBMediumNoIso = true;
            if (Muon.CBMedium .pass[Muon.i]) GenPart().passMuoCBMedium = true;
            if (Muon.Veto  .pass[Muon.i]) GenPart().passLepVeto   = true;
            if (Muon.NoIso .pass[Muon.i]) GenPart().passLepNoIso  = true;
            if (Muon.NonIso.pass[Muon.i]) GenPart().passLepNonIso = true;
          }
        } else if (GenPart.Tau.define(std::abs(GenPart().pdgId)==15)) {
          // Unlike other leptons, the taus decay to leptons or hadrons
          // It is only fair to match to the ID of the corresponding decay product
          // for ele/mu, the matching will be done above for the tau daughters
          // The reco taus are for the hadronic decay, so require Tau_genPartFlav==5
          while (Tau.Loop()) if (Tau().genPartIdx==GenPart.i) {
            Tau().matchGenTau = true;
            if (Tau.Veto.pass[Tau.i]&&Tau().genPartFlav==5) GenPart().passLepVeto = true;
          }
        }
        if (debug>1) std::cout<<"Variables::define_genparticles_: gen 2nd loop "<<GenPart.i<<" lep-reco object matching ok"<<std::endl;
  
        // leptons from the hard process
        if (GenPart.LeptonFromHardProcess.pass[GenPart.i]) {
          if (debug>1) std::cout<<"Variables::define_genparticles_: gen 2nd loop "<<GenPart.i<<" gen lep from hp"<<std::endl;
          if (GenPart.EleFromHardProcess.define(std::abs(GenPart().pdgId)==11))
            while (Electron.Loop()) if (Electron().genPartIdx==GenPart.i) Electron().matchGenEleFromHardProcess = 1;
          if (GenPart.MuFromHardProcess.define(std::abs(GenPart().pdgId)==13))
            while (Muon.Loop()) if (Muon().genPartIdx==GenPart.i) Muon().matchGenMuFromHardProcess = 1;
          GenPart.TauFromHardProcess.define(std::abs(GenPart().pdgId)==15);
        }
        if (debug>1) std::cout<<"Variables::define_genparticles_: gen 2nd loop "<<GenPart.i<<" gen lep from hard process ok"<<std::endl;
  
        // leptons from W decay
        if (GenPart.LeptonFromW.pass[GenPart.i]) {
          if (debug>1) std::cout<<"Variables::define_genparticles_: gen 2nd loop "<<GenPart.i<<" gen lep from W"<<std::endl;
          if (GenPart.EleFromW.define(std::abs(GenPart().pdgId)==11))
            while (Electron.Loop()) if (Electron().genPartIdx==GenPart.i) Electron().matchGenEleFromW = 1;
          if (GenPart.MuFromW.define(std::abs(GenPart().pdgId)==13))
            while (Muon.Loop()) if (Muon().genPartIdx==GenPart.i) Muon().matchGenMuFromW = 1;
          GenPart.TauFromW.define(std::abs(GenPart().pdgId)==15);
        }
        if (debug>1) std::cout<<"Variables::define_genparticles_: gen 2nd loop "<<GenPart.i<<" gen lep from W ok"<<std::endl;
  
        // leptons from Z decay
        if (GenPart.LeptonFromZ.pass[GenPart.i]) {
          if (debug>1) std::cout<<"Variables::define_genparticles_: gen 2nd loop "<<GenPart.i<<" gen lep from Z"<<std::endl;
          if (GenPart.EleFromZ.define(std::abs(GenPart().pdgId)==11))
            while (Electron.Loop()) if (Electron().genPartIdx==GenPart.i) Electron().matchGenEleFromZ = 1;
          if (GenPart.MuFromZ.define(std::abs(GenPart().pdgId)==13))
            while (Muon.Loop()) if (Muon().genPartIdx==GenPart.i) Muon().matchGenMuFromZ = 1;
          GenPart.TauFromZ.define(std::abs(GenPart().pdgId)==15);
        }
        if (debug>1) std::cout<<"Variables::define_genparticles_: gen 2nd loop "<<GenPart.i<<" gen lep from Z ok"<<std::endl;
  
        // leptons from Higgs decay
        if (GenPart.LeptonFromH.pass[GenPart.i]) {
          if (debug>1) std::cout<<"Variables::define_genparticles_: gen 2nd loop "<<GenPart.i<<" gen lep from H"<<std::endl;
          if (GenPart.EleFromH.define(std::abs(GenPart().pdgId)==11))
            while (Electron.Loop()) if (Electron().genPartIdx==GenPart.i) Electron().matchGenEleFromH = 1;
          if (GenPart.MuFromH.define(std::abs(GenPart().pdgId)==13))
            while (Muon.Loop()) if (Muon().genPartIdx==GenPart.i) Muon().matchGenMuFromH = 1;
          GenPart.TauFromH.define(std::abs(GenPart().pdgId)==15);
        }
        if (debug>1) std::cout<<"Variables::define_genparticles_: gen 2nd loop "<<GenPart.i<<" gen lep from Higgs ok"<<std::endl;
  
        // leptons from top decay
        if (GenPart.LeptonFromTop.pass[GenPart.i]) {
          if (debug>1) std::cout<<"Variables::define_genparticles_: gen 2nd loop "<<GenPart.i<<" gen lep from top"<<std::endl;
          if (GenPart.EleFromTop.define(std::abs(GenPart().pdgId)==11))
            while (Electron.Loop()) if (Electron().genPartIdx==GenPart.i) Electron().matchGenEleFromTop = 1;
          if (GenPart.MuFromTop.define(std::abs(GenPart().pdgId)==13))
            while (Muon.Loop()) if (Muon().genPartIdx==GenPart.i) Muon().matchGenMuFromTop = 1;
          GenPart.TauFromTop.define(std::abs(GenPart().pdgId)==15);
        }
        if (debug>1) std::cout<<"Variables::define_genparticles_: gen 2nd loop "<<GenPart.i<<" gen lep from top ok"<<std::endl;
      }
      if (debug>1) std::cout<<"Variables::define_genparticles_: gen 2nd loop "<<GenPart.i<<" gen leptons ok"<<std::endl;
  
      // Generator neutrinos
      if (GenPart.Nu.pass[GenPart.i]) {
        if (debug>1) std::cout<<"Variables::define_genparticles_: gen 2nd loop "<<GenPart.i<<" gen nu"<<std::endl;
        if (GenPart.EleNu.define(std::abs(GenPart().pdgId)==12)) {
          GenPart.EleNuFromW  .define(GenPart.NuFromW  .pass[GenPart.i]);
          GenPart.EleNuFromTop.define(GenPart.NuFromTop.pass[GenPart.i]);
          double dR_min = NOVAL_F;
          size_t iMatch = -1;
          while (Electron.Loop()) {
            double dR = DeltaR(Electron().nu, GenPart.v4());
            if (dR<0.4 && dR<dR_min) { iMatch = Electron.i; dR_min = dR; }
          }
          if (dR_min<0.4) {
            Electron(iMatch).nuMatchGenEleNu = 1;
            if (GenPart.NuFromW  .pass[GenPart.i]) Electron(iMatch).nuMatchGenEleNuFromW   = 1;
            if (GenPart.NuFromTop.pass[GenPart.i]) Electron(iMatch).nuMatchGenEleNuFromTop = 1;
          }
        }
        if (GenPart.MuNu.define(std::abs(GenPart().pdgId)==14)) {
          GenPart.MuNuFromW  .define(GenPart.NuFromW  .pass[GenPart.i]);
          GenPart.MuNuFromTop.define(GenPart.NuFromTop.pass[GenPart.i]);
          double dR_min = NOVAL_F;
          size_t iMatch = -1;
          while (Muon.Loop()) {
            double dR = DeltaR(Muon().nu, GenPart.v4());
            if (dR<0.4 && dR<dR_min) { iMatch = Muon.i; dR_min = dR; }
          }
          if (dR_min<0.4) {
            Muon(iMatch).nuMatchGenMuNu = 1;
            if (GenPart.NuFromW  .pass[GenPart.i]) Muon(iMatch).nuMatchGenMuNuFromW   = 1;
            if (GenPart.NuFromTop.pass[GenPart.i]) Muon(iMatch).nuMatchGenMuNuFromTop = 1;
          }
        }
        if (GenPart.TauNu.define(std::abs(GenPart().pdgId)==16)) {
          GenPart.TauNuFromW  .define(GenPart.NuFromW  .pass[GenPart.i]);
          GenPart.TauNuFromTop.define(GenPart.NuFromTop.pass[GenPart.i]);
        }
        if (debug>1) std::cout<<"Variables::define_genparticles_: gen 2nd loop "<<GenPart.i<<" gen nu ok"<<std::endl;
      }
      if (debug>1) std::cout<<"Variables::define_genparticles_: gen 2nd loop "<<GenPart.i<<" gen neutrinos ok"<<std::endl;
  
      // Prompt photons (from gq mother)
      if(GenPart.PromptPhoton.define(std::abs(GenPart().pdgId) == 22 && 
                                     ( std::abs(GenPart().motherId)==1 ||
                                       std::abs(GenPart().motherId)==2 ||
                                       std::abs(GenPart().motherId)==3 ||
                                       std::abs(GenPart().motherId)==4 ||
                                       std::abs(GenPart().motherId)==5 ||
                                       std::abs(GenPart().motherId)==6 ||
                                       std::abs(GenPart().motherId)==21 ))) {
        // RECO matching
        while (Photon.Loop()) if (Photon().genPartIdx == GenPart.i) {
          Photon().matchGenFake         = 0;
          Photon().matchGenPrompt       = 1;
          GenPart.PromptDirectPhoton.define((Photon().matchGenPromptDirect = !Photon().fromFrag));
          GenPart.PromptFragPhoton  .define((Photon().matchGenPromptFrag   =  Photon().fromFrag));
        }
      }
      if (debug>1) std::cout<<"Variables::define_genparticles_: gen 2nd loop "<<GenPart.i<<" prompt photons ok"<<std::endl;
      
      // Final non-decayed state particles
      if ((GenPart().FinalState = GenPart().status==1 || GenPart().status==2 || (GenPart().LastCopyCand && GenPart().NoSameDaughter))) {
        // gen bs
        GenPart.b.define(std::abs(GenPart().pdgId)==5);
        if (debug>1) std::cout<<"Variables::define_genparticles_: gen 2nd loop "<<GenPart.i<<" final state bs ok"<<std::endl;
        
        // gen Ws
        // Consider hadronically decaying Ws and match them to AK8 jet
        if (GenPart.W.define(std::abs(GenPart().pdgId)==24)) {
          if (GenPart.HadW.define(!GenPart.LepW.define(GenPart.LepW.pass[GenPart.i]))) {
            if (GenPart().matchAK8) {
              FatJet(GenPart().iMatchedAK8).matchGenHadW = true;
              GenPart().passHadWTag = FatJet.HadW.pass[GenPart().iMatchedAK8];
            }
          }
        }
        if (debug>1) std::cout<<"Variables::define_genparticles_: gen 2nd loop "<<GenPart.i<<" final state Ws ok"<<std::endl;
        
        // gen Zs
        // Consider hadronically decaying Zs and match them to AK8 jet
        if (GenPart.Z.define(std::abs(GenPart().pdgId)==23)) {
          if (GenPart.HadZ.define(!GenPart.LepZ.define(GenPart.LepZ.pass[GenPart.i]))) {
            if (GenPart().matchAK8) {
              FatJet(GenPart().iMatchedAK8).matchGenHadZ = true;
              GenPart().passHadZTag = FatJet.HadZ.pass[GenPart().iMatchedAK8];
            }
          }
        }
        if (debug>1) std::cout<<"Variables::define_genparticles_: gen 2nd loop "<<GenPart.i<<" final state Zs ok"<<std::endl;
        
        // gen Higgs
        // Consider hadronically decaying Higgs and match it to AK8 jet
        if (GenPart.H.define(std::abs(GenPart().pdgId)==25)) {
          if (GenPart.HadH.define(!GenPart.LepH.define(GenPart.LepH.pass[GenPart.i]))) {
            if (GenPart().matchAK8) {
              FatJet(GenPart().iMatchedAK8).matchGenHadH = true;
              GenPart().passHadHTag = FatJet.HadH.pass[GenPart().iMatchedAK8];
            }
          }
        }
        if (debug>1) std::cout<<"Variables::define_genparticles_: gen 2nd loop "<<GenPart.i<<" final state Higgs' k"<<std::endl;
        
        // gen tops
        if (GenPart.Top.define(std::abs(GenPart().pdgId)==6)) {
          if (debug>1) std::cout<<"Variables::define_genparticles_: gen 2nd loop "<<GenPart.i<<" gentop start"<<std::endl;
          if (GenPart.LepTop.define(GenPart.LepTop.pass[GenPart.i])) {
            if (debug>1) std::cout<<"Variables::define_genparticles_: gen 2nd loop "<<GenPart.i<<" gentop lep start"<<std::endl;
            if (GenPart().matchAK8) {
              FatJet(GenPart().iMatchedAK8).matchGenLepTop = true;
              GenPart().passLepTopTag = FatJet.LepTop.pass[GenPart().iMatchedAK8];
            }
            if (debug>1) std::cout<<"Variables::define_genparticles_: gen 2nd loop "<<GenPart.i<<" gentop lep ok"<<std::endl;
          } else {
            if (debug>1) std::cout<<"Variables::define_genparticles_: gen 2nd loop "<<GenPart.i<<" gentop had start"<<std::endl;
            GenPart.HadTop.define(!GenPart.LepTop.pass[GenPart.i]);
            GenPart.LepTop.define(1);
            if (GenPart().matchAK8) {
              FatJet(GenPart().iMatchedAK8).matchGenHadTop = true;
              GenPart().passHadTopTag = FatJet.HadTop.pass[GenPart().iMatchedAK8];
            }
            if (debug>1) std::cout<<"Variables::define_genparticles_: gen 2nd loop "<<GenPart.i<<" gentop had ok"<<std::endl;
          }
        }
        if (debug>1) std::cout<<"Variables::define_genparticles_: gen 2nd loop "<<GenPart.i<<" final state tops ok"<<std::endl;
      } // End FinalState

    } // End GenPart 2nd loop
    if (debug) std::cout<<"Variables::define_genparticles_: end gen 2nd loop"<<std::endl;

  }

  std::vector<LorentzVector> saved_megajets;
  std::vector<LorentzVector> saved_megajets_nophoton;
  std::vector<LorentzVector> newmegajets;
  std::vector<LorentzVector> newmegajets_nophoton;
  std::vector<LorentzVector> saved_newmegajets;
  std::vector<LorentzVector> saved_newmegajets_nophoton;
  std::vector<LorentzVector> saved_hemJets;
  
  void define_event_(const unsigned int& syst_index, int debug = 0) {
    resetEventData();

    if (debug) std::cout<<"Variables::define_event_: start"<<std::endl;

    // Lepton selections combined  
    nLepVetoNoIso = Electron.VetoNoIso.n + Muon.VetoNoIso.n;
    nLepVeto      = Electron.Veto.n      + Muon.Veto.n;
    nLepSelect    = Electron.Select.n    + Muon.Select.n;
    nLepTight     = Electron.Tight.n     + Muon.Tight.n;
    nLepNoIso     = Electron.NoIso.n     + Muon.NoIso.n;
    nLepNonIso    = Electron.NonIso.n    + Muon.NonIso.n;

    std::vector<LorentzVector> noniso_leptons;
    while (Electron.NonIso.Loop()) noniso_leptons.push_back(Electron.NonIso.v4());
    while (Muon.NonIso.Loop())     noniso_leptons.push_back(Muon.NonIso.v4());

    if (FatJet.LepTop.n>0) lepNeutrinoDR = FatJet.LepTop(0).lepNonIsoNuDR;

    // MET variables
    double MET_px = MET_pt*std::cos(MET_phi);
    double MET_py = MET_pt*std::sin(MET_phi);
    MET.SetXYZ(MET_px, MET_py, 0);
    if (debug) std::cout<<"Variables::define_event_: end MET variables"<<std::endl;    
    
    // -----------------------------------------------
    //              Lepton/Phton + MET

    // MT (l, MET), MET + 1l
    MET_1vl.SetXYZ(MET_px, MET_py, 0);
    if (Electron.Veto.n>=1&&Muon.Veto.n==0) {
      MET_1vl += Vector3(Electron.Veto.v4(0).Px(),Electron.Veto.v4(0).Py(),0);
      MT_lepVeto = sqrt( 2*Electron.Veto(0).pt*MET_pt * (1 - std::cos(MET_phi-Electron.Veto(0).phi)) );
    } else if (Electron.Veto.n==0&&Muon.Veto.n>=1) {
      MET_1vl += Vector3(Muon.Veto.v4(0).Px(),Muon.Veto.v4(0).Py(),0);
      MT_lepVeto = sqrt( 2*Muon.    Veto(0).pt*MET_pt * (1 - std::cos(MET_phi-Muon.    Veto(0).phi)) );
    }
    MET_1l.SetXYZ(MET_px, MET_py, 0);
    if (Electron.Select.n>=1&&Muon.Select.n==0) {
      MET_1l += Vector3(Electron.Select.v4(0).Px(),Electron.Select.v4(0).Py(),0);
      MT = sqrt( 2*Electron.Select(0).pt*MET_pt * (1 - std::cos(MET_phi-Electron.Select(0).phi)) );
    } else if (Electron.Select.n==0&&Muon.Select.n>=1) {
      MET_1l += Vector3(Muon.Select.v4(0).Px(),Muon.Select.v4(0).Py(),0);
      MT = sqrt( 2*Muon.    Select(0).pt*MET_pt * (1 - std::cos(MET_phi-Muon.    Select(0).phi)) );
    }
    if (Electron.Tight.n>=1&&Muon.Tight.n==0) {
      MT_lepTight = sqrt( 2*Electron.Tight(0).pt*MET_pt * (1 - std::cos(MET_phi-Electron.Tight(0).phi)) );
    } else if (Electron.Tight.n==0&&Muon.Tight.n>=1) {
      MT_lepTight = sqrt( 2*Muon.    Tight(0).pt*MET_pt * (1 - std::cos(MET_phi-Muon.    Tight(0).phi)) );
    }
    if (Electron.NonIso.n>=1&&Muon.NonIso.n==0) {
      MT_lepNonIso = sqrt( 2*Electron.NonIso(0).pt*MET_pt * (1 - std::cos(MET_phi-Electron.NonIso(0).phi)) );
    } else if (Electron.NonIso.n==0&&Muon.NonIso.n>=1) {
      MT_lepNonIso = sqrt( 2*Muon.    NonIso(0).pt*MET_pt * (1 - std::cos(MET_phi-Muon.    NonIso(0).phi)) );
    }
    // M(2l), dPhi(2l, MET), MET + 2l
    MET_2l.SetXYZ(MET_px, MET_py, 0);
    LorentzVector lep_pair;
    if (Electron.Select.n==2&&Muon.Veto.n==0) {
      MET_2l += Vector3(Electron.Select.v4(0).Px(), Electron.Select.v4(0).Py(), 0);
      MET_2l += Vector3(Electron.Select.v4(1).Px(), Electron.Select.v4(1).Py(), 0);
      lep_pair = Electron.Select.v4(0)+Electron.Select.v4(1);
      M_2l = lep_pair.M();
      dPhi_2l_met = std::abs(DeltaPhi(lep_pair.Phi(), MET_phi));
    } else if (Electron.Veto.n==0&&Muon.Select.n==2) {
      MET_2l += Vector3(Muon.Select.v4(0).Px(), Muon.Select.v4(0).Py(), 0);
      MET_2l += Vector3(Muon.Select.v4(1).Px(), Muon.Select.v4(1).Py(), 0);
      lep_pair = Muon.Select.v4(0)+Muon.Select.v4(1);
      M_2l = lep_pair.M();
      dPhi_2l_met = std::abs(DeltaPhi(lep_pair.Phi(), MET_phi));
    }
    // Add randomly one of the two veto leptons to MET
    MET_dilep.SetXYZ(MET_px, MET_py, 0);
    if (Electron.Veto.n==2&&Muon.Veto.n==0) {
      if (rnd_.Rndm()<0.5) {
        MET_dilep += Vector3(Electron.Veto.v4(0).Px(), Electron.Veto.v4(0).Py(), 0);
      } else {
        MET_dilep += Vector3(Electron.Veto.v4(1).Px(), Electron.Veto.v4(1).Py(), 0);
      }
    } else if (Electron.Veto.n==1&&Muon.Veto.n==1) {
      if (rnd_.Rndm()<0.5) {
        MET_dilep += Vector3(Electron.Veto.v4(0).Px(), Electron.Veto.v4(0).Py(), 0);
      } else {
        MET_dilep += Vector3(Muon.Veto.v4(0).Px(), Muon.Veto.v4(0).Py(), 0);
      }
    } else if (Electron.Veto.n==0&&Muon.Veto.n==2) {
      if (rnd_.Rndm()<0.5) {
        MET_dilep += Vector3(Muon.Veto.v4(0).Px(), Muon.Veto.v4(0).Py(), 0);
      } else {
        MET_dilep += Vector3(Muon.Veto.v4(1).Px(), Muon.Veto.v4(1).Py(), 0);
      }
    }

    // Add the (pre-)selected/fake photon to MET
    MET_pho.SetXYZ(MET_px, MET_py, 0);
    MET_fakepho.SetXYZ(MET_px, MET_py, 0);
    if      (Photon.Select.n==1)    MET_pho     += Vector3(Photon.Select   .v4(0).Px(), Photon.Select   .v4(0).Py(), 0);
    else if (Photon.PreSelect.n==1) MET_pho     += Vector3(Photon.PreSelect.v4(0).Px(), Photon.PreSelect.v4(0).Py(), 0);
    if      (Photon.Fake.n==1)      MET_fakepho += Vector3(Photon.Fake     .v4(0).Px(), Photon.Fake     .v4(0).Py(), 0);
    if (debug) std::cout<<"Variables::define_event_variables_: end met + lep/pho variables"<<std::endl;

    // -----------------------------------------------
    //                    Jet/MET

    // Online HT
    while (Jet.Jet.Loop()) {
      if (debug>1) std::cout<<"Variables::define_event_variables_: Jet "<<Jet.Jet.i<<" min dphi start"<<std::endl;
      // AK4 HT
      AK4_Ht += Jet.Jet().pt;

      // minDeltaPhi
      if (Jet.Jet.i<4) {
        double dphi = std::abs(DeltaPhi(MET_phi, Jet.Jet().phi));
        if (dphi<minDeltaPhi) minDeltaPhi = dphi;
        // with added (veto) lepton
        double dphi_met1vl = std::abs(DeltaPhi(MET_1vl.Phi(), Jet.Jet().phi));
        if (dphi_met1vl<minDeltaPhi_1vl) minDeltaPhi_1vl = dphi_met1vl;
        double dphi_met1l = std::abs(DeltaPhi(MET_1l.Phi(), Jet.Jet().phi));
        if (dphi_met1l<minDeltaPhi_1l) minDeltaPhi_1l = dphi_met1l;
        // with added lepton pair
        double dphi_met2l = std::abs(DeltaPhi(MET_2l.Phi(), Jet.Jet().phi));
        if (dphi_met2l<minDeltaPhi_2l) minDeltaPhi_2l = dphi_met2l;
        // with added photon
        double dphi_metpho = std::abs(DeltaPhi(MET_pho.Phi(), Jet.Jet().phi));
        if (dphi_metpho<minDeltaPhi_pho) minDeltaPhi_pho = dphi_metpho;
        // jet lep-pair angle
        if (M_2l!=-NOVAL_F) {
          double dphi_2l = std::abs(DeltaPhi(lep_pair.Phi(), Jet.Jet().phi));
          if (dphi_2l<dPhi_2l_jet) dPhi_2l_jet = dphi_2l;
        }
      }

      // ISR jets counting
      // Taken from:
      // https://github.com/manuelfs/babymaker/blob/0136340602ee28caab14e3f6b064d1db81544a0a/bmaker/plugins/bmaker_full.cc#L1268-L1295
      if (!isData) {
        bool matched = false;
        while (GenPart.Loop()) {
          if (matched) break;
          if (GenPart().status!=23 || std::abs(GenPart().pdgId)>5) continue;
          
          int momid = -NOVAL_I;
          int iMother = GenPart().genPartIdxMother;
          if (iMother!=-1) momid = std::abs(GenPart(iMother).pdgId);
          if(!(momid==6 || momid==23 || momid==24 || momid==25 || momid>1e6)) continue; 
          
          //check against daughter in case of hard initial splitting
          for (const auto& iDau : GenPart().iDaughters) {
            double dR = DeltaR(Jet.Jet.v4(), GenPart.v4(iDau));
            if(dR<0.3){
              matched = true;
              break;
            }
          }
        } // Loop over MC particles
        if(!matched) ++nJetISR;
      }

      if (debug>1) std::cout<<"Variables::define_event_variables_: Jet "<<Jet.Jet.i<<" min dphi ok"<<std::endl;
    }
    
    while (Jet.Loop()) {
      if (Jet.Jet.pass[Jet.i]) {
        if (Jet.MuonJet.define(Jet().pt>=200 && Jet().muEF>=0.5)) {
          double dPhi = std::abs(DeltaPhi(MET_phi, Jet().phi));
          if (dPhi>=dPhiMuonJetMET) dPhiMuonJetMET = dPhi;
        }
      }
    }
    
    if (debug) std::cout<<"Variables::define_event_variables_: end"<<std::endl;
  } // End define_event_variables_


}; // End Variables class definition

#endif // End header guard
