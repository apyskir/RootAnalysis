#ifndef RootAnalysis_HTTAnalyzer_H
#define RootAnalysis_HTTAnalyzer_H

#include <string>
#include <vector>
#include <map>

#include "ObjectMessenger.h"
#include "EventProxyBase.h"
#include "EventProxyHTT.h"

#include "strbitset.h"
#include "TDirectory.h"

//ROOT includes
#include "TTree.h"
#include "TList.h"

#include "Analyzer.h"
#include "ChannelSpecifics.h"
#include "AnalysisEnums.h"

class HTTHistograms;

class TH1F;
class TH2F;
class TH3F;
class TLorentzVector;

class HTTAnalyzer: public Analyzer{

  friend class ChannelSpecifics;
  friend class MuTauSpecifics;
  friend class TauTauSpecifics;
  friend class MuMuSpecifics;

 public:

  HTTAnalyzer(const std::string & aName, const std::string & aDecayMode = "None");

  virtual ~HTTAnalyzer();

  ///Initialize the analyzer
  virtual void initialize(TDirectory* aDir,
			  pat::strbitset *aSelections);

  virtual bool analyze(const EventProxyBase& iEvent);

  virtual bool analyze(const EventProxyBase& iEvent, ObjectMessenger *aMessenger){return analyze(iEvent); }

  virtual void finalize();

  virtual void clear(){;};

  virtual void addBranch(TTree *);
  
  void clearTTreeVariables();

  Analyzer* clone() const;

  bool filter() const{ return filterEvent_;};

  void setAnalysisObjects(const EventProxyHTT & myEventProxy);

  ///Check it the event passes given category selections.
  bool passCategory(unsigned int iCategory);

  ///Return human readable sample name (Data, WJets, etc).
  std::string getSampleName(const EventProxyHTT & myEventProxy);

  ///Return human readable sample name (Data, WJets, etc).
  ///Make the methos static, so other modules can use it.
  ///Method used when sample coding in TTree is not present.
  ///In this case a ROOT file name is used to decode the sample type.
  std::string getSampleNameFromFileName(const EventProxyHTT & myEventProxy);

  ///Return sample name for DY. Name encoded jet bin, and decay mode.
  std::string getDYSampleName(const EventProxyHTT & myEventProxy);

  //Return name sample name suffix for different particles matched to reconstructed tau
  std::string getMatchingName(const EventProxyHTT & myEventProxy);

  ///Return pileup reweighting weight.
  ///Weight is calculatedon fly using the ration of nPU
  ///histograms for data and analyased sample.
  float getPUWeight(const EventProxyHTT & myEventProxy);

  ///Return event weight for systematic effects
  ///implemented by a global event weight.
  float getSystWeight(const HTTAnalysis::sysEffects & aSystEffect=HTTAnalysis::NOMINAL);

  ///Fill pulls between generator and various reco vertices.
  bool fillVertices(const std::string & sysType, float eventWeight);

  ///Return generator weight. Most samples have large values of weights
  ///which are constant up to + or - sign. We normalise those weights to +-1.
  float getGenWeight(const EventProxyHTT & myEventProxy);

  ///Fill histograms for all control plots.
  ///Histogram names will end with hNameSuffix
  void fillControlHistos(const std::string & hNameSuffix, float eventWeight,
			 const HTTAnalysis::sysEffects & aSystEffect=HTTAnalysis::NOMINAL);


  ///Fill histograms with cos(phi), where phi is the decay
  ///between tau decay planes. Method used for reconstructed
  ///mu+tau_h mode
  void fillDecayPlaneAngle(const std::string & hNameSuffix, float eventWeight,
  const HTTAnalysis::sysEffects & aSystEffect=HTTAnalysis::NOMINAL);

  ///Fill histograms with cos(phi), where phi is the decay
  ///between tau decay planes. Method used for
  ///generator level taus for all decay modes.
  void fillGenDecayPlaneAngle(const std::string & hNameSuffix, float eventWeight);

  ///Calculate angle between tau decay planes (first element of pair)
  //and angle betwee decay products (second element of pair)
  std::pair<float,float> angleBetweenPlanes(const TLorentzVector& tau1, const TLorentzVector& tau1Daughter,
					    const TLorentzVector& tau2, const TLorentzVector& tau2Daughter,
					    bool sgn=true);

  ///Get jets separated by deltaR from tau an muon.
  std::vector<HTTParticle> getSeparatedJets(const EventProxyHTT & myEventProxy,
					    float deltaR);

 protected:

  pat::strbitset *mySelections_;

  ///Types of the selection flow
  std::vector<std::string> selectionFlavours_;

 private:

  void getPreselectionEff(const EventProxyHTT & myEventProxy);

  void setHistos(HTTHistograms *histos) { myHistos_ = histos;};

  ///Parts of code specific to give decay channel.
  ///In particular category and object selection.
  ChannelSpecifics *myChannelSpecifics;

  ///Histograms storage.
  HTTHistograms *myHistos_;

  ///ROOT file with PU histogram
  TFile *puDataFile_, *puMCFile_;

  ///ROOT file containing current TTree
  TFile *ntupleFile_;

  ///Histogram with event counts filled during preselection step.
  TH1F *hStatsFromFile;

  ///Vector of PU histograms for MC samples
  std::vector<TH1F*> hPUVec_;

  //should this HTTAnalyzer be able to filter events
  bool filterEvent_;

  ///Map from file name to sample name.
  std::map<std::string, std::string> fileName2sampleName;

  ///Reconstructed objects selected for given event.
  HTTEvent aEvent;
  HTTPair aPair;
  std::string sampleName;

  HTTParticle aLeg2, aLeg1, aMET;
  HTTParticle aGenLeg1, aGenLeg2;
  HTTParticle aJet1, aJet2, aBJet1;
  std::vector<HTTParticle> aSeparatedJets;
  int nJets30;
  int nJetsInGap30;
  int nBJets;
  std::vector<bool> categoryDecisions;
  unsigned int myNumberOfCategories;

  //cut on nPCA
  float nPCAMin_;
  
  //TTree variables
  Bool_t tree_isData;
  //event ID variables
  ULong64_t tree_run;
  ULong64_t tree_lumi;
  ULong64_t tree_evt;
  //Pielup
  ULong64_t tree_npv;
  Float_t tree_npu;
  Float_t tree_rho;
  Float_t tree_puweight;
  //Weights
  Float_t tree_trackingweight_1;
  Float_t tree_trackingweight_2;
  Float_t tree_effweight;
  Float_t tree_weight;
  //Leg 1 (leading tau for tt; electon for et,em; muon for mt)
  Float_t tree_pt_1;
  Float_t tree_phi_1;
  Float_t tree_eta_1;
  Float_t tree_m_1;
  Float_t tree_q_1;
  Float_t tree_d0_1;
  Float_t tree_dZ_1;
  Float_t tree_mt_1;
  Float_t tree_pfmt_1;
  Float_t tree_puppimt_1;
  Float_t tree_iso_1;
  Float_t tree_id_e_mva_nt_loose_1;
  Int_t tree_gen_match_1;
  Float_t tree_againstElectronLooseMVA6_1;
  Float_t tree_againstElectronMediumMVA6_1;
  Float_t tree_againstElectronTightMVA6_1;
  Float_t tree_againstElectronVLooseMVA6_1;
  Float_t tree_againstElectronVTightMVA6_1;
  Float_t tree_againstMuonLoose3_1;
  Float_t tree_againstMuonTight3_1;
  Float_t tree_byCombinedIsolationDeltaBetaCorrRaw3Hits_1;
  Float_t tree_byIsolationMVA3newDMwoLTraw_1;
  Float_t tree_byIsolationMVA3oldDMwoLTraw_1;
  Float_t tree_byIsolationMVA3newDMwLTraw_1;
  Float_t tree_byIsolationMVA3oldDMwLTraw_1;
  Float_t tree_chargedIsoPtSum_1;
  Bool_t tree_decayModeFindingOldDMs_1;
  Float_t tree_neutralIsoPtSum_1;
  Float_t tree_puCorrPtSum_1;
  Float_t tree_trigweight_1;
  Float_t tree_idisoweight_1;
  Int_t tree_tau_decay_mode_1;
  Bool_t tree_mva_olddm_medium_1;
  Bool_t tree_mva_olddm_tight_1;
  Bool_t tree_mva_olddm_vtight_1;
  //Leg 2 (trailing tau for tt, tau for et,mt, muon for em)
  Float_t tree_pt_2;
  Float_t tree_phi_2;
  Float_t tree_eta_2;
  Float_t tree_m_2;
  Float_t tree_q_2;
  Float_t tree_d0_2;
  Float_t tree_dZ_2;
  Float_t tree_mt_2;
  Float_t tree_pfmt_2;
  Float_t tree_puppimt_2;
  Float_t tree_iso_2;
  Float_t tree_id_e_mva_nt_loose_2;
  Int_t tree_gen_match_2;
  Float_t tree_againstElectronLooseMVA6_2;
  Float_t tree_againstElectronMediumMVA6_2;
  Float_t tree_againstElectronTightMVA6_2;
  Float_t tree_againstElectronVLooseMVA6_2;
  Float_t tree_againstElectronVTightMVA6_2;
  Float_t tree_againstMuonLoose3_2;
  Float_t tree_againstMuonTight3_2;
  Float_t tree_byCombinedIsolationDeltaBetaCorrRaw3Hits_2;
  Float_t tree_byIsolationMVA3newDMwoLTraw_2;
  Float_t tree_byIsolationMVA3oldDMwoLTraw_2;
  Float_t tree_byIsolationMVA3newDMwLTraw_2;
  Float_t tree_byIsolationMVA3oldDMwLTraw_2;
  Float_t tree_chargedIsoPtSum_2;
  Bool_t tree_decayModeFindingOldDMs_2;
  Float_t tree_neutralIsoPtSum_2;
  Float_t tree_puCorrPtSum_2;
  Float_t tree_trigweight_2;
  Float_t tree_idisoweight_2;
  Int_t tree_tau_decay_mode_2;
  Bool_t tree_mva_olddm_medium_2;
  Bool_t tree_mva_olddm_tight_2;
  Bool_t tree_mva_olddm_vtight_2;
  //di-tau system
  Float_t tree_pt_tt;
  Float_t tree_mt_tot;
  Float_t tree_m_vis;
  Float_t tree_m_sv;
  Float_t tree_mt_sv;
  Float_t tree_pt_sv;
  Float_t tree_eta_sv;
  Float_t tree_phi_sv;
  Bool_t tree_os;
  //MET
  Float_t tree_met;
  Float_t tree_metphi;
  Float_t tree_pfmet;
  Float_t tree_pfmetphi;
  Float_t tree_puppimet;
  Float_t tree_puppimetphi;
  Float_t tree_mvamet;
  Float_t tree_mvametphi;
  Float_t tree_pzetavis;
  Float_t tree_pzetamiss;
  Float_t tree_pfpzetamiss;
  Float_t tree_puppipzetamiss;
  Float_t tree_mvacov00;
  Float_t tree_mvacov01;
  Float_t tree_mvacov10;
  Float_t tree_mvacov11;
  Float_t tree_metcov00;
  Float_t tree_metcov01;
  Float_t tree_metcov10;
  Float_t tree_metcov11;
  //VBF system
  Float_t tree_mjj;
  Float_t tree_jdeta;
  Float_t tree_njetingap;
  Float_t tree_njetingap20;
  Float_t tree_jdphi;
  //additional jets
  Float_t tree_nbtag;
  Float_t tree_njets;
  Float_t tree_njetspt20;
  //leading jet sorted by pt
  Float_t tree_jpt_1;
  Float_t tree_jeta_1;
  Float_t tree_jphi_1;
  Float_t tree_jrawf_1;
  Float_t tree_jmva_1;
  //trailing jet sorted by pt
  Float_t tree_jpt_2;
  Float_t tree_jeta_2;
  Float_t tree_jphi_2;
  Float_t tree_jrawf_2;
  Float_t tree_jmva_2;
  //leading b-jet sorted by pt
  Float_t tree_bpt_1;
  Float_t tree_beta_1;
  Float_t tree_bphi_1;
  Float_t tree_brawf_1;
  Float_t tree_bmva_1;
  Float_t tree_bcsv_1;
  //trailing b-jet sorted by pt
  Float_t tree_bpt_2;
  Float_t tree_beta_2;
  Float_t tree_bphi_2;
  Float_t tree_brawf_2;
  Float_t tree_bmva_2;
  Float_t tree_bcsv_2;
  //Extra lepton vetos
  Bool_t tree_dilepton_veto = 0;
  Bool_t tree_extraelec_veto = 0;
  Bool_t tree_extramuon_veto = 0;
  //Trigger
  Bool_t tree_trg_singlemuon = 0;
  Bool_t tree_trg_singleelectron = 0;
  Bool_t tree_trg_singletau_1 = 0;
  Bool_t tree_trg_singletau_2 = 0;
  Bool_t tree_trg_doubletau = 0;
  Bool_t tree_trg_muonelectron = 0;

};

#endif
