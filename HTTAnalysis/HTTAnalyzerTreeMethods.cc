#include <sstream>
#include <bitset>

#include "HTTAnalyzer.h"

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::clearTTreeVariables(){

  tree_evt = -999;
  tree_sampleNumber = -999;
  tree_weight = -999;
  //Pielup
  tree_npv = -999;
  tree_npu = -999;
  tree_rho = -999;
  tree_puweight = -999;
  //Leg 1 (leading tau for tt; electon for et,em; muon for mt)
  tree_pt_1 = -999;
  tree_phi_1 = -999;
  tree_eta_1 = -999;
  tree_m_1 = -999;
  tree_q_1 = -999;
  tree_d0_1 = -999;
  tree_dZ_1 = -999;
  tree_mt_1 = -999;
  tree_pfmt_1 = -999;
  tree_puppimt_1 = -999;
  tree_iso_1 = -999;
  tree_id_e_mva_nt_loose_1 = -999;
  tree_gen_match_1 = -999;
  tree_againstElectronLooseMVA6_1 = -999;
  tree_againstElectronMediumMVA6_1 = -999;
  tree_againstElectronTightMVA6_1 = -999;
  tree_againstElectronVLooseMVA6_1 = -999;
  tree_againstElectronVTightMVA6_1 = -999;
  tree_againstMuonLoose3_1 = -999;
  tree_againstMuonTight3_1 = -999;
  tree_byCombinedIsolationDeltaBetaCorrRaw3Hits_1 = -999;
  tree_byIsolationMVA3newDMwoLTraw_1 = -999;
  tree_byIsolationMVA3oldDMwoLTraw_1 = -999;
  tree_byIsolationMVA3newDMwLTraw_1 = -999;
  tree_byIsolationMVA3oldDMwLTraw_1 = -999;
  tree_chargedIsoPtSum_1 = -999;
  tree_decayModeFindingOldDMs_1 = 0;
  tree_neutralIsoPtSum_1 = -999;
  tree_puCorrPtSum_1 = -999;
  tree_trigweight_1 = -999;
  tree_idisoweight_1 = -999;
  tree_tau_decay_mode_1 = -999 ;
  tree_mva_olddm_medium_1 = 0;
  tree_mva_olddm_tight_1 = 0;
  tree_mva_olddm_vtight_1 = 0;
  //Leg 2 (trailing tau for tt; tau for et,mt; muon for em)
  tree_pt_2 = -999;
  tree_phi_2 = -999;
  tree_eta_2 = -999;
  tree_m_2 = -999;
  tree_q_2 = -999;
  tree_d0_2 = -999;
  tree_dZ_2 = -999;
  tree_mt_2 = -999;
  tree_pfmt_2 = -999;
  tree_puppimt_2 = -999;
  tree_iso_2 = -999;
  tree_id_e_mva_nt_loose_2 = -999;
  tree_gen_match_2 = -999;
  tree_againstElectronLooseMVA6_2 = -999;
  tree_againstElectronMediumMVA6_2 = -999;
  tree_againstElectronTightMVA6_2 = -999;
  tree_againstElectronVLooseMVA6_2 = -999;
  tree_againstElectronVTightMVA6_2 = -999;
  tree_againstMuonLoose3_2 = -999;
  tree_againstMuonTight3_2 = -999;
  tree_byCombinedIsolationDeltaBetaCorrRaw3Hits_2 = -999;
  tree_byIsolationMVA3newDMwoLTraw_2 = -999;
  tree_byIsolationMVA3oldDMwoLTraw_2 = -999;
  tree_byIsolationMVA3newDMwLTraw_2 = -999;
  tree_byIsolationMVA3oldDMwLTraw_2 = -999;
  tree_chargedIsoPtSum_2 = -999;
  tree_decayModeFindingOldDMs_2 = 0;
  tree_neutralIsoPtSum_2 = -999;
  tree_puCorrPtSum_2 = -999;
  tree_trigweight_2 = -999;
  tree_idisoweight_2 = -999;
  tree_tau_decay_mode_2 = -999 ;
  tree_mva_olddm_medium_2 = 0;
  tree_mva_olddm_tight_2 = 0;
  tree_mva_olddm_vtight_2 = 0;
  //di-tau system
  tree_pt_tt = -999;
  tree_mt_tot = -999;
  tree_m_vis = -999;
  tree_m_sv = -999;
  tree_mt_sv = -999;
  tree_mt_sv = -999;
  tree_pt_sv = -999;
  tree_eta_sv = -999;
  tree_phi_sv = -999;
  tree_os = 0;
  //MET
  tree_met = -999;
  tree_metphi = -999;
  tree_pfmet = -999;
  tree_pfmetphi = -999;
  tree_puppimet = -999;
  tree_puppimetphi = -999;
  tree_mvamet = -999;
  tree_mvametphi = -999;
  tree_pzetavis = -999;
  tree_pzetamiss = -999;
  tree_pfpzetamiss = -999;
  tree_puppipzetamiss = -999;
  tree_mvacov00 = -999;
  tree_mvacov01 = -999;
  tree_mvacov10 = -999;
  tree_mvacov11 = -999;
  tree_metcov00 = -999;
  tree_metcov01 = -999;
  tree_metcov10 = -999;
  tree_metcov11 = -999;
  //VBF system
  tree_mjj = -999;
  tree_jdeta = -999;
  tree_njetingap = -999;
  tree_njetingap20 = -999;
  tree_jdphi = -999;
  //additional jets
  tree_nbtag = -999;
  tree_njets = -999;
  tree_njetspt20 = -999;
  //leading jet sorted by pt
  tree_jpt_1 = -999;
  tree_jeta_1 = -999;
  tree_jphi_1 = -999;
  tree_jrawf_1 = -999;
  tree_jmva_1 = -999;
  //trailing jet sorted by pt
  tree_jpt_2 = -999;
  tree_jeta_2 = -999;
  tree_jphi_2 = -999;
  tree_jrawf_2 = -999;
  tree_jmva_2 = -999;
  //leading b-jet sorted by pt
  tree_bpt_1 = -999;
  tree_beta_1 = -999;
  tree_bphi_1 = -999;
  tree_brawf_1 = -999;
  tree_bmva_1 = -999;
  tree_bcsv_1 = -999;
  //trailing b-jet sorted by pt
  tree_bpt_2 = -999;
  tree_beta_2 = -999;
  tree_bphi_2 = -999;
  tree_brawf_2 = -999;
  tree_bmva_2 = -999;
  tree_bcsv_2 = -999;
  //Extra lepton vetos
  tree_dilepton_veto = 0;
  tree_extraelec_veto = 0;
  tree_extramuon_veto = 0;
  //Trigger
  tree_trg_singlemuon = 0;
  tree_trg_singleelectron = 0;
  tree_trg_singletau_1 = 0;
  tree_trg_singletau_2 = 0;
  tree_trg_doubletau = 0;
  tree_trg_muonelectron = 0;
  //ML quantities and tree helpful variables
  tree_dR = -999;
  tree_pt_tot = -999;
  tree_pt_tt_vis = -999;
  tree_min_deta = -999;
  tree_leg1_centrality = -999;
  tree_leg2_centrality = -999;
  tree_pt_ratio = -999;
  tree_met_centrality = -999;
  tree_category = -999;

  pt_tot = TLorentzVector();
  prediction = -999;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::addBranch(TTree *tree){

  tree->Branch("evt",&tree_evt);
  tree->Branch("sampleNumber",&tree_sampleNumber);
  tree->Branch("weight",&tree_weight);
  //Pielup
  tree->Branch("npv",&tree_npv);
  tree->Branch("npu",&tree_npu);
  tree->Branch("rho",&tree_rho);
  tree->Branch("puweight",&tree_puweight);
  //Leg 1 (leading tau for tt; electon for et,em; muon for mt)
  tree->Branch("pt_1",&tree_pt_1);
  tree->Branch("phi_1",&tree_phi_1);
  tree->Branch("eta_1",&tree_eta_1);
  tree->Branch("m_1",&tree_m_1);
  tree->Branch("q_1",&tree_q_1);
  tree->Branch("d0_1",&tree_d0_1);
  tree->Branch("dZ_1",&tree_dZ_1);
  tree->Branch("mt_1",&tree_mt_1);
  tree->Branch("pfmt_1",&tree_pfmt_1);
  tree->Branch("puppimt_1",&tree_puppimt_1);
  tree->Branch("iso_1",&tree_iso_1);
  tree->Branch("id_e_mva_nt_loose_1",&tree_id_e_mva_nt_loose_1);
  tree->Branch("gen_match_1",&tree_gen_match_1);
  tree->Branch("againstElectronLooseMVA6_1",&tree_againstElectronLooseMVA6_1);
  tree->Branch("againstElectronMediumMVA6_1",&tree_againstElectronMediumMVA6_1);
  tree->Branch("againstElectronTightMVA6_1",&tree_againstElectronTightMVA6_1);
  tree->Branch("againstElectronVLooseMVA6_1",&tree_againstElectronVLooseMVA6_1);
  tree->Branch("againstElectronVTightMVA6_1",&tree_againstElectronVTightMVA6_1);
  tree->Branch("againstMuonLoose3_1",&tree_againstMuonLoose3_1);
  tree->Branch("againstMuonTight3_1",&tree_againstMuonTight3_1);
  tree->Branch("byCombinedIsolationDeltaBetaCorrRaw3Hits_1",&tree_byCombinedIsolationDeltaBetaCorrRaw3Hits_1);
  tree->Branch("byIsolationMVA3newDMwoLTraw_1",&tree_byIsolationMVA3newDMwoLTraw_1);
  tree->Branch("byIsolationMVA3oldDMwoLTraw_1",&tree_byIsolationMVA3oldDMwoLTraw_1);
  tree->Branch("byIsolationMVA3newDMwLTraw_1",&tree_byIsolationMVA3newDMwLTraw_1);
  tree->Branch("byIsolationMVA3oldDMwLTraw_1",&tree_byIsolationMVA3oldDMwLTraw_1);
  tree->Branch("chargedIsoPtSum_1",&tree_chargedIsoPtSum_1);
  tree->Branch("decayModeFindingOldDMs_1",&tree_decayModeFindingOldDMs_1);
  tree->Branch("neutralIsoPtSum_1",&tree_neutralIsoPtSum_1);
  tree->Branch("puCorrPtSum_1",&tree_puCorrPtSum_1);
  tree->Branch("trigweight_1",&tree_trigweight_1);
  tree->Branch("idisoweight_1",&tree_idisoweight_1);
  tree->Branch("trackingweight_1",&tree_trackingweight_1);
  tree->Branch("tau_decay_mode_1",&tree_tau_decay_mode_1);
  tree->Branch("mva_olddm_medium_1",&tree_mva_olddm_medium_1);
  tree->Branch("mva_olddm_tight_1",&tree_mva_olddm_tight_1);
  tree->Branch("mva_olddm_vtight_1",&tree_mva_olddm_vtight_1);
  //Leg 2 (trailing tau for tt, tau for et,mt, muon for em)
  tree->Branch("pt_2",&tree_pt_2);
  tree->Branch("phi_2",&tree_phi_2);
  tree->Branch("eta_2",&tree_eta_2);
  tree->Branch("m_2",&tree_m_2);
  tree->Branch("q_2",&tree_q_2);
  tree->Branch("d0_2",&tree_d0_2);
  tree->Branch("dZ_2",&tree_dZ_2);
  tree->Branch("mt_2",&tree_mt_2);
  tree->Branch("pfmt_2",&tree_pfmt_2);
  tree->Branch("puppimt_2",&tree_puppimt_2);
  tree->Branch("iso_2",&tree_iso_2);
  tree->Branch("id_e_mva_nt_loose_2",&tree_id_e_mva_nt_loose_2);
  tree->Branch("gen_match_2",&tree_gen_match_2);
  tree->Branch("againstElectronLooseMVA6_2",&tree_againstElectronLooseMVA6_2);
  tree->Branch("againstElectronMediumMVA6_2",&tree_againstElectronMediumMVA6_2);
  tree->Branch("againstElectronTightMVA6_2",&tree_againstElectronTightMVA6_2);
  tree->Branch("againstElectronVLooseMVA6_2",&tree_againstElectronVLooseMVA6_2);
  tree->Branch("againstElectronVTightMVA6_2",&tree_againstElectronVTightMVA6_2);
  tree->Branch("againstMuonLoose3_2",&tree_againstMuonLoose3_2);
  tree->Branch("againstMuonTight3_2",&tree_againstMuonTight3_2);
  tree->Branch("byCombinedIsolationDeltaBetaCorrRaw3Hits_2",&tree_byCombinedIsolationDeltaBetaCorrRaw3Hits_2);
  tree->Branch("byIsolationMVA3newDMwoLTraw_2",&tree_byIsolationMVA3newDMwoLTraw_2);
  tree->Branch("byIsolationMVA3oldDMwoLTraw_2",&tree_byIsolationMVA3oldDMwoLTraw_2);
  tree->Branch("byIsolationMVA3newDMwLTraw_2",&tree_byIsolationMVA3newDMwLTraw_2);
  tree->Branch("byIsolationMVA3oldDMwLTraw_2",&tree_byIsolationMVA3oldDMwLTraw_2);
  tree->Branch("chargedIsoPtSum_2",&tree_chargedIsoPtSum_2);
  tree->Branch("decayModeFindingOldDMs_2",&tree_decayModeFindingOldDMs_2);
  tree->Branch("neutralIsoPtSum_2",&tree_neutralIsoPtSum_2);
  tree->Branch("puCorrPtSum_2",&tree_puCorrPtSum_2);
  tree->Branch("trigweight_2",&tree_trigweight_2);
  tree->Branch("idisoweight_2",&tree_idisoweight_2);
  tree->Branch("trackingweight_2",&tree_trackingweight_2);
  tree->Branch("tau_decay_mode_2",&tree_tau_decay_mode_2);
  tree->Branch("mva_olddm_medium_2",&tree_mva_olddm_medium_2);
  tree->Branch("mva_olddm_tight_2",&tree_mva_olddm_tight_2);
  tree->Branch("mva_olddm_vtight_2",&tree_mva_olddm_vtight_2);
  //di-tau system
  tree->Branch("pt_tt",&tree_pt_tt);
  tree->Branch("mt_tot",&tree_mt_tot);
  tree->Branch("m_vis",&tree_m_vis);
  tree->Branch("m_sv",&tree_m_sv);
  tree->Branch("mt_sv",&tree_mt_sv);
  tree->Branch("pt_sv",&tree_pt_sv);
  tree->Branch("eta_sv",&tree_eta_sv);
  tree->Branch("phi_sv",&tree_phi_sv);
  tree->Branch("os",&tree_os);
  //MET
  tree->Branch("met",&tree_met);
  tree->Branch("metphi",&tree_metphi);
  tree->Branch("pfmet",&tree_pfmet);
  tree->Branch("pfmetphi",&tree_pfmetphi);
  tree->Branch("puppimet",&tree_puppimet);
  tree->Branch("puppimetphi",&tree_puppimetphi);
  tree->Branch("mvamet",&tree_mvamet);
  tree->Branch("mvametphi",&tree_mvametphi);
  tree->Branch("pzetavis",&tree_pzetavis);
  tree->Branch("pzetamiss",&tree_pzetamiss);
  tree->Branch("pfpzetamiss",&tree_pfpzetamiss);
  tree->Branch("puppipzetamiss",&tree_puppipzetamiss);
  tree->Branch("mvacov00",&tree_mvacov00);
  tree->Branch("mvacov01",&tree_mvacov01);
  tree->Branch("mvacov10",&tree_mvacov10);
  tree->Branch("mvacov11",&tree_mvacov11);
  tree->Branch("metcov00",&tree_metcov00);
  tree->Branch("metcov01",&tree_metcov01);
  tree->Branch("metcov10",&tree_metcov10);
  tree->Branch("metcov11",&tree_metcov11);
  //VBF system
  tree->Branch("mjj",&tree_mjj);
  tree->Branch("jdeta",&tree_jdeta);
  tree->Branch("njetingap",&tree_njetingap);
  tree->Branch("njetingap20",&tree_njetingap20);
  tree->Branch("jdphi",&tree_jdphi);
  //additional jets
  tree->Branch("nbtag",&tree_nbtag);
  tree->Branch("njets",&tree_njets);
  tree->Branch("njetspt20",&tree_njetspt20);
  //leading jet sorted by pt
  tree->Branch("jpt_1",&tree_jpt_1);
  tree->Branch("jeta_1",&tree_jeta_1);
  tree->Branch("jphi_1",&tree_jphi_1);
  tree->Branch("jrawf_1",&tree_jrawf_1);
  tree->Branch("jmva_1",&tree_jmva_1);
  //trailing jet sorted by pt
  tree->Branch("jpt_2",&tree_jpt_2);
  tree->Branch("jeta_2",&tree_jeta_2);
  tree->Branch("jphi_2",&tree_jphi_2);
  tree->Branch("jrawf_2",&tree_jrawf_2);
  tree->Branch("jmva_2",&tree_jmva_2);
  //leading b-jet sorted by pt
  tree->Branch("bpt_1",&tree_bpt_1);
  tree->Branch("beta_1",&tree_beta_1);
  tree->Branch("bphi_1",&tree_bphi_1);
  tree->Branch("brawf_1",&tree_brawf_1);
  tree->Branch("bmva_1",&tree_bmva_1);
  tree->Branch("bcsv_1",&tree_bcsv_1);
  //trailing b-jet sorted by pt
  tree->Branch("bpt_2",&tree_bpt_2);
  tree->Branch("beta_2",&tree_beta_2);
  tree->Branch("bphi_2",&tree_bphi_2);
  tree->Branch("brawf_2",&tree_brawf_2);
  tree->Branch("bmva_2",&tree_bmva_2);
  tree->Branch("bcsv_2",&tree_bcsv_2);
  //Extra lepton vetos
  tree->Branch("dilepton_veto",&tree_dilepton_veto);
  tree->Branch("extraelec_veto",&tree_extraelec_veto);
  tree->Branch("extramuon_veto",&tree_extramuon_veto);
  //Trigger
  tree->Branch("trg_singlemuon",&tree_trg_singlemuon);
  tree->Branch("trg_singleelectron",&tree_trg_singleelectron);
  tree->Branch("trg_singletau_1",&tree_trg_singletau_1);
  tree->Branch("trg_singletau_2",&tree_trg_singletau_2);
  tree->Branch("trg_doubletau",&tree_trg_doubletau);
  tree->Branch("trg_muonelectron",&tree_trg_muonelectron);
  //ML quantities and tree helpful variables
  tree->Branch("dR",&tree_dR);
  tree->Branch("pt_tot",&tree_pt_tot);
  tree->Branch("pt_tt_vis",&tree_pt_tt_vis);
  tree->Branch("min_deta",&tree_min_deta);
  tree->Branch("leg1_centrality",&tree_leg1_centrality);
  tree->Branch("leg2_centrality",&tree_leg2_centrality);
  tree->Branch("pt_ratio",&tree_pt_ratio);
  tree->Branch("met_centrality",&tree_met_centrality);
  tree->Branch("category",&tree_category);

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::fillMLtree(const EventProxyHTT & myEventProxy, double eventWeightWithSyst){

  std::vector<HTTParticle> &aJets = *myEventProxy.jets;

  ///Filling TTree
  // Event ID variables
  fillEventID(aEvent, myEventProxy);
  // Fill legs
  fillLegs(aPair.getLeg1(), aPair.getLeg2());
  // Fill di-tau system, but also MET which can be pairwise
  fillPair(aEvent, aPair);
  // Fill jet variables includig VBF system
  fillJets(aJets);
  // Fill extra lepton vetoes info
  fillVetoes(aEvent);
  //fill ML specific quantities
  fillMLquantities(myEventProxy);
  fillCategory();

  if(sampleName=="Data") tree_sampleNumber = 0;
  else if(sampleName.find("DY")!=std::string::npos) {
    if(sampleName.find("MatchT")!=std::string::npos) tree_sampleNumber = 1;
    else if(sampleName.find("MatchJ")!=std::string::npos) tree_sampleNumber = 2;
    else if(sampleName.find("MatchL")!=std::string::npos) tree_sampleNumber = 3; }
  else if(sampleName.find("W")!=std::string::npos && sampleName.find("Jets")!=std::string::npos) tree_sampleNumber = 4;
  else if(sampleName.find("ggHTT125")!=std::string::npos) tree_sampleNumber = 5;
  else if(sampleName.find("qqHTT125")!=std::string::npos) tree_sampleNumber = 6;
  else if(sampleName.find("WplusHTT125")!=std::string::npos) tree_sampleNumber = 7;
  else if(sampleName.find("WminusHTT125")!=std::string::npos) tree_sampleNumber = 8;
  else if(sampleName.find("ZHTT125")!=std::string::npos) tree_sampleNumber = 9;

  /*tree_weight = eventWeightWithSyst;
  //Pielup
  tree_npv = myEventProxy.event->getNPV();
  tree_npu = myEventProxy.event->getNPU();
  tree_rho = -999;
  tree_puweight = -999;
*/
  if(apply_preds){
    pList = PyList_New(featuresCount+1);
    PyList_SetItem(pList, 0, PyFloat_FromDouble(tree_dR));
    PyList_SetItem(pList, 1, PyFloat_FromDouble(tree_pt_tot));
    PyList_SetItem(pList, 2, PyFloat_FromDouble(tree_pt_tt_vis));
    PyList_SetItem(pList, 3, PyFloat_FromDouble(tree_min_deta));
    PyList_SetItem(pList, 4, PyFloat_FromDouble(tree_leg1_centrality));
    PyList_SetItem(pList, 5, PyFloat_FromDouble(tree_leg2_centrality));
    PyList_SetItem(pList, 6, PyFloat_FromDouble(tree_pt_ratio));
    PyList_SetItem(pList, 7, PyFloat_FromDouble(tree_met_centrality));
    PyList_SetItem(pList, 8, PyFloat_FromDouble(tree_m_vis));
    PyList_SetItem(pList, 9, PyFloat_FromDouble(tree_m_sv));
    PyList_SetItem(pList, 10, PyFloat_FromDouble(tree_mt_sv));
    PyList_SetItem(pList, 11, PyFloat_FromDouble(tree_mjj));
    PyList_SetItem(pList, 12, PyInt_FromLong(tree_category));
    pArgs = Py_BuildValue("(O)", pList);
    pValue = PyObject_CallObject(pFunc, pArgs);

    //std::cout<<"XGB przed: "<<prediction<<std::endl;//test
    //std::cout<<"Kategoria: "<<tree_category<<std::endl;//test
    prediction = PyFloat_AsDouble(pValue);

    //std::cout<<"XGB: "<<prediction<<std::endl;//test
  }

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::fillEventID(const HTTEvent &event, const EventProxyHTT & myEventProxy){

  tree_run = event.getRunId();
  tree_lumi = 0; //FIXME, not in ntuples
  tree_evt = event.getEventId();
  if(tree_evt==0) std::cout<<getSampleName(myEventProxy)<<std::endl;
  //Pielup
  tree_npv = event.getNPV();
  tree_npu = event.getNPU();
  tree_rho = 0; //FIXME, not in ntuples
  tree_puweight = getPUWeight(myEventProxy);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::fillLegs(const HTTParticle &leg1, const HTTParticle &leg2){

  //Leg 1 (leading tau for tt; electon for et,em; muon for mt)
  tree_pt_1 = leg1.getP4().Pt();
  tree_phi_1 = leg1.getP4().Phi();
  tree_eta_1 = leg1.getP4().Eta();
  tree_m_1 = leg1.getP4().M();
  if(std::abs(leg1.getPDGid())==13)
    tree_m_1 = 0.1057; //muon mass
  else if(std::abs(leg1.getPDGid())==15 && leg1.getProperty(PropertyEnum::decayMode)==0)
    tree_m_1 = 0.13957; //pi+/- mass
  if(tree_m_1 < 0) tree_m_1 = 0; //protect against numerical instability
  tree_q_1 = leg1.getCharge();
  tree_d0_1 = leg1.getProperty(PropertyEnum::dxy);
  tree_dZ_1 = leg1.getProperty(PropertyEnum::dz);
  tree_gen_match_1 = leg1.getProperty(PropertyEnum::mc_match);

  //Leg 2 (trailing tau for tt, electon for et,em muon for mt)
  tree_pt_2 = leg2.getP4().Pt();
  tree_phi_2 = leg2.getP4().Phi();
  tree_eta_2 = leg2.getP4().Eta();
  tree_m_2 = leg2.getP4().M();
  if(std::abs(leg2.getPDGid())==13)
    tree_m_2 = 0.1057; //muon mass
  else if(std::abs(leg2.getPDGid())==15 && leg2.getProperty(PropertyEnum::decayMode)==0)
    tree_m_2 = 0.13957; //pi+/- mass
  if(tree_m_2 < 0) tree_m_2 = 0; //protect against numerical instability
  tree_q_2 = leg2.getCharge();
  tree_d0_2 = leg2.getProperty(PropertyEnum::dxy);
  tree_dZ_2 = leg2.getProperty(PropertyEnum::dz);
  tree_gen_match_2 = leg2.getProperty(PropertyEnum::mc_match); //according to: https://twiki.cern.ch/twiki/bin/view/CMS/HiggsToTauTauWorking2016#MC_Matching

  tree_os = tree_q_1*tree_q_2<0;

  //Decay channel specific
  fillLegsSpecific(leg1,leg2);
 }
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::fillLegsSpecific(const HTTParticle &leg1, const HTTParticle &leg2){
  return;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::fillPair(const HTTEvent &event, HTTPair &pair){

  //Legs
  TLorentzVector leg1P4 = pair.getLeg1().getP4();
  TLorentzVector leg2P4 = pair.getLeg2().getP4();

  //MET
  TLorentzVector pfmetP4(event.getMET().X(), event.getMET().Y(), 0, event.getMET().Mod());
  tree_pfmet = event.getMET().Mod();
  tree_pfmetphi = event.getMET().Phi_mpi_pi(event.getMET().Phi());
  tree_pfmt_1 = TMath::Sqrt(2.*leg1P4.Pt()*pfmetP4.Pt()*(1.-TMath::Cos(leg1P4.Phi()-pfmetP4.Phi())));
  tree_pfmt_2 = TMath::Sqrt(2.*leg2P4.Pt()*pfmetP4.Pt()*(1.-TMath::Cos(leg2P4.Phi()-pfmetP4.Phi())));
  //Pair met: by default called MVAMET which can be wrong...
  TLorentzVector mvametP4(pair.getMET().X(), pair.getMET().Y(), 0, pair.getMET().Mod());
  tree_mvamet = pair.getMET().Mod();
  tree_mvametphi = pair.getMET().Phi_mpi_pi(pair.getMET().Phi());
  tree_met = tree_mvamet;
  tree_metphi = tree_mvametphi;
  tree_mt_1 = TMath::Sqrt(2.*leg1P4.Pt()*mvametP4.Pt()*(1.-TMath::Cos(leg1P4.Phi()-mvametP4.Phi())));
  tree_mt_2 = TMath::Sqrt(2.*leg2P4.Pt()*mvametP4.Pt()*(1.-TMath::Cos(leg2P4.Phi()-mvametP4.Phi())));
  /*
    puppimet;
    puppimetphi;
    puppimt_1;
    puppimt_2;
  */

  //di-tau system
  tree_pt_tt = (mvametP4 + leg1P4 + leg2P4).Pt();
  tree_mt_tot = TMath::Sqrt( 2. * leg1P4.Pt() * mvametP4.Pt() * (1. - TMath::Cos( leg1P4.Phi() - mvametP4.Phi() ) ) + 2. * leg2P4.Pt() * mvametP4.Pt() * (1. - TMath::Cos( leg2P4.Phi() - mvametP4.Phi() ) ) + 2. * leg1P4.Pt() * leg2P4.Pt() * (1. - TMath::Cos(leg1P4.Phi() - leg2P4.Phi() ) ) );
  tree_m_vis = (leg1P4 + leg2P4).M();
  tree_m_sv = pair.getP4().M();

  float leg1Px = leg1P4.Px(), leg1Py = leg1P4.Py(), leg1Phi = leg1P4.Phi();
  float leg2Px = leg2P4.Px(), leg2Py = leg2P4.Py(), leg2Phi = leg2P4.Phi();
  tree_pzetavis = ( (leg1Px+leg2Px)*(TMath::Cos(leg1Phi) + TMath::Cos(leg2Phi)) + (leg1Py + leg2Py)*(TMath::Sin(leg1Phi) + TMath::Sin(leg2Phi)) ) / ( TMath::Sqrt( TMath::Power(TMath::Cos(leg1Phi) + TMath::Cos(leg2Phi),2) + TMath::Power(TMath::Sin(leg1Phi) + TMath::Sin(leg2Phi),2) ) );
  tree_pzetamiss = ( pair.getMET().X()*(TMath::Cos(leg1Phi) + TMath::Cos(leg2Phi)) +  pair.getMET().Y()*(TMath::Sin(leg1Phi) + TMath::Sin(leg2Phi)) ) / ( TMath::Sqrt( TMath::Power(TMath::Cos(leg1Phi) + TMath::Cos(leg2Phi),2) + TMath::Power(TMath::Sin(leg1Phi) + TMath::Sin(leg2Phi),2) ) );
  tree_pfpzetamiss = ( event.getMET().X()*(TMath::Cos(leg1Phi) + TMath::Cos(leg2Phi)) +  event.getMET().Y()*(TMath::Sin(leg1Phi) + TMath::Sin(leg2Phi)) ) / ( TMath::Sqrt( TMath::Power(TMath::Cos(leg1Phi) + TMath::Cos(leg2Phi),2) + TMath::Power(TMath::Sin(leg1Phi) + TMath::Sin(leg2Phi),2) ) );
  //puppipzetamiss;

  //The following does not like that pair is a const reference
  tree_mvacov00 = pair.getMETMatrix().at(0);
  tree_mvacov01 = pair.getMETMatrix().at(1);
  tree_mvacov10 = pair.getMETMatrix().at(2);
  tree_mvacov11 = pair.getMETMatrix().at(3);
  tree_metcov00 = tree_mvacov00;
  tree_metcov01 = tree_mvacov01;
  tree_metcov10 = tree_mvacov10;
  tree_metcov11 = tree_mvacov11;

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::fillJets(const std::vector<HTTParticle> &jets){

  tree_njetspt20 = jets.size();
  tree_njets = 0;
  tree_nbtag = 0;
  std::vector<HTTParticle> bjets;
  for(unsigned int iJet=0; iJet<jets.size(); ++iJet){
    if(jets.at(iJet).getP4().Pt()>30)
      tree_njets++;
    if(std::abs(jets.at(iJet).getP4().Eta())<2.4 &&
       jets.at(iJet).getP4().Pt()>20 &&
       jets.at(iJet).getProperty(PropertyEnum::bCSVscore)>0.8484 && //FIXME 0.8484 Correct??
       promoteBJet(jets.at(iJet)) &&
       true){
      tree_nbtag++;
      bjets.push_back(jets.at(iJet));
    }
  }
  if(jets.size()==0) return;

  HTTParticle aLeadingJet = jets.at(0);
  //leading jet sorted by pt
  tree_jpt_1 = aLeadingJet.getP4().Pt();
  tree_jeta_1 = aLeadingJet.getP4().Eta();
  tree_jphi_1 = aLeadingJet.getP4().Phi();
  tree_jrawf_1 =  aLeadingJet.getProperty(PropertyEnum::rawPt)/aLeadingJet.getP4().Pt();
  tree_jmva_1 = aLeadingJet.getProperty(PropertyEnum::PUJetID);

  if(jets.size()>1){
    HTTParticle aTrailingJet = jets.at(1);

    //trailing jet sorted by pt
    tree_jpt_2 = aTrailingJet.getP4().Pt();
    tree_jeta_2 = aTrailingJet.getP4().Eta();
    tree_jphi_2 = aTrailingJet.getP4().Phi();
    tree_jrawf_2 = aTrailingJet.getProperty(PropertyEnum::rawPt)/aTrailingJet.getP4().Pt();
    tree_jmva_2 = aTrailingJet.getProperty(PropertyEnum::PUJetID);

    //VBF system
    tree_mjj = (aLeadingJet.getP4()+aTrailingJet.getP4()).M();
    tree_jdeta = std::abs(aLeadingJet.getP4().Eta()-aTrailingJet.getP4().Eta());
    tree_jdphi = aLeadingJet.getP4().Phi()-aTrailingJet.getP4().Phi();
    while(tree_jdphi>TMath::Pi()) tree_jdphi -= 2.*TMath::Pi();
    while(tree_jdphi<=-TMath::Pi()) tree_jdphi += 2.*TMath::Pi();
    tree_jdphi = std::abs(tree_jdphi);
    //jets in eta gap
    for(unsigned int iJet=2; iJet<jets.size(); ++iJet){
      if( (jets.at(iJet).getP4().Eta()>aLeadingJet.getP4().Eta()&&jets.at(iJet).getP4().Eta()<aTrailingJet.getP4().Eta()) ||
	  (jets.at(iJet).getP4().Eta()<aLeadingJet.getP4().Eta()&&jets.at(iJet).getP4().Eta()>aTrailingJet.getP4().Eta()) ){
	tree_njetingap20++;
	if(jets.at(iJet).getP4().Pt()>30)
	  tree_njetingap++;
      }
    }
  }
  //b-jets
  if(bjets.size()>0){
    //leading b-jet sorted by pt
    tree_bpt_1 = bjets.at(0).getP4().Pt();
    tree_beta_1 = bjets.at(0).getP4().Eta();
    tree_bphi_1 = bjets.at(0).getP4().Phi();
    tree_brawf_1 = bjets.at(0).getProperty(PropertyEnum::rawPt)/bjets.at(0).getP4().Pt();
    tree_bmva_1 = bjets.at(0).getProperty(PropertyEnum::PUJetID);
    tree_bcsv_1 = bjets.at(0).getProperty(PropertyEnum::bCSVscore);
    if(bjets.size()>1){
      //trailing b-jet sorted by pt
      tree_bpt_2 = bjets.at(1).getP4().Pt();
      tree_beta_2 = bjets.at(1).getP4().Eta();
      tree_bphi_2 = bjets.at(1).getP4().Phi();
      tree_brawf_2 = bjets.at(1).getProperty(PropertyEnum::rawPt)/bjets.at(1).getP4().Pt();
      tree_bmva_2 = bjets.at(1).getProperty(PropertyEnum::PUJetID);
      tree_bcsv_2 = bjets.at(1).getProperty(PropertyEnum::bCSVscore);
    }
  }
  return;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool HTTAnalyzer::promoteBJet(const HTTParticle &jet){
  return true;//test AP
  //MB: https://twiki.cern.ch/twiki/bin/view/CMS/BTagCalibration#Standalone
  /*bool decision = false;

  if(!reader) initializeBTagCorrections();

  BTagEntry::JetFlavor jetFlavour;
  if(std::abs(jet.getProperty(PropertyEnum::Flavour))==5)//b-quark
    jetFlavour = BTagEntry::FLAV_B;
  else if(std::abs(jet.getProperty(PropertyEnum::Flavour))==4)//c-quark
    jetFlavour = BTagEntry::FLAV_C;
  else //light quark, gluon or undefined
    jetFlavour = BTagEntry::FLAV_UDSG;
  // Note: this is for b jets, for c jets (light jets) use FLAV_C (FLAV_UDSG)
  double btag_SF = reader->eval_auto_bounds("central",
					    jetFlavour,
					    jet.getP4().Eta(),
					    jet.getP4().Pt()
					    //,jet.getProperty(PropertyEnum::bCSVscore) //MB: it is not needed when WP is definied
					    );
  rand_->SetSeed((int)((jet.getP4().Eta()+5)*100000));
  double rand_num = rand_->Rndm();
  //debug
  //std::cout<<"\tbtag_SF(flav,CSVv2): "<<btag_SF
  //	   <<"("<<jetFlavour<<","
  //	   <<jet.getProperty(PropertyEnum::bCSVscore)<<")"<<std::endl;
  //std::cout<<"\tbtag_rand_num: "<<rand_num<<std::endl;
  if(btag_SF>1){
    double tagging_efficiency = 1;
    TH2F *histo_eff = btag_eff_oth_;
    if(jetFlavour == BTagEntry::FLAV_B)
      histo_eff = btag_eff_b_;
    else if(jetFlavour == BTagEntry::FLAV_C)
      histo_eff = btag_eff_c_;
    if( jet.getP4().Pt() > histo_eff->GetXaxis()->GetBinLowEdge(histo_eff->GetNbinsX()+1) ){
      tagging_efficiency = histo_eff->GetBinContent( histo_eff->GetNbinsX(),histo_eff->GetYaxis()->FindBin(std::abs(jet.getP4().Eta())) );
    }
    else{
      tagging_efficiency = histo_eff->GetBinContent( histo_eff->GetXaxis()->FindBin(jet.getP4().Pt()),histo_eff->GetYaxis()->FindBin(std::abs(jet.getP4().Eta())) );
    }
    //debug
    //std::cout<<"\tbtag_eff: "<<tagging_efficiency<<std::endl;
    if(tagging_efficiency < 1e-9)//protection
      decision = false;
    else if(tagging_efficiency > 1.-1e-9)//protection
      decision = true;
    else
      decision = (rand_num < (1. - btag_SF)/(1. - 1./tagging_efficiency) );
  }
  else{
    decision = (rand_num < 1. - btag_SF);
  }
  //debug
  //std::cout<<"\tbtag_decision: "<<decision<<std::endl;
  return !decision;*/
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::fillVetoes(const HTTEvent &event){

  tree_dilepton_veto = event.checkSelectionBit(SelectionBitsEnum::diMuonVeto);
  tree_extraelec_veto = event.checkSelectionBit(SelectionBitsEnum::extraElectronVeto);
  tree_extramuon_veto = event.checkSelectionBit(SelectionBitsEnum::extraMuonVeto);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::fillMLquantities(const EventProxyHTT & myEventProxy){

  std::vector<HTTParticle> &aJets = *myEventProxy.jets;
  TLorentzVector leg1P4 = aPair.getLeg1().getP4();
  TLorentzVector leg2P4 = aPair.getLeg2().getP4();
  TLorentzVector mvametP4(aPair.getMET().X(), aPair.getMET().Y(), 0, aPair.getMET().Mod());

  tree_dR = leg1P4.DeltaR(leg2P4);
  pt_tot = mvametP4 + leg1P4 + leg2P4;
  tree_pt_tt_vis = leg1P4.Pt() + leg2P4.Pt();
  tree_pt_ratio = leg1P4.Eta()/leg2P4.Eta();
  float A = std::sin(mvametP4.Phi() - leg2P4.Phi())/std::sin(leg1P4.Phi()-leg2P4.Phi());
  float B = std::sin(leg1P4.Phi() - mvametP4.Phi())/std::sin(leg1P4.Phi()-leg2P4.Phi());
  tree_met_centrality = (A+B)/std::sqrt(pow(A,2)+pow(B,2));

  HTTParticle aLeadingJet, aTrailingJet;
  if(aJets.size()>0) {
    aLeadingJet = aJets.at(0);
    pt_tot += aLeadingJet.getP4();
    }
  if(aJets.size()>1) {
    aTrailingJet = aJets.at(1);
    pt_tot += aTrailingJet.getP4();
    tree_min_deta = std::min(std::abs((leg1P4+leg2P4).Eta() - aLeadingJet.getP4().Eta()), std::abs((leg1P4+leg2P4).Eta() - aTrailingJet.getP4().Eta()));
    float eta1 = leg1P4.Eta();
    float eta2 = leg2P4.Eta();
    float jeta1 = aLeadingJet.getP4().Eta();
    float jeta2 = aTrailingJet.getP4().Eta();
    tree_leg1_centrality = std::exp(-4/std::pow(jeta1-jeta2,2)*std::pow(eta1-(jeta1+jeta2)/2,2));
    tree_leg2_centrality = std::exp(-4/std::pow(jeta1-jeta2,2)*std::pow(eta2-(jeta1+jeta2)/2,2));
    }

  tree_pt_tot = pt_tot.Pt();
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::fillCategory(){

    std::vector<std::string> mainCategories = {"0jet", "boosted", "vbf"};
    for(unsigned int iCategory = 0; iCategory<myNumberOfCategories; ++iCategory) {
      bool isMain = false;
      //std::cout<<myChannelSpecifics->getCategoryRejester().at(iCategory)->name()<<": "<<iCategory<<std::endl;
      for(auto s: mainCategories) if(s==myChannelSpecifics->getCategoryRejester().at(iCategory)->name()) isMain=true;
      if(!isMain) continue;
      if(categoryDecisions.at(iCategory)==true) tree_category=iCategory;
    }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
