#include <sstream>
#include <bitset>

#include "HTTAnalyzer.h"
#include "HTTHistograms.h"
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::getPreselectionEff(const EventProxyHTT & myEventProxy){

        if(true || ntupleFile_!=myEventProxy.getTTree()->GetCurrentFile()) {
                ntupleFile_ = myEventProxy.getTTree()->GetCurrentFile();
                if(hStatsFromFile) delete hStatsFromFile;
                hStatsFromFile = (TH1F*)ntupleFile_->Get("hStats");

                std::string hName = "h1DStats"+getSampleName(myEventProxy);
                TH1F *hStats = myHistos_->get1DHistogram(hName,true);

                float genWeight = getGenWeight(myEventProxy);

                hStats->SetBinContent(2,std::abs(hStatsFromFile->GetBinContent(hStatsFromFile->FindBin(1))*genWeight));
                hStats->SetBinContent(3,std::abs(hStatsFromFile->GetBinContent(hStatsFromFile->FindBin(3))*genWeight));
        }
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::string HTTAnalyzer::getSampleName(const EventProxyHTT & myEventProxy){

        return getSampleNameFromFileName(myEventProxy);

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
float HTTAnalyzer::getSystWeight(const HTTAnalysis::sysEffects & aSystEffect){

        if(aSystEffect==HTTAnalysis::NOMINAL) return 1.0;

        if(aSystEffect==HTTAnalysis::ZPtUp && sampleName.find("DYJets")!=std::string::npos) return aEvent.getPtReWeight();
        else if(aSystEffect==HTTAnalysis::ZPtDown && sampleName.find("DYJets")!=std::string::npos) return 1.0/aEvent.getPtReWeight();
        else if(aSystEffect==HTTAnalysis::TTUp && sampleName.find("TTbar")!=std::string::npos) return aEvent.getPtReWeight();
        else if(aSystEffect==HTTAnalysis::TTDown && sampleName.find("TTbar")!=std::string::npos) return 1.0/aEvent.getPtReWeight();
        else if((aSystEffect==HTTAnalysis::J2TUp || aSystEffect==HTTAnalysis::J2TDown)
                && aLeg2.getProperty(PropertyEnum::mc_match)==6) {
                float delta = 0.2*aLeg2.getP4().Perp()/100.0;
                if(aLeg2.getP4().Perp()>200) delta = 0.4;
                if(aSystEffect==HTTAnalysis::J2TDown) delta*=-1;
                return 1-delta;
        }
        else return 1.0;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::string HTTAnalyzer::getSampleNameFromFileName(const EventProxyHTT & myEventProxy){

        std::string fileName = myEventProxy.getTTree()->GetCurrentFile()->GetName();

        //std::map<std::string, std::string>::const_iterator sampleNameIt = fileName2sampleName.find(fileName);
        //if(sampleNameIt!=fileName2sampleName.end()) return sampleNameIt->second;

        std::string sampleName = "Unknown";

        if(fileName.find("W1JetsToLNu")!=std::string::npos) sampleName = "W1Jets";
        else if(fileName.find("W2JetsToLNu")!=std::string::npos) sampleName = "W2Jets";
        else if(fileName.find("W3JetsToLNu")!=std::string::npos) sampleName = "W3Jets";
        else if(fileName.find("W4JetsToLNu")!=std::string::npos) sampleName = "W4Jets";
        else if(fileName.find("WJetsToLNu")!=std::string::npos && myEventProxy.event->getLHEnOutPartons()==0) sampleName = "W0Jets";
        else if(fileName.find("WJetsToLNu")!=std::string::npos && myEventProxy.event->getLHEnOutPartons()==1) sampleName = "W1JetsIncl";
        else if(fileName.find("WJetsToLNu")!=std::string::npos && myEventProxy.event->getLHEnOutPartons()==2) sampleName = "W2JetsIncl";
        else if(fileName.find("WJetsToLNu")!=std::string::npos && myEventProxy.event->getLHEnOutPartons()==3) sampleName = "W3JetsIncl";
        else if(fileName.find("WJetsToLNu")!=std::string::npos && myEventProxy.event->getLHEnOutPartons()==4) sampleName = "W4JetsIncl";
        //else if(fileName.find("WJetsToLNu")!=std::string::npos && myEventProxy.event->getLHEnOutPartons()>0) sampleName = "WAllJets";

        else if(fileName.find("Run201")!=std::string::npos) sampleName =  "Data";
        else if(fileName.find("SUSYGluGluToHToTauTau")!=std::string::npos) sampleName =  "ATT";

        else if(fileName.find("GluGluHToTauTauM110")!=std::string::npos) sampleName =  "ggHTT110";
        else if(fileName.find("GluGluHToTauTauM120")!=std::string::npos) sampleName =  "ggHTT120";
        else if(fileName.find("GluGluHToTauTauM125")!=std::string::npos) sampleName =  "ggHTT125";
        else if(fileName.find("GluGluHToTauTauM130")!=std::string::npos) sampleName =  "ggHTT130";
        else if(fileName.find("GluGluHToTauTauM140")!=std::string::npos) sampleName =  "ggHTT140";

        else if(fileName.find("VBFHToTauTauM110")!=std::string::npos) sampleName =  "qqHTT110";
        else if(fileName.find("VBFHToTauTauM120")!=std::string::npos) sampleName =  "qqHTT120";
        else if(fileName.find("VBFHToTauTauM125")!=std::string::npos) sampleName =  "qqHTT125";
        else if(fileName.find("VBFHToTauTauM130")!=std::string::npos) sampleName =  "qqHTT130";
        else if(fileName.find("VBFHToTauTauM140")!=std::string::npos) sampleName =  "qqHTT140";

        else if(fileName.find("WplusHToTauTauM110")!=std::string::npos) sampleName =  "WplusHTT110";
        else if(fileName.find("WplusHToTauTauM120")!=std::string::npos) sampleName =  "WplusHTT120";
        else if(fileName.find("WplusHToTauTauM125")!=std::string::npos) sampleName =  "WplusHTT125";
        else if(fileName.find("WplusHToTauTauM130")!=std::string::npos) sampleName =  "WplusHTT130";
        else if(fileName.find("WplusHToTauTauM140")!=std::string::npos) sampleName =  "WplusHTT140";

        else if(fileName.find("WminusHToTauTauM110")!=std::string::npos) sampleName =  "WminusHTT110";
        else if(fileName.find("WminusHToTauTauM120")!=std::string::npos) sampleName =  "WminusHTT120";
        else if(fileName.find("WminusHToTauTauM125")!=std::string::npos) sampleName =  "WminusHTT125";
        else if(fileName.find("WminusHToTauTauM130")!=std::string::npos) sampleName =  "WminusHTT130";
        else if(fileName.find("WminusHToTauTauM140")!=std::string::npos) sampleName =  "WminusHTT140";

        else if(fileName.find("ZHToTauTauM110")!=std::string::npos) sampleName =  "ZHTT110";
        else if(fileName.find("ZHToTauTauM120")!=std::string::npos) sampleName =  "ZHTT120";
        else if(fileName.find("ZHToTauTauM125")!=std::string::npos) sampleName =  "ZHTT125";
        else if(fileName.find("ZHToTauTauM130")!=std::string::npos) sampleName =  "ZHTT130";
        else if(fileName.find("ZHToTauTauM140")!=std::string::npos) sampleName =  "ZHTT140";

        else if(fileName.find("STtWantitop")!=std::string::npos) sampleName =  "Wantitop";
        else if(fileName.find("STtWtop")!=std::string::npos) sampleName =  "Wtop";
        else if(fileName.find("STtchannel__antitop")!=std::string::npos) sampleName =  "t-channel_antitop";
        else if(fileName.find("STtchannel__top")!=std::string::npos) sampleName =  "t-channel_top";
        else if(fileName.find("ZZTo2L2Q")!=std::string::npos) sampleName =  "ZZTo2L2Q";
        else if(fileName.find("ZZTo4L")!=std::string::npos) sampleName =  "ZZTo4L";
        else if(fileName.find("WZTo1L3Nu")!=std::string::npos) sampleName =  "WZTo1L3Nu";
        else if(fileName.find("WZJToLLLNu")!=std::string::npos) sampleName =  "WZJToLLLNu";
        else if(fileName.find("WWTo1L1Nu2Q")!=std::string::npos) sampleName =  "WWTo1L1Nu2Q";
        else if(fileName.find("WWToLNuQQ")!=std::string::npos) sampleName =  "WWTo1L1Nu2Q";
        else if(fileName.find("WZTo1L1Nu2Q")!=std::string::npos) sampleName =  "WZTo1L1Nu2Q";
        else if(fileName.find("VVTo2L2Nu")!=std::string::npos) sampleName =  "VVTo2L2Nu";
        else if(fileName.find("WZTo2L2Q")!=std::string::npos) sampleName =  "WZTo2L2Q";
        else if(fileName.find("EWKWMinus")!=std::string::npos) sampleName =  "EWKWMinus";
        else if(fileName.find("EWKWPlus")!=std::string::npos) sampleName =  "EWKWPlus";
        else if(fileName.find("EWKZ2JetsZToLL")!=std::string::npos) sampleName =  "EWKZ2JetsZToLL";
        else if(fileName.find("EWKZ2JetsZToNuNu")!=std::string::npos) sampleName =  "EWKZ2JetsZToNuNu";
        else if(fileName.find("QCD")!=std::string::npos) sampleName =  "QCD_MC";
        else if(fileName.find("DY")!=std::string::npos) sampleName =  getDYSampleName(myEventProxy);
        else if(fileName.find("TTTune")!=std::string::npos) sampleName =  "TTbar";
        std::string matchingMode = getMatchingName(myEventProxy);

        if(sampleName=="Unknown") std::cout<<"Unkwown sample type. "<<fileName<<std::endl;

        sampleName += matchingMode;
        fileName2sampleName[fileName] = sampleName;

        return sampleName;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::string HTTAnalyzer::getDYSampleName(const EventProxyHTT & myEventProxy){

        std::string fileName = myEventProxy.getTTree()->GetCurrentFile()->GetName();

        std::string jetsName = "";
        if(fileName.find("DY1JetsToLLM50")!=std::string::npos) jetsName ="1Jets";
        else if(fileName.find("DY2JetsToLLM50")!=std::string::npos) jetsName = "2Jets";
        else if(fileName.find("DY3JetsToLLM50")!=std::string::npos) jetsName = "3Jets";
        else if(fileName.find("DY4JetsToLLM50")!=std::string::npos) jetsName = "4Jets";
        else if(fileName.find("DYJetsToLLM10to50")!=std::string::npos) jetsName =  "LowM";
        else if(fileName.find("DYJetsToLLM50")!=std::string::npos && myEventProxy.event->getLHEnOutPartons()==0) jetsName =  "0Jets";
        else if(fileName.find("DYJetsToLLM50")!=std::string::npos && myEventProxy.event->getLHEnOutPartons()==1) jetsName =  "1JetsIncl";
        else if(fileName.find("DYJetsToLLM50")!=std::string::npos && myEventProxy.event->getLHEnOutPartons()==2) jetsName =  "2JetsIncl";
        else if(fileName.find("DYJetsToLLM50")!=std::string::npos && myEventProxy.event->getLHEnOutPartons()==3) jetsName =  "3JetsIncl";
        else if(fileName.find("DYJetsToLLM50")!=std::string::npos && myEventProxy.event->getLHEnOutPartons()==4) jetsName =  "4JetsIncl";
        //else if(fileName.find("DYJetsToLLM50")!=std::string::npos && myEventProxy.event->getLHEnOutPartons()>0) jetsName =  "AllJets";

        int decayModeBoson = myEventProxy.event->getDecayModeBoson();
        int leg1MCMatch = 6, leg2MCMatch = 6;
        if(myEventProxy.pairs->size()) {
                HTTPair aPair = (*myEventProxy.pairs)[0];
                HTTParticle aLeg1 = aPair.getLeg1();
                HTTParticle aLeg2 = aPair.getLeg2();
                leg1MCMatch = aLeg1.getProperty(PropertyEnum::mc_match);
                leg2MCMatch = aLeg2.getProperty(PropertyEnum::mc_match);
        }
        std::string decayName = "Unknown";
        if(fileName.find("MT_")!=std::string::npos) {
                if(leg2MCMatch<5) decayName = "L";
                else if(leg2MCMatch==5) decayName = "T";
                else decayName = "J";
        }
        if(fileName.find("TT_")!=std::string::npos) {
                if(leg1MCMatch==5 && leg2MCMatch==5) decayName = "T";
                else if(leg1MCMatch<6 && leg2MCMatch<6) decayName = "L";
                else decayName = "J";
        }
        return "DY"+jetsName+"Match"+decayName;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::string HTTAnalyzer::getMatchingName(const EventProxyHTT & myEventProxy){

        std::string fileName = myEventProxy.getTTree()->GetCurrentFile()->GetName();
        std::vector<std::string> sampleNames = {"TTTune", "ZZTo2L2Q", "ZZTo4L","WZTo1L3Nu", "WZJToLLLNu", "WWTo1L1Nu2Q", "WZTo1L1Nu2Q", "VVTo2L2Nu", "WZTo2L2Q"};

        bool sampleToAnalyze = false;
        for(auto sampleName : sampleNames) {
                if(fileName.find(sampleName)!=std::string::npos) sampleToAnalyze = true;
        }
        if(!sampleToAnalyze) return "";
        if(!(fileName.find("TT_")!=std::string::npos || fileName.find("MT_")!=std::string::npos)) return "";

        int tauMCMatch_1 = 6, tauMCMatch_2 = 6;
        if(myEventProxy.pairs->size() && fileName.find("MT_")!=std::string::npos) {
                HTTPair aPair = (*myEventProxy.pairs)[0];
                HTTParticle aLeg2 = aPair.getTau();
                tauMCMatch_1 = aLeg2.getProperty(PropertyEnum::mc_match);
                if(tauMCMatch_1 == 5) return "MatchT"; else return "MatchJ";
        }

        if(myEventProxy.pairs->size() && fileName.find("TT_")!=std::string::npos) {
                HTTPair aPair = (*myEventProxy.pairs)[0];
                HTTParticle aLeg1 = aPair.getLeg1(), aLeg2 = aPair.getLeg2();
                tauMCMatch_1 = aLeg1.getProperty(PropertyEnum::mc_match);
                tauMCMatch_2 = aLeg2.getProperty(PropertyEnum::mc_match);
                if(tauMCMatch_1 == 5 && tauMCMatch_2 == 5) return "MatchT"; else return "MatchJ";
        }

        return "";
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
float HTTAnalyzer::getPUWeight(const EventProxyHTT & myEventProxy){

        if(getSampleName(myEventProxy)=="Data") return 1.0;

        if(!puDataFile_ || !puMCFile_ ||
           puDataFile_->IsZombie() ||
           puMCFile_->IsZombie()) { return 1.0; }

        if(!hPUVec_.size()) hPUVec_.resize(1);

        if(!hPUVec_[0]) {
                std::string hName = "pileup";
                TH1F *hPUData = (TH1F*)puDataFile_->Get(hName.c_str());
                TH1F *hPUSample = (TH1F*)puMCFile_->Get(hName.c_str());
                ///Normalise both histograms.
                //TEST hPUData->Scale(1.0/hPUData->Integral(0,hPUData->GetNbinsX()+1));
                hPUSample->Scale(1.0/hPUSample->Integral(0,hPUSample->GetNbinsX()+1));
                ///
                hPUData->SetDirectory(0);
                hPUSample->SetDirectory(0);
                //TEST hPUData->Divide(hPUSample);
                hPUData->SetName(("h1DPUWeight"+getSampleName(myEventProxy)).c_str());
                hPUVec_[0] =  hPUData;
        }

        int iBinPU = hPUVec_[0]->FindBin(myEventProxy.event->getNPU());
        return hPUVec_[0]->GetBinContent(iBinPU);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
float HTTAnalyzer::getGenWeight(const EventProxyHTT & myEventProxy){

        ///MC weights cab be quite large, but are always +-const.
        ///to avoid counter overflow we keep only sign.
        return myEventProxy.event->getMCWeight()/fabs(myEventProxy.event->getMCWeight());
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::clearTTreeVariables(){

  tree_isData = 0;
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
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTAnalyzer::addBranch(TTree *tree){

  tree->Branch("isData",&tree_isData);
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

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
