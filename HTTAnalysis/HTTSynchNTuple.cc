#include <sstream>

#include "HTTSynchNTuple.h"
#include "HTTHistograms.h"
#include "EventProxyHTT.h"

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
HTTSynchNTuple::HTTSynchNTuple(const std::string & aName):Analyzer(aName){ tmpName = "h1DXSignal";}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
HTTSynchNTuple::~HTTSynchNTuple(){ if(myHistos_) delete myHistos_; }
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
Analyzer* HTTSynchNTuple::clone() const{

  HTTSynchNTuple* clone = new HTTSynchNTuple(name());
  clone->setHistos(myHistos_);
  return clone;

};
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTSynchNTuple::initialize(TDirectory* aDir,
			      pat::strbitset *aSelections){

  mySelections_ = aSelections;
  
  ///The histograms for this analyzer will be saved into "HTTSynchNTuple"
  ///directory of the ROOT file
  ///NOTE: due to a bug hists land in the Summary directory
  myHistos_ = new HTTHistograms(aDir, selectionFlavours_);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTSynchNTuple::finalize(){ 

  myHistos_->finalizeHistograms(0,1.0);
 
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void HTTSynchNTuple::addBranch(TTree *tree){

	//event ID variables
	tree->Branch("run",&run);
	tree->Branch("lumi",&lumi);
	tree->Branch("evt",&evt);
	//Pielup
	tree->Branch("npv",&npv);
	tree->Branch("npu",&npu);
	tree->Branch("rho",&rho);
	//Leg 1 (leading tau for tt, electon for et,em muon for mt)
	tree->Branch("pt_1",&pt_1);
	tree->Branch("phi_1",&phi_1);
	tree->Branch("eta_1",&eta_1);
	tree->Branch("m_1",&m_1);
	tree->Branch("q_1",&q_1);
	tree->Branch("d0_1",&d0_1);
	tree->Branch("dZ_1",&dZ_1);
	tree->Branch("mt_1",&mt_1);
	tree->Branch("pfmt_1",&pfmt_1);
	tree->Branch("puppimt_1",&puppimt_1);
	tree->Branch("iso_1",&iso_1);
	tree->Branch("id_e_mva_nt_loose_1",&id_e_mva_nt_loose_1);
	tree->Branch("gen_match_1",&gen_match_1);
	tree->Branch("againstElectronLooseMVA6_1",&againstElectronLooseMVA6_1);
	tree->Branch("againstElectronMediumMVA6_1",&againstElectronMediumMVA6_1);
	tree->Branch("againstElectronTightMVA6_1",&againstElectronTightMVA6_1);
	tree->Branch("againstElectronVLooseMVA6_1",&againstElectronVLooseMVA6_1);
	tree->Branch("againstElectronVTightMVA6_1",&againstElectronVTightMVA6_1);
	tree->Branch("againstMuonLoose3_1",&againstMuonLoose3_1);
	tree->Branch("againstMuonTight3_1",&againstMuonTight3_1);
	tree->Branch("byCombinedIsolationDeltaBetaCorrRaw3Hits_1",&byCombinedIsolationDeltaBetaCorrRaw3Hits_1);
	tree->Branch("byIsolationMVA3newDMwoLTraw_1",&byIsolationMVA3newDMwoLTraw_1);
	tree->Branch("byIsolationMVA3oldDMwoLTraw_1",&byIsolationMVA3oldDMwoLTraw_1);
	tree->Branch("byIsolationMVA3newDMwLTraw_1",&byIsolationMVA3newDMwLTraw_1);
	tree->Branch("byIsolationMVA3oldDMwLTraw_1",&byIsolationMVA3oldDMwLTraw_1);
	tree->Branch("chargedIsoPtSum_1",&chargedIsoPtSum_1);
	tree->Branch("decayModeFindingOldDMs_1",&decayModeFindingOldDMs_1);
	tree->Branch("neutralIsoPtSum_1",&neutralIsoPtSum_1);
	tree->Branch("puCorrPtSum_1",&puCorrPtSum_1);
	tree->Branch("trigweight_1",&trigweight_1);
	tree->Branch("idisoweight_1",&idisoweight_1);
	//Leg 2 (trailing tau for tt, electon for et,em muon for mt)
	tree->Branch("pt_2",&pt_2);
	tree->Branch("phi_2",&phi_2);
	tree->Branch("eta_2",&eta_2);
	tree->Branch("m_2",&m_2);
	tree->Branch("q_2",&q_2);
	tree->Branch("d0_2",&d0_2);
	tree->Branch("dZ_2",&dZ_2);
	tree->Branch("mt_2",&mt_2);
	tree->Branch("pfmt_2",&pfmt_2);
	tree->Branch("puppimt_2",&puppimt_2);
	tree->Branch("iso_2",&iso_2);
	tree->Branch("id_e_mva_nt_loose_2",&id_e_mva_nt_loose_2);
	tree->Branch("gen_match_2",&gen_match_2);
	tree->Branch("againstElectronLooseMVA6_2",&againstElectronLooseMVA6_2);
	tree->Branch("againstElectronMediumMVA6_2",&againstElectronMediumMVA6_2);
	tree->Branch("againstElectronTightMVA6_2",&againstElectronTightMVA6_2);
	tree->Branch("againstElectronVLooseMVA6_2",&againstElectronVLooseMVA6_2);
	tree->Branch("againstElectronVTightMVA6_2",&againstElectronVTightMVA6_2);
	tree->Branch("againstMuonLoose3_2",&againstMuonLoose3_2);
	tree->Branch("againstMuonTight3_2",&againstMuonTight3_2);
	tree->Branch("byCombinedIsolationDeltaBetaCorrRaw3Hits_2",&byCombinedIsolationDeltaBetaCorrRaw3Hits_2);
	tree->Branch("byIsolationMVA3newDMwoLTraw_2",&byIsolationMVA3newDMwoLTraw_2);
	tree->Branch("byIsolationMVA3oldDMwoLTraw_2",&byIsolationMVA3oldDMwoLTraw_2);
	tree->Branch("byIsolationMVA3newDMwLTraw_2",&byIsolationMVA3newDMwLTraw_2);
	tree->Branch("byIsolationMVA3oldDMwLTraw_2",&byIsolationMVA3oldDMwLTraw_2);
	tree->Branch("chargedIsoPtSum_2",&chargedIsoPtSum_2);
	tree->Branch("decayModeFindingOldDMs_2",&decayModeFindingOldDMs_2);
	tree->Branch("neutralIsoPtSum_2",&neutralIsoPtSum_2);
	tree->Branch("puCorrPtSum_2",&puCorrPtSum_2);
	tree->Branch("trigweight_2",&trigweight_2);
	tree->Branch("idisoweight_2",&idisoweight_2);
	//di-tau system
	tree->Branch("pt_tt",&pt_tt);
	tree->Branch("mt_tot",&mt_tot);
	tree->Branch("m_vis",&m_vis);
	tree->Branch("m_sv",&m_sv);
	tree->Branch("mt_sv",&mt_sv);
	//MET
	tree->Branch("met",&met);
	tree->Branch("metphi",&metphi);
	tree->Branch("puppimet",&puppimet);
	tree->Branch("puppimetphi",&puppimetphi);
	tree->Branch("mvamet",&mvamet);
	tree->Branch("mvametphi",&mvametphi);
	tree->Branch("pzetavis",&pzetavis);
	tree->Branch("pzetamiss",&pzetamiss);
	tree->Branch("pfpzetamiss",&pfpzetamiss);
	tree->Branch("puppipzetamiss",&puppipzetamiss);
	tree->Branch("mvacov00",&mvacov00);
	tree->Branch("mvacov01",&mvacov01);
	tree->Branch("mvacov10",&mvacov10);
	tree->Branch("mvacov11",&mvacov11);
	tree->Branch("metcov00",&metcov00);
	tree->Branch("metcov01",&metcov01);
	tree->Branch("metcov10",&metcov10);
	tree->Branch("metcov11",&metcov11);
	//VBF system
	tree->Branch("mjj",&mjj);
	tree->Branch("jdeta",&jdeta);
	tree->Branch("njetingap",&njetingap);
	tree->Branch("njetingap20",&njetingap20);
	tree->Branch("jdphi",&jdphi);
	//additional jets
	tree->Branch("nbtag",&nbtag);
	tree->Branch("njets",&njets);
	tree->Branch("njetspt20",&njetspt20);
	//leading jet sorted by pt
	tree->Branch("jpt_1",&jpt_1);
	tree->Branch("jeta_1",&jeta_1);
	tree->Branch("jphi_1",&jphi_1);
	tree->Branch("jrawf_1",&jrawf_1);
	tree->Branch("jmva_1",&jmva_1);
	//trailing jet sorted by pt
	tree->Branch("jpt_2",&jpt_2);
	tree->Branch("jeta_2",&jeta_2);
	tree->Branch("jphi_2",&jphi_2);
	tree->Branch("jrawf_2",&jrawf_2);
	tree->Branch("jmva_2",&jmva_2);
	//leading b-jet sorted by pt
	tree->Branch("bpt_1",&bpt_1);
	tree->Branch("beta_1",&beta_1);
	tree->Branch("bphi_1",&bphi_1);
	tree->Branch("brawf_1",&brawf_1);
	tree->Branch("bmva_1",&bmva_1);
	tree->Branch("bcsv_1",&bcsv_1);
	//trailing b-jet sorted by pt
	tree->Branch("bpt_2",&bpt_2);
	tree->Branch("beta_2",&beta_2);
	tree->Branch("bphi_2",&bphi_2);
	tree->Branch("brawf_2",&brawf_2);
	tree->Branch("bmva_2",&bmva_2);
	tree->Branch("bcsv_2",&bcsv_2);
	//Extra lepton vetos
	tree->Branch("dilepton_veto",&dilepton_veto);
	tree->Branch("extraelec_veto",&extraelec_veto);
	tree->Branch("extramuon_veto",&extramuon_veto);
	tree->Branch("puweight",&puweight);

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
bool HTTSynchNTuple::analyze(const EventProxyBase& iEvent){

  const EventProxyHTT & myEventProxy = static_cast<const EventProxyHTT&>(iEvent);

  aEvent = *myEventProxy.event;  
  aPair = (*myEventProxy.pairs)[0];
  aTau = aPair.getTau();
  aMuon = aPair.getMuon();

	i_++;
	//muonPt = aMuon.getP4().Vect().Perp();
	if(i_%1000==0){
		//std::cout<<i_<<std::endl;
		//std::cout<<muonPt<<std::endl;
		}
	
	//Filling TTree
	//event ID variables
	run = aEvent.getRunId();
	lumi = 0;						//NEED TO FIX, not in ntuples
	evt = aEvent.getEventId();
	//Pielup
	npv = aEvent.getNPV();
	npu = aEvent.getNPU();
	rho = 0;						//NEED TO FIX, not in ntuples
	//Leg 1 (leading tau for tt, electon for et,em muon for mt) - here Leg 1 is always muon!
	leg1 = aMuon;
	pt_1 = leg1.getP4().Perp();
	phi_1 = leg1.getP4().Phi();
	eta_1 = leg1.getP4().Eta();
	m_1 = leg1.getP4().M();
/*
//////////////////////////////////////////////////////////////////////////
	if(m_1 < 0) {
		std::cout<<"Masa: "<<m_1<<", nr przypadka: "<<evt<<std::endl<<"Px : "<<leg1.getP4().Px()<<" Py: "<<leg1.getP4().Py()<<" Pz: "<<leg1.getP4().Pz()<<" E: "<<leg1.getP4().E()<<"\n"; 
		}
		myHistos_->fillProfile("hProfMuon_PtVsMass",m_1,leg1.getP4().Pt());
		
//////////////////////////////////////////////////////////////////////////
*/
	q_1 = leg1.getCharge();
	d0_1 = leg1.getProperty(PropertyEnum::dxy);
	dZ_1 = leg1.getProperty(PropertyEnum::dz);
	mt_1 = aPair.getMTMuon();/*
	pfmt_1;
	puppimt_1;*/
	iso_1 = leg1.getProperty(PropertyEnum::combreliso);/*
	id_e_mva_nt_loose_1;
	gen_match_1;		
	againstElectronLooseMVA6_1 = leg1.getProperty(PropertyEnum::againstElectronLooseMVA6);				//FIX: following 7 properites need to be added
	againstElectronMediumMVA6_1 = leg1.getProperty(PropertyEnum::againstElectronMediumMVA6);
	againstElectronTightMVA6_1 = leg1.getProperty(PropertyEnum::againstElectronTightMVA6);
	againstElectronVLooseMVA6_1 = leg1.getProperty(PropertyEnum::againstElectronVLooseMVA6);
	againstElectronVTightMVA6_1 = leg1.getProperty(PropertyEnum::againstElectronVTightMVA6);
	againstMuonLoose3_1 = leg1.getProperty(PropertyEnum::againstMuonLoose3);
	againstMuonTight3_1 = leg1.getProperty(PropertyEnum::againstMuonLoose3);*/
	byCombinedIsolationDeltaBetaCorrRaw3Hits_1 = leg1.getProperty(PropertyEnum::byCombinedIsolationDeltaBetaCorrRaw3Hits);/*
	byIsolationMVA3newDMwoLTraw_1;
	byIsolationMVA3oldDMwoLTraw_1;
	byIsolationMVA3newDMwLTraw_1;
	byIsolationMVA3oldDMwLTraw_1;
	chargedIsoPtSum_1;
	decayModeFindingOldDMs_1;
	neutralIsoPtSum_1;
	puCorrPtSum_1;
	trigweight_1;
	idisoweight_1;*/
	//Leg 2 (trailing tau for tt, electon for et,em muon for mt)
	leg2 = aTau;
	pt_2 = leg2.getP4().Perp();
	phi_2 = leg2.getP4().Phi();
	eta_2 = leg2.getP4().Eta();
	m_2 = leg2.getP4().M();
	q_2 = leg2.getCharge();
	d0_2 = leg2.getProperty(PropertyEnum::dxy);
	dZ_2 = leg2.getProperty(PropertyEnum::dz);
	mt_2 = abs(leg2.getPDGid())==15 ? aPair.getMTLeg2() : aPair.getMTLeg1();/*
	pfmt_2;
	puppimt_2;*/
	iso_2 = leg2.getProperty(PropertyEnum::combreliso);/*
	id_e_mva_nt_loose_2;*/
	gen_match_2 = 5;		/*//according to: https://twiki.cern.ch/twiki/bin/view/CMS/HiggsToTauTauWorking2016#MC_Matching
	againstElectronLooseMVA6_2 = leg2.getProperty(PropertyEnum::againstElectronLooseMVA6);				//FIX: following 7 properites need to be added
	againstElectronMediumMVA6_2 = leg2.getProperty(PropertyEnum::againstElectronMediumMVA6);
	againstElectronTightMVA6_2 = leg2.getProperty(PropertyEnum::againstElectronTightMVA6);
	againstElectronVLooseMVA6_2 = leg2.getProperty(PropertyEnum::againstElectronVLooseMVA6);
	againstElectronVTightMVA6_2 = leg2.getProperty(PropertyEnum::againstElectronVTightMVA6);
	againstMuonLoose3_2 = leg2.getProperty(PropertyEnum::againstMuonLoose3);
	againstMuonTight3_2 = leg2.getProperty(PropertyEnum::againstMuonLoose3);*/
	byCombinedIsolationDeltaBetaCorrRaw3Hits_2 = leg2.getProperty(PropertyEnum::byCombinedIsolationDeltaBetaCorrRaw3Hits);/*
	byIsolationMVA3newDMwoLTraw_2;
	byIsolationMVA3oldDMwoLTraw_2;
	byIsolationMVA3newDMwLTraw_2;
	byIsolationMVA3oldDMwLTraw_2;
	chargedIsoPtSum_2;
	decayModeFindingOldDMs_2;
	neutralIsoPtSum_2;
	puCorrPtSum_2;
	trigweight_2;
	idisoweight_2;
	//di-tau system
	pt_tt;
	mt_tot;
	m_vis;
	m_sv;
	mt_sv;
	//MET
*/
	met = aPair.getMET().Mod();
/*
	metphi = aPair.getProperty(PropertyEnum::metphi);
	puppimet;
	puppimetphi;
	mvamet;
	mvametphi;
	pzetavis;
	pzetamiss;
	pfpzetamiss;
	puppipzetamiss;
	mvacov00;
	mvacov01;
	mvacov10;
	mvacov11;
	metcov00;
	metcov01;
	metcov10;
	metcov11;
	//VBF system
	mjj;
	jdeta;
	njetingap;
	njetingap20;
	jdphi;
	//additional jets
	nbtag;
	njets;
	njetspt20;
	//leading jet sorted by pt
	jpt_1;
	jeta_1;
	jphi_1;
	jrawf_1;
	jmva_1;
	//trailing jet sorted by pt
	jpt_2;
	jeta_2;
	jphi_2;
	jrawf_2;
	jmva_2;
	//leading b-jet sorted by pt
	bpt_1;
	beta_1;
	bphi_1;
	brawf_1;
	bmva_1;
	bcsv_1;
	//trailing b-jet sorted by pt
	bpt_2;
	beta_2;
	bphi_2;
	brawf_2;
	bmva_2;
	bcsv_2;
	//Extra lepton vetos
	dilepton_veto = aEvent.checkSelectionBit(SelectionBitsEnum::diMuonVeto);
/*	extraelec_veto;
	extramuon_veto;
	puweight;*/

  return true;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

