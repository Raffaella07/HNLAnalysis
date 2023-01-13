#include "setStyle.C"


int BkgComposition(std::string var , int bin, float min, float max){


	const int nQCD =9; 
	std::string QCDsamples[nQCD] = {"QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8/crab_QCD_Pt-15to20/220114_163748/0000/",
                                        "QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8/crab_QCD_Pt-20to30/220114_163711/0000/",
                                        "QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8/crab_QCD_Pt-30to50/220114_164140/0000/",
                                        "QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8/crab_QCD_Pt-50to80/220114_163904/0000/",
                                        "QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8/crab_QCD_Pt-80to120/220114_164220/0000/",
                                        "QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8/crab_QCD_Pt-80to120_ext/220114_164025/0000/",
                                        "QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8/crab_QCD_Pt-120to170_ext/220114_163827/0000/",
                                        "QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8/crab_QCD_Pt-120to170/220114_163944/0000/",
                                        "QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/crab_QCD_Pt-170to300/220114_164102/0000/"
					};

	TChain * c[nQCD];
	double weights[nQCD] = { 475.836, 64.5889,35.5222,14.0434, 2.6714, 2.6714, 0.797302, 0.797302, 0.151216};

//	TH1D* h_piUnmatched[nQCD];
//	TH1D* h_eleUnmatched[nQCD];
	TH2D* h_SVmatched[nQCD];

//	TH1D* hSum_piUnmatched = new TH1D("pi_unMatched","pi_unMatched_",bin,min,max);
//	TH1D* hSum_eleUnmatched = new TH1D("ele_unMatched","ele_unMatched",bin,min,max);
	TH2D* hSum_SVmatched= new TH2D("sv_Matched","sv_unMatched",50,0,600,50,0,600);      
 
	std::string TriggerFiring = "&&   Muon_isTriggeringBPark[BToMuEPi_trg_mu_idx]==1 ";//&&   (Muon_fired_HLT_Mu8_IP4[BToMuEPi_trg_mu_idx]==1 || Muon_fired_HLT_Mu8_IP3[BToMuEPi_trg_mu_idx]==1 || Muon_fired_HLT_Mu8_IP5[BToMuEPi_trg_mu_idx]==1 || Muon_fired_HLT_Mu8_IP6[BToMuEPi_trg_mu_idx]==1 || Muon_fired_HLT_Mu8p5_IP3p5[BToMuEPi_trg_mu_idx]==1 || Muon_fired_HLT_Mu9_IP4[BToMuEPi_trg_mu_idx]==1 || Muon_fired_HLT_Mu9_IP5[BToMuEPi_trg_mu_idx]==1 || Muon_fired_HLT_Mu9_IP6[BToMuEPi_trg_mu_idx]==1 || Muon_fired_HLT_Mu10p5_IP3p5[BToMuEPi_trg_mu_idx]==1 || Muon_fired_HLT_Mu12_IP6[BToMuEPi_trg_mu_idx]==1 )";
	std::string pf_sel = " && Muon_pt[BToMuEPi_trg_mu_idx]>7  && Muon_softId[BToMuEPi_trg_mu_idx]==1 && fabs(Muon_eta[BToMuEPi_trg_mu_idx])<1.5 && Electron_pfmvaId[BToMuEPi_sel_ele_idx]> -3  && ProbeTracks_pt[BToMuEPi_pi_idx]>0.7  &&  fabs(ProbeTracks_eta[BToMuEPi_pi_idx])<2  &&  fabs(ProbeTracks_dz[BToMuEPi_pi_idx])>0.005  &&  fabs(ProbeTracks_dxy[BToMuEPi_pi_idx])>0.005  &&  fabs(ProbeTracks_dxyS[BToMuEPi_pi_idx])>3  &&  fabs(ProbeTracks_dzS[BToMuEPi_pi_idx])>1.5  &&  fabs(ProbeTracks_DCASig[BToMuEPi_pi_idx])>5  &&  BToMuEPi_sv_prob>0.001 &&  BToMuEPi_hnl_cos2D>0.99  &&  BToMuEPi_sv_lxy/BToMuEPi_sv_lxye>20  &&  fabs(Electron_dxy[BToMuEPi_sel_ele_idx])>0.008 &&  fabs(Electron_dz[BToMuEPi_sel_ele_idx])>0.008 &&  fabs(Electron_dxy[BToMuEPi_sel_ele_idx]/Electron_dxyErr[BToMuEPi_sel_ele_idx])>1.5 &&  fabs(Electron_dz[BToMuEPi_sel_ele_idx]/Electron_dzErr[BToMuEPi_sel_ele_idx])>1 ";

	//study performed on events with at least trigger muon matched 
	std::string SV_matched = "Muon_genPartIdx[BToMuEPi_trg_mu_idx]!=-1 && Electron_genPartIdx[BToMuEPi_sel_ele_idx]!=-1 && fabs(GenPart_pdgId[Electron_genPartIdx[BToMuEPi_sel_ele_idx]])==11 && ProbeTracks_genPartIdx[BToMuEPi_pi_idx]!=-1 &&  fabs(GenPart_pdgId[GenPart_genPartIdxMother[ProbeTracks_genPartIdx[BToMuEPi_pi_idx]]])>250 &&  fabs(GenPart_pdgId[GenPart_genPartIdxMother[ProbeTracks_genPartIdx[BToMuEPi_pi_idx]]])<350"; //HNL matched
	std::string pion_unmatched = "Muon_genPartIdx[BToMuEPi_trg_mu_idx]!=-1 && Electron_genPartIdx[BToMuEPi_sel_ele_idx]!=-1 && fabs(GenPart_pdgId[Electron_genPartIdx[BToMuEPi_sel_ele_idx]])==11  && ProbeTracks_genPartIdx[BToMuEPi_pi_idx]==-1"; //pion unmatched
	std::string ele_unmatched ="Muon_genPartIdx[BToMuEPi_trg_mu_idx]!=-1 && Electron_genPartIdx[BToMuEPi_sel_ele_idx]==-1 && ProbeTracks_genPartIdx[BToMuEPi_pi_idx]!=-1"; //ele unmatched
	std::string pdgId_eleMother = "GenPart_pdgId[GenPart_genPartIdxMother[Electron_genPartIdx[BToMuEPi_sel_ele_idx]]]";
	std::string pdgId_piMother = "GenPart_pdgId[GenPart_genPartIdxMother[ProbeTracks_genPartIdx[BToMuEPi_pi_idx]]]";
	std::string pdgId_svMother = "GenPart_pdgId[GenPart_genPartIdxMother[ProbeTracks_genPartIdx[BToMuEPi_pi_idx]]]";


	for (int i=0; i < nQCD; i++){
	

//	h_piUnmatched[i] = new TH1D(("pi_unMatched_"+std::to_string(i)).c_str(),("pi_unMatched_"+std::to_string(i)).c_str(),bin,min,max);
//	h_eleUnmatched[i]= new TH1D(("ele_unMatched_"+std::to_string(i)).c_str(),("ele_unMatched_"+std::to_string(i)).c_str(),bin,min,max);
        h_SVmatched[i]= new TH2D(("sv_Matched_"+std::to_string(i)).c_str(),("sv_unMatched_"+std::to_string(i)).c_str(),50,0,600,50,0,600);
	c[i] = new TChain("Events");

	c[i]->Add(("/pnfs/roma1.infn.it/data/cms/store/user/ratramon/HNLGen_ntuples/"+QCDsamples[i]+"*.root").c_str());
	std::cout << c[i]->GetEntries() << std::endl;
//	c[i]->Draw((var+">>"+h_piUnmatched[i]->GetName()).c_str(),(pion_unmatched + TriggerFiring + pf_sel ).c_str());
//	c[i]->Draw((var+">>"+h_eleUnmatched[i]->GetName()).c_str(),(ele_unmatched + TriggerFiring + pf_sel ).c_str());
//	c[i]->Draw((" GenPart_pdgId[GenPart_genPartIdxMother[ProbeTracks_genPartIdx[BToMuEPi_pi_idx]]]:GenPart_pdgId[GenPart_genPartIdxMother[Electron_genPartIdx[BToMuEPi_sel_ele_idx]]]>>"+h_SVmatched[i]->GetName()).c_str(),(SV_matched + TriggerFiring + pf_sel ).c_str());
	c[i]->Draw(("fabs(GenPart_pdgId[GenPart_genPartIdxMother[Electron_genPartIdx[BToMuEPi_sel_ele_idx]]]):fabs(GenPart_pdgId[GenPart_genPartIdxMother[ProbeTracks_genPartIdx[BToMuEPi_pi_idx]]])>>"+std::string(h_SVmatched[i]->GetName())).c_str(),"Electron_genPartIdx[BToMuEPi_sel_ele_idx]!=-1 && ProbeTracks_genPartIdx[BToMuEPi_pi_idx]!=-1 && BToMuEPi_sv_lxy<5 && fabs(GenPart_pdgId[GenPart_genPartIdxMother[ProbeTracks_genPartIdx[BToMuEPi_pi_idx]]])<600 && fabs(GenPart_pdgId[GenPart_genPartIdxMother[Electron_genPartIdx[BToMuEPi_sel_ele_idx]]])<600");
	
//	std::cout << h_piUnmatched[i]->GetEntries() << std::endl;
//	h_piUnmatched[i]->Scale(weights[i]);
//	h_eleUnmatched[i]->Scale(weights[i]);
	h_SVmatched[i]->Scale(weights[i]);
//	std::cout << h_piUnmatched[i]->GetEntries() << std::endl;
//	hSum_piUnmatched->Add(h_piUnmatched[i]);
//	hSum_eleUnmatched->Add(h_eleUnmatched[i]);
	hSum_SVmatched->Add(h_SVmatched[i]);

	}

//	hSum_piUnmatched->Scale(1/hSum_piUnmatched->Integral());
//	hSum_eleUnmatched->Scale(1/hSum_eleUnmatched->Integral());
//	hSum_eleUnmatched->Scale(1/hSum_eleUnmatched->Integral());
	hSum_SVmatched->Scale(1/hSum_SVmatched->Integral());
//	hSum_piUnmatched->SaveAs(("Bkg_421/Pi_Unmatched_"+var+".root").c_str());
//	hSum_eleUnmatched->SaveAs(("Bkg_421/Ele_Unmatched_"+var+".root").c_str());
	hSum_SVmatched->SaveAs(("Bkg_421/SV_matched_"+var+".root").c_str());
	


	return 0;	

	}






