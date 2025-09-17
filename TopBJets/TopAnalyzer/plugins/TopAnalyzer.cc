//// -*- C++ -*-
//
// Package:    TopBJets/TopAnalyzer
// Class:      TopAnalyzer
// 
/**\class TopAnalyzer TopAnalyzer.cc TopBJets/TopAnalyzer/plugins/TopAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Ridhi Chawla
//         Created:  Thu, 19 Sep 2019 18:16:53 GMT
//
//

// system include files
#include <memory>
#include <regex>

// user include files
#include <vector>
#include <iostream>
#include <TLorentzVector.h>
#include <TMath.h>
#include "TTree.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/JetReco/interface/TrackJet.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include "DataFormats/PatCandidates/interface/MET.h"

using namespace edm;
using namespace std;
using namespace reco;
using namespace TMath;

// Typedefs
typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LV;

class TopAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit TopAnalyzer(const edm::ParameterSet&);
      ~TopAnalyzer();
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      void GenAnalysis(const edm::Event&,const edm::EventSetup&,
                       edm::Handle<edm::View<reco::GenParticle> > &,
                       edm::Handle<edm::View<pat::PackedGenParticle> > &);
      void RecoAnalysis(const edm::Event&,const edm::EventSetup&,
                        edm::Handle<edm::View<pat::PackedCandidate> > &,
                        edm::Handle<edm::View<pat::Jet> > &);
      void DataAnalysis(const edm::Event&,
                        const edm::EventSetup&,
                        edm::Handle<edm::View<pat::Jet> > &);

      template<typename TYPE> LV getLV(const TYPE&) const;
      void setbit(UShort_t&, UShort_t);
      bool isAncestor(const reco::Candidate*,const reco::Candidate*);
      bool isAncestor(const reco::Candidate*,double);
      std::vector<const pat::Jet*> isTightJet(edm::Handle<edm::View<pat::Jet> > &);
      bool matchingTrack(const reco::Candidate*,
                         const pat::PackedCandidate&);
      std::vector<const pat::PackedCandidate*> matchingTrack(const reco::Candidate*,
                                                             edm::Handle<edm::View<pat::PackedCandidate> > &);
      std::vector<std::vector<double> > unpackTrackParameters(const pat::PackedCandidate&,
                                                              std::string);
      std::vector<std::vector<double> > unpackTrackCovariance(const pat::PackedCandidate&,
                                                              std::string);
      edm::EDGetTokenT<edm::View<reco::GenParticle> > genToken_;
      edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_;
      edm::EDGetTokenT<edm::View<reco::Vertex> > vertexToken_;
      edm::EDGetTokenT<edm::View<pat::PackedCandidate> > tracksToken_;
      edm::EDGetTokenT<edm::View<pat::Jet> > jetsToken_;
      edm::EDGetTokenT<edm::View<pat::Muon> > muonsToken_;
      edm::EDGetTokenT<edm::View<pat::Electron> > electronsToken_;
      edm::EDGetTokenT<edm::View<pat::MET> > metToken_;
      edm::EDGetTokenT<edm::TriggerResults> triggerToken_;

      // Global variables
      bool misMC;
      const double mk = 0.493677;
      const double mpi = 0.13957039;
      const double mp = 0.9382720813;
      const double me = 0.00051099895;
      const double mmu = 0.1056583755;
      unsigned nkpi, nkpipi, nkkpi, ndstar, npkpi;
      const reco::Vertex *_pv;
      int nlepton;
      const pat::PackedCandidate nulltrack;
      std::vector<const reco::Candidate*> top_quarks;
      std::vector<const reco::Candidate*> bottom_hadrons;
      std::vector<const reco::Candidate*> charm_hadrons;
      std::vector<const reco::Candidate*> leptons;
      std::vector<const reco::Candidate*> d0_hadron, lepton_d0, kaon_d0, pion_d0;
      std::vector<const reco::Candidate*> dt_hadron, lepton_dt, kaon_dt, pion_dt, pionsoft_dt;
      std::vector<const reco::Candidate*> dp_hadron, lepton_dp, kaon_dp, pion1_dp, pion2_dp;
      std::vector<const reco::Candidate*> ds_hadron, lepton_ds, kaon1_ds, kaon2_ds, pion_ds;
      std::vector<const reco::Candidate*> lambda_hadron, lepton_lambda, kaon_lambda, pion_lambda, proton_lambda;
      std::vector<const pat::Muon*> muons;
      std::vector<const pat::Electron*> electrons;
      std::vector<double> para, covmat;
      std::vector<std::vector<double> > par, covm;
      std::vector<std::vector<double> > kd0_par, pid0_par, kd0_cov, pid0_cov;
      std::vector<std::vector<double> > kdt_par, pidt_par, pisdt_par, kdt_cov, pidt_cov, pisdt_cov;
      std::vector<std::vector<double> > kdp_par, pi1dp_par, pi2dp_par, kdp_cov, pi1dp_cov, pi2dp_cov;
      std::vector<std::vector<double> > k1ds_par, k2ds_par, pids_par, k1ds_cov, k2ds_cov, pids_cov;
      std::vector<std::vector<double> > klambda_par, pilambda_par, plambda_par, klambda_cov, pilambda_cov, plambda_cov;
      std::regex Mu17_Mu8_VVL,
                 Mu17_TkMu8_VVL,
                 Mu17_Mu8_DZ,
                 Mu17_TkMu8_DZ,
                 IsoMu24,
                 IsoTkMu24,
                 Ele23_Ele12_DZ,
                 Double_Ele33_MW,
                 Double_Ele33_VL,
                 Ele27_WPTight,
                 Mu8_Ele23_IsoVL,
                 Mu23_Ele12_IsoVL,
                 Mu8_Ele23_DZ,
                 Mu23_Ele12_DZ;
      std::regex Mu17_Mu8_DZ_Mass8,
                 IsoMu27,
                 Ele23_Ele12_IsoVL,
                 Ele35_WPTight;
      std::regex Mu17_Mu8_DZ_Mass3p8,
                 Ele32_WPTight,
                 Double_Ele25_MW;

      // Ntuple and branch variables
      TTree* ntuple_;
      edm::Service<TFileService> fileService;
      UInt_t run_,
             event_,
             lumi_;
      double nPVertex_,
             mVertex1_,
             mVertex2_,
             zPVertex_,
             vxTrack_,
             vyTrack_,
             dz_trackPV_,
             pfMET_,
             nLepGen_,
             nJet_;
      // 2016 triggers
      bool hlt_Mu17_Mu8_VVL_,
           hlt_Mu17_TkMu8_VVL_,
           hlt_Mu17_Mu8_DZ_,
           hlt_Mu17_TkMu8_DZ_,
           hlt_IsoMu24_,
           hlt_IsoTkMu24_,
           hlt_Ele23_Ele12_DZ_,
           hlt_Double_Ele33_MW_,
           hlt_Double_Ele33_VL_,
           hlt_Ele27_WPTight_,
           hlt_Mu8_Ele23_IsoVL_,
           hlt_Mu23_Ele12_IsoVL_,
           hlt_Mu8_Ele23_DZ_,
           hlt_Mu23_Ele12_DZ_;
      // 2017 triggers
      bool hlt_Mu17_Mu8_DZ_Mass8_,     
           hlt_IsoMu27_,
           hlt_Ele23_Ele12_IsoVL_,
           hlt_Ele35_WPTight_;
      // 2018 triggers
      bool hlt_Mu17_Mu8_DZ_Mass3p8_,
           hlt_Ele32_WPTight_,
           hlt_Double_Ele25_MW_;
      std::vector<LV> p4_muon_;
      std::vector<LV> p4_electron_;
      std::vector<double> charge_muon_;
      std::vector<double> charge_electron_;
      std::vector<double> pfIso_muon_;
      std::vector<bool> passLooseId_muon_;
      std::vector<bool> passMediumId_muon_;
      std::vector<bool> passTightId_muon_;
      std::vector<bool> passVetoId_electron_;
      std::vector<bool> passLooseId_electron_;
      std::vector<bool> passMediumId_electron_;
      std::vector<bool> passTightId_electron_;
      //std::vector<UShort_t> passIdbits_electron_;
      //std::vector<std::vector<int> > passId_electron_;
      bool foundGenBlepton_;
      bool foundGenBPlusdecay_;
      bool foundGenB0Dstardecay_;
      bool foundGenB0decay_;
      bool foundGenBs0decay_;
      bool foundGenLambdab0decay_;
      std::vector<LV> p4_gen_D0_lepton_;
      std::vector<LV> p4_gen_Dstar_lepton_;
      std::vector<LV> p4_gen_DPlus_lepton_;
      std::vector<LV> p4_gen_DsPlus_lepton_;
      std::vector<LV> p4_gen_LambdacPlus_lepton_;
      std::vector<LV> p4_jet_;
      std::vector<double> lepton_perJet_;
      std::vector<std::vector<double> > pdgId_lepton_;
      std::vector<std::vector<double> > pt_lepton_;
      std::vector<std::vector<double> > eta_lepton_;
      std::vector<std::vector<double> > phi_lepton_;
      std::vector<std::vector<double> > charge_lepton_;
      std::vector<std::vector<double> > mass_lepton_;
      std::vector<LV> p4_D0_;
      std::vector<LV> p4_D0_kaon_;
      std::vector<LV> p4_D0_pion_;
      std::vector<double> mdiff_Dstar_D0_;
      std::vector<LV> p4_Dstar_;
      std::vector<LV> p4_Dstar_D0_;
      std::vector<LV> p4_Dstar_D0_kaon_;
      std::vector<LV> p4_Dstar_D0_pion_;
      std::vector<LV> p4_Dstar_pionsoft_;
      std::vector<LV> p4_DPlus_;
      std::vector<LV> p4_DPlus_kaon_;
      std::vector<LV> p4_DPlus_pion1_;
      std::vector<LV> p4_DPlus_pion2_;
      std::vector<LV> p4_DsPlus_;
      std::vector<LV> p4_DsPlus_kaon1_;
      std::vector<LV> p4_DsPlus_kaon2_;
      std::vector<LV> p4_DsPlus_pion_;
      std::vector<LV> p4_LambdacPlus_;
      std::vector<LV> p4_LambdacPlus_kaon_;
      std::vector<LV> p4_LambdacPlus_pion_;
      std::vector<LV> p4_LambdacPlus_proton_;
      std::vector<double> charge_D0_kaon_;
      std::vector<double> charge_D0_pion_;
      std::vector<double> charge_Dstar_D0_kaon_;
      std::vector<double> charge_Dstar_D0_pion_;
      std::vector<double> charge_Dstar_pionsoft_;
      std::vector<double> charge_DPlus_kaon_;
      std::vector<double> charge_DPlus_pion1_;
      std::vector<double> charge_DPlus_pion2_;
      std::vector<double> charge_DsPlus_kaon1_;
      std::vector<double> charge_DsPlus_kaon2_;
      std::vector<double> charge_DsPlus_pion_;
      std::vector<double> charge_LambdacPlus_kaon_;
      std::vector<double> charge_LambdacPlus_pion_;
      std::vector<double> charge_LambdacPlus_proton_;
      std::vector<std::vector<double> > trkpara_D0_kaon_;
      std::vector<std::vector<double> > trkpara_D0_pion_;
      std::vector<std::vector<double> > trkcovm_D0_kaon_;
      std::vector<std::vector<double> > trkcovm_D0_pion_;
      std::vector<std::vector<double> > trkpara_Dstar_D0_kaon_;
      std::vector<std::vector<double> > trkpara_Dstar_D0_pion_;
      std::vector<std::vector<double> > trkpara_Dstar_pionsoft_;
      std::vector<std::vector<double> > trkcovm_Dstar_D0_kaon_;
      std::vector<std::vector<double> > trkcovm_Dstar_D0_pion_;
      std::vector<std::vector<double> > trkcovm_Dstar_pionsoft_;
      std::vector<std::vector<double> > trkpara_DPlus_kaon_;
      std::vector<std::vector<double> > trkpara_DPlus_pion1_;
      std::vector<std::vector<double> > trkpara_DPlus_pion2_;
      std::vector<std::vector<double> > trkcovm_DPlus_kaon_;
      std::vector<std::vector<double> > trkcovm_DPlus_pion1_;
      std::vector<std::vector<double> > trkcovm_DPlus_pion2_;
      std::vector<std::vector<double> > trkpara_DsPlus_kaon1_;
      std::vector<std::vector<double> > trkpara_DsPlus_kaon2_;
      std::vector<std::vector<double> > trkpara_DsPlus_pion_;
      std::vector<std::vector<double> > trkcovm_DsPlus_kaon1_;
      std::vector<std::vector<double> > trkcovm_DsPlus_kaon2_;
      std::vector<std::vector<double> > trkcovm_DsPlus_pion_;
      std::vector<std::vector<double> > trkpara_LambdacPlus_kaon_;
      std::vector<std::vector<double> > trkpara_LambdacPlus_pion_;
      std::vector<std::vector<double> > trkpara_LambdacPlus_proton_; 
      std::vector<std::vector<double> > trkcovm_LambdacPlus_kaon_;
      std::vector<std::vector<double> > trkcovm_LambdacPlus_pion_;
      std::vector<std::vector<double> > trkcovm_LambdacPlus_proton_;
      std::vector<double> pv_lxy_D0_;
      std::vector<double> pv_lxy_Dstar_D0_;
      std::vector<double> pv_lxy_DPlus_;
      std::vector<double> pv_lxy_DsPlus_;
      std::vector<double> pv_lxy_LambdacPlus_;
};

TopAnalyzer::TopAnalyzer(const edm::ParameterSet& iConfig) :
  genToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("prunedgenparticles"))),
  packedGenToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packedgenparticles"))),
  vertexToken_(consumes<edm::View<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("offlinePrimaryVertices"))),
  tracksToken_(consumes<edm::View<pat::PackedCandidate> >(iConfig.getParameter<edm::InputTag>("pfCandidates"))),
  jetsToken_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jets"))),
  muonsToken_(consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("muons"))),
  electronsToken_(consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electrons"))),
  metToken_(consumes<edm::View<pat::MET> >(iConfig.getParameter<edm::InputTag>("met"))),
  triggerToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("trigger")))
{
  misMC = iConfig.getUntrackedParameter<bool>("isMC");

  Mu17_Mu8_VVL = "(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v)(.*)";
  Mu17_TkMu8_VVL = "(HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v)(.*)";
  Mu17_Mu8_DZ = "(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v)(.*)";
  Mu17_TkMu8_DZ = "(HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v)(.*)";
  IsoMu24 = "(HLT_IsoMu24_v)(.*)";
  IsoTkMu24 = "(HLT_IsoTkMu24_v)(.*)";
  Ele23_Ele12_DZ = "(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v)(.*)";
  Double_Ele33_MW = "(HLT_DoubleEle33_CaloIdL_MW_v)(.*)";
  Double_Ele33_VL = "(HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v)(.*)";
  Ele27_WPTight = "(HLT_Ele27_WPTight_Gsf_v)(.*)";
  Mu8_Ele23_IsoVL = "(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v)(.*)";
  Mu23_Ele12_IsoVL = "(HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v)(.*)";
  Mu8_Ele23_DZ = "(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v)(.*)";
  Mu23_Ele12_DZ = "(HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v)(.*)";

  Mu17_Mu8_DZ_Mass8 = "(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v)(.*)";
  IsoMu27 = "(HLT_IsoMu27_v)(.*)";
  Ele23_Ele12_IsoVL = "(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v)(.*)";
  Ele35_WPTight = "(HLT_Ele35_WPTight_Gsf_v)(.*)";

  Mu17_Mu8_DZ_Mass3p8 = "(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v)(.*)";
  Ele32_WPTight = "(HLT_Ele32_WPTight_Gsf_v)(.*)";
  Double_Ele25_MW = "(HLT_DoubleEle25_CaloIdL_MW_v)(.*)";

  ntuple_ = fileService->make<TTree>("Ntuple", "Ntuple");
  ntuple_->Branch("event", &event_);
  ntuple_->Branch("nPVertex", &nPVertex_);
  ntuple_->Branch("mVertex1", &mVertex1_);
  ntuple_->Branch("mVertex2", &mVertex2_);
  ntuple_->Branch("zPVertex", &zPVertex_);
  ntuple_->Branch("vxTrack", &vxTrack_);
  ntuple_->Branch("vyTrack", &vyTrack_);
  ntuple_->Branch("dz_trackPV", &dz_trackPV_);
  ntuple_->Branch("pfMET", &pfMET_);
  ntuple_->Branch("nLepGen", &nLepGen_);
  ntuple_->Branch("nJet", &nJet_);
  ntuple_->Branch("hlt_Mu17_Mu8_VVL", &hlt_Mu17_Mu8_VVL_);
  ntuple_->Branch("hlt_Mu17_TkMu8_VVL", &hlt_Mu17_TkMu8_VVL_);
  ntuple_->Branch("hlt_Mu17_Mu8_DZ", &hlt_Mu17_Mu8_DZ_);
  ntuple_->Branch("hlt_Mu17_TkMu8_DZ", &hlt_Mu17_TkMu8_DZ_);
  ntuple_->Branch("hlt_IsoMu24", &hlt_IsoMu24_);
  ntuple_->Branch("hlt_IsoTkMu24", &hlt_IsoTkMu24_);
  ntuple_->Branch("hlt_Ele23_Ele12_DZ", &hlt_Ele23_Ele12_DZ_);
  ntuple_->Branch("hlt_Double_Ele33_MW", &hlt_Double_Ele33_MW_);
  ntuple_->Branch("hlt_Double_Ele33_VL", &hlt_Double_Ele33_VL_);
  ntuple_->Branch("hlt_Ele27_WPTight", &hlt_Ele27_WPTight_);
  ntuple_->Branch("hlt_Mu8_Ele23_IsoVL", &hlt_Mu8_Ele23_IsoVL_);
  ntuple_->Branch("hlt_Mu23_Ele12_IsoVL", &hlt_Mu23_Ele12_IsoVL_);
  ntuple_->Branch("hlt_Mu8_Ele23_DZ", &hlt_Mu8_Ele23_DZ_);
  ntuple_->Branch("hlt_Mu23_Ele12_DZ", &hlt_Mu23_Ele12_DZ_);
  ntuple_->Branch("hlt_Mu17_Mu8_DZ_Mass8", &hlt_Mu17_Mu8_DZ_Mass8_);
  ntuple_->Branch("hlt_IsoMu27", &hlt_IsoMu27_);
  ntuple_->Branch("hlt_Ele23_Ele12_IsoVL", &hlt_Ele23_Ele12_IsoVL_);
  ntuple_->Branch("hlt_Ele35_WPTight", &hlt_Ele35_WPTight_);
  ntuple_->Branch("hlt_Mu17_Mu8_DZ_Mass3p8", &hlt_Mu17_Mu8_DZ_Mass3p8_);
  ntuple_->Branch("hlt_Ele32_WPTight", &hlt_Ele32_WPTight_);
  ntuple_->Branch("hlt_Double_Ele25_MW", &hlt_Double_Ele25_MW_);
  ntuple_->Branch("p4_muon", &p4_muon_);
  ntuple_->Branch("p4_electron", &p4_electron_);
  ntuple_->Branch("charge_muon", &charge_muon_);
  ntuple_->Branch("charge_electron", &charge_electron_);
  ntuple_->Branch("pfIso_muon", &pfIso_muon_);
  ntuple_->Branch("passLooseId_muon", &passLooseId_muon_);
  ntuple_->Branch("passMediumId_muon", &passMediumId_muon_);
  ntuple_->Branch("passTightId_muon", &passTightId_muon_);
  ntuple_->Branch("passVetoId_electron", &passVetoId_electron_);
  ntuple_->Branch("passLooseId_electron", &passLooseId_electron_);
  ntuple_->Branch("passMediumId_electron", &passMediumId_electron_);
  ntuple_->Branch("passTightId_electron", &passTightId_electron_);
  //ntuple_->Branch("passIdbits_electron", &passIdbits_electron_);
  //ntuple_->Branch("passId_electron", &passId_electron_);
  ntuple_->Branch("foundGenBlepton", &foundGenBlepton_);
  ntuple_->Branch("foundGenBPlusdecay", &foundGenBPlusdecay_);
  ntuple_->Branch("foundGenB0Dstardecay", &foundGenB0Dstardecay_);
  ntuple_->Branch("foundGenB0decay", &foundGenB0decay_);
  ntuple_->Branch("foundGenBs0decay", &foundGenBs0decay_);
  ntuple_->Branch("foundGenLambdab0decay", &foundGenLambdab0decay_);
  ntuple_->Branch("p4_gen_D0_lepton", &p4_gen_D0_lepton_);
  ntuple_->Branch("p4_gen_Dstar_lepton", &p4_gen_Dstar_lepton_);
  ntuple_->Branch("p4_gen_DPlus_lepton", &p4_gen_DPlus_lepton_);
  ntuple_->Branch("p4_gen_DsPlus_lepton", &p4_gen_DsPlus_lepton_);
  ntuple_->Branch("p4_gen_LambdacPlus_lepton", &p4_gen_LambdacPlus_lepton_);
  ntuple_->Branch("p4_jet", &p4_jet_);
  ntuple_->Branch("lepton_perJet", &lepton_perJet_);
  ntuple_->Branch("pdgId_lepton", &pdgId_lepton_);
  ntuple_->Branch("pt_lepton", &pt_lepton_);
  ntuple_->Branch("eta_lepton", &eta_lepton_);
  ntuple_->Branch("phi_lepton", &phi_lepton_);
  ntuple_->Branch("charge_lepton", &charge_lepton_);
  ntuple_->Branch("mass_lepton", &mass_lepton_);
  ntuple_->Branch("p4_D0", &p4_D0_);
  ntuple_->Branch("p4_D0_kaon", &p4_D0_kaon_);
  ntuple_->Branch("p4_D0_pion", &p4_D0_pion_);
  ntuple_->Branch("charge_D0_kaon", &charge_D0_kaon_);
  ntuple_->Branch("charge_D0_pion", &charge_D0_pion_);
  ntuple_->Branch("trkpara_D0_kaon", &trkpara_D0_kaon_);
  ntuple_->Branch("trkpara_D0_pion", &trkpara_D0_pion_);
  ntuple_->Branch("trkcovm_D0_kaon", &trkcovm_D0_kaon_);
  ntuple_->Branch("trkcovm_D0_pion", &trkcovm_D0_pion_);
  ntuple_->Branch("mdiff_Dstar_D0", &mdiff_Dstar_D0_);
  ntuple_->Branch("p4_Dstar", &p4_Dstar_);
  ntuple_->Branch("p4_Dstar_D0", &p4_Dstar_D0_);
  ntuple_->Branch("p4_Dstar_D0_kaon", &p4_Dstar_D0_kaon_);
  ntuple_->Branch("p4_Dstar_D0_pion", &p4_Dstar_D0_pion_);
  ntuple_->Branch("p4_Dstar_pionsoft", &p4_Dstar_pionsoft_);
  ntuple_->Branch("charge_Dstar_D0_kaon", &charge_Dstar_D0_kaon_);
  ntuple_->Branch("charge_Dstar_D0_pion", &charge_Dstar_D0_pion_);
  ntuple_->Branch("charge_Dstar_pionsoft", &charge_Dstar_pionsoft_);
  ntuple_->Branch("trkpara_Dstar_D0_kaon", &trkpara_Dstar_D0_kaon_);
  ntuple_->Branch("trkpara_Dstar_D0_pion", &trkpara_Dstar_D0_pion_);
  ntuple_->Branch("trkpara_Dstar_pionsoft", &trkpara_Dstar_pionsoft_);
  ntuple_->Branch("trkcovm_Dstar_D0_kaon", &trkcovm_Dstar_D0_kaon_);
  ntuple_->Branch("trkcovm_Dstar_D0_pion", &trkcovm_Dstar_D0_pion_);
  ntuple_->Branch("trkcovm_Dstar_pionsoft", &trkcovm_Dstar_pionsoft_);
  ntuple_->Branch("p4_DPlus", &p4_DPlus_);
  ntuple_->Branch("p4_DPlus_kaon", &p4_DPlus_kaon_);
  ntuple_->Branch("p4_DPlus_pion1", &p4_DPlus_pion1_);
  ntuple_->Branch("p4_DPlus_pion2", &p4_DPlus_pion2_);
  ntuple_->Branch("charge_DPlus_kaon", &charge_DPlus_kaon_);
  ntuple_->Branch("charge_DPlus_pion1", &charge_DPlus_pion1_);
  ntuple_->Branch("charge_DPlus_pion2", &charge_DPlus_pion2_);
  ntuple_->Branch("trkpara_DPlus_kaon", &trkpara_DPlus_kaon_);
  ntuple_->Branch("trkpara_DPlus_pion1", &trkpara_DPlus_pion1_);
  ntuple_->Branch("trkpara_DPlus_pion2", &trkpara_DPlus_pion2_);
  ntuple_->Branch("trkcovm_DPlus_kaon", &trkcovm_DPlus_kaon_);
  ntuple_->Branch("trkcovm_DPlus_pion1", &trkcovm_DPlus_pion1_);
  ntuple_->Branch("trkcovm_DPlus_pion2", &trkcovm_DPlus_pion2_);
  ntuple_->Branch("p4_DsPlus", &p4_DsPlus_);
  ntuple_->Branch("p4_DsPlus_kaon1", &p4_DsPlus_kaon1_);
  ntuple_->Branch("p4_DsPlus_kaon2", &p4_DsPlus_kaon2_);
  ntuple_->Branch("p4_DsPlus_pion", &p4_DsPlus_pion_);
  ntuple_->Branch("charge_DsPlus_kaon1", &charge_DsPlus_kaon1_);
  ntuple_->Branch("charge_DsPlus_kaon2", &charge_DsPlus_kaon2_);
  ntuple_->Branch("charge_DsPlus_pion", &charge_DsPlus_pion_);
  ntuple_->Branch("trkpara_DsPlus_kaon1", &trkpara_DsPlus_kaon1_);
  ntuple_->Branch("trkpara_DsPlus_kaon2", &trkpara_DsPlus_kaon2_);
  ntuple_->Branch("trkpara_DsPlus_pion", &trkpara_DsPlus_pion_);
  ntuple_->Branch("trkcovm_DsPlus_kaon1", &trkcovm_DsPlus_kaon1_);
  ntuple_->Branch("trkcovm_DsPlus_kaon2", &trkcovm_DsPlus_kaon2_);
  ntuple_->Branch("trkcovm_DsPlus_pion", &trkcovm_DsPlus_pion_);
  ntuple_->Branch("p4_LambdacPlus", &p4_LambdacPlus_);
  ntuple_->Branch("p4_LambdacPlus_kaon", &p4_LambdacPlus_kaon_);
  ntuple_->Branch("p4_LambdacPlus_pion", &p4_LambdacPlus_pion_);
  ntuple_->Branch("p4_LambdacPlus_proton", &p4_LambdacPlus_proton_);
  ntuple_->Branch("charge_LambdacPlus_kaon", &charge_LambdacPlus_kaon_);
  ntuple_->Branch("charge_LambdacPlus_pion", &charge_LambdacPlus_pion_);
  ntuple_->Branch("charge_LambdacPlus_proton", &charge_LambdacPlus_proton_);
  ntuple_->Branch("trkpara_LambdacPlus_kaon", &trkpara_LambdacPlus_kaon_);
  ntuple_->Branch("trkpara_LambdacPlus_pion", &trkpara_LambdacPlus_pion_);
  ntuple_->Branch("trkpara_LambdacPlus_proton", &trkpara_LambdacPlus_proton_);
  ntuple_->Branch("trkcovm_LambdacPlus_kaon", &trkcovm_LambdacPlus_kaon_);
  ntuple_->Branch("trkcovm_LambdacPlus_pion", &trkcovm_LambdacPlus_pion_);
  ntuple_->Branch("trkcovm_LambdacPlus_proton", &trkcovm_LambdacPlus_proton_);
  ntuple_->Branch("pv_lxy_D0", &pv_lxy_D0_);
  ntuple_->Branch("pv_lxy_Dstar_D0", &pv_lxy_Dstar_D0_);
  ntuple_->Branch("pv_lxy_DPlus", &pv_lxy_DPlus_);
  ntuple_->Branch("pv_lxy_DsPlus", &pv_lxy_DsPlus_);
  ntuple_->Branch("pv_lxy_LambdacPlus", &pv_lxy_LambdacPlus_);
}

TopAnalyzer::~TopAnalyzer()
{
}

// ------------ method called for each event  ------------
void
TopAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  run_ = iEvent.id().run();
  event_ = iEvent.id().event();
  lumi_ = iEvent.id().luminosityBlock();

  //if(lumi_==21 || lumi_==22) {
  //cout<<"TopAnalyzer::analyze()"<<endl;
  //cout<<"Run "<<run_<<"  Event "<<event_<<"  Lumi "<<lumi_<<endl;

  edm::Handle<edm::View<reco::Vertex> > vertexHandle_;
  iEvent.getByToken(vertexToken_, vertexHandle_);
  //cout<<"Vertex collection size = "<<vertexHandle_->size()<<endl;
  nPVertex_ = vertexHandle_->size();
  if(nPVertex_ == 0) return;
  int iv = 0;
  double zpv = 0;
  _pv = NULL;
  View<reco::Vertex>::const_iterator iVertex;
  for(iVertex = vertexHandle_->begin(); iVertex != vertexHandle_->end(); iVertex++) {
    if(iv == 2) break;
    if(iv == 0) {
      mVertex1_ = iVertex->p4().mass();
      zpv = iVertex->z();
      _pv = &(*iVertex);
    }
    if(iv == 1) {
      mVertex2_ = iVertex->p4().mass();
    }
    iv += 1;
  }
  zPVertex_ = zpv;

  edm::Handle<edm::View<pat::PackedCandidate>> trackHandle_;
  iEvent.getByToken(tracksToken_,trackHandle_);
  //cout<<"Track collection size = "<<trackHandle_->size()<<endl;
  for(View<pat::PackedCandidate>::const_iterator iTrack = trackHandle_->begin();
        iTrack != trackHandle_->end(); iTrack++) {
    vxTrack_ = iTrack->vx();
    vyTrack_ = iTrack->vy();
    dz_trackPV_ = iTrack->vz()-zpv;
  }

  edm::Handle<edm::View<pat::MET> > metHandle_;
  iEvent.getByToken(metToken_,metHandle_);
  //cout<<"MET collection size = "<<metHandle_->size()<<endl;
  const pat::MET &met = metHandle_->front();
  pfMET_ = met.pt();

  edm::Handle<edm::TriggerResults > triggerHandle_;
  iEvent.getByToken(triggerToken_,triggerHandle_);
  const edm::TriggerNames &trigNames = iEvent.triggerNames(*triggerHandle_);
  for(unsigned int i=0; i<triggerHandle_->size(); i++)
  {
    if(std::regex_match(trigNames.triggerName(i),Mu17_Mu8_VVL)) hlt_Mu17_Mu8_VVL_ = triggerHandle_->accept(i);
    if(std::regex_match(trigNames.triggerName(i),Mu17_TkMu8_VVL)) hlt_Mu17_TkMu8_VVL_ = triggerHandle_->accept(i);
    if(std::regex_match(trigNames.triggerName(i),Mu17_Mu8_DZ)) hlt_Mu17_Mu8_DZ_ = triggerHandle_->accept(i);
    if(std::regex_match(trigNames.triggerName(i),Mu17_TkMu8_DZ)) hlt_Mu17_TkMu8_DZ_ = triggerHandle_->accept(i);
    if(std::regex_match(trigNames.triggerName(i),IsoMu24)) hlt_IsoMu24_ = triggerHandle_->accept(i);
    if(std::regex_match(trigNames.triggerName(i),IsoTkMu24)) hlt_IsoTkMu24_ = triggerHandle_->accept(i);
    if(std::regex_match(trigNames.triggerName(i),Ele23_Ele12_DZ)) hlt_Ele23_Ele12_DZ_ = triggerHandle_->accept(i);
    if(std::regex_match(trigNames.triggerName(i),Double_Ele33_MW)) hlt_Double_Ele33_MW_ = triggerHandle_->accept(i);
    if(std::regex_match(trigNames.triggerName(i),Double_Ele33_VL)) hlt_Double_Ele33_VL_ = triggerHandle_->accept(i);
    if(std::regex_match(trigNames.triggerName(i),Ele27_WPTight)) hlt_Ele27_WPTight_ = triggerHandle_->accept(i);
    if(std::regex_match(trigNames.triggerName(i),Mu8_Ele23_IsoVL)) hlt_Mu8_Ele23_IsoVL_ = triggerHandle_->accept(i);
    if(std::regex_match(trigNames.triggerName(i),Mu23_Ele12_IsoVL)) hlt_Mu23_Ele12_IsoVL_ = triggerHandle_->accept(i);
    if(std::regex_match(trigNames.triggerName(i),Mu8_Ele23_DZ)) hlt_Mu8_Ele23_DZ_ = triggerHandle_->accept(i);
    if(std::regex_match(trigNames.triggerName(i),Mu23_Ele12_DZ)) hlt_Mu23_Ele12_DZ_ = triggerHandle_->accept(i);
    if(std::regex_match(trigNames.triggerName(i),Mu17_Mu8_DZ_Mass8)) hlt_Mu17_Mu8_DZ_Mass8_ = triggerHandle_->accept(i);
    if(std::regex_match(trigNames.triggerName(i),IsoMu27)) hlt_IsoMu27_ = triggerHandle_->accept(i);
    if(std::regex_match(trigNames.triggerName(i),Ele23_Ele12_IsoVL)) hlt_Ele23_Ele12_IsoVL_ = triggerHandle_->accept(i);
    if(std::regex_match(trigNames.triggerName(i),Ele35_WPTight)) hlt_Ele35_WPTight_ = triggerHandle_->accept(i);
    if(std::regex_match(trigNames.triggerName(i),Mu17_Mu8_DZ_Mass3p8)) hlt_Mu17_Mu8_DZ_Mass3p8_ = triggerHandle_->accept(i);
    if(std::regex_match(trigNames.triggerName(i),Ele32_WPTight)) hlt_Ele32_WPTight_ = triggerHandle_->accept(i);
    if(std::regex_match(trigNames.triggerName(i),Double_Ele25_MW)) hlt_Double_Ele25_MW_ = triggerHandle_->accept(i);
  }

  edm::Handle<edm::View<pat::Muon> > muonsHandle_;
  iEvent.getByToken(muonsToken_,muonsHandle_);
  //cout<<"Muon collection size = "<<muonsHandle_->size()<<endl;
  for(size_t i=0; i<muonsHandle_->size(); ++i) {
    const auto mu = muonsHandle_->ptrAt(i);
    if(mu->pt() < 20.0 ) continue;
    if(fabs(mu->eta()) > 2.4) continue;
    double isomu = (mu->pfIsolationR04().sumChargedHadronPt + max(0.,
                    mu->pfIsolationR04().sumNeutralHadronEt +
                    mu->pfIsolationR04().sumPhotonEt - 0.5 *
                    mu->pfIsolationR04().sumPUPt))/mu->pt();
    if(isomu > 0.25) continue;
    //cout<<"  muon: "<<mu->pt()<<"  "<<mu->charge()<<", ";
    p4_muon_.push_back(this->getLV(mu->p4()));
    charge_muon_.push_back(mu->charge());
    passLooseId_muon_.push_back(muon::isLooseMuon(*mu));
    passMediumId_muon_.push_back(muon::isMediumMuon(*mu));
    passTightId_muon_.push_back(muon::isTightMuon(*mu,*_pv));
    pfIso_muon_.push_back(isomu);
    muons.push_back(&(*mu));
  }
  //cout<<endl;
  edm::Handle<edm::View<pat::Electron> > electronsHandle_;
  iEvent.getByToken(electronsToken_,electronsHandle_);
  //cout<<"Electron collection size = "<<electronsHandle_->size()<<endl;
  bool barrel, endcap = false;
  for(size_t i=0; i<electronsHandle_->size(); ++i) {
    const auto el = electronsHandle_->ptrAt(i);
    if(el->pt() < 20.0) continue;
    if(fabs(el->eta()) > 2.4) continue;
    if(fabs(el->superCluster()->eta()) > 1.4442 &&
       fabs(el->superCluster()->eta()) < 1.5660) continue;
    if(fabs(el->superCluster()->eta()) < 1.4442 &&
       fabs(el->gsfTrack()->dxy(_pv->position())) < 0.05 &&
       fabs(el->gsfTrack()->dz(_pv->position())) < 0.10) barrel = true;
    if(fabs(el->superCluster()->eta()) > 1.5660 &&
       fabs(el->gsfTrack()->dxy(_pv->position())) < 0.10 &&
       fabs(el->gsfTrack()->dz(_pv->position())) < 0.20) endcap = true;
    if(!barrel) continue;
    if(!endcap) continue;
    p4_electron_.push_back(this->getLV(el->p4()));
    charge_electron_.push_back(el->charge());
    UShort_t tmpeleIdbit = 0;   
    bool isPassVeto = el->electronID("cutBasedElectronID-Fall17-94X-V2-veto");
    if(isPassVeto) setbit(tmpeleIdbit, 0);    
    bool isPassLoose = el->electronID("cutBasedElectronID-Fall17-94X-V2-loose");
    if(isPassLoose) setbit(tmpeleIdbit, 1);   
    bool isPassMedium = el->electronID("cutBasedElectronID-Fall17-94X-V2-medium");
    if(isPassMedium) setbit(tmpeleIdbit, 2);    
    bool isPassTight = el->electronID("cutBasedElectronID-Fall17-94X-V2-tight");
    if(isPassTight) setbit(tmpeleIdbit, 3);
    //cout<<"  electron: "<<el->pt()<<"  "<<el->charge()<<", ";
    passVetoId_electron_.push_back(isPassVeto);
    passLooseId_electron_.push_back(isPassLoose);
    passMediumId_electron_.push_back(isPassMedium);
    passTightId_electron_.push_back(isPassTight);
    //passIdbits_electron_.push_back(tmpeleIdbit);
    //passId_electron_.push_back({
      //el->userInt("cutBasedElectronID-Fall17-94X-V2-veto"),
      //el->userInt("cutBasedElectronID-Fall17-94X-V2-loose"),
      //el->userInt("cutBasedElectronID-Fall17-94X-V2-medium"),
      //el->userInt("cutBasedElectronID-Fall17-94X-V2-tight")});
    electrons.push_back(&(*el));
  }
  //cout<<endl;

  edm::Handle<edm::View<pat::Jet> > jetsHandle_;
  iEvent.getByToken(jetsToken_,jetsHandle_);
  edm::Handle<edm::View<reco::GenParticle>> genHandle_;
  iEvent.getByToken(genToken_,genHandle_);
  edm::Handle<edm::View<pat::PackedGenParticle>> packedGenHandle_;
  iEvent.getByToken(packedGenToken_,packedGenHandle_);
  if(misMC) {
    GenAnalysis(iEvent,iSetup,genHandle_,packedGenHandle_);
    RecoAnalysis(iEvent,iSetup,trackHandle_,jetsHandle_);
  }
  else {
    DataAnalysis(iEvent,iSetup,jetsHandle_);
  }
  //cout<<""<<endl;
  ntuple_->Fill();

  // Clear the variables
  nkpi = 0;
  nkpipi = 0;
  nkkpi = 0;
  ndstar = 0;
  npkpi = 0;
  top_quarks.clear();
  bottom_hadrons.clear();
  charm_hadrons.clear();
  leptons.clear();
  d0_hadron.clear();
  lepton_d0.clear();
  kaon_d0.clear();
  pion_d0.clear();
  dt_hadron.clear();
  lepton_dt.clear();
  kaon_dt.clear();
  pion_dt.clear();
  pionsoft_dt.clear();
  dp_hadron.clear();
  lepton_dp.clear();
  kaon_dp.clear();
  pion1_dp.clear();
  pion2_dp.clear();
  ds_hadron.clear();
  lepton_ds.clear();
  kaon1_ds.clear();
  kaon2_ds.clear();
  pion_ds.clear();
  lambda_hadron.clear();
  lepton_lambda.clear();
  kaon_lambda.clear();
  pion_lambda.clear();
  proton_lambda.clear();
  muons.clear();
  electrons.clear();
  par.clear();
  covm.clear();
  kd0_par.clear();
  pid0_par.clear();
  kdt_par.clear();
  pidt_par.clear();
  pisdt_par.clear();
  kdp_par.clear();
  pi1dp_par.clear();
  pi2dp_par.clear();
  k1ds_par.clear();
  k2ds_par.clear();
  pids_par.clear();
  klambda_par.clear();
  pilambda_par.clear();
  plambda_par.clear();
  kd0_cov.clear();
  pid0_cov.clear();
  kdt_cov.clear();
  pidt_cov.clear();
  pisdt_cov.clear();
  kdp_cov.clear();
  pi1dp_cov.clear();
  pi2dp_cov.clear();
  k1ds_cov.clear();
  k2ds_cov.clear();
  pids_cov.clear();
  klambda_cov.clear();
  pilambda_cov.clear();
  plambda_cov.clear();

  // Branch variables
  run_ = 0;
  lumi_ = 0;
  event_ = 0;
  nPVertex_ = 0;
  mVertex1_ = 0;
  mVertex2_ = 0;
  zPVertex_ = 0;
  vxTrack_ = 0;
  vyTrack_ = 0;
  dz_trackPV_ = 0;
  pfMET_ = 0;
  nLepGen_ = 0;
  nJet_ = 0;
  p4_muon_.clear();
  p4_electron_.clear();
  charge_muon_.clear();
  charge_electron_.clear();
  passLooseId_muon_.clear();
  passMediumId_muon_.clear();
  passTightId_muon_.clear();
  passVetoId_electron_.clear();
  passLooseId_electron_.clear();
  passMediumId_electron_.clear();
  passTightId_electron_.clear();
  pfIso_muon_.clear();
  //passIdbits_electron_.clear();
  //passId_electron_.clear();
  foundGenBlepton_ = false;
  foundGenBPlusdecay_ = false;
  foundGenB0Dstardecay_ = false;
  foundGenB0decay_ = false;
  foundGenBs0decay_ = false;
  foundGenLambdab0decay_ = false;
  p4_gen_D0_lepton_.clear();
  p4_gen_Dstar_lepton_.clear();
  p4_gen_DPlus_lepton_.clear();
  p4_gen_DsPlus_lepton_.clear();
  p4_gen_LambdacPlus_lepton_.clear();
  p4_jet_.clear();
  lepton_perJet_.clear();
  pdgId_lepton_.clear();
  pt_lepton_.clear();
  eta_lepton_.clear();
  phi_lepton_.clear();
  charge_lepton_.clear();
  mass_lepton_.clear();
  p4_D0_.clear();
  p4_D0_kaon_.clear();
  p4_D0_pion_.clear();
  charge_D0_kaon_.clear();
  charge_D0_pion_.clear();
  mdiff_Dstar_D0_.clear();
  p4_Dstar_.clear();
  p4_Dstar_D0_.clear();
  p4_Dstar_D0_kaon_.clear();
  p4_Dstar_D0_pion_.clear();
  p4_Dstar_pionsoft_.clear();
  charge_Dstar_D0_kaon_.clear();
  charge_Dstar_D0_pion_.clear();
  charge_Dstar_pionsoft_.clear();
  p4_DPlus_.clear();
  p4_DPlus_kaon_.clear();
  p4_DPlus_pion1_.clear();
  p4_DPlus_pion2_.clear();
  charge_DPlus_kaon_.clear();
  charge_DPlus_pion1_.clear();
  charge_DPlus_pion2_.clear();
  p4_DsPlus_.clear();
  p4_DsPlus_kaon1_.clear();
  p4_DsPlus_kaon2_.clear();
  p4_DsPlus_pion_.clear();
  charge_DsPlus_kaon1_.clear();
  charge_DsPlus_kaon2_.clear();
  charge_DsPlus_pion_.clear();
  p4_LambdacPlus_.clear();
  p4_LambdacPlus_kaon_.clear();
  p4_LambdacPlus_pion_.clear();
  p4_LambdacPlus_proton_.clear();
  charge_LambdacPlus_kaon_.clear();
  charge_LambdacPlus_pion_.clear();
  charge_LambdacPlus_proton_.clear();
  trkpara_D0_kaon_.clear();
  trkpara_D0_pion_.clear();
  trkcovm_D0_kaon_.clear();
  trkcovm_D0_pion_.clear();
  trkpara_Dstar_D0_kaon_.clear();
  trkpara_Dstar_D0_pion_.clear();
  trkpara_Dstar_pionsoft_.clear();
  trkcovm_Dstar_D0_kaon_.clear();
  trkcovm_Dstar_D0_pion_.clear();
  trkcovm_Dstar_pionsoft_.clear();
  trkpara_DPlus_kaon_.clear();
  trkpara_DPlus_pion1_.clear();
  trkpara_DPlus_pion2_.clear();
  trkcovm_DPlus_kaon_.clear();
  trkcovm_DPlus_pion1_.clear();
  trkcovm_DPlus_pion2_.clear();
  trkpara_DsPlus_kaon1_.clear();
  trkpara_DsPlus_kaon2_.clear();
  trkpara_DsPlus_pion_.clear();
  trkcovm_DsPlus_kaon1_.clear();
  trkcovm_DsPlus_kaon2_.clear();
  trkcovm_DsPlus_pion_.clear();
  trkpara_LambdacPlus_kaon_.clear();
  trkpara_LambdacPlus_pion_.clear();
  trkpara_LambdacPlus_proton_.clear();
  trkcovm_LambdacPlus_kaon_.clear();
  trkcovm_LambdacPlus_pion_.clear();
  trkcovm_LambdacPlus_proton_.clear();
  pv_lxy_D0_.clear();
  pv_lxy_Dstar_D0_.clear();
  pv_lxy_DPlus_.clear();
  pv_lxy_DsPlus_.clear();
  pv_lxy_LambdacPlus_.clear();
//}
}

void TopAnalyzer::GenAnalysis(const edm::Event& iEvent,
                              const edm::EventSetup& iSetup,
                              edm::Handle<edm::View<reco::GenParticle> > &handle,
                              edm::Handle<edm::View<pat::PackedGenParticle> > &handle1) {
  //cout<<"Gen collection size = "<<handle->size()<<endl;
  int ij = 0;
  for(size_t igen=0; igen<handle->size(); igen++) {
    const Candidate *lepton = &(*handle)[igen];
    if(((abs(lepton->pdgId()) == 11 ||
         abs(lepton->pdgId()) == 13) &&
         lepton->status() == 1) ||
        (abs(lepton->pdgId()) == 15 && 
        (lepton->status() == 2 || 
         lepton->status() == 1))) {
      if(!isAncestor(lepton,6)) continue;
      if(!isAncestor(lepton,24)) continue;
      ij += 1;
      //cout<<"  lepton: "<<lepton->pdgId()<<"  "<<lepton->pt()<<"  "<<lepton->eta()<<endl;
    }
    const Candidate *bhadron = &(*handle)[igen];
    if(abs(bhadron->pdgId()) == 511 ||
       abs(bhadron->pdgId()) == 521 ||
       abs(bhadron->pdgId()) == 531 ||
       abs(bhadron->pdgId()) == 5122) {
      if(!isAncestor(bhadron,6)) continue;
      //if(fabs(bhadron->eta()) > 2.5) continue;
      bottom_hadrons.push_back(bhadron);
      for(size_t ilep=0; ilep<bhadron->numberOfDaughters(); ilep++) {
        const Candidate *lepton = bhadron->daughter(ilep);
        if(abs(lepton->pdgId()) == 11 ||
           abs(lepton->pdgId()) == 13) {
          if(lepton->pt() < 3.) continue;
          if(fabs(lepton->eta()) > 2.5) continue;
          foundGenBlepton_ = true;
          leptons.push_back(lepton);
          for(size_t ichad=0; ichad<bhadron->numberOfDaughters(); ichad++) {
            if(ilep == ichad) continue;
            const Candidate *chadron = bhadron->daughter(ichad);
            if(abs(chadron->pdgId()) == 411 ||
               abs(chadron->pdgId()) == 413 ||
               abs(chadron->pdgId()) == 421 ||
               abs(chadron->pdgId()) == 431 ||
               abs(chadron->pdgId()) == 4122) {
              charm_hadrons.push_back(chadron);
              if(abs(chadron->pdgId()) == 421) {
                for(size_t idau=0; idau<chadron->numberOfDaughters(); idau++) {
                  const Candidate *k = chadron->daughter(idau);
                  if(lepton->pdgId()*k->pdgId() == -11*321 ||
                     lepton->pdgId()*k->pdgId() == -13*321) {
                    for(size_t jdau=0; jdau<chadron->numberOfDaughters(); jdau++) {
                      if(idau == jdau) continue;
                      const Candidate *pi = chadron->daughter(jdau);
                      if(k->pdgId()*pi->pdgId() != -321*211) continue;
                      if(k->pt() < 1. || pi->pt() < 1.) continue;
                      if(fabs(k->eta()) > 2.5 || fabs(pi->eta()) > 2.5) continue;
                      foundGenBPlusdecay_ =  true;
                      d0_hadron.push_back(chadron);
                      lepton_d0.push_back(lepton);
                      kaon_d0.push_back(k);
                      pion_d0.push_back(pi);
                      p4_gen_D0_lepton_.push_back(this->getLV(lepton->p4()));
                    } // pion candidate
                  } // lepton and kaon have same charge
                } // kaon candidate
              } // d0 candidate
              if(abs(chadron->pdgId()) == 413) {
                for(size_t idau=0; idau<chadron->numberOfDaughters(); idau++) {
                  const Candidate *pis = chadron->daughter(idau);
                  if(abs(pis->pdgId()) == 211) {
                    for(size_t cdau = 0; cdau<chadron->numberOfDaughters(); cdau++) {
                      if(idau == cdau) continue;
                      if(abs(chadron->daughter(cdau)->pdgId()) == 421) {
                        for(size_t jdau=0; jdau<chadron->daughter(cdau)->numberOfDaughters(); jdau++) {
                          const Candidate *k = chadron->daughter(cdau)->daughter(jdau);
                          if(lepton->pdgId()*k->pdgId() == -11*321 ||
                             lepton->pdgId()*k->pdgId() == -13*321) {
                            for(size_t kdau=0; kdau<chadron->daughter(cdau)->numberOfDaughters(); kdau++) {
                              if(jdau == kdau) continue;
                              const Candidate *pi = chadron->daughter(cdau)->daughter(kdau);
                              if(k->pdgId()*pi->pdgId() != -321*211) continue;
                              if(k->pt() < 1. || pi->pt() < 1. || pis->pt() < 0.5) continue;
                              if(fabs(k->eta()) > 2.5 || fabs(pi->eta()) > 2.5 || fabs(pis->eta()) > 2.5) continue;
                              foundGenB0Dstardecay_ = true;
                              dt_hadron.push_back(chadron);
                              lepton_dt.push_back(lepton);
                              kaon_dt.push_back(k);
                              pion_dt.push_back(pi);
                              pionsoft_dt.push_back(pis);
                              p4_gen_Dstar_lepton_.push_back(this->getLV(lepton->p4()));
                            } // pion candidate
                          } // lepton and kaon have same charge
                        } // kaon candidate
                      } // d0 candidate
                    } // d*+ daughters
                  } // soft pion
                } // d*+ daughters
              } // d*+ candidate
              if(abs(chadron->pdgId()) == 411) {
                for(size_t idau=0; idau<chadron->numberOfDaughters(); idau++) {
                  const Candidate *k = chadron->daughter(idau);
                  if(lepton->pdgId()*k->pdgId() == -11*321 ||
                     lepton->pdgId()*k->pdgId() == -13*321) {
                    for(size_t jdau=0; jdau<chadron->numberOfDaughters(); jdau++) {
                      if(idau == jdau) continue;
                      const Candidate *pi1 = chadron->daughter(jdau);
                      if(k->pdgId()*pi1->pdgId() != -321*211) continue;
                      for(size_t kdau=0; kdau<jdau; kdau++) {
                        if(jdau == kdau) continue;
                        const Candidate *pi2 = chadron->daughter(kdau);
                        if(k->pdgId()*pi2->pdgId() != -321*211) continue;
                        if(k->pt() < 1. || pi1->pt() < 1. || pi2->pt() < 1.0) continue;
                        if(fabs(k->eta()) > 2.5 || fabs(pi1->eta()) > 2.5 || fabs(pi2->eta()) > 2.5) continue;
                        foundGenB0decay_  = true;
                        dp_hadron.push_back(chadron);
                        lepton_dp.push_back(lepton);
                        kaon_dp.push_back(k);
                        pion1_dp.push_back(pi1);
                        pion2_dp.push_back(pi2);
                        p4_gen_DPlus_lepton_.push_back(this->getLV(lepton->p4()));
                      } // second pion candidate
                    } // first pion candidate
                  } // lepton and kaon have same charge
                } // kaon candidate
              } // d+ chandidate
              if(abs(chadron->pdgId()) == 431) {
                for(size_t idau=0; idau<chadron->numberOfDaughters(); idau++) {
                  const Candidate *pi = chadron->daughter(idau);
                  if(lepton->pdgId()*pi->pdgId() == 11*211 ||
                     lepton->pdgId()*pi->pdgId() == 13*211) {
                    for(size_t jdau=0; jdau<chadron->numberOfDaughters(); jdau++) {
                      if(idau == jdau) continue;
                      const Candidate *k1 = chadron->daughter(jdau);
                      for(size_t kdau=0; kdau<jdau; kdau++) {
                        if(jdau == kdau) continue;
                        const Candidate *k2 = chadron->daughter(kdau);
                        if(k1->pdgId()*k2->pdgId() != -321*321) continue;
                        if(k1->pt() < 1. || k2->pt() < 1. || pi->pt() < 1.0) continue;
                        if(fabs(k1->eta()) > 2.5 || fabs(k2->eta()) > 2.5 || fabs(pi->eta()) > 2.5) continue;
                        foundGenBs0decay_ = true;
                        ds_hadron.push_back(chadron);
                        lepton_ds.push_back(lepton);
                        kaon1_ds.push_back(k1);
                        kaon2_ds.push_back(k2);
                        pion_ds.push_back(pi);
                        p4_gen_DsPlus_lepton_.push_back(this->getLV(lepton->p4()));
                      } // second kaon candidate
                    } // first kaon candidate
                  } // lepton and pion have opposite charge
                } // pion candidate 
              } // ds+ candidate
              if(abs(chadron->pdgId()) == 4122) {
                for(size_t idau=0; idau<chadron->numberOfDaughters(); idau++) {
                  const Candidate *proton = chadron->daughter(idau);
                  if(lepton->pdgId()*proton->pdgId() == 11*2212 ||
                     lepton->pdgId()*proton->pdgId() == 13*2212) {
                    for(size_t ksdau=0; ksdau<chadron->numberOfDaughters(); ksdau++) {
                      if(idau == ksdau) continue;
                      if(abs(chadron->daughter(ksdau)->pdgId()) == 313) {
                        for(size_t jdau=0; jdau<handle1->size(); jdau++) {
                          const Candidate *ks = (*handle1)[jdau].mother(0);
                          if(ks!=nullptr && isAncestor(chadron->daughter(ksdau),ks)) {
                            const Candidate *k = &(*handle1)[jdau];
                            if(lepton->pdgId()*k->pdgId() == -11*321 ||
                               lepton->pdgId()*k->pdgId() == -13*321) {
                              for(size_t kdau=0; kdau<handle1->size(); kdau++) {
                                if(jdau == kdau) continue;
                                const Candidate *ks = (*handle1)[kdau].mother(0);
                                if(ks!=nullptr && isAncestor(chadron->daughter(ksdau),ks)) {
                                  const Candidate *pi = &(*handle1)[kdau];
                                  if(k->pdgId()*pi->pdgId() != -321*211) continue;
                                  if(k->pt() < 1. || pi->pt() < 1. || proton->pt() < 1.0) continue;
                                  if(fabs(k->eta()) > 2.5 || fabs(pi->eta()) > 2.5 || fabs(proton->eta()) > 2.5) continue;
                                  foundGenLambdab0decay_ = true;
                                  lambda_hadron.push_back(chadron);
                                  lepton_lambda.push_back(lepton);
                                  kaon_lambda.push_back(k);
                                  pion_lambda.push_back(pi);
                                  proton_lambda.push_back(proton);
                                  p4_gen_LambdacPlus_lepton_.push_back(this->getLV(lepton->p4()));
                                } // k*0 is the ancestor
                              } // pion candidate
                            } // kaon and lepton have same charge
                          } // k*0 is the ancestor
                        } // kaon candidate
                      } // k*0 candidate
                    } // lambda_c+ daughters
                  } // lepton and proton have opposite charge
                } // proton candidate
              } // lambda_c+ candidate
            } // charm hadrons
          } // b daughters
        } // leptons
      } // b daughters
    } // b0 or b+ or bs+ or lambda_b+ hadron
  } // gen particles
  nLepGen_ = ij;
  //cout<<"  n-lep: "<<ij<<endl;
  //cout<<"Number of top quarks: "<<top_quarks.size()<<endl;
  //cout<<"Number of bottom hadrons: "<<bottom_hadrons.size()<<"  "<<endl;
}

void TopAnalyzer::RecoAnalysis(const edm::Event& iEvent,
                               const edm::EventSetup& iSetup,
                               edm::Handle<edm::View<pat::PackedCandidate> > &handle,
                               edm::Handle<edm::View<pat::Jet> > &jhandle) {
  //cout<<"Jet collection size = "<<jhandle->size()<<endl;
  if(!handle.isValid()) return;
  int ij = 0;
  std::vector<const pat::Jet*> jets = isTightJet(jhandle);
  for(unsigned int iJet=0; iJet<jets.size(); iJet++) {
    if(jets[iJet]->hadronFlavour() != 5) continue;
    ij += 1;
    p4_jet_.push_back(this->getLV(jets[iJet]->p4()));
    //cout<<"  jet: "<<jets[iJet]->pt()<<"  "<<jets[iJet]->eta()<<endl;
    unsigned int tracks = jets[iJet]->numberOfDaughters();
    int nlepton = 0;
    std::vector<double> pdgId_lep, pt_lep, eta_lep, phi_lep, charge_lep, mass_lep;
    for(unsigned int i = 0; i<tracks; i++) {
      const pat::PackedCandidate &itrack = dynamic_cast<const pat::PackedCandidate &>(*jets[iJet]->daughter(i));
      if(abs(itrack.pdgId()) != 11 && abs(itrack.pdgId()) != 13) continue;
      if(abs(itrack.pdgId()) == 11 && itrack.lostInnerHits() > 1) continue;
      if(abs(itrack.pdgId()) == 13 && !itrack.isTrackerMuon() && !itrack.isGlobalMuon()) continue;
      if(itrack.pt() < 3.0) continue;
      if(fabs(itrack.eta()) > 2.5) continue;
      if(itrack.numberOfPixelHits() == 0) continue;
      if(!itrack.trackHighPurity()) continue;
      nlepton += 1;
      pdgId_lep.push_back(itrack.pdgId());
      pt_lep.push_back(itrack.pt());
      eta_lep.push_back(itrack.eta());
      phi_lep.push_back(itrack.phi());
      charge_lep.push_back(itrack.charge());
      if(abs(itrack.pdgId()) == 11) mass_lep.push_back(me);
      else mass_lep.push_back(mmu);
      //cout<<"    lepton in b-jet: "<<itrack.pdgId()<<"  "<<itrack.pt()<<"  "<<itrack.eta()<<endl;
    }
    lepton_perJet_.push_back(nlepton);
    pdgId_lepton_.push_back(pdgId_lep);
    pt_lepton_.push_back(pt_lep);
    eta_lepton_.push_back(eta_lep);
    phi_lepton_.push_back(phi_lep);
    charge_lepton_.push_back(charge_lep);
    mass_lepton_.push_back(mass_lep);
    //cout<<"    n-lep per jet: "<<nlepton<<endl;
    if(nlepton == 0) continue;
    for(unsigned int i = 0; i<tracks; i++) {
      const pat::PackedCandidate &itrack = dynamic_cast<const pat::PackedCandidate &>(*jets[iJet]->daughter(i));
      if(abs(itrack.pdgId()) != 211) continue;
      if(itrack.pt() < 1.) continue;
      if(fabs(itrack.eta()) > 2.5) continue;
      if(itrack.numberOfPixelHits() == 0) continue;
      if(!itrack.trackHighPurity()) continue;
      for(unsigned int j = 0; j<tracks; j++) {
        if(i == j) continue;
        const pat::PackedCandidate &jtrack = dynamic_cast<const pat::PackedCandidate &>(*jets[iJet]->daughter(j));
        if(abs(jtrack.pdgId()) != 211) continue;
        if(jtrack.pt() < 1.) continue;
        if(fabs(jtrack.eta()) > 2.5) continue; 
        if(jtrack.numberOfPixelHits() == 0) continue;
        if(!jtrack.trackHighPurity()) continue;
        if(itrack.charge()*jtrack.charge() > 0) continue;
        for(unsigned int igen=0; igen<d0_hadron.size(); igen++) {
          //cout<<"found D0 hadron: "<<d0_hadron.at(igen)->pdgId()<<endl;
          if(!matchingTrack(kaon_d0.at(igen),itrack)) continue;
          if(!matchingTrack(pion_d0.at(igen),jtrack)) continue;
          double ek = sqrt(itrack.px()*itrack.px()+
                           itrack.py()*itrack.py()+
                           itrack.pz()*itrack.pz()+mk*mk);
          double epi = sqrt(jtrack.px()*jtrack.px()+
                            jtrack.py()*jtrack.py()+
                            jtrack.pz()*jtrack.pz()+mpi*mpi);
          LV lvk, lvpi, d0;
          lvk.SetPxPyPzE(itrack.px(),itrack.py(),itrack.pz(),ek);
          lvpi.SetPxPyPzE(jtrack.px(),jtrack.py(),jtrack.pz(),epi);
          d0 = lvk+lvpi;
          if(d0.Pt() < 1.) continue;
          if(d0.M() < 1.7 || d0.M() > 2.0) continue;
          p4_D0_.push_back(d0);
          p4_D0_kaon_.push_back(lvk);
          p4_D0_pion_.push_back(lvpi);
          charge_D0_kaon_.push_back(itrack.charge());
          charge_D0_pion_.push_back(jtrack.charge());
          para.clear(); covmat.clear();
          trkpara_D0_kaon_ = unpackTrackParameters(itrack,"kaon_d0");
          trkpara_D0_pion_ = unpackTrackParameters(jtrack,"pion_d0");
          trkcovm_D0_kaon_ = unpackTrackCovariance(itrack,"kaon_d0");
          trkcovm_D0_pion_ = unpackTrackCovariance(jtrack,"pion_d0");
          double pv_lxy = (_pv->x()*d0.Px()+
                           _pv->y()*d0.Py())/d0.Pt();
          pv_lxy_D0_.push_back(pv_lxy);
        }
      }
    }
    for(unsigned int i = 0; i<tracks; i++) {
      const pat::PackedCandidate &itrack = dynamic_cast<const pat::PackedCandidate &>(*jets[iJet]->daughter(i));
      if(abs(itrack.pdgId()) != 211) continue;
      if(itrack.pt() < 1.) continue;
      if(fabs(itrack.eta()) > 2.5) continue;
      if(itrack.numberOfPixelHits() == 0) continue;
      if(!itrack.trackHighPurity()) continue;
      for(unsigned int j = 0; j<tracks; j++) {
        const pat::PackedCandidate &jtrack = dynamic_cast<const pat::PackedCandidate &>(*jets[iJet]->daughter(j));
        if(abs(jtrack.pdgId()) != 211) continue;
        if(jtrack.pt() < 1.) continue;
        if(fabs(jtrack.eta()) > 2.5) continue;
        if(jtrack.numberOfPixelHits() == 0) continue;
        if(!jtrack.trackHighPurity()) continue;
        for(unsigned int k = 0; k<tracks; k++) {
          if(i == j) continue;
          if(i == k) continue;
          if(j == k) continue;
          const pat::PackedCandidate &ktrack = dynamic_cast<const pat::PackedCandidate &>(*jets[iJet]->daughter(k));
          if(abs(ktrack.pdgId()) != 211) continue;
          if(ktrack.pt() < 0.5) continue;
          if(fabs(ktrack.eta()) > 2.5) continue;
          if(ktrack.numberOfPixelHits() == 0) continue;
          if(!ktrack.trackHighPurity()) continue;
          if(itrack.charge()*jtrack.charge() > 0) continue;
          if(itrack.charge()*ktrack.charge() > 0) continue;
          for(unsigned int igen=0; igen<dt_hadron.size(); igen++) {
            //cout<<"found D* hadron: "<<dt_hadron.at(igen)->pdgId()<<endl;
            if(!matchingTrack(kaon_dt.at(igen),itrack)) continue;
            if(!matchingTrack(pion_dt.at(igen),jtrack)) continue;
            if(!matchingTrack(pionsoft_dt.at(igen),ktrack)) continue;
            double ek = sqrt(itrack.px()*itrack.px()+
                             itrack.py()*itrack.py()+
                             itrack.pz()*itrack.pz()+mk*mk);
            double epi = sqrt(jtrack.px()*jtrack.px()+
                              jtrack.py()*jtrack.py()+
                              jtrack.pz()*jtrack.pz()+mpi*mpi);
            double epis = sqrt(ktrack.px()*ktrack.px()+
                               ktrack.py()*ktrack.py()+
                               ktrack.pz()*ktrack.pz()+mpi*mpi);
            LV lvk, lvpi, lvpis, d0, dt;
            lvk.SetPxPyPzE(itrack.px(),itrack.py(),itrack.pz(),ek);
            lvpi.SetPxPyPzE(jtrack.px(),jtrack.py(),jtrack.pz(),epi);
            lvpis.SetPxPyPzE(ktrack.px(),ktrack.py(),ktrack.pz(),epis);
            d0 = lvk+lvpi;
            dt = lvk+lvpi+lvpis;
            if(d0.Pt() < 1.) continue;
            if(d0.M() < 1.7 || d0.M() > 2.0) continue;
            mdiff_Dstar_D0_.push_back(dt.M()-d0.M());
            p4_Dstar_.push_back(dt);
            p4_Dstar_D0_.push_back(d0);
            p4_Dstar_D0_kaon_.push_back(lvk);
            p4_Dstar_D0_pion_.push_back(lvpi);
            p4_Dstar_pionsoft_.push_back(lvpis);
            charge_Dstar_D0_kaon_.push_back(itrack.charge());
            charge_Dstar_D0_pion_.push_back(jtrack.charge());
            charge_Dstar_pionsoft_.push_back(ktrack.charge());
            para.clear(); covmat.clear();
            trkpara_Dstar_D0_kaon_ = unpackTrackParameters(itrack,"kaon_dt");
            trkpara_Dstar_D0_pion_ = unpackTrackParameters(jtrack,"pion_dt");
            trkpara_Dstar_pionsoft_ = unpackTrackParameters(ktrack,"pionsoft_dt");
            trkcovm_Dstar_D0_kaon_ = unpackTrackCovariance(itrack,"kaon_dt");
            trkcovm_Dstar_D0_pion_ = unpackTrackCovariance(jtrack,"pion_dt");
            trkcovm_Dstar_pionsoft_ = unpackTrackCovariance(ktrack,"pionsoft_dt");
            double pv_lxy = (_pv->x()*d0.Px()+
                             _pv->y()*d0.Py())/d0.Pt();
            pv_lxy_Dstar_D0_.push_back(pv_lxy);
          }
        }
      }
    }
    for(unsigned int i = 0; i<tracks; i++) {
      const pat::PackedCandidate &itrack = dynamic_cast<const pat::PackedCandidate &>(*jets[iJet]->daughter(i));
      if(abs(itrack.pdgId()) != 211) continue;
      if(itrack.pt() < 1.) continue;
      if(fabs(itrack.eta()) > 2.5) continue;
      if(itrack.numberOfPixelHits() == 0) continue;
      if(!itrack.trackHighPurity()) continue;
      for(unsigned int j = 0; j<tracks; j++) {
        const pat::PackedCandidate &jtrack = dynamic_cast<const pat::PackedCandidate &>(*jets[iJet]->daughter(j));
        if(abs(jtrack.pdgId()) != 211) continue;
        if(jtrack.pt() < 1.) continue;
        if(fabs(jtrack.eta()) > 2.5) continue;
        if(jtrack.numberOfPixelHits() == 0) continue;
        if(!jtrack.trackHighPurity()) continue;
        for(unsigned int k = 0; k<tracks; k++) {
          if(i == j) continue;
          if(i == k) continue;
          if(j == k) continue;
          const pat::PackedCandidate &ktrack = dynamic_cast<const pat::PackedCandidate &>(*jets[iJet]->daughter(k));
          if(abs(ktrack.pdgId()) != 211) continue;
          if(ktrack.pt() < 1.) continue;
          if(fabs(ktrack.eta()) > 2.5) continue;
          if(ktrack.numberOfPixelHits() == 0) continue;
          if(!ktrack.trackHighPurity()) continue;
          if(itrack.charge()*jtrack.charge() > 0) continue;
          if(itrack.charge()*ktrack.charge() > 0) continue;
          for(unsigned int igen=0; igen<dp_hadron.size(); igen++) {
            //cout<<"  found D+ hadron: "<<dp_hadron.at(igen)->pdgId()<<endl;
            if(!matchingTrack(kaon_dp.at(igen),itrack)) continue;
            if(!matchingTrack(pion1_dp.at(igen),jtrack)) continue;
            if(!matchingTrack(pion2_dp.at(igen),ktrack)) continue;
            double ek = sqrt(itrack.px()*itrack.px()+
                             itrack.py()*itrack.py()+
                             itrack.pz()*itrack.pz()+mk*mk);
            double epi1 = sqrt(jtrack.px()*jtrack.px()+
                               jtrack.py()*jtrack.py()+
                               jtrack.pz()*jtrack.pz()+mpi*mpi);
            double epi2 = sqrt(ktrack.px()*ktrack.px()+
                               ktrack.py()*ktrack.py()+
                               ktrack.pz()*ktrack.pz()+mpi*mpi);
            LV lvk, lvpi1, lvpi2, dp;
            lvk.SetPxPyPzE(itrack.px(),itrack.py(),itrack.pz(),ek);
            lvpi1.SetPxPyPzE(jtrack.px(),jtrack.py(),jtrack.pz(),epi1);
            lvpi2.SetPxPyPzE(ktrack.px(),ktrack.py(),ktrack.pz(),epi2);
            dp = lvk+lvpi1+lvpi2;
            if(dp.Pt() < 1.) continue;
            if(dp.M() < 1.7 || dp.M() > 2.0) continue;
            p4_DPlus_.push_back(dp);
            p4_DPlus_kaon_.push_back(lvk);
            p4_DPlus_pion1_.push_back(lvpi1);
            p4_DPlus_pion2_.push_back(lvpi2);
            charge_DPlus_kaon_.push_back(itrack.charge());
            charge_DPlus_pion1_.push_back(jtrack.charge());
            charge_DPlus_pion2_.push_back(ktrack.charge());
            para.clear(); covmat.clear();
            trkpara_DPlus_kaon_ = unpackTrackParameters(itrack,"kaon_dp");
            trkpara_DPlus_pion1_ = unpackTrackParameters(jtrack,"pion1_dp");
            trkpara_DPlus_pion2_ = unpackTrackParameters(ktrack,"pion2_dp");
            trkcovm_DPlus_kaon_ = unpackTrackCovariance(itrack,"kaon_dp");
            trkcovm_DPlus_pion1_ = unpackTrackCovariance(jtrack,"pion1_dp");
            trkcovm_DPlus_pion2_ = unpackTrackCovariance(ktrack,"pion2_dp");
            double pv_lxy = (_pv->x()*dp.Px()+
                             _pv->y()*dp.Py())/dp.Pt();
            pv_lxy_DPlus_.push_back(pv_lxy);
          }
        }
      }
    }
    for(unsigned int i = 0; i<tracks; i++) {
      const pat::PackedCandidate &itrack = dynamic_cast<const pat::PackedCandidate &>(*jets[iJet]->daughter(i));
      if(abs(itrack.pdgId()) != 211) continue;
      if(itrack.pt() < 1.) continue;
      if(fabs(itrack.eta()) > 2.5) continue;
      if(itrack.numberOfPixelHits() == 0) continue;
      if(!itrack.trackHighPurity()) continue;
      for(unsigned int j = 0; j<tracks; j++) {
        const pat::PackedCandidate &jtrack = dynamic_cast<const pat::PackedCandidate &>(*jets[iJet]->daughter(j));
        if(abs(jtrack.pdgId()) != 211) continue;
        if(jtrack.pt() < 1.) continue;
        if(fabs(jtrack.eta()) > 2.5) continue;
        if(jtrack.numberOfPixelHits() == 0) continue;
        if(!jtrack.trackHighPurity()) continue;
        for(unsigned int k = 0; k<tracks; k++) {
          if(i == j) continue;
          if(i == k) continue;
          if(j == k) continue;
          const pat::PackedCandidate &ktrack = dynamic_cast<const pat::PackedCandidate &>(*jets[iJet]->daughter(k));
          if(abs(ktrack.pdgId()) != 211) continue;
          if(ktrack.pt() < 1.) continue;
          if(fabs(ktrack.eta()) > 2.5) continue;
          if(ktrack.numberOfPixelHits() == 0) continue;
          if(!ktrack.trackHighPurity()) continue;
          if(itrack.charge()*jtrack.charge() > 0) continue;
          for(unsigned int igen=0; igen<ds_hadron.size(); igen++) {
            //cout<<"found Ds hadron: "<<ds_hadron.at(igen)->pdgId()<<endl;
            if(!matchingTrack(kaon1_ds.at(igen),itrack)) continue;
            if(!matchingTrack(kaon2_ds.at(igen),jtrack)) continue;
            if(!matchingTrack(pion_ds.at(igen),ktrack)) continue;
            double ek1 = sqrt(itrack.px()*itrack.px()+
                              itrack.py()*itrack.py()+
                              itrack.pz()*itrack.pz()+mk*mk);
            double ek2 = sqrt(jtrack.px()*jtrack.px()+
                              jtrack.py()*jtrack.py()+
                              jtrack.pz()*jtrack.pz()+mk*mk);
            double epi = sqrt(ktrack.px()*ktrack.px()+
                              ktrack.py()*ktrack.py()+
                              ktrack.pz()*ktrack.pz()+mpi*mpi);
            LV lvk1, lvk2, lvpi, ds;
            lvk1.SetPxPyPzE(itrack.px(),itrack.py(),itrack.pz(),ek1);
            lvk2.SetPxPyPzE(jtrack.px(),jtrack.py(),jtrack.pz(),ek2);
            lvpi.SetPxPyPzE(ktrack.px(),ktrack.py(),ktrack.pz(),epi);
            ds = lvk1+lvk2+lvpi;
            if(ds.Pt() < 1.) continue;
            if(ds.M() < 1.8 || ds.M() > 2.1) continue;
            p4_DsPlus_.push_back(ds);
            p4_DsPlus_kaon1_.push_back(lvk1);
            p4_DsPlus_kaon2_.push_back(lvk2);
            p4_DsPlus_pion_.push_back(lvpi);
            charge_DsPlus_kaon1_.push_back(itrack.charge());
            charge_DsPlus_kaon2_.push_back(jtrack.charge());
            charge_DsPlus_pion_.push_back(ktrack.charge());
            para.clear(); covmat.clear();
            trkpara_DsPlus_kaon1_ = unpackTrackParameters(itrack,"kaon1_ds");
            trkpara_DsPlus_kaon2_ = unpackTrackParameters(jtrack,"kaon2_ds");
            trkpara_DsPlus_pion_ = unpackTrackParameters(ktrack,"pion_ds");
            trkcovm_DsPlus_kaon1_ = unpackTrackCovariance(itrack,"kaon1_ds");
            trkcovm_DsPlus_kaon2_ = unpackTrackCovariance(jtrack,"kaon2_ds");
            trkcovm_DsPlus_pion_ = unpackTrackCovariance(ktrack,"pion_ds");
            double pv_lxy = (_pv->x()*ds.Px()+
                           _pv->y()*ds.Py())/ds.Pt();
            pv_lxy_DsPlus_.push_back(pv_lxy);
          }
        }
      }
    }
    for(unsigned int i = 0; i<tracks; i++) {
      const pat::PackedCandidate &itrack = dynamic_cast<const pat::PackedCandidate &>(*jets[iJet]->daughter(i));
      if(abs(itrack.pdgId()) != 211) continue;
      if(itrack.pt() < 1.) continue;
      if(fabs(itrack.eta()) > 2.5) continue;
      if(itrack.numberOfPixelHits() == 0) continue;
      if(!itrack.trackHighPurity()) continue;
      for(unsigned int j = 0; j<tracks; j++) {
        const pat::PackedCandidate &jtrack = dynamic_cast<const pat::PackedCandidate &>(*jets[iJet]->daughter(j));
        if(abs(jtrack.pdgId()) != 211) continue; 
        if(jtrack.pt() < 1.) continue;
        if(fabs(jtrack.eta()) > 2.5) continue;
        if(jtrack.numberOfPixelHits() == 0) continue;
        if(!jtrack.trackHighPurity()) continue;
        for(unsigned int k = 0; k<tracks; k++) {
          if(i == j) continue;
          if(i == k) continue;
          if(j == k) continue;
          const pat::PackedCandidate &ktrack = dynamic_cast<const pat::PackedCandidate &>(*jets[iJet]->daughter(k));
          if(abs(ktrack.pdgId()) != 211) continue;
          if(ktrack.pt() < 1.) continue;
          if(fabs(ktrack.eta()) > 2.5) continue;
          if(ktrack.numberOfPixelHits() == 0) continue;
          if(!ktrack.trackHighPurity()) continue;
          if(itrack.charge()*jtrack.charge() > 0) continue;
          if(itrack.charge()*ktrack.charge() > 0) continue;
          for(unsigned int igen=0; igen<lambda_hadron.size(); igen++) {
            //cout<<"found Lambda b0 hadron: "<<lambda_hadron.at(igen)->pdgId()<<endl;
            if(!matchingTrack(kaon_lambda.at(igen),itrack)) continue;
            if(!matchingTrack(pion_lambda.at(igen),jtrack)) continue;
            if(!matchingTrack(proton_lambda.at(igen),ktrack)) continue;
            double ek = sqrt(itrack.px()*itrack.px()+
                             itrack.py()*itrack.py()+
                             itrack.pz()*itrack.pz()+mk*mk);
            double epi = sqrt(jtrack.px()*jtrack.px()+
                              jtrack.py()*jtrack.py()+
                              jtrack.pz()*jtrack.pz()+mpi*mpi);
            double ep = sqrt(ktrack.px()*ktrack.px()+
                             ktrack.py()*ktrack.py()+
                             ktrack.pz()*ktrack.pz()+mp*mp);
            LV lvk, lvpi, lvp, lambdac;
            lvk.SetPxPyPzE(itrack.px(),itrack.py(),itrack.pz(),ek);
            lvpi.SetPxPyPzE(jtrack.px(),jtrack.py(),jtrack.pz(),epi);
            lvp.SetPxPyPzE(ktrack.px(),ktrack.py(),ktrack.pz(),ep);
            lambdac = lvk+lvpi+lvp;
            if(lambdac.Pt() < 1.) continue;
            if(lambdac.M() < 2.1 || lambdac.M() > 2.4) continue;
            p4_LambdacPlus_.push_back(lambdac);
            p4_LambdacPlus_kaon_.push_back(lvk);
            p4_LambdacPlus_pion_.push_back(lvpi);
            p4_LambdacPlus_proton_.push_back(lvp);
            charge_LambdacPlus_kaon_.push_back(itrack.charge());
            charge_LambdacPlus_pion_.push_back(jtrack.charge());
            charge_LambdacPlus_proton_.push_back(ktrack.charge());
            para.clear(); covmat.clear();
            trkpara_LambdacPlus_kaon_ = unpackTrackParameters(itrack,"kaon_lambda");
            trkpara_LambdacPlus_pion_ = unpackTrackParameters(jtrack,"pion_lambda");
            trkpara_LambdacPlus_proton_ = unpackTrackParameters(ktrack,"proton_lambda");
            trkcovm_LambdacPlus_kaon_ = unpackTrackCovariance(itrack,"kaon_lambda");
            trkcovm_LambdacPlus_pion_ = unpackTrackCovariance(jtrack,"pion_lambda");
            trkcovm_LambdacPlus_proton_ = unpackTrackCovariance(ktrack,"proton_lambda");
            double pv_lxy = (_pv->x()*lambdac.Px()+
                             _pv->y()*lambdac.Py())/lambdac.Pt();
            pv_lxy_LambdacPlus_.push_back(pv_lxy);
          }
        }
      }
    }
  }
  nJet_ = ij;
  //cout<<"  n-jet: "<<ij<<endl;
}

void TopAnalyzer::DataAnalysis(const edm::Event& iEvent,
                               const edm::EventSetup& iSetup,
                               edm::Handle<edm::View<pat::Jet> >& jhandle) {
  //cout<<"Jet collection size = "<<jhandle->size()<<endl;
  if(!jhandle.isValid()) return;
  //cout<<"muons.size(): "<<muons.size()<<"  electrons.size(): "<<electrons.size()<<endl;
  //if(muons.size() < 1 && electrons.size() < 1) return;
  int ij = 0;
  std::vector<const pat::Jet*> jets = isTightJet(jhandle);
  for(unsigned int iJet=0; iJet<jets.size(); iJet++) {
    ij += 1;
    p4_jet_.push_back(this->getLV(jets[iJet]->p4()));
    //cout<<"  jets: "<<jets[iJet]->pt()<<"  "<<jets[iJet]->eta()<<endl;
    unsigned int tracks = jets[iJet]->numberOfDaughters();
    int nlepton = 0;
    std::vector<double> pdgId_lep, pt_lep, eta_lep, phi_lep, charge_lep, mass_lep;
    for(unsigned int i = 0; i<tracks; i++) {
      const pat::PackedCandidate &itrack = dynamic_cast<const pat::PackedCandidate &>(*jets[iJet]->daughter(i));
      if(abs(itrack.pdgId()) != 11 && abs(itrack.pdgId()) != 13) continue;
      if(abs(itrack.pdgId()) == 11 && itrack.lostInnerHits() > 1) continue;
      if(abs(itrack.pdgId()) == 13 && !itrack.isTrackerMuon() && !itrack.isGlobalMuon()) continue;
      if(itrack.pt() < 3.) continue;
      if(fabs(itrack.eta()) > 2.5) continue;
      if(itrack.numberOfPixelHits() == 0) continue;
      if(!itrack.trackHighPurity()) continue;
      nlepton += 1;
      pdgId_lep.push_back(itrack.pdgId());
      pt_lep.push_back(itrack.pt());
      eta_lep.push_back(itrack.eta());
      phi_lep.push_back(itrack.phi());
      charge_lep.push_back(itrack.charge());
      if(abs(itrack.pdgId()) == 11) mass_lep.push_back(me);
      else mass_lep.push_back(mmu);
      //cout<<"    lepton in b-jet: "<<itrack.pdgId()<<"  "<<itrack.pt()<<"  "<<itrack.eta()<<endl;
    }
    lepton_perJet_.push_back(nlepton);
    pdgId_lepton_.push_back(pdgId_lep);
    pt_lepton_.push_back(pt_lep);
    eta_lepton_.push_back(eta_lep);
    phi_lepton_.push_back(phi_lep);
    charge_lepton_.push_back(charge_lep);
    mass_lepton_.push_back(mass_lep);
    //cout<<"    n-lep per jet: "<<nlepton<<endl;
    if(nlepton == 0) continue;
    for(unsigned int i = 0; i<tracks; i++) {
      const pat::PackedCandidate &itrack = dynamic_cast<const pat::PackedCandidate &>(*jets[iJet]->daughter(i));
      if(abs(itrack.pdgId()) != 211) continue;
      if(itrack.pt() < 1.) continue;
      if(fabs(itrack.eta()) > 2.5) continue;
      if(itrack.numberOfPixelHits() == 0) continue;
      if(!itrack.trackHighPurity()) continue;
      for(unsigned int j = 0; j<tracks; j++) {
        if(i == j) continue;
        const pat::PackedCandidate &jtrack = dynamic_cast<const pat::PackedCandidate &>(*jets[iJet]->daughter(j));
        if(abs(jtrack.pdgId()) != 211) continue;
        if(jtrack.pt() < 1.) continue;
        if(fabs(jtrack.eta()) > 2.5) continue;
        if(jtrack.numberOfPixelHits() == 0) continue;
        if(!jtrack.trackHighPurity()) continue; 
        if(itrack.charge()*jtrack.charge() > 0) continue;
        double ek = sqrt(itrack.px()*itrack.px()+
                         itrack.py()*itrack.py()+
                         itrack.pz()*itrack.pz()+mk*mk);
        double epi = sqrt(jtrack.px()*jtrack.px()+
                          jtrack.py()*jtrack.py()+
                          jtrack.pz()*jtrack.pz()+mpi*mpi);
        LV lvk, lvpi, d0;
        lvk.SetPxPyPzE(itrack.px(),itrack.py(),itrack.pz(),ek);
        lvpi.SetPxPyPzE(jtrack.px(),jtrack.py(),jtrack.pz(),epi);
        d0 = lvk+lvpi;
        if(d0.Pt() < 1.) continue;
        if(d0.M() < 1.7 || d0.M() > 2.0) continue;
        //cout<<"    D0: "<<d0<<endl;
        p4_D0_.push_back(d0);
        p4_D0_kaon_.push_back(lvk);
        p4_D0_pion_.push_back(lvpi);
        charge_D0_kaon_.push_back(itrack.charge());
        charge_D0_pion_.push_back(jtrack.charge());
        para.clear(); covmat.clear();
        trkpara_D0_kaon_ = unpackTrackParameters(itrack,"kaon_d0");
        trkpara_D0_pion_ = unpackTrackParameters(jtrack,"pion_d0");
        trkcovm_D0_kaon_ = unpackTrackCovariance(itrack,"kaon_d0");
        trkcovm_D0_pion_ = unpackTrackCovariance(jtrack,"pion_d0");
        double pv_lxy = (_pv->x()*d0.Px()+
                         _pv->y()*d0.Py())/d0.Pt();
        pv_lxy_D0_.push_back(pv_lxy);
      }
    }
    for(unsigned int i = 0; i<tracks; i++) {
      const pat::PackedCandidate &itrack = dynamic_cast<const pat::PackedCandidate &>(*jets[iJet]->daughter(i));
      if(abs(itrack.pdgId()) != 211) continue;
      if(itrack.pt() < 1.) continue;
      if(fabs(itrack.eta()) > 2.5) continue;
      if(itrack.numberOfPixelHits() == 0) continue;
      if(!itrack.trackHighPurity()) continue;
      for(unsigned int j = 0; j<tracks; j++) {
        const pat::PackedCandidate &jtrack = dynamic_cast<const pat::PackedCandidate &>(*jets[iJet]->daughter(j));
        if(abs(jtrack.pdgId()) != 211) continue;
        if(jtrack.pt() < 1.) continue;
        if(fabs(jtrack.eta()) > 2.5) continue;
        if(jtrack.numberOfPixelHits() == 0) continue;
        if(!jtrack.trackHighPurity()) continue;
        for(unsigned int k = 0; k<tracks; k++) {
          if(i == j) continue;
          if(i == k) continue;
          if(j == k) continue;
          const pat::PackedCandidate &ktrack = dynamic_cast<const pat::PackedCandidate &>(*jets[iJet]->daughter(k));
          if(abs(ktrack.pdgId()) != 211) continue;
          if(ktrack.pt() < 0.5) continue;
          if(fabs(ktrack.eta()) > 2.5) continue;
          if(ktrack.numberOfPixelHits() == 0) continue;
          if(!ktrack.trackHighPurity()) continue;
          //if(itrack.charge()*jtrack.charge() > 0) continue;
          if(itrack.charge()*ktrack.charge() > 0) continue;
          double ek = sqrt(itrack.px()*itrack.px()+
                           itrack.py()*itrack.py()+
                           itrack.pz()*itrack.pz()+mk*mk);
          double epi = sqrt(jtrack.px()*jtrack.px()+
                            jtrack.py()*jtrack.py()+
                            jtrack.pz()*jtrack.pz()+mpi*mpi);
          double epis = sqrt(ktrack.px()*ktrack.px()+
                             ktrack.py()*ktrack.py()+
                             ktrack.pz()*ktrack.pz()+mpi*mpi);
          LV lvk, lvpi, lvpis, d0, dt;
          lvk.SetPxPyPzE(itrack.px(),itrack.py(),itrack.pz(),ek);
          lvpi.SetPxPyPzE(jtrack.px(),jtrack.py(),jtrack.pz(),epi);
          lvpis.SetPxPyPzE(ktrack.px(),ktrack.py(),ktrack.pz(),epis);
          d0 = lvk+lvpi;
          dt = lvk+lvpi+lvpis;
          if(d0.Pt() < 1.) continue;
          if(d0.M() < 1.7 || d0.M() > 2.0) continue;
          //cout<<"    D*-D0: "<<d0<<endl;
          mdiff_Dstar_D0_.push_back(dt.M()-d0.M());
          p4_Dstar_.push_back(dt);
          p4_Dstar_D0_.push_back(d0);
          p4_Dstar_D0_kaon_.push_back(lvk);
          p4_Dstar_D0_pion_.push_back(lvpi);
          p4_Dstar_pionsoft_.push_back(lvpis);
          charge_Dstar_D0_kaon_.push_back(itrack.charge());
          charge_Dstar_D0_pion_.push_back(jtrack.charge());
          charge_Dstar_pionsoft_.push_back(ktrack.charge());
          para.clear(); covmat.clear();
          trkpara_Dstar_D0_kaon_ = unpackTrackParameters(itrack,"kaon_dt");
          trkpara_Dstar_D0_pion_ = unpackTrackParameters(jtrack,"pion_dt");
          trkpara_Dstar_pionsoft_ = unpackTrackParameters(ktrack,"pionsoft_dt");
          trkcovm_Dstar_D0_kaon_ = unpackTrackCovariance(itrack,"kaon_dt");
          trkcovm_Dstar_D0_pion_ = unpackTrackCovariance(jtrack,"pion_dt");
          trkcovm_Dstar_pionsoft_ = unpackTrackCovariance(ktrack,"pionsoft_dt");
          double pv_lxy = (_pv->x()*d0.Px()+
                           _pv->y()*d0.Py())/d0.Pt();
          pv_lxy_Dstar_D0_.push_back(pv_lxy);
        }
      }
    }
    for(unsigned int i = 0; i<tracks; i++) {
      const pat::PackedCandidate &itrack = dynamic_cast<const pat::PackedCandidate &>(*jets[iJet]->daughter(i));
      if(abs(itrack.pdgId()) != 211) continue;
      if(itrack.pt() < 1.) continue;
      if(fabs(itrack.eta()) > 2.5) continue;
      if(itrack.numberOfPixelHits() == 0) continue;
      if(!itrack.trackHighPurity()) continue;
      for(unsigned int j = 0; j<tracks; j++) {
        const pat::PackedCandidate &jtrack = dynamic_cast<const pat::PackedCandidate &>(*jets[iJet]->daughter(j));
        if(abs(jtrack.pdgId()) != 211) continue;
        if(jtrack.pt() < 1.) continue;
        if(fabs(jtrack.eta()) > 2.5) continue;
        if(jtrack.numberOfPixelHits() == 0) continue;
        if(!jtrack.trackHighPurity()) continue;
        for(unsigned int k = 0; k<j; k++) {
          if(i == j) continue;
          if(i == k) continue;
          if(j == k) continue;
          const pat::PackedCandidate &ktrack = dynamic_cast<const pat::PackedCandidate &>(*jets[iJet]->daughter(k));
          if(abs(ktrack.pdgId()) != 211) continue;
          if(ktrack.pt() < 1.) continue;
          if(fabs(ktrack.eta()) > 2.5) continue;
          if(ktrack.numberOfPixelHits() == 0) continue;
          if(!ktrack.trackHighPurity()) continue;
          if(itrack.charge()*jtrack.charge() > 0) continue;
          if(itrack.charge()*ktrack.charge() > 0) continue;
          double ek = sqrt(itrack.px()*itrack.px()+
                           itrack.py()*itrack.py()+
                           itrack.pz()*itrack.pz()+mk*mk);
          double epi1 = sqrt(jtrack.px()*jtrack.px()+
                             jtrack.py()*jtrack.py()+
                             jtrack.pz()*jtrack.pz()+mpi*mpi);
          double epi2 = sqrt(ktrack.px()*ktrack.px()+
                             ktrack.py()*ktrack.py()+
                             ktrack.pz()*ktrack.pz()+mpi*mpi);
          LV lvk, lvpi1, lvpi2, dp;
          lvk.SetPxPyPzE(itrack.px(),itrack.py(),itrack.pz(),ek);
          lvpi1.SetPxPyPzE(jtrack.px(),jtrack.py(),jtrack.pz(),epi1);
          lvpi2.SetPxPyPzE(ktrack.px(),ktrack.py(),ktrack.pz(),epi2);
          dp = lvk+lvpi1+lvpi2;
          if(dp.Pt() < 1.) continue;
          if(dp.M() < 1.7 || dp.M() > 2.0) continue;
          //cout<<"    D+: "<<dp<<endl;
          p4_DPlus_.push_back(dp);
          p4_DPlus_kaon_.push_back(lvk);
          p4_DPlus_pion1_.push_back(lvpi1);
          p4_DPlus_pion2_.push_back(lvpi2);
          charge_DPlus_kaon_.push_back(itrack.charge());
          charge_DPlus_pion1_.push_back(jtrack.charge());
          charge_DPlus_pion2_.push_back(ktrack.charge());
          para.clear(); covmat.clear();
          trkpara_DPlus_kaon_ = unpackTrackParameters(itrack,"kaon_dp");
          trkpara_DPlus_pion1_ = unpackTrackParameters(jtrack,"pion1_dp");
          trkpara_DPlus_pion2_ = unpackTrackParameters(ktrack,"pion2_dp");
          trkcovm_DPlus_kaon_ = unpackTrackCovariance(itrack,"kaon_dp");
          trkcovm_DPlus_pion1_ = unpackTrackCovariance(jtrack,"pion1_dp");
          trkcovm_DPlus_pion2_ = unpackTrackCovariance(ktrack,"pion2_dp");
          double pv_lxy = (_pv->x()*dp.Px()+
                           _pv->y()*dp.Py())/dp.Pt();
          pv_lxy_DPlus_.push_back(pv_lxy);
        }
      }
    }
    for(unsigned int i = 0; i<tracks; i++) {
      const pat::PackedCandidate &itrack = dynamic_cast<const pat::PackedCandidate &>(*jets[iJet]->daughter(i));
      if(abs(itrack.pdgId()) != 211) continue;
      if(itrack.pt() < 1.) continue;
      if(fabs(itrack.eta()) > 2.5) continue;
      if(itrack.numberOfPixelHits() == 0) continue;
      if(!itrack.trackHighPurity()) continue;
      for(unsigned int j = 0; j<i; j++) {
        const pat::PackedCandidate &jtrack = dynamic_cast<const pat::PackedCandidate &>(*jets[iJet]->daughter(j));
        if(abs(jtrack.pdgId()) != 211) continue;
        if(jtrack.pt() < 1.) continue;
        if(fabs(jtrack.eta()) > 2.5) continue;
        if(jtrack.numberOfPixelHits() == 0) continue;
        if(!jtrack.trackHighPurity()) continue;
        for(unsigned int k = 0; k<tracks; k++) {
          if(i == j) continue;
          if(i == k) continue;
          if(j == k) continue;
          const pat::PackedCandidate &ktrack = dynamic_cast<const pat::PackedCandidate &>(*jets[iJet]->daughter(k));
          if(abs(ktrack.pdgId()) != 211) continue;
          if(ktrack.pt() < 1.) continue;
          if(fabs(ktrack.eta()) > 2.5) continue;
          if(ktrack.numberOfPixelHits() == 0) continue;
          if(!ktrack.trackHighPurity()) continue;
          if(itrack.charge()*jtrack.charge() > 0) continue;
          double ek1 = sqrt(itrack.px()*itrack.px()+
                            itrack.py()*itrack.py()+
                            itrack.pz()*itrack.pz()+mk*mk);
          double ek2 = sqrt(jtrack.px()*jtrack.px()+
                            jtrack.py()*jtrack.py()+
                            jtrack.pz()*jtrack.pz()+mk*mk);
          double epi = sqrt(ktrack.px()*ktrack.px()+
                            ktrack.py()*ktrack.py()+
                            ktrack.pz()*ktrack.pz()+mpi*mpi);
          LV lvk1, lvk2, lvpi, ds;
          lvk1.SetPxPyPzE(itrack.px(),itrack.py(),itrack.pz(),ek1);
          lvk2.SetPxPyPzE(jtrack.px(),jtrack.py(),jtrack.pz(),ek2);
          lvpi.SetPxPyPzE(ktrack.px(),ktrack.py(),ktrack.pz(),epi);
          ds = lvk1+lvk2+lvpi;
          if(ds.Pt() < 1.) continue;
          if(ds.M() < 1.8 || ds.M() > 2.1) continue;
          //cout<<"    Ds+: "<<ds<<endl;
          p4_DsPlus_.push_back(ds);
          p4_DsPlus_kaon1_.push_back(lvk1);
          p4_DsPlus_kaon2_.push_back(lvk2);
          p4_DsPlus_pion_.push_back(lvpi);
          charge_DsPlus_kaon1_.push_back(itrack.charge());
          charge_DsPlus_kaon2_.push_back(jtrack.charge());
          charge_DsPlus_pion_.push_back(ktrack.charge());
          para.clear(); covmat.clear();
          trkpara_DsPlus_kaon1_ = unpackTrackParameters(itrack,"kaon1_ds");
          trkpara_DsPlus_kaon2_ = unpackTrackParameters(jtrack,"kaon2_ds");
          trkpara_DsPlus_pion_ = unpackTrackParameters(ktrack,"pion_ds");
          trkcovm_DsPlus_kaon1_ = unpackTrackCovariance(itrack,"kaon1_ds");
          trkcovm_DsPlus_kaon2_ = unpackTrackCovariance(jtrack,"kaon2_ds");
          trkcovm_DsPlus_pion_ = unpackTrackCovariance(ktrack,"pion_ds");
          double pv_lxy = (_pv->x()*ds.Px()+
                           _pv->y()*ds.Py())/ds.Pt();
          pv_lxy_DsPlus_.push_back(pv_lxy);
        }
      }
    }
    for(unsigned int i = 0; i<tracks; i++) {
      const pat::PackedCandidate &itrack = dynamic_cast<const pat::PackedCandidate &>(*jets[iJet]->daughter(i));
      if(abs(itrack.pdgId()) != 211) continue;
      if(itrack.pt() < 1.) continue;
      if(fabs(itrack.eta()) > 2.5) continue;
      if(itrack.numberOfPixelHits() == 0) continue;
      if(!itrack.trackHighPurity()) continue;
      for(unsigned int j = 0; j<tracks; j++) {
        const pat::PackedCandidate &jtrack = dynamic_cast<const pat::PackedCandidate &>(*jets[iJet]->daughter(j));
        if(abs(jtrack.pdgId()) != 211) continue;
        if(jtrack.pt() < 1.) continue;
        if(fabs(jtrack.eta()) > 2.5) continue;
        if(jtrack.numberOfPixelHits() == 0) continue;
        if(!jtrack.trackHighPurity()) continue;
        for(unsigned int k = 0; k<tracks; k++) {
          if(i == j) continue;
          if(i == k) continue;
          if(j == k) continue;
          const pat::PackedCandidate &ktrack = dynamic_cast<const pat::PackedCandidate &>(*jets[iJet]->daughter(k));
          if(abs(ktrack.pdgId()) != 211) continue;
          if(ktrack.pt() < 1.) continue;
          if(fabs(ktrack.eta()) > 2.5) continue;
          if(ktrack.numberOfPixelHits() == 0) continue;
          if(!ktrack.trackHighPurity()) continue;
          if(itrack.charge()*jtrack.charge() > 0) continue;
          if(itrack.charge()*ktrack.charge() > 0) continue;
          double ek = sqrt(itrack.px()*itrack.px()+
                           itrack.py()*itrack.py()+
                           itrack.pz()*itrack.pz()+mk*mk);
          double epi = sqrt(jtrack.px()*jtrack.px()+
                            jtrack.py()*jtrack.py()+
                            jtrack.pz()*jtrack.pz()+mpi*mpi);
          double ep = sqrt(ktrack.px()*ktrack.px()+
                           ktrack.py()*ktrack.py()+
                           ktrack.pz()*ktrack.pz()+mp*mp);
          LV lvk, lvpi, lvp, lambdac;
          lvk.SetPxPyPzE(itrack.px(),itrack.py(),itrack.pz(),ek);
          lvpi.SetPxPyPzE(jtrack.px(),jtrack.py(),jtrack.pz(),epi);
          lvp.SetPxPyPzE(ktrack.px(),ktrack.py(),ktrack.pz(),ep);
          lambdac = lvk+lvpi+lvp;
          if(lambdac.Pt() < 1.) continue;
          if(lambdac.M() < 2.1 || lambdac.M() > 2.4) continue;
          //cout<<"    Lambdac+: "<<lambdac<<endl;
          p4_LambdacPlus_.push_back(lambdac);
          p4_LambdacPlus_kaon_.push_back(lvk);
          p4_LambdacPlus_pion_.push_back(lvpi);
          p4_LambdacPlus_proton_.push_back(lvp);
          charge_LambdacPlus_kaon_.push_back(itrack.charge());
          charge_LambdacPlus_pion_.push_back(jtrack.charge());
          charge_LambdacPlus_proton_.push_back(ktrack.charge());
          para.clear(); covmat.clear();
          trkpara_LambdacPlus_kaon_ = unpackTrackParameters(itrack,"kaon_lambda");
          trkpara_LambdacPlus_pion_ = unpackTrackParameters(jtrack,"pion_lambda");
          trkpara_LambdacPlus_proton_ = unpackTrackParameters(ktrack,"proton_lambda");
          trkcovm_LambdacPlus_kaon_ = unpackTrackCovariance(itrack,"kaon_lambda");
          trkcovm_LambdacPlus_pion_ = unpackTrackCovariance(jtrack,"pion_lambda");
          trkcovm_LambdacPlus_proton_ = unpackTrackCovariance(ktrack,"proton_lambda");
          double pv_lxy = (_pv->x()*lambdac.Px()+
                           _pv->y()*lambdac.Py())/lambdac.Pt();
          pv_lxy_LambdacPlus_.push_back(pv_lxy);
        }
      }
    }
  }
  nJet_ = ij;
  //cout<<"  n-jet: "<<ij<<endl;
}

template<typename TYPE>
LV TopAnalyzer::getLV(const TYPE& p4) const {
  return LV(p4.pt(), p4.eta(), p4.phi(), p4.mass());
}

void TopAnalyzer::setbit(UShort_t& x, UShort_t bit) {
}

bool TopAnalyzer::isAncestor(const reco::Candidate* ancestor,
                             const reco::Candidate * particle) {
  if(ancestor == particle) return true;
  
  for(size_t i=0; i<particle->numberOfMothers(); i++) {
    if(isAncestor(ancestor,particle->mother(i))) {
      return true;
    }
  }
  return false;
}

bool TopAnalyzer::isAncestor(const reco::Candidate* particle,
                             double id) {
  if(abs(particle->pdgId()) == id) {
    //top_quarks.push_back(particle);
    return true;
  }
  for(size_t i=0; i<particle->numberOfMothers(); i++) {
    if(isAncestor(particle->mother(i),id)) {
      return true;
    }
  }
  return false;
}

std::vector<const pat::Jet*> TopAnalyzer::isTightJet(edm::Handle<edm::View<pat::Jet> > &handle) {
  std::vector<const pat::Jet *> list;
  for(edm::View<pat::Jet>::const_iterator i = handle->begin();
                                          i != handle->end(); i++) {
    if(i->pt() < 20.) continue;
    if(fabs(i->eta() > 2.4)) continue;
    //if(i->hadronFlavour() != 5) continue;
    double NHF = i->neutralHadronEnergyFraction();
    double NEMF = i->neutralEmEnergyFraction();
    double NumConst = i->chargedMultiplicity()+i->neutralMultiplicity();
    double MUF = i->muonEnergyFraction();
    double CHF = i->chargedHadronEnergyFraction();
    double CHM = i->chargedMultiplicity();
    double CEMF = i->chargedEmEnergyFraction();
    bool jetId = (NHF<0.9 && NEMF<0.9 && NumConst>1 && MUF<0.8 && CHF>0 && CHM>0 && CEMF<0.8);
    if(!jetId) continue;
    list.push_back(&(*i));
  }
  return list;
}

bool TopAnalyzer::matchingTrack(const reco::Candidate* cand,
                                const pat::PackedCandidate &track) {
  double deta = track.eta() - cand->eta();
  double dphi = reco::deltaPhi(track.phi(), cand->phi());
  double dR = sqrt(deta*deta + dphi*dphi);
  if(fabs(dR) < 0.02) {
    return true;
  }
  return false;
}

std::vector<const pat::PackedCandidate*> TopAnalyzer::matchingTrack(const reco::Candidate* cand, 
                                                 edm::Handle<edm::View<pat::PackedCandidate> > &handle) {
  std::vector<const pat::PackedCandidate *> list;
  for(edm::View<pat::PackedCandidate>::const_iterator i = handle->begin();
                                                      i != handle->end(); i++) {
    if(!(i->trackHighPurity()) ||
         i->numberOfPixelHits() == 0) continue;
    if(i->pt() < 1.0) continue;
    if(fabs(i->eta()) > 2.5) continue;
    double deta = cand->eta() - i->eta();
    double dphi = reco::deltaPhi(cand->phi(), i->phi());
    double dR = sqrt(deta*deta + dphi*dphi);
    if(fabs(dR) < 0.02) {
      list.push_back(&(*i));
    }
  }
  return list;
}

std::vector<std::vector<double> > TopAnalyzer::unpackTrackParameters(const pat::PackedCandidate& track,
                                                                           std::string name) {
  para.clear();
  double curvature, lambda, phi, dxy, dsz;
  curvature = track.charge()/track.p();
  lambda = Pi()/2 - track.theta();
  phi = track.phi();
  dxy = -track.vx()*Sin(track.phi()) + track.vy()*Cos(track.phi());
  dsz = track.vz()*Cos(lambda) - (track.vx()*Cos(track.phi()) + track.vy()*Sin(track.phi()))*Sin(lambda);
  //cout<<"      "<<name<<" parameters: "<<curvature<<"  "<<lambda<<"  "<<phi<<"  "<<dxy<<"  "<<dsz<<endl;
  para.push_back(curvature); para.push_back(lambda); para.push_back(phi); para.push_back(dxy); para.push_back(dsz);
  if(name=="kaon_dp") { kdp_par.push_back(para); return kdp_par; }
  if(name=="pion1_dp") { pi1dp_par.push_back(para); return pi1dp_par; }
  if(name=="pion2_dp") { pi2dp_par.push_back(para); return pi2dp_par; }
  if(name=="kaon_d0") { kd0_par.push_back(para); return kd0_par; }
  if(name=="pion_d0") { pid0_par.push_back(para); return pid0_par; }
  if(name=="kaon_dt") { kdt_par.push_back(para); return kdt_par; }
  if(name=="pion_dt") { pidt_par.push_back(para); return pidt_par; }
  if(name=="pionsoft_dt") { pisdt_par.push_back(para); return pisdt_par; }
  if(name=="kaon1_ds") { k1ds_par.push_back(para); return k1ds_par; }
  if(name=="kaon2_ds") { k2ds_par.push_back(para); return k2ds_par; }
  if(name=="pion_ds") { pids_par.push_back(para); return pids_par; }
  if(name=="kaon_lambda") { klambda_par.push_back(para); return klambda_par; }
  if(name=="pion_lambda") { pilambda_par.push_back(para); return pilambda_par; }
  if(name=="proton_lambda") { plambda_par.push_back(para); return plambda_par; }
  else return par;
}

std::vector<std::vector<double> > TopAnalyzer::unpackTrackCovariance(const pat::PackedCandidate& track,
                                                                           std::string name) {
  covmat.clear();
  reco::TrackBase::CovarianceMatrix cov;
  cov = track.pseudoTrack().covariance();
  for(int i=0; i<cov.kRows; i++)
    for(int j=0; j<cov.kRows; j++)
      covmat.push_back(cov(i,j));
  //cout<<"    covariance matrix: "<<cov<<endl;
  if(name=="kaon_dp") { kdp_cov.push_back(covmat); return kdp_cov; }
  if(name=="pion1_dp") { pi1dp_cov.push_back(covmat); return pi1dp_cov; }
  if(name=="pion2_dp") { pi2dp_cov.push_back(covmat); return pi2dp_cov; }
  if(name=="kaon_d0") { kd0_cov.push_back(covmat); return kd0_cov; }
  if(name=="pion_d0") { pid0_cov.push_back(covmat); return pid0_cov; }
  if(name=="kaon_dt") { kdt_cov.push_back(covmat); return kdt_cov; }
  if(name=="pion_dt") { pidt_cov.push_back(covmat); return pidt_cov; }
  if(name=="pionsoft_dt") { pisdt_cov.push_back(covmat); return pisdt_cov; }
  if(name=="kaon1_ds") { k1ds_cov.push_back(covmat); return k1ds_cov; }
  if(name=="kaon2_ds") { k2ds_cov.push_back(covmat); return k2ds_cov; }
  if(name=="pion_ds") { pids_cov.push_back(covmat); return pids_cov; }
  if(name=="kaon_lambda") { klambda_cov.push_back(covmat); return klambda_cov; }
  if(name=="pion_lambda") { pilambda_cov.push_back(covmat); return pilambda_cov; }
  if(name=="proton_lambda") { plambda_cov.push_back(covmat); return plambda_cov; }
  else return par;
}

// ------------ method called once each job just before starting event loop  ------------
void 
TopAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TopAnalyzer::endJob() 
{
  //cout<<"N(K,pi) = "<<nkpi<<endl;
  //cout<<"N(K,pi,pi) = "<<nkpipi<<endl;
  //cout<<"N(K,K,pi) = "<<nkkpi<<endl;
  //cout<<"N(p,K,pi) = "<<npkpi<<endl;
  //cout<<"N(D*+) = "<<ndstar<<endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TopAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TopAnalyzer);
