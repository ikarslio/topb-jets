#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include "Math/LorentzVector.h"
#include "Math/Vector4D.h"
#include "TMatrixD.h"
#include "TRandom.h"
#include "ConstrainedFit.hh"
#define me 0.00051099895
#define mmu 0.1056583755
#define mk 0.493677
#define mpi 0.13957039 
#define mp 0.9382720813
#define mphi 1.019461
#define M_DPlus 1.86966
#define M_D0 1.86484
#define M_DsPlus 1.96835
#define M_LambdacPlus 2.28646

using namespace std;
using namespace ROOT;

typedef Math::LorentzVector<Math::PtEtaPhiM4D<double> > LV;
TRandom *_rnd;

Double_t deltaPhi(Double_t phi1, Double_t phi2)
{
  Double_t pi = 3.1415927;
  Double_t dphi = fabs(phi1 - phi2);
  if(dphi >= pi) dphi = 2. * pi - dphi;
  return dphi;
}

class SignalFit : public ConstrainedFit {
public:
  SignalFit() : ConstrainedFit(4,10,16) { }
  ~SignalFit() { }
};

void SignalYield() {
  cout<<"Calculating Signal Yield from Monte-Carlo"<<endl;
}

int main() {
  SignalYield();
  //TFile* file = TFile::Open("/depot/cms/top/chawla19/Topbjets/CMSSW_10_6_32/src/TopBJets/TopAnalyzer/ttbarbdecays_ul2018_15may.root");
  TFile* file = TFile::Open("ttbarbdecays_ul2018.root");

  TTree* T1 = (TTree*)file->Get("demo/Ntuple");

  UInt_t event;
  double pfMET;
  double nLepGen;
  double nJet;
  bool foundGenBlepton;
  bool foundGenB0decay;
  bool foundGenBPlusdecay;
  bool foundGenB0Dstardecay;
  bool foundGenBs0decay;
  bool foundGenLambdab0decay;
  std::vector<LV> * p4_muon = 0;
  std::vector<LV> * p4_electron = 0;
  std::vector<double> * charge_muon = 0;
  std::vector<double> * charge_electron = 0;
  std::vector<bool> * passLooseId_muon = 0;
  std::vector<bool> * passMediumId_muon = 0;
  std::vector<bool> * passTightId_muon = 0;
  std::vector<double> * pfIso_muon = 0;
  std::vector<bool> * passVetoId_electron = 0;
  std::vector<bool> * passLooseId_electron = 0;
  std::vector<bool> * passMediumId_electron = 0;
  std::vector<bool> * passTightId_electron = 0;
  std::vector<LV> * p4_gen_DPlus_lepton = 0;
  std::vector<LV> * p4_gen_D0_lepton = 0;
  std::vector<LV> * p4_gen_Dstar_lepton = 0;
  std::vector<LV> * p4_gen_DsPlus_lepton = 0;
  std::vector<LV> * p4_gen_LambdacPlus_lepton = 0;
  std::vector<LV> * p4_jet = 0;
  std::vector<double> * lepton_perJet = 0;
  std::vector<std::vector<double> > * pdgId_lepton = 0;
  std::vector<std::vector<double> > * pt_lepton = 0;
  std::vector<std::vector<double> > * eta_lepton = 0;
  std::vector<std::vector<double> > * phi_lepton = 0;
  std::vector<std::vector<double> > * charge_lepton = 0;
  std::vector<std::vector<double> > * mass_lepton = 0;
  std::vector<LV> * p4_DPlus = 0;
  std::vector<LV> * p4_DPlus_kaon = 0;
  std::vector<LV> * p4_DPlus_pion1 = 0;
  std::vector<LV> * p4_DPlus_pion2 = 0;
  std::vector<double> * charge_DPlus_kaon = 0;
  std::vector<std::vector<double>> * trkpara_DPlus_kaon = 0;
  std::vector<std::vector<double>> * trkpara_DPlus_pion1 = 0;
  std::vector<std::vector<double>> * trkpara_DPlus_pion2 = 0;
  std::vector<std::vector<double>> * trkcovm_DPlus_kaon = 0;
  std::vector<std::vector<double>> * trkcovm_DPlus_pion1 = 0;
  std::vector<std::vector<double>> * trkcovm_DPlus_pion2 = 0;
  std::vector<double> * pv_lxy_DPlus = 0;
  std::vector<LV> * p4_D0 = 0;
  std::vector<LV> * p4_D0_kaon = 0;
  std::vector<LV> * p4_D0_pion = 0;
  std::vector<double> * charge_D0_kaon = 0;
  std::vector<std::vector<double> > * trkpara_D0_kaon = 0;
  std::vector<std::vector<double> > * trkpara_D0_pion = 0;
  std::vector<std::vector<double> > * trkcovm_D0_kaon = 0;
  std::vector<std::vector<double> > * trkcovm_D0_pion = 0;
  std::vector<double> * pv_lxy_D0 = 0;
  std::vector<LV> * p4_Dstar_D0 = 0;
  std::vector<LV> * p4_Dstar_D0_kaon = 0;
  std::vector<LV> * p4_Dstar_D0_pion = 0;
  std::vector<LV> * p4_Dstar_pionsoft = 0;
  std::vector<double> * charge_Dstar_D0_kaon = 0;
  std::vector<std::vector<double>> * trkpara_Dstar_D0_kaon = 0;
  std::vector<std::vector<double>> * trkpara_Dstar_D0_pion = 0;
  std::vector<std::vector<double>> * trkcovm_Dstar_D0_kaon = 0;
  std::vector<std::vector<double>> * trkcovm_Dstar_D0_pion = 0;
  std::vector<double> * pv_lxy_Dstar_D0 = 0;
  std::vector<double> * mdiff_Dstar_D0 = 0;
  std::vector<LV> * p4_DsPlus = 0;
  std::vector<LV> * p4_DsPlus_kaon1 = 0;
  std::vector<LV> * p4_DsPlus_kaon2 = 0;
  std::vector<LV> * p4_DsPlus_pion = 0;
  std::vector<double> * charge_DsPlus_pion = 0;
  std::vector<std::vector<double>> * trkpara_DsPlus_kaon1 = 0;
  std::vector<std::vector<double>> * trkpara_DsPlus_kaon2 = 0;
  std::vector<std::vector<double>> * trkpara_DsPlus_pion = 0;
  std::vector<std::vector<double>> * trkcovm_DsPlus_kaon1 = 0;
  std::vector<std::vector<double>> * trkcovm_DsPlus_kaon2 = 0;
  std::vector<std::vector<double>> * trkcovm_DsPlus_pion = 0;
  std::vector<double> * pv_lxy_DsPlus = 0;
  std::vector<LV> * p4_LambdacPlus = 0;
  std::vector<LV> * p4_LambdacPlus_kaon = 0;
  std::vector<LV> * p4_LambdacPlus_pion = 0;
  std::vector<LV> * p4_LambdacPlus_proton = 0;
  std::vector<double> * charge_LambdacPlus_kaon = 0;
  std::vector<std::vector<double>> * trkpara_LambdacPlus_kaon = 0;
  std::vector<std::vector<double>> * trkpara_LambdacPlus_pion = 0;
  std::vector<std::vector<double>> * trkpara_LambdacPlus_proton = 0;
  std::vector<std::vector<double>> * trkcovm_LambdacPlus_kaon = 0;
  std::vector<std::vector<double>> * trkcovm_LambdacPlus_pion = 0;
  std::vector<std::vector<double>> * trkcovm_LambdacPlus_proton = 0;
  std::vector<double> * pv_lxy_LambdacPlus = 0;
  T1->SetBranchAddress("event", &event);
  T1->SetBranchAddress("pfMET", &pfMET);
  T1->SetBranchAddress("nLepGen", &nLepGen);
  T1->SetBranchAddress("nJet", &nJet);
  T1->SetBranchAddress("foundGenBlepton", &foundGenBlepton);
  T1->SetBranchAddress("foundGenB0decay", &foundGenB0decay);
  T1->SetBranchAddress("foundGenBPlusdecay", &foundGenBPlusdecay);
  T1->SetBranchAddress("foundGenB0Dstardecay", &foundGenB0Dstardecay);
  T1->SetBranchAddress("foundGenBs0decay", &foundGenBs0decay);
  T1->SetBranchAddress("foundGenLambdab0decay", &foundGenLambdab0decay);
  T1->SetBranchAddress("p4_muon", &p4_muon);
  T1->SetBranchAddress("p4_electron", &p4_electron);
  T1->SetBranchAddress("charge_muon", &charge_muon);
  T1->SetBranchAddress("charge_electron", &charge_electron);
  T1->SetBranchAddress("passLooseId_muon", &passLooseId_muon);
  T1->SetBranchAddress("passMediumId_muon", &passMediumId_muon);
  T1->SetBranchAddress("passTightId_muon", &passTightId_muon);
  T1->SetBranchAddress("pfIso_muon", &pfIso_muon);
  T1->SetBranchAddress("passVetoId_electron", &passVetoId_electron);
  T1->SetBranchAddress("passLooseId_electron", &passLooseId_electron);
  T1->SetBranchAddress("passMediumId_electron", &passMediumId_electron);
  T1->SetBranchAddress("passTightId_electron", &passTightId_electron);
  T1->SetBranchAddress("p4_gen_DPlus_lepton", &p4_gen_DPlus_lepton);
  T1->SetBranchAddress("p4_gen_D0_lepton", &p4_gen_D0_lepton);
  T1->SetBranchAddress("p4_gen_Dstar_lepton", &p4_gen_Dstar_lepton);
  T1->SetBranchAddress("p4_gen_DsPlus_lepton", &p4_gen_DsPlus_lepton);
  T1->SetBranchAddress("p4_gen_LambdacPlus_lepton", &p4_gen_LambdacPlus_lepton);
  T1->SetBranchAddress("p4_jet", &p4_jet);
  T1->SetBranchAddress("lepton_perJet", &lepton_perJet);
  T1->SetBranchAddress("pdgId_lepton", &pdgId_lepton);
  T1->SetBranchAddress("pt_lepton", &pt_lepton);
  T1->SetBranchAddress("eta_lepton", &eta_lepton);
  T1->SetBranchAddress("phi_lepton", &phi_lepton);
  T1->SetBranchAddress("charge_lepton", &charge_lepton);
  T1->SetBranchAddress("mass_lepton", &mass_lepton);
  T1->SetBranchAddress("p4_DPlus", &p4_DPlus);
  T1->SetBranchAddress("p4_DPlus_kaon", &p4_DPlus_kaon);
  T1->SetBranchAddress("p4_DPlus_pion1", &p4_DPlus_pion1);
  T1->SetBranchAddress("p4_DPlus_pion2", &p4_DPlus_pion2);
  T1->SetBranchAddress("charge_DPlus_kaon", &charge_DPlus_kaon);
  T1->SetBranchAddress("trkpara_DPlus_kaon", &trkpara_DPlus_kaon);
  T1->SetBranchAddress("trkpara_DPlus_pion1", &trkpara_DPlus_pion1);
  T1->SetBranchAddress("trkpara_DPlus_pion2", &trkpara_DPlus_pion2);
  T1->SetBranchAddress("trkcovm_DPlus_kaon", &trkcovm_DPlus_kaon);
  T1->SetBranchAddress("trkcovm_DPlus_pion1", &trkcovm_DPlus_pion1);
  T1->SetBranchAddress("trkcovm_DPlus_pion2", &trkcovm_DPlus_pion2);
  T1->SetBranchAddress("pv_lxy_DPlus", &pv_lxy_DPlus);
  T1->SetBranchAddress("p4_D0", &p4_D0);
  T1->SetBranchAddress("p4_D0_kaon", &p4_D0_kaon);
  T1->SetBranchAddress("p4_D0_pion", &p4_D0_pion);
  T1->SetBranchAddress("charge_D0_kaon", &charge_D0_kaon);
  T1->SetBranchAddress("trkpara_D0_kaon", &trkpara_D0_kaon);
  T1->SetBranchAddress("trkpara_D0_pion", &trkpara_D0_pion);
  T1->SetBranchAddress("trkcovm_D0_kaon", &trkcovm_D0_kaon);
  T1->SetBranchAddress("trkcovm_D0_pion", &trkcovm_D0_pion);
  T1->SetBranchAddress("pv_lxy_D0", &pv_lxy_D0);
  T1->SetBranchAddress("p4_Dstar_D0", &p4_Dstar_D0);
  T1->SetBranchAddress("p4_Dstar_D0_kaon", &p4_Dstar_D0_kaon);
  T1->SetBranchAddress("p4_Dstar_D0_pion", &p4_Dstar_D0_pion);
  T1->SetBranchAddress("p4_Dstar_pionsoft", &p4_Dstar_pionsoft);
  T1->SetBranchAddress("charge_Dstar_D0_kaon", &charge_Dstar_D0_kaon);
  T1->SetBranchAddress("trkpara_Dstar_D0_kaon", &trkpara_Dstar_D0_kaon);
  T1->SetBranchAddress("trkpara_Dstar_D0_pion", &trkpara_Dstar_D0_pion);
  T1->SetBranchAddress("trkcovm_Dstar_D0_kaon", &trkcovm_Dstar_D0_kaon);
  T1->SetBranchAddress("trkcovm_Dstar_D0_pion", &trkcovm_Dstar_D0_pion);
  T1->SetBranchAddress("pv_lxy_Dstar_D0", &pv_lxy_Dstar_D0);
  T1->SetBranchAddress("mdiff_Dstar_D0", &mdiff_Dstar_D0);
  T1->SetBranchAddress("p4_DsPlus", &p4_DsPlus);
  T1->SetBranchAddress("p4_DsPlus_kaon1", &p4_DsPlus_kaon1);
  T1->SetBranchAddress("p4_DsPlus_kaon2", &p4_DsPlus_kaon2);
  T1->SetBranchAddress("p4_DsPlus_pion", &p4_DsPlus_pion);
  T1->SetBranchAddress("charge_DsPlus_pion", &charge_DsPlus_pion);
  T1->SetBranchAddress("trkpara_DsPlus_kaon1", &trkpara_DsPlus_kaon1);
  T1->SetBranchAddress("trkpara_DsPlus_kaon2", &trkpara_DsPlus_kaon2);
  T1->SetBranchAddress("trkpara_DsPlus_pion", &trkpara_DsPlus_pion);
  T1->SetBranchAddress("trkcovm_DsPlus_kaon1", &trkcovm_DsPlus_kaon1);
  T1->SetBranchAddress("trkcovm_DsPlus_kaon2", &trkcovm_DsPlus_kaon2);
  T1->SetBranchAddress("trkcovm_DsPlus_pion", &trkcovm_DsPlus_pion);
  T1->SetBranchAddress("pv_lxy_DsPlus", &pv_lxy_DsPlus);
  T1->SetBranchAddress("p4_LambdacPlus", &p4_LambdacPlus);
  T1->SetBranchAddress("p4_LambdacPlus_kaon", &p4_LambdacPlus_kaon);
  T1->SetBranchAddress("p4_LambdacPlus_pion", &p4_LambdacPlus_pion);
  T1->SetBranchAddress("p4_LambdacPlus_proton", &p4_LambdacPlus_proton);
  T1->SetBranchAddress("charge_LambdacPlus_kaon", &charge_LambdacPlus_kaon);
  T1->SetBranchAddress("trkpara_LambdacPlus_kaon", &trkpara_LambdacPlus_kaon);
  T1->SetBranchAddress("trkpara_LambdacPlus_pion", &trkpara_LambdacPlus_pion);
  T1->SetBranchAddress("trkpara_LambdacPlus_proton", &trkpara_LambdacPlus_proton);
  T1->SetBranchAddress("trkcovm_LambdacPlus_kaon", &trkcovm_LambdacPlus_kaon);
  T1->SetBranchAddress("trkcovm_LambdacPlus_pion", &trkcovm_LambdacPlus_pion);
  T1->SetBranchAddress("trkcovm_LambdacPlus_proton", &trkcovm_LambdacPlus_proton);
  T1->SetBranchAddress("pv_lxy_LambdacPlus", &pv_lxy_LambdacPlus);

  TFile* f1 = new TFile("ul2018.root", "RECREATE");
  TTree *tree = new TTree("tree", "ml-tree");
  // Branch variables
  std::vector<double> fPt, fLxy, fLxysig, fCt, fChi2, fMass, fMass_lepD;
  tree->Branch("fPt", &fPt);
  tree->Branch("fLxy", &fLxy);
  tree->Branch("fLxysig", &fLxysig);
  tree->Branch("fCt", &fCt);
  tree->Branch("fChi2", &fChi2);
  tree->Branch("fMass", &fMass);
  tree->Branch("fMass_lepD", &fMass_lepD);

  // Histograms
  TH1D *h_D0_mass = new TH1D("h_D0_mass","m(K,#pi)",50,1.73,1.98);
  TH1D *h_D0_pik_mass = new TH1D("h_D0_pik_mass","m(K,#pi)",50,1.73,1.98);
  TH1D *h_DstarD0_mass = new TH1D("h_DstarD0_mass","m(D^{*})-m(D^{0})",70,0.135,0.17);
  TH1D *h_DPlus_mass = new TH1D("h_DPlus_mass","m(K,#pi,#pi)",50,1.73,1.98);
  TH1D *h_DsPlus_DPlus_mass = new TH1D("h_DsPlus_DPlus_mass","m(K,#pi,#pi)",50,1.73,1.98);
  TH1D *h_LambdacPlus_DPlus_mass = new TH1D("h_LambdacPlus_DPlus_mass","m(K,#pi,#pi)",50,1.73,1.98);
  TH1D *h_DsPlus_phipi_mass = new TH1D("h_DsPlus_phipi_mass", "m(#phi,#pi)",50,1.85,2.1);
  TH1D *h_DsPlus_mass = new TH1D("h_DsPlus_mass","m(K,K,#pi)",50,1.85,2.1);
  TH1D *h_Dstar_DsPlus_mass = new TH1D("h_Dstar_DsPlus_mass","m(K,K,#pi)",50,1.85,2.1);
  TH1D *h_DPlus_DsPlus_mass = new TH1D("h_DPlus_DsPlus_mass","m(K,K,#pi)",50,1.85,2.1);
  TH1D *h_LambdacPlus_DsPlus_mass = new TH1D("h_LambdacPlus_DsPlus_mass","m(K,K,#pi)",50,1.85,2.1);
  TH1D *h_LambdacPlus_mass = new TH1D("h_LambdacPlus_mass","m(K,#pi,p)",50,2.15,2.4);
  TH1D *h_Dstar_LambdacPlus_mass = new TH1D("h_Dstar_LambdacPlus_mass","m(K,#pi,p)",50,2.15,2.4);
  TH1D *h_DPlus_LambdacPlus_mass = new TH1D("h_DPlus_LambdacPlus_mass","m(K,#pi,p)",50,2.15,2.4);
  TH1D *h_DsPlus_LambdacPlus_mass = new TH1D("h_DsPlus_LambdacPlus_mass","m(K,#pi,p)",50,2.15,2.4);

  int nentries = T1->GetEntries();
  cout<<"entries: "<<nentries<<endl;
  int gen = 0; int genlep = 0; int dilep = 0; int reco = 0; int bjets = 0; int lep1 = 0; int lep2 = 0;
  int genb0 = 0; int genbplus = 0; int genb0dstar = 0; int genbs0 = 0; int genlambdab0 = 0;
  bool mumu, mue, ee, emu, zmm, zee;
  std::vector<double> newlep1_pt, newlep2_pt,
                      newlep1_eta, newlep2_eta,
                      newlep1_phi, newlep2_phi;
  std::vector<double> pt_DPlus_lepton, eta_DPlus_lepton, phi_DPlus_lepton, charge_DPlus_lepton, mass_DPlus_lepton;
  std::vector<double> pt_D0_lepton, eta_D0_lepton, phi_D0_lepton, charge_D0_lepton, mass_D0_lepton;
  std::vector<double> pt_Dstar_lepton, eta_Dstar_lepton, phi_Dstar_lepton, charge_Dstar_lepton, mass_Dstar_lepton;
  std::vector<double> pt_DsPlus_lepton, eta_DsPlus_lepton, phi_DsPlus_lepton, charge_DsPlus_lepton, mass_DsPlus_lepton;
  std::vector<double> pt_LambdacPlus_lepton, eta_LambdacPlus_lepton, phi_LambdacPlus_lepton,
                      charge_LambdacPlus_lepton, mass_LambdacPlus_lepton;
  double i_par[5], j_par[5], k_par[5];
  double i_para[5], j_para[5], k_para[5], i_covm[25], j_covm[25], k_covm[25];
  _rnd = new TRandom();
  for(unsigned int jentry=0; jentry<nentries; jentry++) {
    T1->GetEntry(jentry);
    if(jentry%50000 == 0) cout<<"Events Processed :  "<<jentry<<endl;
    fPt.clear();
    fLxy.clear();
    fLxysig.clear();
    fCt.clear();
    fChi2.clear();
    fMass.clear();
    fMass_lepD.clear();

    cout<<"Event "<<event<<endl;
    newlep1_pt.clear();
    newlep1_eta.clear();
    newlep1_phi.clear();
    newlep2_pt.clear();
    newlep2_eta.clear();
    newlep2_phi.clear();
    pt_DPlus_lepton.clear();
    eta_DPlus_lepton.clear();
    phi_DPlus_lepton.clear();
    charge_DPlus_lepton.clear();
    mass_DPlus_lepton.clear();
    pt_D0_lepton.clear();
    eta_D0_lepton.clear();
    phi_D0_lepton.clear();
    charge_D0_lepton.clear();
    mass_D0_lepton.clear();
    pt_Dstar_lepton.clear();
    eta_Dstar_lepton.clear();
    phi_Dstar_lepton.clear();
    charge_Dstar_lepton.clear();
    mass_Dstar_lepton.clear();
    pt_DsPlus_lepton.clear();
    eta_DsPlus_lepton.clear();
    phi_DsPlus_lepton.clear();
    charge_DsPlus_lepton.clear();
    mass_DsPlus_lepton.clear();
    pt_LambdacPlus_lepton.clear();
    eta_LambdacPlus_lepton.clear();
    phi_LambdacPlus_lepton.clear();
    charge_LambdacPlus_lepton.clear();
    mass_LambdacPlus_lepton.clear();

    if(nLepGen > 1) gen++;
    if(foundGenBlepton) genlep++;
    if(foundGenB0decay) genb0++;
    if(foundGenBPlusdecay) genbplus++;
    if(foundGenB0Dstardecay) genb0dstar++;
    if(foundGenBs0decay) genbs0++;
    if(foundGenLambdab0decay) genlambdab0++;
    int indx[p4_muon->size()];float ptmu[p4_muon->size()];
    for(unsigned int i=0; i<p4_muon->size(); i++) {
      ptmu[i] = p4_muon->at(i).Pt();
    }
    int indy[p4_electron->size()];float ptel[p4_electron->size()];
    for(unsigned int i=0; i<p4_electron->size(); i++) {
      ptel[i] = p4_electron->at(i).Pt();
    }
    TMath::Sort(int(sizeof(ptmu)/sizeof(ptmu[0])),ptmu,indx,true);
    TMath::Sort(int(sizeof(ptel)/sizeof(ptel[0])),ptel,indy,true);

    mumu = false; mue = false; ee = false; emu = false;
    zmm = false; zee = false;
    for(unsigned int i=0; i<p4_muon->size(); i++) {
      for(unsigned int j=0; j<p4_muon->size(); j++) {
        if(i == j) continue;
        //cout<<"  muon-muon: "<<p4_muon->at(indx[i]).Pt()<<", "<<p4_muon->at(indx[j]).Pt()<<", "<<charge_muon->at(indx[i])<<", "<<charge_muon->at(indx[j])<<", "<<passTightId_muon->at(indx[i])<<", "<<passTightId_muon->at(indx[j])<<", "<<pfIso_muon->at(indx[i])<<", "<<pfIso_muon->at(indx[j])<<", ";
        if(p4_muon->at(indx[i]).Pt() < 25.0) continue;
        if(p4_muon->at(indx[i]).Pt() < p4_muon->at(indx[j]).Pt()) continue;
        if(charge_muon->at(indx[i])*charge_muon->at(indx[j]) > 0) continue;
        if(!passTightId_muon->at(indx[i])) continue;
        if(!passTightId_muon->at(indx[j])) continue;
        if(pfIso_muon->at(indx[i]) > 0.15) continue;
        if(pfIso_muon->at(indx[j]) > 0.15) continue;
        TLorentzVector tvm1, tvm2, tvmm;
        tvm1.SetPtEtaPhiM(p4_muon->at(indx[i]).Pt(),p4_muon->at(indx[i]).Eta(),p4_muon->at(indx[i]).Phi(),mmu);
        tvm2.SetPtEtaPhiM(p4_muon->at(indx[j]).Pt(),p4_muon->at(indx[j]).Eta(),p4_muon->at(indx[j]).Phi(),mmu);
        tvmm = tvm1+tvm2;
        if(tvmm.M() < 20.) continue;
        if(tvmm.M() > 76. && tvmm.M() < 106.) zmm = true;
        mumu = true;
        cout<<"  mu-mu: "<<p4_muon->at(indx[i]).Pt()<<"  "<<p4_muon->at(indx[j]).Pt()<<endl;
        newlep1_pt.push_back(p4_muon->at(indx[i]).Pt());
        newlep1_eta.push_back(p4_muon->at(indx[i]).Eta());
        newlep1_phi.push_back(p4_muon->at(indx[i]).Phi());
        newlep2_pt.push_back(p4_muon->at(indx[j]).Pt());
        newlep2_eta.push_back(p4_muon->at(indx[j]).Eta());
        newlep2_phi.push_back(p4_muon->at(indx[j]).Phi());        
      }
    }
    //cout<<endl;
    for(unsigned int i=0; i<p4_muon->size(); i++) {
      for(unsigned int j=0; j<p4_electron->size(); j++) {
        //cout<<"  muon-electron: "<<p4_muon->at(indx[i]).Pt()<<", "<<p4_electron->at(indy[j]).Pt()<<", "<<charge_muon->at(indx[i])<<", "<<charge_electron->at(indy[j])<<", "<<passTightId_muon->at(indx[i])<<", "<<pfIso_muon->at(indx[i])<<", "<<passTightId_electron->at(indy[j])<<", ";
        if(p4_muon->at(indx[i]).Pt() < 25.0) continue;
        if(p4_muon->at(indx[i]).Pt() < p4_electron->at(indy[j]).Pt()) continue;
        if(charge_muon->at(indx[i])*charge_electron->at(indy[j]) > 0) continue;
        if(!passTightId_muon->at(indx[i])) continue;
        if(pfIso_muon->at(indx[i]) > 0.15) continue;
        if(!passTightId_electron->at(indy[j])) continue;
        TLorentzVector tvm1, tve2, tvme;
        tvm1.SetPtEtaPhiM(p4_muon->at(indx[i]).Pt(),p4_muon->at(indx[i]).Eta(),p4_muon->at(indx[i]).Phi(),mmu);
        tve2.SetPtEtaPhiM(p4_electron->at(indy[j]).Pt(),p4_electron->at(indy[j]).Eta(),p4_electron->at(indy[j]).Phi(),me);
        tvme = tvm1+tve2;
        if(tvme.M() < 20.) continue;
        mue = true;
        cout<<"  mu-e: "<<p4_muon->at(indx[i]).Pt()<<"  "<<p4_electron->at(indy[j]).Pt()<<endl;
        newlep1_pt.push_back(p4_muon->at(indx[i]).Pt());
        newlep1_eta.push_back(p4_muon->at(indx[i]).Eta());
        newlep1_phi.push_back(p4_muon->at(indx[i]).Phi());
        newlep2_pt.push_back(p4_electron->at(indy[j]).Pt());
        newlep2_eta.push_back(p4_electron->at(indy[j]).Eta());
        newlep2_phi.push_back(p4_electron->at(indy[j]).Phi());
      }
    }
    //cout<<endl;
    for(unsigned int i=0; i<p4_electron->size(); i++) {
      for(unsigned int j=0; j<p4_electron->size(); j++) {
        if(i == j) continue;
        //cout<<"  electron-electron: "<<p4_electron->at(indy[i]).Pt()<<", "<<p4_electron->at(indy[j]).Pt()<<", "<<charge_electron->at(indy[i])<<", "<<charge_electron->at(indy[j])<<", "<<passTightId_electron->at(indy[i])<<", "<<passTightId_electron->at(indy[j])<<", ";
        if(p4_electron->at(indy[i]).Pt() < 25.0) continue;
        if(p4_electron->at(indy[i]).Pt() < p4_electron->at(indy[j]).Pt()) continue;
        if(charge_electron->at(indy[i])*charge_electron->at(indy[j]) > 0) continue;
        if(!passTightId_electron->at(indy[i])) continue;
        if(!passTightId_electron->at(indy[j])) continue;
        TLorentzVector tve1, tve2, tvee;
        tve1.SetPtEtaPhiM(p4_electron->at(indy[i]).Pt(),p4_electron->at(indy[i]).Eta(),p4_electron->at(indy[i]).Phi(),me);
        tve2.SetPtEtaPhiM(p4_electron->at(indy[j]).Pt(),p4_electron->at(indy[j]).Eta(),p4_electron->at(indy[j]).Phi(),me);
        tvee = tve1+tve2;
        if(tvee.M() < 20.) continue;
        if(tvee.M() > 76. && tvee.M() < 106.) zee = true;
        ee = true;
        cout<<"  e-e: "<<p4_electron->at(indy[i]).Pt()<<"  "<<p4_electron->at(indy[j]).Pt()<<endl;
        newlep1_pt.push_back(p4_electron->at(indy[i]).Pt());
        newlep1_eta.push_back(p4_electron->at(indy[i]).Eta());
        newlep1_phi.push_back(p4_electron->at(indy[i]).Phi());
        newlep2_pt.push_back(p4_electron->at(indy[j]).Pt());
        newlep2_eta.push_back(p4_electron->at(indy[j]).Eta());
        newlep2_phi.push_back(p4_electron->at(indy[j]).Phi());
      }
    }
    //cout<<endl;
    for(unsigned int i=0; i<p4_electron->size(); i++) {
      for(unsigned int j=0; j<p4_muon->size(); j++) {
        //cout<<"  electron-muon: "<<p4_electron->at(indy[i]).Pt()<<", "<<p4_muon->at(indx[j]).Pt()<<", "<<charge_electron->at(indy[i])<<", "<<charge_muon->at(indx[j])<<", "<<passTightId_electron->at(indy[i])<<", "<<passTightId_muon->at(indx[j])<<", "<<pfIso_muon->at(indx[j])<<", ";
        if(p4_electron->at(indy[i]).Pt() < 25.0) continue;
        if(p4_electron->at(indy[i]).Pt() < p4_muon->at(indx[j]).Pt()) continue;
        if(charge_electron->at(indy[i])*charge_muon->at(indx[j]) > 0) continue;
        if(!passTightId_electron->at(indy[i])) continue;
        if(!passTightId_muon->at(indx[j])) continue;
        if(pfIso_muon->at(indx[j]) > 0.15) continue;
        TLorentzVector tve1, tvm2, tvem;
        tve1.SetPtEtaPhiM(p4_electron->at(indy[i]).Pt(),p4_electron->at(indy[i]).Eta(),p4_electron->at(indy[i]).Phi(),me);
        tvm2.SetPtEtaPhiM(p4_muon->at(indx[j]).Pt(),p4_muon->at(indx[j]).Eta(),p4_muon->at(indx[j]).Phi(),mmu);
        tvem = tve1+tvm2;
        if(tvem.M() < 20.) continue;
        emu = true;
        cout<<"  e-mu: "<<p4_electron->at(indy[i]).Pt()<<"  "<<p4_muon->at(indx[j]).Pt()<<endl;
        newlep1_pt.push_back(p4_electron->at(indy[i]).Pt());
        newlep1_eta.push_back(p4_electron->at(indy[i]).Eta());
        newlep1_phi.push_back(p4_electron->at(indy[i]).Phi());
        newlep2_pt.push_back(p4_muon->at(indx[j]).Pt());
        newlep2_eta.push_back(p4_muon->at(indx[j]).Eta());
        newlep2_phi.push_back(p4_muon->at(indx[j]).Phi());
      }
    }
    //cout<<endl;
    //if(!mumu && !mue && !ee && !emu) continue;
    //if(zmm || zee) continue;
    //if(pfMET < 40.) continue;
    //dilep++;
    int ij = 0;
    if((mumu || mue || ee || emu) && (!zmm && !zee) && pfMET > 40.) {
      dilep++;
      for(unsigned int i=0; i<p4_jet->size(); i++) {
        for(unsigned int j=0; j<newlep1_pt.size(); j++) {
          cout<<"  i: "<<i<<"  j: "<<j<<endl;
          double deta1 = p4_jet->at(i).Eta() - newlep1_eta.at(j);
          double dphi1 = deltaPhi(p4_jet->at(i).Phi(), newlep1_phi.at(j));
          double dR1 = sqrt(deta1*deta1+dphi1*dphi1);
          double deta2 = p4_jet->at(i).Eta() - newlep2_eta.at(j);
          double dphi2 = deltaPhi(p4_jet->at(i).Phi(), newlep2_phi.at(j));
          double dR2 = sqrt(deta2*deta2+dphi2*dphi2);
          cout<<"  jet: "<<p4_jet->at(i).Pt()<<endl;
          cout<<"  lep-lep: "<<newlep1_pt.at(j)<<"  "<<newlep2_pt.at(j)<<endl;
          cout<<"    dR between lepton and jet: "<<dR1<<"  "<<dR2<<endl;
          if(dR1 < 0.4 || dR2 < 0.4) continue;
          ij += 1;
          cout<<"    ij: "<<ij<<endl;
        }
        for(unsigned int k=0; k<pt_lepton->at(i).size(); k++) {
          cout<<"    lepton in b-jet: "<<pt_lepton->at(i).at(k)<<"  "<<eta_lepton->at(i).at(k)<<endl;
        }
        cout<<"    n-lep per jet: "<<lepton_perJet->at(i)<<endl;
      }
    }
    bool foundlep = false;
    bool genmatch = false;
    if(ij > 0) reco++;
    bjets = bjets+ij;
    cout<<"  jets: "<<ij<<endl;
    double deta, dphi, dR;
    for(int i=0; i<p4_jet->size(); i++) {
      for(unsigned int j=0; j<pt_lepton->at(i).size(); j++) {
        cout<<"    reco: "<<pt_lepton->at(i).at(j)<<"  "
                          <<eta_lepton->at(i).at(j)<<"  "
                          <<phi_lepton->at(i).at(j)<<endl;
        if(!foundGenB0decay &&
           !foundGenBPlusdecay &&
           !foundGenB0Dstardecay &&
           !foundGenBs0decay &&
           !foundGenLambdab0decay) continue;
        if(foundGenB0decay) {
          for(int k=0; k<p4_gen_DPlus_lepton->size(); k++) {
            cout<<"    gen (b0): "<<p4_gen_DPlus_lepton->at(k).Pt()<<"  "
                                  <<p4_gen_DPlus_lepton->at(k).Eta()<<"  "
                                  <<p4_gen_DPlus_lepton->at(k).Phi()<<endl;
            deta = eta_lepton->at(i).at(j) - p4_gen_DPlus_lepton->at(k).Eta();
            dphi = deltaPhi(phi_lepton->at(i).at(j), p4_gen_DPlus_lepton->at(k).Phi());
            dR = sqrt(deta*deta+dphi*dphi);
            cout<<"    dR: "<<dR<<endl;
            if(dR < 0.02) {
              genmatch = true;
              pt_DPlus_lepton.push_back(pt_lepton->at(i).at(j));
              eta_DPlus_lepton.push_back(eta_lepton->at(i).at(j));
              phi_DPlus_lepton.push_back(phi_lepton->at(i).at(j));
              charge_DPlus_lepton.push_back(charge_lepton->at(i).at(j));
              mass_DPlus_lepton.push_back(mass_lepton->at(i).at(j));
            }
          }
        }
        if(foundGenBPlusdecay) {
          for(int k=0; k<p4_gen_D0_lepton->size(); k++) {
            cout<<"    gen (b+): "<<p4_gen_D0_lepton->at(k).Pt()<<"  "
                                  <<p4_gen_D0_lepton->at(k).Eta()<<"  "
                                  <<p4_gen_D0_lepton->at(k).Phi()<<endl;
            deta = eta_lepton->at(i).at(j) - p4_gen_D0_lepton->at(k).Eta();
            dphi = deltaPhi(phi_lepton->at(i).at(j), p4_gen_D0_lepton->at(k).Phi());
            dR = sqrt(deta*deta+dphi*dphi);
            cout<<"    dR: "<<dR<<endl;
            if(dR < 0.02) {
              genmatch = true;
              pt_D0_lepton.push_back(pt_lepton->at(i).at(j));
              eta_D0_lepton.push_back(eta_lepton->at(i).at(j));
              phi_D0_lepton.push_back(phi_lepton->at(i).at(j));
              charge_D0_lepton.push_back(charge_lepton->at(i).at(j));
              mass_D0_lepton.push_back(mass_lepton->at(i).at(j));
            }
          }
        }
        if(foundGenB0Dstardecay) {
          for(int k=0; k<p4_gen_Dstar_lepton->size(); k++) {
            cout<<"    gen (b0 d*+): "<<p4_gen_Dstar_lepton->at(k).Pt()<<"  "
                                      <<p4_gen_Dstar_lepton->at(k).Eta()<<"  "
                                      <<p4_gen_Dstar_lepton->at(k).Phi()<<endl;
            deta = eta_lepton->at(i).at(j) - p4_gen_Dstar_lepton->at(k).Eta();
            dphi = deltaPhi(phi_lepton->at(i).at(j), p4_gen_Dstar_lepton->at(k).Phi());
            dR = sqrt(deta*deta+dphi*dphi);
            cout<<"    dR: "<<dR<<endl;
            if(dR < 0.02) {
              genmatch = true;
              pt_Dstar_lepton.push_back(pt_lepton->at(i).at(j));
              eta_Dstar_lepton.push_back(eta_lepton->at(i).at(j));
              phi_Dstar_lepton.push_back(phi_lepton->at(i).at(j));
              charge_Dstar_lepton.push_back(charge_lepton->at(i).at(j));
              mass_Dstar_lepton.push_back(mass_lepton->at(i).at(j));
            }
          }
        }
        if(foundGenBs0decay) {
          for(int k=0; k<p4_gen_DsPlus_lepton->size(); k++) {
            cout<<"    gen (bs0): "<<p4_gen_DsPlus_lepton->at(k).Pt()<<"  "
                                   <<p4_gen_DsPlus_lepton->at(k).Eta()<<"  "
                                   <<p4_gen_DsPlus_lepton->at(k).Phi()<<endl;
            deta = eta_lepton->at(i).at(j) - p4_gen_DsPlus_lepton->at(k).Eta();
            dphi = deltaPhi(phi_lepton->at(i).at(j), p4_gen_DsPlus_lepton->at(k).Phi());
            dR = sqrt(deta*deta+dphi*dphi);
            cout<<"    dR: "<<dR<<endl;
            if(dR < 0.02) {
              genmatch = true;
              pt_DsPlus_lepton.push_back(pt_lepton->at(i).at(j));
              eta_DsPlus_lepton.push_back(eta_lepton->at(i).at(j));
              phi_DsPlus_lepton.push_back(phi_lepton->at(i).at(j));
              charge_DsPlus_lepton.push_back(charge_lepton->at(i).at(j));
              mass_DsPlus_lepton.push_back(mass_lepton->at(i).at(j));
            }
          }
        }
        if(foundGenLambdab0decay) {
          for(int k=0; k<p4_gen_LambdacPlus_lepton->size(); k++) {
            cout<<"    gen (lambda_b0): "<<p4_gen_LambdacPlus_lepton->at(k).Pt()<<"  "
                                         <<p4_gen_LambdacPlus_lepton->at(k).Eta()<<"  "
                                         <<p4_gen_LambdacPlus_lepton->at(k).Phi()<<endl;
            deta = eta_lepton->at(i).at(j) - p4_gen_LambdacPlus_lepton->at(k).Eta();
            dphi = deltaPhi(phi_lepton->at(i).at(j), p4_gen_LambdacPlus_lepton->at(k).Phi());
            dR = sqrt(deta*deta+dphi*dphi);
            cout<<"    dR: "<<dR<<endl;
            if(dR < 0.02) {
              genmatch = true;
              pt_LambdacPlus_lepton.push_back(pt_lepton->at(i).at(j));
              eta_LambdacPlus_lepton.push_back(eta_lepton->at(i).at(j));
              phi_LambdacPlus_lepton.push_back(phi_lepton->at(i).at(j));
              charge_LambdacPlus_lepton.push_back(charge_lepton->at(i).at(j));
              mass_LambdacPlus_lepton.push_back(mass_lepton->at(i).at(j));
            }
          }
        }
      }
      if(lepton_perJet->at(i) > 0) foundlep = true;
    }
    cout<<"  is gen-reco matched lepton: "<<genmatch<<endl;
    cout<<"  found lepton: "<<foundlep<<endl;
    if(ij > 0 && foundlep) lep1++;
    if(ij > 0 && genmatch) lep2++;
    //if(!genmatch) continue;
    // D0 Geometric Fitting - j and k correspond to kaon and pion tracks
    for(unsigned int i=0; i<p4_D0->size(); i++) {
      cout<<"  reco d0 mass: "<<p4_D0->at(i).M()<<endl;
      TLorentzVector k, pi, kpi;
      k.SetPtEtaPhiM(p4_D0_kaon->at(i).Pt(),
                     p4_D0_kaon->at(i).Eta(),
                     p4_D0_kaon->at(i).Phi(),mpi);
      pi.SetPtEtaPhiM(p4_D0_pion->at(i).Pt(),
                      p4_D0_pion->at(i).Eta(),
                      p4_D0_pion->at(i).Phi(),mk);
      kpi = k+pi;
      for(unsigned int j=0; j<trkpara_D0_kaon->at(i).size(); j++) {
        j_par[j] = { trkpara_D0_kaon->at(i).at(j) };
        k_par[j] = { trkpara_D0_pion->at(i).at(j) };
        const int index[5] = {0, 6, 12, 18, 24};
        j_para[j] = j_par[j] + _rnd->Gaus()*sqrt(trkcovm_D0_kaon->at(i).at(index[j]));
        k_para[j] = k_par[j] + _rnd->Gaus()*sqrt(trkcovm_D0_pion->at(i).at(index[j]));
      }
      for(unsigned int j=0; j<trkcovm_D0_kaon->at(i).size(); j++) {
        j_covm[j] = { trkcovm_D0_kaon->at(i).at(j) };
        k_covm[j] = { trkcovm_D0_pion->at(i).at(j) };
      }
      CmsFit3D fit(2);
      fit.SetVerbose(2);
      fit.AddTrack(j_para,j_covm,mk);
      fit.AddTrack(k_para,k_covm,mpi);
      fit.Fit();
      double pv_lxy = pv_lxy_D0->at(i);
      double pt = fit.Pt();
      double lxy = fit.S();
      double decay_length = lxy - pv_lxy;
      double signif = decay_length/fit.ErrS();
      double proper_time = (decay_length/pt)*M_D0;
      cout<<"  fitted mass: "<<fit.VertexMass()<<endl;
      //h_D0_mass->Fill(p4_D0->at(i).M());
      //h_D0_mass->Fill(fit.VertexMass());
      TLorentzVector lep, d0, bp;
      d0.SetPtEtaPhiM(pt,p4_D0->at(i).Eta(),p4_D0->at(i).Phi(),fit.VertexMass());
      for(unsigned int j = 0; j<pt_D0_lepton.size(); j++) {
        lep.SetPtEtaPhiM(pt_D0_lepton.at(j),eta_D0_lepton.at(j),phi_D0_lepton.at(j),mass_D0_lepton.at(j));
        bp = d0+lep;
        if(bp.M() < 2.0 || bp.M() > 6.0) continue;
        cout<<"    lep: "<<pt_D0_lepton.at(j)<<endl;
        cout<<"    B+: "<<bp.Pt()<<"  "<<bp.M()<<endl;
        double charge = charge_D0_lepton.at(j)*charge_D0_kaon->at(i);
        if(charge > 0) {
          fPt.push_back(pt);
          fLxy.push_back(decay_length);
          fLxysig.push_back(signif);
          fCt.push_back(proper_time);
          fChi2.push_back(fit.Chi2());
          fMass.push_back(p4_D0->at(i).M());
          fMass_lepD.push_back(bp.M());
          h_D0_mass->Fill(p4_D0->at(i).M());
        }
        h_D0_pik_mass->Fill(kpi.M());
      }
    }
    // D0-D* Geometric Fitting - i and j correspond to kaon and pion tracks
    for(unsigned int i=0; i<p4_Dstar_D0->size(); i++) {
      cout<<"  reco d*-d0 mass: "<<p4_Dstar_D0->at(i).M()<<endl;
      if(p4_Dstar_pionsoft->at(i).Pt() < 1.) continue;
      TLorentzVector k, pi1, pis, k1, k2, kkpi, kpik, p1, p2, kppi, kpip;
      k.SetPtEtaPhiM(p4_Dstar_D0_kaon->at(i).Pt(),
                     p4_Dstar_D0_kaon->at(i).Eta(),
                     p4_Dstar_D0_kaon->at(i).Phi(),
                     p4_Dstar_D0_kaon->at(i).M());
      pi1.SetPtEtaPhiM(p4_Dstar_D0_pion->at(i).Pt(),
                       p4_Dstar_D0_pion->at(i).Eta(),
                       p4_Dstar_D0_pion->at(i).Phi(),
                       p4_Dstar_D0_pion->at(i).M());
      pis.SetPtEtaPhiM(p4_Dstar_pionsoft->at(i).Pt(),
                       p4_Dstar_pionsoft->at(i).Eta(),
                       p4_Dstar_pionsoft->at(i).Phi(),
                       p4_Dstar_pionsoft->at(i).M());
      k1.SetPtEtaPhiM(p4_Dstar_D0_pion->at(i).Pt(),
                      p4_Dstar_D0_pion->at(i).Eta(),
                      p4_Dstar_D0_pion->at(i).Phi(),mk);
      k2.SetPtEtaPhiM(p4_Dstar_pionsoft->at(i).Pt(),
                      p4_Dstar_pionsoft->at(i).Eta(),
                      p4_Dstar_pionsoft->at(i).Phi(),mk);
      p1.SetPtEtaPhiM(p4_Dstar_D0_pion->at(i).Pt(),
                      p4_Dstar_D0_pion->at(i).Eta(),
                      p4_Dstar_D0_pion->at(i).Phi(),mp);
      p2.SetPtEtaPhiM(p4_Dstar_pionsoft->at(i).Pt(),
                      p4_Dstar_pionsoft->at(i).Eta(),
                      p4_Dstar_pionsoft->at(i).Phi(),mp);
      kkpi = k+k1+pis; // Ds+
      kpik = k+pi1+k2; // Ds+
      kppi = k+p1+pis; // Lambdac+
      kpip = k+pi1+p2; // Lambdac+
      for(unsigned int j=0; j<trkpara_Dstar_D0_kaon->at(i).size(); j++) {
        j_par[j] = { trkpara_Dstar_D0_kaon->at(i).at(j) };
        k_par[j] = { trkpara_Dstar_D0_pion->at(i).at(j) };
        const int index[5] = {0, 6, 12, 18, 24};
        j_para[j] = j_par[j] + _rnd->Gaus()*sqrt(trkcovm_Dstar_D0_kaon->at(i).at(index[j]));
        k_para[j] = k_par[j] + _rnd->Gaus()*sqrt(trkcovm_Dstar_D0_pion->at(i).at(index[j]));
      }
      for(unsigned int j=0; j<trkcovm_Dstar_D0_kaon->at(i).size(); j++) {
        j_covm[j] = { trkcovm_Dstar_D0_kaon->at(i).at(j) };
        k_covm[j] = { trkcovm_Dstar_D0_pion->at(i).at(j) };
      }
      CmsFit3D fit(2);
      fit.SetVerbose(0);
      fit.AddTrack(j_para,j_covm,mk);
      fit.AddTrack(k_para,k_covm,mpi);
      fit.Fit();
      LV pi, d0, dt;
      double epi = sqrt(p4_Dstar_pionsoft->at(i).Px()*p4_Dstar_pionsoft->at(i).Px()+
                        p4_Dstar_pionsoft->at(i).Py()*p4_Dstar_pionsoft->at(i).Py()+
                        p4_Dstar_pionsoft->at(i).Pz()*p4_Dstar_pionsoft->at(i).Pz()+mpi*mpi);
      double ed0 = sqrt(fit.Px()*fit.Px()+
                        fit.Py()*fit.Py()+
                        fit.Pz()*fit.Pz()+
                        fit.VertexMass()*fit.VertexMass());
      pi.SetPxPyPzE(p4_Dstar_pionsoft->at(i).Px(),p4_Dstar_pionsoft->at(i).Py(),p4_Dstar_pionsoft->at(i).Pz(),epi);
      d0.SetPxPyPzE(fit.Px(),fit.Py(),fit.Pz(),ed0);
      dt = pi+d0;
      double pv_lxy = pv_lxy_Dstar_D0->at(i);
      double pt = fit.Pt();
      double lxy = fit.S();
      double decay_length = lxy - pv_lxy;
      double signif = decay_length/fit.ErrS();
      double proper_time = (decay_length/pt)*M_D0;
      cout<<"  fitted mass: "<<fit.VertexMass()<<endl;
      //h_DstarD0_mass->Fill(mdiff_Dstar_D0->at(i));
      //h_DstarD0_mass->Fill(dt.M()-fit.VertexMass());
      TLorentzVector lep, tvdt, b0;
      tvdt.SetPtEtaPhiM(dt.Pt(),dt.Eta(),dt.Phi(),dt.M());
      for(unsigned int j = 0; j<pt_Dstar_lepton.size(); j++) {
        lep.SetPtEtaPhiM(pt_Dstar_lepton.at(j),eta_Dstar_lepton.at(j),phi_Dstar_lepton.at(j),mass_Dstar_lepton.at(j));
        b0 = tvdt+lep;
        if(b0.M() < 2.0 || b0.M() > 6.0) continue;
        cout<<"    lep: "<<pt_Dstar_lepton.at(j)<<endl;
        cout<<"    B0: "<<b0.Pt()<<"  "<<b0.M()<<endl;
        double charge = charge_Dstar_lepton.at(j)*charge_Dstar_D0_kaon->at(i);
        if(charge > 0) h_DstarD0_mass->Fill(mdiff_Dstar_D0->at(i));
        h_Dstar_DsPlus_mass->Fill(kkpi.M());
        h_Dstar_DsPlus_mass->Fill(kpik.M());
        h_Dstar_LambdacPlus_mass->Fill(kppi.M());
        h_Dstar_LambdacPlus_mass->Fill(kpip.M());
      }
    }
    // D+ Geometric Fitting - i, j and k correspond to kaon, pion1 and pion2 tracks
    for(unsigned int i=0; i<p4_DPlus->size(); i++) {
      cout<<"  reco d+ mass: "<<p4_DPlus->at(i).M()<<endl;
      TLorentzVector k, pi1, pi2, k1, k2, kkpi, kpik, p1, p2, kppi, kpip;
      k.SetPtEtaPhiM(p4_DPlus_kaon->at(i).Pt(),
                     p4_DPlus_kaon->at(i).Eta(),
                     p4_DPlus_kaon->at(i).Phi(),
                     p4_DPlus_kaon->at(i).M());
      pi1.SetPtEtaPhiM(p4_DPlus_pion1->at(i).Pt(),
                       p4_DPlus_pion1->at(i).Eta(),
                       p4_DPlus_pion1->at(i).Phi(),
                       p4_DPlus_pion1->at(i).M());
      pi2.SetPtEtaPhiM(p4_DPlus_pion2->at(i).Pt(),
                       p4_DPlus_pion2->at(i).Eta(),
                       p4_DPlus_pion2->at(i).Phi(),
                       p4_DPlus_pion2->at(i).M());
      k1.SetPtEtaPhiM(p4_DPlus_pion1->at(i).Pt(),
                      p4_DPlus_pion1->at(i).Eta(),
                      p4_DPlus_pion1->at(i).Phi(),mk);
      k2.SetPtEtaPhiM(p4_DPlus_pion2->at(i).Pt(),
                      p4_DPlus_pion2->at(i).Eta(),
                      p4_DPlus_pion2->at(i).Phi(),mk);
      p1.SetPtEtaPhiM(p4_DPlus_pion1->at(i).Pt(),
                      p4_DPlus_pion1->at(i).Eta(),
                      p4_DPlus_pion1->at(i).Phi(),mp);
      p2.SetPtEtaPhiM(p4_DPlus_pion2->at(i).Pt(),
                      p4_DPlus_pion2->at(i).Eta(),
                      p4_DPlus_pion2->at(i).Phi(),mp);
      kkpi = k+k1+pi2; // Ds+
      kpik = k+pi1+k2; // Ds+
      kppi = k+p1+pi2; // Lambdac+
      kpip = k+pi1+p2; // Lambdac+
      for(unsigned int j=0; j<trkpara_DPlus_kaon->at(i).size(); j++) {
        i_par[j] = { trkpara_DPlus_kaon->at(i).at(j) };
        j_par[j] = { trkpara_DPlus_pion1->at(i).at(j) };
        k_par[j] = { trkpara_DPlus_pion2->at(i).at(j) };
        const int index[5] = {0, 6, 12, 18, 24};
        i_para[j] = i_par[j] + _rnd->Gaus()*sqrt(trkcovm_DPlus_kaon->at(i).at(index[j]));
        j_para[j] = j_par[j] + _rnd->Gaus()*sqrt(trkcovm_DPlus_pion1->at(i).at(index[j]));
        k_para[j] = k_par[j] + _rnd->Gaus()*sqrt(trkcovm_DPlus_pion2->at(i).at(index[j]));
      }
      for(unsigned int j=0; j<trkcovm_DPlus_kaon->at(i).size(); j++) {
        i_covm[j] = { trkcovm_DPlus_kaon->at(i).at(j) };
        j_covm[j] = { trkcovm_DPlus_pion1->at(i).at(j) };
        k_covm[j] = { trkcovm_DPlus_pion2->at(i).at(j) };
      }
      CmsFit3D fit(3);
      fit.SetVerbose(0);
      fit.AddTrack(i_para,i_covm,mk);
      fit.AddTrack(j_para,j_covm,mpi);
      fit.AddTrack(k_para,k_covm,mpi);
      fit.Fit();
      double pv_lxy = pv_lxy_DPlus->at(i);
      double pt = fit.Pt();
      double lxy = fit.S();
      double decay_length = lxy - pv_lxy;
      double signif = decay_length/fit.ErrS();
      double proper_time = (decay_length/pt)*M_DPlus;
      cout<<"  fitted mass: "<<fit.VertexMass()<<endl;
      //h_DPlus_mass->Fill(p4_DPlus->at(i).M());
      //h_DPlus_mass->Fill(fit.VertexMass());
      TLorentzVector lep, dp, b0;
      dp.SetPtEtaPhiM(fit.Pt(),p4_DPlus->at(i).Eta(),p4_DPlus->at(i).Phi(),fit.VertexMass());
      for(unsigned int j = 0; j<pt_DPlus_lepton.size(); j++) {
        lep.SetPtEtaPhiM(pt_DPlus_lepton.at(j),eta_DPlus_lepton.at(j),phi_DPlus_lepton.at(j),mass_DPlus_lepton.at(j));
        b0 = dp+lep;
        if(b0.M() < 2.0 || b0.M() > 6.0) continue;
        cout<<"    lep: "<<pt_DPlus_lepton.at(j)<<"  "<<charge_DPlus_lepton.at(j)<<endl;
        cout<<"    B0: "<<b0.Pt()<<"  "<<b0.M()<<endl;
        double charge = charge_DPlus_lepton.at(j)*charge_DPlus_kaon->at(i);
        if(charge > 0) h_DPlus_mass->Fill(p4_DPlus->at(i).M());
        h_DPlus_DsPlus_mass->Fill(kkpi.M());
        h_DPlus_DsPlus_mass->Fill(kpik.M());
        h_DPlus_LambdacPlus_mass->Fill(kppi.M());
        h_DPlus_LambdacPlus_mass->Fill(kpip.M());
      }
    }
    // Ds+ Geometric Fitting - i, j and k correspond to kaon1, kaon2 and pion tracks
    for(unsigned int i=0; i<p4_DsPlus->size(); i++) {
      cout<<"  reco ds+ mass: "<<p4_DsPlus->at(i).M()<<endl;
      TLorentzVector k1, k2, pi, pi1, pi2, pikpi, kpipi, dplus, p1, p2, p1kpi, kp2pi;
      k1.SetPtEtaPhiM(p4_DsPlus_kaon1->at(i).Pt(),
                      p4_DsPlus_kaon1->at(i).Eta(),
                      p4_DsPlus_kaon1->at(i).Phi(),
                      p4_DsPlus_kaon1->at(i).M());
      k2.SetPtEtaPhiM(p4_DsPlus_kaon2->at(i).Pt(),
                      p4_DsPlus_kaon2->at(i).Eta(),
                      p4_DsPlus_kaon2->at(i).Phi(),
                      p4_DsPlus_kaon2->at(i).M());
      pi.SetPtEtaPhiM(p4_DsPlus_pion->at(i).Pt(),
                      p4_DsPlus_pion->at(i).Eta(),
                      p4_DsPlus_pion->at(i).Phi(),
                      p4_DsPlus_pion->at(i).M());
      pi1.SetPtEtaPhiM(p4_DsPlus_kaon1->at(i).Pt(),
                       p4_DsPlus_kaon1->at(i).Eta(),
                       p4_DsPlus_kaon1->at(i).Phi(),mpi);
      pi2.SetPtEtaPhiM(p4_DsPlus_kaon2->at(i).Pt(),
                       p4_DsPlus_kaon2->at(i).Eta(),
                       p4_DsPlus_kaon2->at(i).Phi(),mpi);
      p1.SetPtEtaPhiM(p4_DsPlus_kaon1->at(i).Pt(),
                      p4_DsPlus_kaon1->at(i).Eta(),
                      p4_DsPlus_kaon1->at(i).Phi(),mp);
      p2.SetPtEtaPhiM(p4_DsPlus_kaon2->at(i).Pt(),
                      p4_DsPlus_kaon2->at(i).Eta(),
                      p4_DsPlus_kaon2->at(i).Phi(),mp);
      pikpi = pi1+k2+pi; // Ds+ reflection into D+
      kpipi = k1+pi2+pi; // Ds+ reflection into D+
      dplus = pikpi+kpipi;
      TLorentzVector dsplus, kk;
      dsplus = k1+k2+pi;
      kk = k1+k2;
      p1kpi = p1+k2+pi; // Lambdac+
      kp2pi = k1+p2+pi; // Lambdac+
      for(unsigned int j=0; j<trkpara_DsPlus_pion->at(i).size(); j++) {
        i_par[j] = { trkpara_DsPlus_kaon1->at(i).at(j) };
        j_par[j] = { trkpara_DsPlus_kaon2->at(i).at(j) };
        k_par[j] = { trkpara_DsPlus_pion->at(i).at(j) };
        const int index[5] = {0, 6, 12, 18, 24};
        i_para[j] = i_par[j] + _rnd->Gaus()*sqrt(trkcovm_DsPlus_kaon1->at(i).at(index[j]));
        j_para[j] = j_par[j] + _rnd->Gaus()*sqrt(trkcovm_DsPlus_kaon2->at(i).at(index[j]));
        k_para[j] = k_par[j] + _rnd->Gaus()*sqrt(trkcovm_DsPlus_pion->at(i).at(index[j]));
      }
      for(unsigned int j=0; j<trkcovm_DsPlus_pion->at(i).size(); j++) {
        i_covm[j] = { trkcovm_DsPlus_kaon1->at(i).at(j) };
        j_covm[j] = { trkcovm_DsPlus_kaon2->at(i).at(j) };
        k_covm[j] = { trkcovm_DsPlus_pion->at(i).at(j) };
      }
      CmsFit3D fit(3);
      fit.SetVerbose(0);
      fit.AddTrack(i_para,i_covm,mk);
      fit.AddTrack(j_para,j_covm,mk);
      fit.AddTrack(k_para,k_covm,mpi);
      fit.Fit();
      double pv_lxy = pv_lxy_DsPlus->at(i);
      double pt = fit.Pt();
      double lxy = fit.S();
      double decay_length = lxy - pv_lxy;
      double signif = decay_length/fit.ErrS();
      double proper_time = (decay_length/pt)*M_DsPlus;
      cout<<"  fitted mass: "<<fit.VertexMass()<<endl;
      //h_DsPlus_mass->Fill(p4_DsPlus->at(i).M());
      //h_DsPlus_mass->Fill(fit.VertexMass());
      TLorentzVector lep, ds, bs;
      ds.SetPtEtaPhiM(fit.Pt(),p4_DsPlus->at(i).Eta(),p4_DsPlus->at(i).Phi(),fit.VertexMass());
      for(unsigned int j = 0; j<pt_DsPlus_lepton.size(); j++) {
        lep.SetPtEtaPhiM(pt_DsPlus_lepton.at(j),
                         eta_DsPlus_lepton.at(j),
                         phi_DsPlus_lepton.at(j),
                         mass_DsPlus_lepton.at(j));
        bs = ds+lep;
        if(bs.M() < 2.1 || bs.M() > 6.0) continue;
        cout<<"    lep: "<<pt_DsPlus_lepton.at(j)<<endl;
        cout<<"    Bs0: "<<bs.Pt()<<"  "<<bs.M()<<endl;
        double charge = charge_DsPlus_lepton.at(j)*charge_DsPlus_pion->at(i);
        if(charge < 0) h_DsPlus_mass->Fill(p4_DsPlus->at(i).M());
        h_DsPlus_DPlus_mass->Fill(pikpi.M());
        h_DsPlus_DPlus_mass->Fill(kpipi.M());
        h_DsPlus_LambdacPlus_mass->Fill(p1kpi.M());
        h_DsPlus_LambdacPlus_mass->Fill(kp2pi.M());
        if((pikpi.M() > 1.78 || kpipi.M() > 1.78) && (pikpi.M() < 1.95 || kpipi.M() < 1.95)) {
          if(abs(kk.M()-mphi) < 0.01) {
            cout<<"    Ds+ reflection into D+: "<<pikpi.M()<<"  "<<kpipi.M()<<endl;
            cout<<"    m(KK) - 1.019: "<<kk.M()-mphi<<endl;
            h_DsPlus_phipi_mass->Fill(dsplus.M());
          }
        }
      }
    }
    // Lambdac+ Geometric Fitting - i, j and k correspond to kaon, pion and proton tracks
    for(unsigned int i=0; i<p4_LambdacPlus->size(); i++) {
      cout<<"  reco lambdac+ mass: "<<p4_LambdacPlus->at(i).M()<<endl;
     TLorentzVector k, pi, p1, p2, kpipi, kpik;
      k.SetPtEtaPhiM(p4_LambdacPlus_kaon->at(i).Pt(),
                     p4_LambdacPlus_kaon->at(i).Eta(),
                     p4_LambdacPlus_kaon->at(i).Phi(),
                     p4_LambdacPlus_kaon->at(i).M());
      pi.SetPtEtaPhiM(p4_LambdacPlus_pion->at(i).Pt(),
                      p4_LambdacPlus_pion->at(i).Eta(),
                      p4_LambdacPlus_pion->at(i).Phi(),
                      p4_LambdacPlus_pion->at(i).M());
      p1.SetPtEtaPhiM(p4_LambdacPlus_proton->at(i).Pt(),
                      p4_LambdacPlus_proton->at(i).Eta(),
                      p4_LambdacPlus_proton->at(i).Phi(),mpi);
      p2.SetPtEtaPhiM(p4_LambdacPlus_proton->at(i).Pt(),
                      p4_LambdacPlus_proton->at(i).Eta(),
                      p4_LambdacPlus_proton->at(i).Phi(),mk);
      kpipi = k+pi+p1; // D+
      kpik = k+pi+p2; // Ds+
      for(unsigned int j=0; j<trkpara_LambdacPlus_kaon->at(i).size(); j++) {
        i_par[j] = { trkpara_LambdacPlus_kaon->at(i).at(j) };
        j_par[j] = { trkpara_LambdacPlus_pion->at(i).at(j) };
        k_par[j] = { trkpara_LambdacPlus_proton->at(i).at(j) };
        const int index[5] = {0, 6, 12, 18, 24};
        i_para[j] = i_par[j] + _rnd->Gaus()*sqrt(trkcovm_LambdacPlus_kaon->at(i).at(index[j]));
        j_para[j] = j_par[j] + _rnd->Gaus()*sqrt(trkcovm_LambdacPlus_pion->at(i).at(index[j]));
        k_para[j] = k_par[j] + _rnd->Gaus()*sqrt(trkcovm_LambdacPlus_proton->at(i).at(index[j]));
      }
      for(unsigned int j=0; j<trkcovm_LambdacPlus_kaon->at(i).size(); j++) {
        i_covm[j] = { trkcovm_LambdacPlus_kaon->at(i).at(j) };
        j_covm[j] = { trkcovm_LambdacPlus_pion->at(i).at(j) };
        k_covm[j] = { trkcovm_LambdacPlus_proton->at(i).at(j) };
      }
      CmsFit3D fit(3);
      fit.SetVerbose(0);
      fit.AddTrack(i_para,i_covm,mk);
      fit.AddTrack(j_para,j_covm,mpi);
      fit.AddTrack(k_par,k_covm,mp);
      fit.Fit();
      double pv_lxy = pv_lxy_LambdacPlus->at(i);
      double pt = fit.Pt();
      double lxy = fit.S();
      double decay_length = lxy - pv_lxy;
      double signif = decay_length/fit.ErrS();
      double proper_time = (decay_length/pt)*M_LambdacPlus;
      cout<<"  fitted mass: "<<fit.VertexMass()<<endl;
      //h_LambdacPlus_mass->Fill(p4_LambdacPlus->at(i).M());
      //h_LambdacPlus_mass->Fill(fit.VertexMass());
      TLorentzVector lep, lambdac, lambdab;
      lambdac.SetPtEtaPhiM(fit.Pt(),p4_LambdacPlus->at(i).Eta(),p4_LambdacPlus->at(i).Phi(),fit.VertexMass());
      for(unsigned int j = 0; j<pt_LambdacPlus_lepton.size(); j++) {
        lep.SetPtEtaPhiM(pt_LambdacPlus_lepton.at(j),
                         eta_LambdacPlus_lepton.at(j),
                         phi_LambdacPlus_lepton.at(j),
                         mass_LambdacPlus_lepton.at(j));
        lambdab = lambdac+lep;
        if(lambdab.M() < 2.4 || lambdab.M() > 6.0) continue;
        cout<<"    lep: "<<pt_LambdacPlus_lepton.at(j)<<endl;
        cout<<"    Lambda_b0: "<<lambdab.Pt()<<"  "<<lambdab.M()<<endl;
        double charge = charge_LambdacPlus_lepton.at(j)*charge_LambdacPlus_kaon->at(i);
        if(charge > 0) h_LambdacPlus_mass->Fill(p4_LambdacPlus->at(i).M());
        h_LambdacPlus_DPlus_mass->Fill(kpipi.M());
        h_LambdacPlus_DsPlus_mass->Fill(kpik.M());
      } 
    }
    cout<<endl;
    tree->Fill();
  }
  cout<<"gen: "<<gen<<"  gen b-lep: "<<genlep<<endl;
  cout<<"gen: "<<genbplus<<"  "<<genb0dstar<<"  "<<genb0<<"  "<<genbs0<<"  "<<genlambdab0<<endl;
  cout<<"di-lep: "<<dilep<<"  reco: "<<reco<<"  bjets: "<<bjets<<endl;
  cout<<"lep-bjet: "<<lep1<<"  "<<lep2<<endl;
  f1->Write();
  f1->Close();
}
