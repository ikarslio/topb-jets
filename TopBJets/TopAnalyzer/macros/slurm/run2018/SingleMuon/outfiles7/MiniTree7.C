#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TTree.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include "Math/LorentzVector.h"
#include "Math/Vector4D.h"
#include "TMatrixD.h"
#include "TRandom.h"
#define me 0.00051099895
#define mmu 0.1056583755

using namespace std;
using namespace ROOT;
typedef Math::LorentzVector<Math::PtEtaPhiM4D<double> > LV;
Double_t deltaPhi(Double_t phi1, Double_t phi2)
{
  Double_t pi = 3.1415927;
  Double_t dphi = fabs(phi1 - phi2);
  if(dphi >= pi) dphi = 2. * pi - dphi;
  return dphi;
}

void MiniTree7() {
  TFile *f1;
  TTree* T1;
  std::ifstream file("infile7.txt");
  std::string fname;
  std::string path = "root://eos.cms.rcac.purdue.edu//";
  //std::string base = "outfiles1/";
  while(!file.eof()) {
    file>>fname;
    if(file.eof()) break;
    if(fname.substr(fname.size() - 5,5) == ".root") {
      //cout<<fname.c_str()<<endl;
      std::size_t f = fname.find("store");
      std::string name = fname.substr(f);
      f1 = TFile::Open((path+name).c_str());
      T1 = (TTree*)f1->Get("demo/Ntuple");
      std::size_t opos = fname.find("data");
      std::string outfile = fname.substr(opos);
      std::size_t epos = fname.find("run2018");
      std::string ext = fname.substr(epos,8);
      //cout<<outfile<<endl;
      //cout<<ext<<endl;

      // Branch variables
      UInt_t event;
      double pfMET;
      double nJet;
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
      std::vector<double> * charge_DPlus_pion1 = 0;
      std::vector<double> * charge_DPlus_pion2 = 0;
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
      std::vector<double> * charge_D0_pion = 0;
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
      std::vector<double> * charge_Dstar_D0_pion = 0;
      std::vector<double> * charge_Dstar_pionsoft = 0;
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
      std::vector<double> * charge_DsPlus_kaon1 = 0;
      std::vector<double> * charge_DsPlus_kaon2 = 0;
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
      std::vector<double> * charge_LambdacPlus_pion = 0;
      std::vector<double> * charge_LambdacPlus_proton = 0;
      std::vector<std::vector<double>> * trkpara_LambdacPlus_kaon = 0;
      std::vector<std::vector<double>> * trkpara_LambdacPlus_pion = 0;
      std::vector<std::vector<double>> * trkpara_LambdacPlus_proton = 0;
      std::vector<std::vector<double>> * trkcovm_LambdacPlus_kaon = 0;
      std::vector<std::vector<double>> * trkcovm_LambdacPlus_pion = 0;
      std::vector<std::vector<double>> * trkcovm_LambdacPlus_proton = 0;
      std::vector<double> * pv_lxy_LambdacPlus = 0;
      T1->SetBranchAddress("event", &event);
      T1->SetBranchAddress("pfMET", &pfMET);
      T1->SetBranchAddress("nJet", &nJet);
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
      T1->SetBranchAddress("charge_DPlus_pion1", &charge_DPlus_pion1);
      T1->SetBranchAddress("charge_DPlus_pion2", &charge_DPlus_pion2);
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
      T1->SetBranchAddress("charge_D0_pion", &charge_D0_pion);
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
      T1->SetBranchAddress("charge_Dstar_D0_pion", &charge_Dstar_D0_pion);
      T1->SetBranchAddress("charge_Dstar_pionsoft", &charge_Dstar_pionsoft);
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
      T1->SetBranchAddress("charge_DsPlus_kaon1", &charge_DsPlus_kaon1);
      T1->SetBranchAddress("charge_DsPlus_kaon2", &charge_DsPlus_kaon2);
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
      T1->SetBranchAddress("charge_LambdacPlus_pion", &charge_LambdacPlus_pion);
      T1->SetBranchAddress("charge_LambdacPlus_proton", &charge_LambdacPlus_proton);
      T1->SetBranchAddress("trkpara_LambdacPlus_kaon", &trkpara_LambdacPlus_kaon);
      T1->SetBranchAddress("trkpara_LambdacPlus_pion", &trkpara_LambdacPlus_pion);
      T1->SetBranchAddress("trkpara_LambdacPlus_proton", &trkpara_LambdacPlus_proton);
      T1->SetBranchAddress("trkcovm_LambdacPlus_kaon", &trkcovm_LambdacPlus_kaon);
      T1->SetBranchAddress("trkcovm_LambdacPlus_pion", &trkcovm_LambdacPlus_pion);
      T1->SetBranchAddress("trkcovm_LambdacPlus_proton", &trkcovm_LambdacPlus_proton);
      T1->SetBranchAddress("pv_lxy_LambdacPlus", &pv_lxy_LambdacPlus);

      TFile* f2 = new TFile((ext+"_"+outfile).c_str(), "RECREATE");
      TTree *tree = new TTree("tree","minitree");
      // Branch variables
      std::vector<double> pt_ttbar_lepton1;
      std::vector<double> pt_ttbar_lepton2;
      std::vector<double> eta_ttbar_lepton1;
      std::vector<double> eta_ttbar_lepton2;
      std::vector<double> phi_ttbar_lepton1;
      std::vector<double> phi_ttbar_lepton2;
      tree->Branch("event", &event);
      tree->Branch("pt_ttbar_lepton1", &pt_ttbar_lepton1);
      tree->Branch("eta_ttbar_lepton1", &eta_ttbar_lepton1);
      tree->Branch("phi_ttbar_lepton1", &phi_ttbar_lepton1);
      tree->Branch("pt_ttbar_lepton2", &pt_ttbar_lepton2);
      tree->Branch("eta_ttbar_lepton2", &eta_ttbar_lepton2);
      tree->Branch("phi_ttbar_lepton2", &phi_ttbar_lepton2);
      tree->Branch("nJet", &nJet);
      tree->Branch("p4_jet", &p4_jet);
      tree->Branch("lepton_perJet", &lepton_perJet);
      tree->Branch("pdgId_lepton", &pdgId_lepton);
      tree->Branch("pt_lepton", &pt_lepton);
      tree->Branch("eta_lepton", &eta_lepton);
      tree->Branch("phi_lepton", &phi_lepton);
      tree->Branch("charge_lepton", &charge_lepton);
      tree->Branch("mass_lepton", &mass_lepton);
      tree->Branch("p4_DPlus", &p4_DPlus);
      tree->Branch("p4_DPlus_kaon", &p4_DPlus_kaon);
      tree->Branch("p4_DPlus_pion1", &p4_DPlus_pion1);
      tree->Branch("p4_DPlus_pion2", &p4_DPlus_pion2);
      tree->Branch("charge_DPlus_kaon", &charge_DPlus_kaon);
      tree->Branch("charge_DPlus_pion1", &charge_DPlus_pion1);
      tree->Branch("charge_DPlus_pion2", &charge_DPlus_pion2);
      tree->Branch("trkpara_DPlus_kaon", &trkpara_DPlus_kaon);
      tree->Branch("trkpara_DPlus_pion1", &trkpara_DPlus_pion1);
      tree->Branch("trkpara_DPlus_pion2", &trkpara_DPlus_pion2);
      tree->Branch("trkcovm_DPlus_kaon", &trkcovm_DPlus_kaon);
      tree->Branch("trkcovm_DPlus_pion1", &trkcovm_DPlus_pion1);
      tree->Branch("trkcovm_DPlus_pion2", &trkcovm_DPlus_pion2);
      tree->Branch("pv_lxy_DPlus", &pv_lxy_DPlus);
      tree->Branch("p4_D0", &p4_D0);
      tree->Branch("p4_D0_kaon", &p4_D0_kaon);
      tree->Branch("p4_D0_pion", &p4_D0_pion);
      tree->Branch("charge_D0_kaon", &charge_D0_kaon);
      tree->Branch("charge_D0_pion", &charge_D0_pion);
      tree->Branch("trkpara_D0_kaon", &trkpara_D0_kaon);
      tree->Branch("trkpara_D0_pion", &trkpara_D0_pion);
      tree->Branch("trkcovm_D0_kaon", &trkcovm_D0_kaon);
      tree->Branch("trkcovm_D0_pion", &trkcovm_D0_pion);
      tree->Branch("pv_lxy_D0", &pv_lxy_D0);
      tree->Branch("p4_Dstar_D0", &p4_Dstar_D0);
      tree->Branch("p4_Dstar_D0_kaon", &p4_Dstar_D0_kaon);
      tree->Branch("p4_Dstar_D0_pion", &p4_Dstar_D0_pion);
      tree->Branch("p4_Dstar_pionsoft", &p4_Dstar_pionsoft);
      tree->Branch("charge_Dstar_D0_kaon", &charge_Dstar_D0_kaon);
      tree->Branch("charge_Dstar_D0_pion", &charge_Dstar_D0_pion);
      tree->Branch("charge_Dstar_pionsoft", &charge_Dstar_pionsoft);
      tree->Branch("trkpara_Dstar_D0_kaon", &trkpara_Dstar_D0_kaon);
      tree->Branch("trkpara_Dstar_D0_pion", &trkpara_Dstar_D0_pion);
      tree->Branch("trkcovm_Dstar_D0_kaon", &trkcovm_Dstar_D0_kaon);
      tree->Branch("trkcovm_Dstar_D0_pion", &trkcovm_Dstar_D0_pion);
      tree->Branch("pv_lxy_Dstar_D0", &pv_lxy_Dstar_D0);
      tree->Branch("mdiff_Dstar_D0", &mdiff_Dstar_D0);
      tree->Branch("p4_DsPlus", &p4_DsPlus);
      tree->Branch("p4_DsPlus_kaon1", &p4_DsPlus_kaon1);
      tree->Branch("p4_DsPlus_kaon2", &p4_DsPlus_kaon2);
      tree->Branch("p4_DsPlus_pion", &p4_DsPlus_pion);
      tree->Branch("charge_DsPlus_kaon1", &charge_DsPlus_kaon1);
      tree->Branch("charge_DsPlus_kaon2", &charge_DsPlus_kaon2);
      tree->Branch("charge_DsPlus_pion", &charge_DsPlus_pion);
      tree->Branch("trkpara_DsPlus_kaon1", &trkpara_DsPlus_kaon1);
      tree->Branch("trkpara_DsPlus_kaon2", &trkpara_DsPlus_kaon2);
      tree->Branch("trkpara_DsPlus_pion", &trkpara_DsPlus_pion);
      tree->Branch("trkcovm_DsPlus_kaon1", &trkcovm_DsPlus_kaon1);
      tree->Branch("trkcovm_DsPlus_kaon2", &trkcovm_DsPlus_kaon2);
      tree->Branch("trkcovm_DsPlus_pion", &trkcovm_DsPlus_pion);
      tree->Branch("pv_lxy_DsPlus", &pv_lxy_DsPlus);
      tree->Branch("p4_LambdacPlus", &p4_LambdacPlus);
      tree->Branch("p4_LambdacPlus_kaon", &p4_LambdacPlus_kaon);
      tree->Branch("p4_LambdacPlus_pion", &p4_LambdacPlus_pion);
      tree->Branch("p4_LambdacPlus_proton", &p4_LambdacPlus_proton);
      tree->Branch("charge_LambdacPlus_kaon", &charge_LambdacPlus_kaon);
      tree->Branch("charge_LambdacPlus_pion", &charge_LambdacPlus_pion);
      tree->Branch("charge_LambdacPlus_proton", &charge_LambdacPlus_proton);
      tree->Branch("trkpara_LambdacPlus_kaon", &trkpara_LambdacPlus_kaon);
      tree->Branch("trkpara_LambdacPlus_pion", &trkpara_LambdacPlus_pion);
      tree->Branch("trkpara_LambdacPlus_proton", &trkpara_LambdacPlus_proton);
      tree->Branch("trkcovm_LambdacPlus_kaon", &trkcovm_LambdacPlus_kaon);
      tree->Branch("trkcovm_LambdacPlus_pion", &trkcovm_LambdacPlus_pion);
      tree->Branch("trkcovm_LambdacPlus_proton", &trkcovm_LambdacPlus_proton);
      tree->Branch("pv_lxy_LambdacPlus", &pv_lxy_LambdacPlus);
      unsigned int nentries = T1->GetEntries();
      cout<<"entries: "<<nentries<<endl;
      int j = 0;
      bool mumu, mue, ee, emu, zmm, zee;
      for(unsigned int jentry=0; jentry<nentries; jentry++) {
        T1->GetEntry(jentry);
        if(jentry%100000 == 0) cout<<"Events Processed :  "<<jentry<<endl;
        //cout<<"Event "<<event<<endl;

        pt_ttbar_lepton1.clear();
        eta_ttbar_lepton1.clear();
        phi_ttbar_lepton1.clear();
        pt_ttbar_lepton2.clear();
        eta_ttbar_lepton2.clear();
        phi_ttbar_lepton2.clear();
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
            //cout<<"  mu-mu: "<<p4_muon->at(indx[i]).Pt()<<"  "<<p4_muon->at(indx[j]).Pt()<<endl;
            pt_ttbar_lepton1.push_back(p4_muon->at(indx[i]).Pt());
            eta_ttbar_lepton1.push_back(p4_muon->at(indx[i]).Eta());
            phi_ttbar_lepton1.push_back(p4_muon->at(indx[i]).Phi());
            pt_ttbar_lepton2.push_back(p4_muon->at(indx[j]).Pt());
            eta_ttbar_lepton2.push_back(p4_muon->at(indx[j]).Eta());
            phi_ttbar_lepton2.push_back(p4_muon->at(indx[j]).Phi());
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
            //cout<<"  mu-e: "<<p4_muon->at(indx[i]).Pt()<<"  "<<p4_electron->at(indy[j]).Pt()<<endl;
            pt_ttbar_lepton1.push_back(p4_muon->at(indx[i]).Pt());
            eta_ttbar_lepton1.push_back(p4_muon->at(indx[i]).Eta());
            phi_ttbar_lepton1.push_back(p4_muon->at(indx[i]).Phi());
            pt_ttbar_lepton2.push_back(p4_electron->at(indy[j]).Pt());
            eta_ttbar_lepton2.push_back(p4_electron->at(indy[j]).Eta());
            phi_ttbar_lepton2.push_back(p4_electron->at(indy[j]).Phi());
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
            //cout<<"  e-e: "<<p4_electron->at(indy[i]).Pt()<<"  "<<p4_electron->at(indy[j]).Pt()<<endl;
            pt_ttbar_lepton1.push_back(p4_electron->at(indy[i]).Pt());
            eta_ttbar_lepton1.push_back(p4_electron->at(indy[i]).Eta());
            phi_ttbar_lepton1.push_back(p4_electron->at(indy[i]).Phi());
            pt_ttbar_lepton2.push_back(p4_electron->at(indy[j]).Pt());
            eta_ttbar_lepton2.push_back(p4_electron->at(indy[j]).Eta());
            phi_ttbar_lepton2.push_back(p4_electron->at(indy[j]).Phi());
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
            //cout<<"  e-mu: "<<p4_electron->at(indy[i]).Pt()<<"  "<<p4_muon->at(indx[j]).Pt()<<endl;
            pt_ttbar_lepton1.push_back(p4_electron->at(indy[i]).Pt());
            eta_ttbar_lepton1.push_back(p4_electron->at(indy[i]).Eta());
            phi_ttbar_lepton1.push_back(p4_electron->at(indy[i]).Phi());
            pt_ttbar_lepton2.push_back(p4_muon->at(indx[j]).Pt());
            eta_ttbar_lepton2.push_back(p4_muon->at(indx[j]).Eta());
            phi_ttbar_lepton2.push_back(p4_muon->at(indx[j]).Phi());
          }
        }
        //cout<<endl;
        if(!mumu && !mue && !ee && !emu) continue;
        if(zmm || zee) continue;
        if(pfMET < 40.) continue;
        tree->Fill();
        //cout<<endl;
      }
      f2->Write();
      f2->Close();
    }
    else {
      cout<<"** ERROR: Can't open "<<fname.c_str()<<" for input **"<<endl;
    }
  }
}
