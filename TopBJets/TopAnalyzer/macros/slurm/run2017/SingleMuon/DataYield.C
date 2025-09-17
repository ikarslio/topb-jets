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
#include "ConstrainedFit.hh"
#define me 0.00051099895
#define mmu 0.1056583755
#define mk 0.493677
#define mpi 0.13957039 
#define mp 0.9382720813
#define mphi 1.019461
#define mkstar 0.89167
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

double coshel(TLorentzVector particle, TLorentzVector parent, TLorentzVector grandparent) {
  TVector3 boosttoparent = -(parent.BoostVector());
  particle.Boost(boosttoparent);
  grandparent.Boost(boosttoparent);
  TVector3 particle3 = particle.Vect();
  TVector3 grandparent3 = grandparent.Vect();
  double numerator = particle3.Dot(grandparent3);
  double denominator = (particle3.Mag())*(grandparent3.Mag());
  double temp = numerator/denominator;
  return temp;
}

class DataFit : public ConstrainedFit {
public:
  DataFit() : ConstrainedFit(4,10,16) { }
  ~DataFit() { }
};

void DataYield() {
  cout<<"Calculating Data Yield"<<endl;
}

int main() {
  DataYield();
  TFile *f1;
  TTree* T1;
  std::ifstream file("outputFileList.txt");
  //std::ifstream file("outputFileList");
  std::string fname;
  //std::string path = "root://eos.cms.rcac.purdue.edu//";
  TFile* f2 = new TFile("smu_run2017.root", "RECREATE");
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
  TH1D *h_DPlus_pt = new TH1D("h_DPlus_pt","p_{T} D^{+}",100,0,100);
  TH1D *h_DPlus_lxy = new TH1D("h_DPlus_lxy","L_{xy}(K,#pi)",100,-1,1);
  TH1D *h_DPlus_lxysig = new TH1D("h_DPlus_lxysig","L_{xy}/#sigma_{L_{xy}}(K,#pi)",200,-10,40);
  TH1D *h_DPlus_ct = new TH1D("h_DPlus_ct","ct(K,#pi)",100,-0.05,0.2);
  TH1D *h_DPlus_chi2 = new TH1D("h_DPlus_chi2","#chi^{2}(K,#pi)",100,0,25);
  TH1D *h_DPlus_SameSign_mass = new TH1D("h_DPlus_SameSign_mass","m(K,#pi,#pi)",60,1.7,2.0);
  TH1D *h_DPlus_OppSign_mass = new TH1D("h_DPlus_OppSign_mass","m(K,#pi,#pi)",60,1.7,2.0);
  TH1D *h_D0_pt = new TH1D("h_D0_pt","p_{T} D^{0}",100,0,100);
  TH1D *h_D0_lxy = new TH1D("h_D0_lxy","L_{xy}(K,#pi)",100,-1,1);
  TH1D *h_D0_lxysig = new TH1D("h_D0_lxysig","L_{xy}/#sigma_{L_{xy}}(K,#pi)",200,-10,40);
  TH1D *h_D0_ct = new TH1D("h_D0_ct","ct(K,#pi)",100,-0.05,0.2);
  TH1D *h_D0_chi2 = new TH1D("h_D0_chi2","#chi^{2}(K,#pi)",100,0,25);
  TH1D *h_D0_SameSign_mass = new TH1D("h_D0_SameSign_mass","m(K,#pi)",50,1.73,1.98);
  TH1D *h_D0_OppSign_mass = new TH1D("h_D0_OppSign_mass","m(K,#pi)",50,1.73,1.98);
  TH1D *h_DstarD0_pt = new TH1D("h_DstarD0_pt","p_{T} D^{*}-D^{0}",100,0,100);
  TH1D *h_DstarD0_lxy = new TH1D("h_DstarD0_lxy","L_{xy}(K,#pi)",100,-1,1);
  TH1D *h_DstarD0_lxysig = new TH1D("h_DstarD0_lxysig","L_{xy}/#sigma_{L_{xy}}(K,#pi)",200,-10,40);
  TH1D *h_DstarD0_ct = new TH1D("h_DstarD0_ct","ct(K,#pi)",100,-0.05,0.2);
  TH1D *h_DstarD0_chi2 = new TH1D("h_DstarD0_chi2","#chi^{2}(K,#pi)",100,0,25);
  TH1D *h_DstarD0_SameSign_mass = new TH1D("h_DstarD0_SameSign_mass","m(D^{*})-m(D^{0})",70,0.135,0.17);
  TH1D *h_DstarD0_OppSign_mass = new TH1D("h_DstarD0_OppSign_mass","m(D^{*})-m(D^{0})",70,0.135,0.17);
  TH1D *h_DsPlus_pt = new TH1D("h_DsPlus_pt","p_{T} D_{s}^{+}",100,0,100);
  TH1D *h_DsPlus_lxy = new TH1D("h_DsPlus_lxy","L_{xy}(K,#pi)",100,-1,1);
  TH1D *h_DsPlus_lxysig = new TH1D("h_DsPlus_lxysig","L_{xy}/#sigma_{L_{xy}}(K,#pi)",200,-10,40);
  TH1D *h_DsPlus_ct = new TH1D("h_DsPlus_ct","ct(K,#pi)",100,-0.05,0.2);
  TH1D *h_DsPlus_chi2 = new TH1D("h_DsPlus_chi2","#chi^{2}(K,#pi)",100,0,25);
  TH1D *h_DsPlus_K2Phi_mass = new TH1D("h_DsPlus_K2Phi_mass","m(K^{+}K^{-} - #phi)",105,-0.05,1);
  TH1D *h_DsPlus_helicity = new TH1D("h_DsPlus_helicity","cos#psi",101,-1.01,1.01);
  TH1D *h_DsPlus_SameSign_mass = new TH1D("h_DsPlus_SameSign_mass","m(K,K,#pi)",60,1.8,2.1);
  TH1D *h_DsPlus_OppSign_mass = new TH1D("h_DsPlus_OppSign_mass","m(K,K,#pi)",60,1.8,2.1);
  TH1D *h_LambdacPlus_pt = new TH1D("h_LambdacPlus_pt","p_{T} #Lambda_{c}^{+}",100,0,100);
  TH1D *h_LambdacPlus_lxy = new TH1D("h_LambdacPlus_lxy","L_{xy}(K,#pi)",100,-1,1);
  TH1D *h_LambdacPlus_lxysig = new TH1D("h_LambdacPlus_lxysig","L_{xy}/#sigma_{L_{xy}}(K,#pi)",200,-10,40);
  TH1D *h_LambdacPlus_ct = new TH1D("h_LambdacPlus_ct","ct(K,#pi)",100,-0.05,0.2);
  TH1D *h_LambdacPlus_chi2 = new TH1D("h_LambdacPlus_chi2","#chi^{2}(K,#pi)",100,0,25);
  TH1D *h_LambdacPlus_angle = new TH1D("h_LambdacPlus_angle","cos#theta",67,-1,1.01);
  TH1D *h_LambdacPlus_Kpipi = new TH1D("h_LambdacPlus_Kpipi", "m(K,#pi,#pi)",60,1.7,2.0);
  TH1D *h_LambdacPlus_KKpi = new TH1D("h_LambdacPlus_KKpi", "m(K,K,#pi)",60,1.8,2.1);
  TH1D *h_LambdacPlus_Kpipi_Kpi = new TH1D("h_LambdacPlus_Kpipi_Kpi", "m(K,#pi,p)",70,0.135,0.17);
  TH1D *h_LambdacPlus_SameSign_mass = new TH1D("h_LambdacPlus_SameSign_mass","m(K,#pi,p)",60,2.1,2.4);
  TH1D *h_LambdacPlus_OppSign_mass = new TH1D("h_LambdacPlus_OppSign_mass","m(K,#pi,p)",60,2.1,2.4);
  while(!file.eof()) {
    file>>fname;
    if(file.eof()) break;
    if(fname.substr(fname.size() - 5,5) == ".root") {
      //cout<<fname.c_str()<<endl;
      //std::size_t f = fname.find("store");
      //std::string name = fname.substr(f);
      //f1 = TFile::Open((path+name).c_str());
      f1 = TFile::Open(fname.c_str());
      T1 = (TTree*)f1->Get("tree");
      //std::size_t pos = fname.find("data");
      //std::string outfile = fname.substr(pos);
      //cout<<outfile<<endl;

      // Branch variables
      UInt_t event;
      double nJet;
      std::vector<double> * pt_ttbar_lepton1 = 0;
      std::vector<double> * pt_ttbar_lepton2 = 0;
      std::vector<double> * eta_ttbar_lepton1 = 0;
      std::vector<double> * eta_ttbar_lepton2 = 0;
      std::vector<double> * phi_ttbar_lepton1 = 0;
      std::vector<double> * phi_ttbar_lepton2 = 0;
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
      T1->SetBranchAddress("pt_ttbar_lepton1", &pt_ttbar_lepton1);
      T1->SetBranchAddress("eta_ttbar_lepton1", &eta_ttbar_lepton1);
      T1->SetBranchAddress("phi_ttbar_lepton1", &phi_ttbar_lepton1);
      T1->SetBranchAddress("pt_ttbar_lepton2", &pt_ttbar_lepton2);
      T1->SetBranchAddress("eta_ttbar_lepton2", &eta_ttbar_lepton2);
      T1->SetBranchAddress("phi_ttbar_lepton2", &phi_ttbar_lepton2);
      T1->SetBranchAddress("nJet", &nJet);
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
      unsigned int nentries = T1->GetEntries();
      cout<<"entries: "<<nentries<<endl;

      std::vector<double> pt_jet, eta_jet;
      std::vector<std::vector<double> > newlepton_pdgId, newlepton_pt, newlepton_eta,
                                        newlepton_phi, newlepton_charge, newlepton_mass;
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

        pt_jet.clear();
        eta_jet.clear();
        newlepton_pdgId.clear();
        newlepton_pt.clear();
        newlepton_eta.clear();
        newlepton_phi.clear();
        newlepton_charge.clear();
        newlepton_mass.clear();
        cout<<"Event "<<event<<endl;
        int ij = 0;
        int foundlep = 0;
        for(unsigned int i=0; i<p4_jet->size(); i++) {
          cout<<"  jet: "<<p4_jet->at(i).Pt()<<endl;
          cout<<"  n-lep per jet: "<<lepton_perJet->at(i)<<endl;
          if(lepton_perJet->at(i) < 1) continue; 
          for(unsigned int j=0; j<pt_ttbar_lepton1->size(); j++) {
            cout<<"  i: "<<i<<"  j: "<<j<<endl;
            double deta1 = p4_jet->at(i).Eta() - eta_ttbar_lepton1->at(j);
            double dphi1 = deltaPhi(p4_jet->at(i).Phi(), phi_ttbar_lepton1->at(j));
            double dR1 = sqrt(deta1*deta1+dphi1*dphi1);
            double deta2 = p4_jet->at(i).Eta() - eta_ttbar_lepton2->at(j);
            double dphi2 = deltaPhi(p4_jet->at(i).Phi(), phi_ttbar_lepton2->at(j));
            double dR2 = sqrt(deta2*deta2+dphi2*dphi2);
            //cout<<"  lep-lep: "<<pt_ttbar_lepton1->at(j)<<"  "<<pt_ttbar_lepton2->at(j)<<endl;
            cout<<"    dR between lepton and jet: "<<dR1<<"  "<<dR2<<endl;
            if(dR1 > 0.4 && dR2 > 0.4) {
              ij += 1;
              cout<<"    ij: "<<ij<<endl; 
              pt_jet.push_back(p4_jet->at(i).Pt());
              eta_jet.push_back(p4_jet->at(i).Eta());
              std::vector<double> newlep_pdgId, newlep_pt, newlep_eta,
                                  newlep_phi, newlep_charge, newlep_mass;
              for(unsigned int k=0; k<pt_lepton->at(i).size(); k++) {
                foundlep += 1;
                cout<<"    lepton in b-jet: "<<pdgId_lepton->at(i).at(k)<<"  "<<pt_lepton->at(i).at(k)<<"  "<<eta_lepton->at(i).at(k)<<endl;
                newlep_pdgId.push_back(pdgId_lepton->at(i).at(k));
                newlep_pt.push_back(pt_lepton->at(i).at(k));
                newlep_eta.push_back(eta_lepton->at(i).at(k));
                newlep_phi.push_back(phi_lepton->at(i).at(k));
                newlep_charge.push_back(charge_lepton->at(i).at(k));
                newlep_mass.push_back(mass_lepton->at(i).at(k));
              }
              newlepton_pdgId.push_back(newlep_pdgId);
              newlepton_pt.push_back(newlep_pt);
              newlepton_eta.push_back(newlep_eta);
              newlepton_phi.push_back(newlep_phi);
              newlepton_charge.push_back(newlep_charge);
              newlepton_mass.push_back(newlep_mass);
            }
          }
        }
        if(ij < 1) continue;
        if(foundlep < 1) continue;
        cout<<"  jets: "<<ij<<endl;
        cout<<"  found lepton: "<<foundlep<<endl;
        for(unsigned int ijet = 0; ijet<pt_jet.size(); ijet++) {
          cout<<"  jet: "<<pt_jet.at(ijet)<<"  "<<eta_jet.at(ijet)<<endl;
          for(unsigned int ilep = 0; ilep<newlepton_pt.at(ijet).size(); ilep++) {
            TLorentzVector lep;
            lep.SetPtEtaPhiM(newlepton_pt.at(ijet).at(ilep),newlepton_eta.at(ijet).at(ilep),
                             newlepton_phi.at(ijet).at(ilep),newlepton_mass.at(ijet).at(ilep));
            // D0 Geometric Fitting - j and k correspond to kaon and pion tracks
            for(unsigned int i=0; i<p4_D0->size(); i++) {
              //cout<<"  reco d0 mass: "<<p4_D0->at(i).M()<<endl;
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
              fit.SetVerbose(0);
              fit.AddTrack(j_para,j_covm,mk);
              fit.AddTrack(k_para,k_covm,mpi);
              fit.Fit();
              TLorentzVector d0, bp;
              d0.SetPtEtaPhiM(fit.Pt(),p4_D0->at(i).Eta(),p4_D0->at(i).Phi(),fit.VertexMass());
              bp = d0+lep;
              double pv_lxy = pv_lxy_D0->at(i);
              double pt = fit.Pt();
              double lxy = fit.S();
              double decay_length = lxy - pv_lxy;
              double signif = decay_length/fit.ErrS();
              double proper_time = (decay_length/pt)*M_D0;
              double chi2 = fit.Chi2();
              h_D0_pt->Fill(pt);
              h_D0_lxy->Fill(decay_length);
              h_D0_lxysig->Fill(signif);
              h_D0_ct->Fill(proper_time);
              h_D0_chi2->Fill(fit.Chi2());
              fPt.push_back(pt);
              fLxy.push_back(decay_length);
              fLxysig.push_back(signif);
              fCt.push_back(proper_time);
              fChi2.push_back(fit.Chi2());
              fMass.push_back(p4_D0->at(i).M());
              fMass_lepD.push_back(bp.M());
              if(abs(decay_length) < 0.5) continue;
              if(abs(signif) < 10) continue;
              if(abs(proper_time) < 0.05) continue;
              if(chi2 > 1.) continue;
              if(bp.M() < 3. || bp.M() > 5.5) continue;
              double charge = newlepton_charge.at(ijet).at(ilep)*charge_D0_kaon->at(i);
              // lepton and kaon have same charge
              if(charge > 0) {
                //cout<<"  D0: "<<fit.Pt()<<"  "<<fit.VertexMass()<<"  "<<charge_D0_kaon->at(i)<<endl;
                //cout<<"  lepton: "<<newlepton_pdgId.at(ijet).at(ilep)<<"  "<<newlepton_pt.at(ijet).at(ilep)<<"  "<<newlepton_charge.at(ijet).at(ilep)<<endl;
                //cout<<"    B+: "<<bp.Pt()<<"  "<<bp.M()<<endl;
                h_D0_SameSign_mass->Fill(fit.VertexMass());
              }
              if(charge < 0) h_D0_OppSign_mass->Fill(fit.VertexMass());
            }
            // D0-D* Geometric Fitting - i and j correspond to kaon and pion tracks
            for(unsigned int i=0; i<p4_Dstar_D0->size(); i++) {
              //cout<<"  reco d*-d0 mass: "<<p4_Dstar_D0->at(i).M()<<endl;
              //if(p4_Dstar_D0->at(i).Pt() < 5.) continue;
              //if(p4_Dstar_pionsoft->at(i).Pt() < 1.) continue;
              for(unsigned int j=0; j<trkpara_Dstar_D0_kaon->at(i).size(); j++) {
                if(p4_Dstar_pionsoft->at(i).Pt() < 0.8) continue;
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
              TLorentzVector lep, tvdt, b0;
              tvdt.SetPtEtaPhiM(dt.Pt(),dt.Eta(),dt.Phi(),dt.M());
              b0 = tvdt+lep;
              double pv_lxy = pv_lxy_Dstar_D0->at(i);
              double pt = fit.Pt();
              double lxy = fit.S();
              double decay_length = lxy - pv_lxy;
              double signif = decay_length/fit.ErrS();
              double proper_time = (decay_length/pt)*M_D0;
              double chi2 = fit.Chi2();
              h_DstarD0_pt->Fill(pt);
              h_DstarD0_lxy->Fill(decay_length);
              h_DstarD0_lxysig->Fill(signif);
              h_DstarD0_ct->Fill(proper_time);
              h_DstarD0_chi2->Fill(fit.Chi2());
              if(abs(decay_length) < 0.6) continue;
              if(abs(signif) < 10) continue;
              if(abs(proper_time) < 0.05) continue;
              if(chi2 > 4.) continue;
              if(fit.VertexMass() < 1.8 || fit.VertexMass() > 1.95) continue;
              if(b0.M() < 2. || b0.M() > 5.5) continue;
              //double charge = newlepton_charge.at(ijet).at(ilep)*charge_Dstar_D0_kaon->at(i);
              double charge = charge_Dstar_D0_pion->at(i)*charge_Dstar_pionsoft->at(i);
              // lepton and kaon have same charge
              if(charge > 0) {
                //cout<<"  D0-D*: "<<fit.Pt()<<"  "<<dt.M()-fit.VertexMass()<<"  "<<charge_Dstar_D0_kaon->at(i)<<endl;
                //cout<<"  lepton: "<<newlepton_pdgId.at(ijet).at(ilep)<<"  "<<newlepton_pt.at(ijet).at(ilep)<<"  "<<newlepton_charge.at(ijet).at(ilep)<<endl;
                //cout<<"    B0-D*: "<<b0.Pt()<<"  "<<b0.M()<<endl;
                h_DstarD0_SameSign_mass->Fill(dt.M()-fit.VertexMass());
              }
              if(charge < 0) h_DstarD0_OppSign_mass->Fill(dt.M()-fit.VertexMass());
            }
            // D+ Geometric Fitting - i, j and k correspond to kaon, pion1 and pion2 tracks
            for(unsigned int i=0; i<p4_DPlus->size(); i++) {
              //cout<<"  reco d+ mass: "<<p4_DPlus->at(i).M()<<endl;
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
              TLorentzVector dp, b0;
              dp.SetPtEtaPhiM(fit.Pt(),p4_DPlus->at(i).Eta(),p4_DPlus->at(i).Phi(),fit.VertexMass());
              b0 = dp+lep;
              double pv_lxy = pv_lxy_DPlus->at(i);
              double pt = fit.Pt();
              double lxy = fit.S();
              double decay_length = lxy - pv_lxy;
              double signif = decay_length/fit.ErrS();
              double proper_time = (decay_length/pt)*M_DPlus;
              double chi2 = fit.Chi2();
              h_DPlus_pt->Fill(pt);
              h_DPlus_lxy->Fill(decay_length);
              h_DPlus_lxysig->Fill(signif);
              h_DPlus_ct->Fill(proper_time);
              h_DPlus_chi2->Fill(fit.Chi2());
              if(abs(decay_length) < 0.5) continue;
              if(abs(signif) < 10) continue;
              if(abs(proper_time) < 0.05) continue;
              if(chi2 > 10) continue;
              if(b0.M() < 3. || b0.M() > 5.) continue;
              double charge = newlepton_charge.at(ijet).at(ilep)*charge_DPlus_kaon->at(i);
              // lepton and kaon have same charge
              if(charge > 0) {
                //cout<<"  D+: "<<fit.Pt()<<"  "<<fit.VertexMass()<<"  "<<charge_DPlus_kaon->at(i)<<endl;
                //cout<<"  lepton: "<<newlepton_pdgId.at(ijet).at(ilep)<<"  "<<newlepton_pt.at(ijet).at(ilep)<<"  "<<newlepton_charge.at(ijet).at(ilep)<<endl;
                //cout<<"    B0: "<<b0.Pt()<<"  "<<b0.M()<<endl;
                h_DPlus_SameSign_mass->Fill(fit.VertexMass());
              }
              if(charge < 0) h_DPlus_OppSign_mass->Fill(fit.VertexMass());
            }
            // Ds+ Geometric Fitting - i, j and k correspond to kaon1, kaon2 and pion tracks
            for(unsigned int i=0; i<p4_DsPlus->size(); i++) {
              //cout<<"  reco ds+ mass: "<<p4_DsPlus->at(i).M()<<endl;
              //if(p4_DsPlus->at(i).Pt() < 5.) continue;
              TLorentzVector k1, k2, kk;
              k1.SetPtEtaPhiM(p4_DsPlus_kaon1->at(i).Pt(),
                              p4_DsPlus_kaon1->at(i).Eta(),
                              p4_DsPlus_kaon1->at(i).Phi(),
                              p4_DsPlus_kaon1->at(i).M());
              k2.SetPtEtaPhiM(p4_DsPlus_kaon2->at(i).Pt(),
                              p4_DsPlus_kaon2->at(i).Eta(),
                              p4_DsPlus_kaon2->at(i).Phi(),
                              p4_DsPlus_kaon2->at(i).M());
              kk = k1+k2;
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
              TLorentzVector ds, bs;
              ds.SetPtEtaPhiM(fit.Pt(),p4_DsPlus->at(i).Eta(),p4_DsPlus->at(i).Phi(),fit.VertexMass());
              bs = ds+lep;
              double pv_lxy = pv_lxy_DsPlus->at(i);
              double pt = fit.Pt();
              double lxy = fit.S();
              double decay_length = lxy - pv_lxy;
              double signif = decay_length/fit.ErrS();
              double proper_time = (decay_length/pt)*M_DsPlus;
              double chi2 = fit.Chi2();
              double hel = coshel(k1,kk,ds);
              h_DsPlus_pt->Fill(pt);
              h_DsPlus_lxy->Fill(decay_length);
              h_DsPlus_lxysig->Fill(signif);
              h_DsPlus_ct->Fill(proper_time);
              h_DsPlus_chi2->Fill(fit.Chi2());
              h_DsPlus_K2Phi_mass->Fill(kk.M()-mphi);
              h_DsPlus_helicity->Fill(hel);
              if(p4_DsPlus_kaon1->at(i).Pt() < 2) continue;
              if(p4_DsPlus_kaon2->at(i).Pt() < 2) continue;
              if(p4_DsPlus_pion->at(i).Pt() < 2) continue;
              if(abs(decay_length) < 0.2) continue;
              if(abs(signif) < 5) continue;
              if(abs(proper_time) < 0.02) continue;
              if(chi2 > 15) continue;
              if(abs(kk.M()-mphi) > 0.01) continue;
              if(abs(hel) < 0.4) continue;
              if(bs.M() < 2.5 || bs.M() > 5.5) continue;
              double charge = newlepton_charge.at(ijet).at(ilep)*charge_DsPlus_pion->at(i);
              // lepton and pion have opposite charge
              if(charge < 0) {
                //cout<<"  Ds+: "<<fit.Pt()<<"  mass: "<<fit.VertexMass()<<"  "<<charge_DsPlus_pion->at(i)<<endl;
                //cout<<"  lepton: "<<newlepton_pdgId.at(ijet).at(ilep)<<"  "<<newlepton_pt.at(ijet).at(ilep)<<"  "<<newlepton_charge.at(ijet).at(ilep)<<endl;
                //cout<<"    Bs0: "<<bs.Pt()<<"  "<<bs.M()<<endl;
                h_DsPlus_SameSign_mass->Fill(fit.VertexMass());
              }
              if(charge > 0) h_DsPlus_OppSign_mass->Fill(fit.VertexMass());
            }
            // Lambdac+ Geometric Fitting - i, j and k correspond to kaon, pion and proton tracks
            for(unsigned int i=0; i<p4_LambdacPlus->size(); i++) {
              //cout<<"  reco lambdac+ mass: "<<p4_LambdacPlus->at(i).M()<<endl;
              TLorentzVector k, pi, m1p, m2p, kpipi, kkpi, kpi;
              k.SetPtEtaPhiM(p4_LambdacPlus_kaon->at(i).Pt(),
                             p4_LambdacPlus_kaon->at(i).Eta(),
                             p4_LambdacPlus_kaon->at(i).Phi(),
                             p4_LambdacPlus_kaon->at(i).M());
              pi.SetPtEtaPhiM(p4_LambdacPlus_pion->at(i).Pt(),
                              p4_LambdacPlus_pion->at(i).Eta(),
                              p4_LambdacPlus_pion->at(i).Phi(),
                              p4_LambdacPlus_pion->at(i).M());
              m1p.SetPtEtaPhiM(p4_LambdacPlus_proton->at(i).Pt(),
                               p4_LambdacPlus_proton->at(i).Eta(),
                               p4_LambdacPlus_proton->at(i).Phi(),mpi);
              m2p.SetPtEtaPhiM(p4_LambdacPlus_proton->at(i).Pt(),
                               p4_LambdacPlus_proton->at(i).Eta(),
                               p4_LambdacPlus_proton->at(i).Phi(),mk);
              kpipi = k+pi+m1p;
              kkpi = k+pi+m2p;
              kpi = k+pi;
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
              TLorentzVector lambdac, lambdab;
              lambdac.SetPtEtaPhiM(fit.Pt(),p4_LambdacPlus->at(i).Eta(),p4_LambdacPlus->at(i).Phi(),fit.VertexMass());
              lambdab = lambdac+lep;
              double pv_lxy = pv_lxy_LambdacPlus->at(i);
              double pt = fit.Pt();
              double lxy = fit.S();
              double decay_length = lxy - pv_lxy;
              double signif = decay_length/fit.ErrS();
              double proper_time = (decay_length/pt)*M_LambdacPlus;
              double chi2 = fit.Chi2();
              TVector3 lepton = lep.Vect();
              TVector3 lambda = lambdac.Vect();
              double numerator = lambda.Dot(lepton);
              double denominator = (lambda.Mag())*(lepton.Mag());
              double angle = numerator/denominator;
              h_LambdacPlus_pt->Fill(pt);
              h_LambdacPlus_lxy->Fill(decay_length);
              h_LambdacPlus_lxysig->Fill(signif);
              h_LambdacPlus_ct->Fill(proper_time);
              h_LambdacPlus_chi2->Fill(fit.Chi2());
              h_LambdacPlus_angle->Fill(angle);
              if(p4_LambdacPlus_kaon->at(i).Pt() < 2) continue;
              if(p4_LambdacPlus_proton->at(i).Pt() < 4) continue;
              if(abs(decay_length) < 0.8) continue;
              if(abs(signif) < 12) continue;
              if(abs(proper_time) < 0.08) continue;
              if(chi2 > 5) continue;
              //if(angle < 0.707) continue; //angle < 45deg;
              cout<<"  Kpipi: "<<kpipi.M()<<"  KKpi: "<<kkpi.M()<<"  Kpipi-Kpi: "<<kpipi.M()<<" - "<<kpi.M()<<endl;
              h_LambdacPlus_Kpipi->Fill(kpipi.M());
              h_LambdacPlus_KKpi->Fill(kkpi.M());
              h_LambdacPlus_Kpipi_Kpi->Fill(kpipi.M()-kpi.M());
              if(abs(kpipi.M()-M_DPlus) < 0.020) continue;
              if(abs(kkpi.M()-M_DsPlus) < 0.020) continue;
              if((kpipi.M()-kpi.M()) > 0.143 && (kpipi.M()-kpi.M()) < 0.149) continue;
              //if(pt < 5) continue;
              if(lambdab.Pt() < 10) continue;
              if(lambdab.M() < 4. || lambdab.M() > 5.5) continue;
              double charge = newlepton_charge.at(ijet).at(ilep)*charge_LambdacPlus_kaon->at(i);
              // lepton and kaon have same charge
              if(charge > 0) {
                cout<<"  Lambdac+: "<<fit.Pt()<<"  "<<fit.VertexMass()<<"  "<<charge_LambdacPlus_kaon->at(i)<<endl;
                cout<<"  lepton: "<<newlepton_pdgId.at(ijet).at(ilep)<<"  "<<newlepton_pt.at(ijet).at(ilep)<<"  "<<newlepton_charge.at(ijet).at(ilep)<<endl;
                cout<<"    Lambda_b0: "<<lambdab.Pt()<<"  "<<lambdab.M()<<endl;
                h_LambdacPlus_SameSign_mass->Fill(fit.VertexMass());
              }
              if(charge < 0) h_LambdacPlus_OppSign_mass->Fill(fit.VertexMass());
            }
          }
        }
        //cout<<endl;
        tree->Fill();
      }
    }
  }
  f2->Write();
  f2->Close();
}
