#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
//#include <gsl/gsl_linalg.h>

// Helper functions, forward declared
std::string ftoa1(double i); std::string ftoa2(double i); // Converts double to string
void beautifyHist(TH1D *h);
void beautifyFunc(TF1 *f);

Double_t background(Double_t *x, Double_t *par);
Double_t doubleGaussianPeak(Double_t *x, Double_t *par);
Double_t fitFunction(Double_t *x, Double_t *par);

void MassFitsSignal() {

  TFile* f1 = TFile::Open("ul2018.root");

  TH1D *h1_mass = (TH1D*)f1->Get("h_D0_mass");
  TH1D *h2_mass = (TH1D*)f1->Get("h_DstarD0_mass");
  TH1D *h3_mass = (TH1D*)f1->Get("h_DPlus_mass");
  TH1D *h4_mass = (TH1D*)f1->Get("h_DsPlus_mass");
  TH1D *h5_mass = (TH1D*)f1->Get("h_LambdacPlus_mass");

  // ***************************************************************** c_MassFit_D0 ***************************************************************** //

  double binWidth1 = h1_mass->GetBinWidth(1);
  std::cout<<"Bin Width = "<<binWidth1<<std::endl;

  h1_mass->SetTitle(("; M(D^{0} #rightarrow K^{-} #pi^{+}) GeV; Entries / "+ftoa1(binWidth1*1000)+" MeV").c_str());
  TCanvas *c_MassFit_D0 = new TCanvas("c_MassFit_D0", "",125,111,700,600);
  TF1 *f_function_D0 = new TF1("f_function_D0", "fitFunction", 1.73, 1.98, 8);
  f_function_D0->SetParameters(10, 0, 300, 1.865, 0.005, 300, 1.875, 0.005);
  h1_mass->Fit("f_function_D0", "L", "ep");
  Double_t parameters_1[8];
  f_function_D0->GetParameters(parameters_1);

  TF1 *f_background_D0 = new TF1("f_background_D0", background, 1.73, 1.98, 2);
 f_background_D0->SetParameters(parameters_1);
  TF1 *f_signal_D0 = new TF1("f_signal_D0", doubleGaussianPeak, 1.73, 1.98, 6);
  f_signal_D0->SetParameters(parameters_1+2);

  // Beautify everything before putting on canvas
  beautifyHist(h1_mass);
  //h1_mass->SetLineColor(kBlue);
  beautifyFunc(f_function_D0);
  beautifyFunc(f_background_D0);
  f_function_D0->SetLineColor(kBlue);
  f_background_D0->SetLineColor(kRed);

  h1_mass->SetMaximum(8000);
  h1_mass->Draw("HIST pe");
  f_function_D0->Draw("same");
  f_background_D0->Draw("same");
  TLegend *legend_1=new TLegend(0.14,0.69,0.44,0.86);
  legend_1->SetLineColor(0);
  legend_1->SetTextSize(0.045);
  legend_1->AddEntry(h1_mass,"D^{0}","lpe");
  legend_1->AddEntry(f_function_D0,"Signal fit","l");
  legend_1->AddEntry(f_background_D0,"Background fit","l");
  legend_1->Draw();

  cout<<"Integral = "<<f_signal_D0->Integral(1.73, 1.98)<<endl;
  cout<<"Integral = "<<f_signal_D0->Integral(1.8, 1.9)<<endl;

  c_MassFit_D0->SaveAs("Plots/Signal/c_MassFit_D0.pdf");

  // *************************************************************** c_MassFit_D*-D0 ***************************************************************** //

  double binWidth2 = h2_mass->GetBinWidth(1);
  std::cout<<"Bin Width = "<<binWidth2<<std::endl;

  h2_mass->SetTitle(("; M(D^{*+}) - M(D^{0}) GeV; Entries / "+ftoa2(binWidth2*1000)+" MeV").c_str());
  TCanvas *c_MassFit_Dstar = new TCanvas("c_MassFit_Dstar", "",125,111,700,600);
  TF1 *f_function_Dstar = new TF1("f_function_Dstar", "fitFunction", 0.14, 0.17,8);
  f_function_Dstar->SetParameters(10, 1, 6000, 0.145, 0.001, 6000, 0.146, 0.001);
  h2_mass->Fit("f_function_Dstar", "L", "ep");
  Double_t parameters_2[8];
  f_function_Dstar->GetParameters(parameters_2);

  TF1 *f_background_Dstar = new TF1("f_background_Dstar", background, 0.14, 0.17, 2);
  f_background_Dstar->SetParameters(parameters_2);
  TF1 *f_signal_Dstar = new TF1("f_signal_Dstar", doubleGaussianPeak, 0.135, 0.17, 6);
  f_signal_Dstar->SetParameters(parameters_2+2);

  // Beautify everything before putting on canvas
  beautifyHist(h2_mass);
  //h2_mass->SetLineColor(kBlue);
  beautifyFunc(f_function_Dstar);
  beautifyFunc(f_background_Dstar);
  f_function_Dstar->SetLineColor(kBlue);
  f_background_Dstar->SetLineColor(kRed);

  h2_mass->SetMaximum(9000);
  h2_mass->Draw("HIST pe");
  f_function_Dstar->Draw("same");
  f_background_Dstar->Draw("same");
  TLegend *legend_2=new TLegend(0.14,0.69,0.44,0.86);
  legend_2->SetLineColor(0);
  legend_2->SetTextSize(0.045);
  legend_2->AddEntry(h2_mass,"#DeltaM(D^{*+}- D^{0})","lpe");
  legend_2->AddEntry(f_function_Dstar,"Signal fit","l");
  legend_2->AddEntry(f_background_Dstar,"Background fit","l");
  legend_2->Draw();

  cout<<"Integral = "<<f_signal_Dstar->Integral(0.135, 0.17)<<endl;

  c_MassFit_Dstar->SaveAs("Plots/Signal/c_MassFit_Dstar.pdf");

  // *************************************************************** c_MassFit_DPlus ***************************************************************** //

  double binWidth3 = h3_mass->GetBinWidth(1);
  std::cout<<"Bin Width = "<<binWidth1<<std::endl;

  h3_mass->SetTitle(("; M(D^{+} #rightarrow K^{-} #pi^{+} #pi^{+}) GeV; Entries / "+ftoa1(binWidth3*1000)+" MeV").c_str());
  TCanvas *c_MassFit_DPlus = new TCanvas("c_MassFit_DPlus", "",125,111,700,600);
  TF1 *f_function_DPlus = new TF1("f_function_DPlus", "fitFunction", 1.73, 1.98, 8);
  f_function_DPlus->SetParameters(1, 0, 1200, 1.865, 0.005, 1200, 1.875, 0.005);
  h3_mass->Fit("f_function_DPlus", "L", "ep");
  Double_t parameters_3[8];
  f_function_DPlus->GetParameters(parameters_3);

  TF1 *f_background_DPlus = new TF1("f_background_DPlus", background, 1.73, 1.98, 2);
 f_background_DPlus->SetParameters(parameters_3);
  TF1 *f_signal_DPlus = new TF1("f_signal_DPlus", doubleGaussianPeak, 1.73, 1.98, 6);
  f_signal_DPlus->SetParameters(parameters_3+2);

  // Beautify everything before putting on canvas
  beautifyHist(h3_mass);
  //h3_mass->SetLineColor(kBlue);
  beautifyFunc(f_function_DPlus);
  beautifyFunc(f_background_DPlus);
  f_function_DPlus->SetLineColor(kBlue);
  f_background_DPlus->SetLineColor(kRed);

  h3_mass->SetMaximum(1800);
  h3_mass->Draw("HIST pe");
  f_function_DPlus->Draw("same");
  f_background_DPlus->Draw("same");
  TLegend *legend_3=new TLegend(0.14,0.69,0.44,0.86);
  legend_3->SetLineColor(0);
  legend_3->SetTextSize(0.045);
  legend_3->AddEntry(h3_mass,"D^{+}","lpe");
  legend_3->AddEntry(f_function_DPlus,"Signal fit","l");
  legend_3->AddEntry(f_background_DPlus,"Background fit","l");
  legend_3->Draw();

  cout<<"Integral = "<<f_signal_DPlus->Integral(1.73, 1.98)<<endl;

  c_MassFit_DPlus->SaveAs("Plots/Signal/c_MassFit_DPlus.pdf");

  // *************************************************************** c_MassFit_DsPlus **************************************************************** //

  double binWidth4 = h4_mass->GetBinWidth(1);
  std::cout<<"Bin Width = "<<binWidth4<<std::endl;

  h4_mass->SetTitle(("; M(D_{s}^{+} #rightarrow K^{-} K^{+} #pi^{+}) GeV; Entries / "+ftoa1(binWidth4*1000)+" MeV").c_str());
  TCanvas *c_MassFit_DsPlus = new TCanvas("c_MassFit_DsPlus", "",125,111,700,600);
  TF1 *f_function_DsPlus = new TF1("f_function_DsPlus", "fitFunction", 1.85, 2.1, 8);
  f_function_DsPlus->SetParameters(1, 0, 1200, 1.96, 0.005, 1200, 1.97, 0.005);
  h4_mass->Fit("f_function_DsPlus", "L", "ep");
  Double_t parameters_4[8];
  f_function_DsPlus->GetParameters(parameters_4);

  TF1 *f_background_DsPlus = new TF1("f_background_DsPlus", background, 1.85, 2.1, 2);
  f_background_DsPlus->SetParameters(parameters_4);
  TF1 *f_signal_DsPlus = new TF1("f_signal_DsPlus", doubleGaussianPeak, 1.85, 2.1, 6);
  f_signal_DsPlus->SetParameters(parameters_4+2);

  // Beautify everything before putting on canvas
  beautifyHist(h4_mass);
  //h4_mass->SetLineColor(kBlue);
  beautifyFunc(f_function_DsPlus);
  beautifyFunc(f_background_DsPlus);
  f_function_DsPlus->SetLineColor(kBlue);
  f_background_DsPlus->SetLineColor(kRed);

  h4_mass->SetMaximum(2500);
  h4_mass->Draw("HIST pe");
  f_function_DsPlus->Draw("same");
  f_background_DsPlus->Draw("same");
  TLegend *legend_4=new TLegend(0.14,0.69,0.44,0.86);
  legend_4->SetLineColor(0);
  legend_4->SetTextSize(0.045);
  legend_4->AddEntry(h4_mass,"D_{s}^{+}","lpe");
  legend_4->AddEntry(f_function_DsPlus,"Signal fit","l");
  legend_4->AddEntry(f_background_DsPlus,"Background fit","l");
  legend_4->Draw();

  cout<<"Integral = "<<f_signal_DsPlus->Integral(1.85, 2.1)<<endl;

  c_MassFit_DsPlus->SaveAs("Plots/Signal/c_MassFit_DsPlus.pdf");

  // *********************************************************** c_MassFit_LambdacPlus ************************************************************** //

  double binWidth5 = h5_mass->GetBinWidth(1);
  std::cout<<"Bin Width = "<<binWidth5<<std::endl;

  h5_mass->SetTitle(("; M(#Lambda_{c}^{+} #rightarrow p K^{-} #pi^{+}) GeV; Entries / "+ftoa1(binWidth1*1000)+" MeV").c_str());
  TCanvas *c_MassFit_LambdacPlus = new TCanvas("c_MassFit_LambdacPlus", "",125,111,700,600);
  TF1 *f_function_LambdacPlus = new TF1("f_function_LambdacPlus", "fitFunction", 2.15, 2.4, 8);
  f_function_LambdacPlus->SetParameters(20, 1, 600, 2.285, 0.005, 500, 2.285, 0.005);
  h5_mass->Fit("f_function_LambdacPlus", "L", "ep");
  Double_t parameters_5[8];
  f_function_LambdacPlus->GetParameters(parameters_5);

  TF1 *f_background_LambdacPlus = new TF1("f_background_LambdacPlus", background, 2.15, 2.4, 2);
  f_background_LambdacPlus->SetParameters(parameters_5);
  TF1 *f_signal_LambdacPlus = new TF1("f_signal_LambdacPlus", doubleGaussianPeak, 2.15, 2.4, 6);
  f_signal_LambdacPlus->SetParameters(parameters_5+2);

  // Beautify everything before putting on canvas
  beautifyHist(h5_mass);
  //h5_mass->SetLineColor(kBlue);
  beautifyFunc(f_function_LambdacPlus);
  beautifyFunc(f_background_LambdacPlus);
  f_function_LambdacPlus->SetLineColor(kBlue);
  f_background_LambdacPlus->SetLineColor(kRed);

  h5_mass->SetMaximum(750);
  h5_mass->Draw("HIST pe");
  f_function_LambdacPlus->Draw("same");
  f_background_LambdacPlus->Draw("same");
  TLegend *legend_5=new TLegend(0.14,0.69,0.44,0.86);
  legend_5->SetLineColor(0);
  legend_5->SetTextSize(0.045);
  legend_5->AddEntry(h5_mass,"#Lambda_{c}^{+}","lpe");
  legend_5->AddEntry(f_function_LambdacPlus,"Signal fit","l");
  legend_5->AddEntry(f_background_LambdacPlus,"Background fit","l");
  legend_5->Draw();

  cout<<"Integral = "<<f_signal_LambdacPlus->Integral(2.15, 2.4)<<endl;

  c_MassFit_LambdacPlus->SaveAs("Plots/Signal/c_MassFit_LambdacPlus.pdf");
}

Double_t background(Double_t *x, Double_t *par) {
   return par[0] + par[1]*x[0];
}

Double_t doubleGaussianPeak(Double_t *x, Double_t *par) {
  return par[0]*TMath::Gaus(x[0],par[1],par[2]) + par[3]*TMath::Gaus(x[0],par[4],par[5]);
}

// Sum of background and peak function
Double_t fitFunction(Double_t *x, Double_t *par) {
  return background(x,par) + doubleGaussianPeak(x,&par[2]);
}

// Beautification functions
void beautifyHist(TH1D *h)
{
  h->SetMarkerStyle(20);
  h->SetLineColor(1);
  h->SetMarkerColor(1);
  h->SetMarkerSize(0.8);
  h->SetStats(0);
  h->SetMinimum(0);
  h->GetYaxis()->SetTitleOffset(1.55);
  h->GetXaxis()->SetTitleOffset(0.9);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.045);
  h->GetYaxis()->SetMaxDigits(2);
  h->GetYaxis()->SetTitleOffset(1);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetLabelSize(0.045);
}

void beautifyFunc(TF1 *f)
{
  f->SetNpx(1000);
  f->SetLineWidth(2);
}

std::string ftoa1(double i)
{
  char res[10];
  sprintf(res, "%.0f", i);
  std::string ret(res);
  return ret;
}

std::string ftoa2(double i)
{
  char res[10];
  sprintf(res, "%.1f", i);
  std::string ret(res);
  return ret;
}
