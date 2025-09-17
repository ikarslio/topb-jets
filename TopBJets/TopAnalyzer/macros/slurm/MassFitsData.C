#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include <TAttFill.h>

// Helper functions, forward declared
std::string ftoa1(double i); std::string ftoa2(double i); // Converts double to string
void beautifyHist(TH1D *h);
void beautifyFunc(TF1 *f);

Double_t exponential(Double_t *x, Double_t *par);
Double_t background(Double_t *x, Double_t *par);
Double_t gaussianPeak(Double_t *x, Double_t *par);
Double_t doubleGaussianPeak(Double_t *x, Double_t *par);
Double_t fitFunction1(Double_t *x, Double_t *par);
Double_t fitFunction2(Double_t *x, Double_t *par);
Double_t fitFunction3(Double_t *x, Double_t *par);

void MassFitsData() {

  gROOT->ForceStyle();

  TFile* f1 = TFile::Open("data.root");

  TH1D *h1_mass1	= (TH1D*)f1->Get("h_D0_SameSign_mass");
  TH1D *h1_mass2 	= (TH1D*)f1->Get("h_D0_OppSign_mass");
  TH1D *h2_mass1	= (TH1D*)f1->Get("h_DstarD0_SameSign_mass");
  TH1D *h2_mass2 	= (TH1D*)f1->Get("h_DstarD0_OppSign_mass");
  TH1D *h3_mass1	= (TH1D*)f1->Get("h_DPlus_SameSign_mass");
  TH1D *h3_mass2 	= (TH1D*)f1->Get("h_DPlus_OppSign_mass");
  TH1D *h4_mass1	= (TH1D*)f1->Get("h_DsPlus_SameSign_mass");
  TH1D *h4_mass2 	= (TH1D*)f1->Get("h_DsPlus_OppSign_mass");
  TH1D *h5_mass1	= (TH1D*)f1->Get("h_LambdacPlus_SameSign_mass");
  TH1D *h5_mass2 	= (TH1D*)f1->Get("h_LambdacPlus_OppSign_mass");

  // ***************************************************************** c_MassFit_D0***************************************************************** //

  double binWidth1 = h1_mass1->GetBinWidth(1);
  std::cout<<"Bin Width = "<<binWidth1<<std::endl;

  h1_mass1->SetTitle(("; M(D^{0} #rightarrow K^{-} #pi^{+}) [GeV]; Entries / "+ftoa1(binWidth1*1000)+" MeV").c_str());
  TCanvas *c_MassFit_D0 = new TCanvas("c_MassFit_D0", "",125,111,700,600);
  TF1 *f_function_D0 = new TF1("f_function_D0", "fitFunction2", 1.73,1.98, 8);
  //f_function_D0->SetParameters(10, 10, 24, 1.861, 0.005, 19, 1.866, 0.005);
  f_function_D0->SetParameters(20, 20, 55, 1.862, 0.005, 46, 1.8695, 0.005);
  h1_mass1->Fit("f_function_D0", "L", "ep");
  Double_t parameters_1[8];
  f_function_D0->GetParameters(parameters_1);

  TF1 *f_background_D0 = new TF1("f_background_D0", exponential, 1.73, 1.98, 2);
  f_background_D0->SetParameters(parameters_1);
  TF1 *f_signal_D0 = new TF1("f_signal_D0", doubleGaussianPeak, 1.73, 1.98, 6);
  f_signal_D0->SetParameters(parameters_1+2);

  // Beautify everything before putting on canvas
  beautifyHist(h1_mass1);
  h1_mass1->SetLineColor(kBlack);
  h1_mass2->SetFillStyle(3002);
  h1_mass2->SetFillColor(1);
  h1_mass2->SetLineColor(0);
  beautifyFunc(f_function_D0);
  beautifyFunc(f_background_D0);
  f_function_D0->SetLineColor(kBlue);
  f_background_D0->SetLineColor(kRed);

  h1_mass1->SetMaximum(130);
  h1_mass1->GetYaxis()->SetNdivisions(509);
  h1_mass1->Draw("HIST pe");
  h1_mass2->Scale(0.5);
  h1_mass2->Draw("HIST same");
  f_function_D0->Draw("same");
  f_background_D0->Draw("same");
  TLegend *legend_1=new TLegend(0.14,0.72,0.47,0.86);
  legend_1->SetLineColor(0);
  legend_1->SetTextSize(0.03);
  //legend_1->SetHeader("Reconstructed Mass Combinations","C");
  legend_1->AddEntry(h1_mass1,"RS","lpe");
  legend_1->AddEntry(h1_mass2,"WS","f");
  legend_1->Draw();

  cout<<"Integral = "<<f_signal_D0->Integral(1.73, 1.98)<<endl;
  cout<<"Signal = "<<f_signal_D0->Integral(1.8, 1.9)<<endl;
  cout<<"Bcakground = "<<f_background_D0->Integral(1.73, 1.8)<<"   "<<f_background_D0->Integral(1.9, 1.98)<<endl;

  c_MassFit_D0->SaveAs("Plots/Data/c_MassFit_D0.pdf");

  // *************************************************************** c_MassFit_D*-D0 ***************************************************************** //

  double binWidth2 = h2_mass1->GetBinWidth(1);
  std::cout<<"Bin Width = "<<binWidth2<<std::endl;

  h2_mass1->SetTitle(("; M(D^{*+}) - M(D^{0}) [GeV]; Entries / "+ftoa2(binWidth2*1000)+" MeV").c_str());
  TCanvas *c_MassFit_Dstar = new TCanvas("c_MassFit_Dstar", "",125,111,700,600);
  TF1 *f_function_Dstar = new TF1("f_function_Dstar", "fitFunction1", 0.135, 0.17, 5);
  f_function_Dstar->SetParameters(1, 1, 25, 0.1466, 0.001);
  h2_mass1->Fit("f_function_Dstar", "L", "ep");
  Double_t parameters_2[5];
  f_function_Dstar->GetParameters(parameters_2);

  TF1 *f_background_Dstar = new TF1("f_background_Dstar", background, 0.135, 0.17, 2);
  f_background_Dstar->SetParameters(parameters_2);
  TF1 *f_signal_Dstar = new TF1("f_signal_Dstar", gaussianPeak, 0.135, 0.17, 3);
  f_signal_Dstar->SetParameters(parameters_2+2);

  // Beautify everything before putting on canvas
  beautifyHist(h2_mass1);
  h2_mass1->SetLineColor(kBlack);
  h2_mass2->SetFillStyle(3002);
  h2_mass2->SetFillColor(1);
  h2_mass2->SetLineColor(0);
  beautifyFunc(f_function_Dstar);
  beautifyFunc(f_background_Dstar);
  f_function_Dstar->SetLineColor(kBlue);
  f_background_Dstar->SetLineColor(kRed);

  h2_mass1->SetMaximum(40);
  h2_mass1->GetYaxis()->SetNdivisions(509);
  h2_mass1->Draw("HIST pe");
  h2_mass2->Scale(0.5);
  h2_mass2->Draw("HIST same");
  f_function_Dstar->Draw("same");
  f_background_Dstar->Draw("same");
  TLegend *legend_2=new TLegend(0.14,0.72,0.47,0.86);
  legend_2->SetLineColor(0);
  legend_2->SetTextSize(0.03);
  //legend_2->SetHeader("Reconstructed Mass Combinations","C");
  legend_2->AddEntry(h2_mass1,"RS","lpe");
  legend_2->AddEntry(h2_mass2,"WS","f");
  legend_2->Draw();

  cout<<"Integral = "<<f_signal_Dstar->Integral(0.135, 0.17)<<endl;

  c_MassFit_Dstar->SaveAs("Plots/Data/c_MassFit_Dstar.pdf");

  // ***************************************************************** c_MassFit_DPlus***************************************************************** //

  double binWidth3 = h3_mass1->GetBinWidth(1);
  std::cout<<"Bin Width = "<<binWidth3<<std::endl;

  h3_mass1->SetTitle(("; M(D^{+} #rightarrow K^{-} #pi^{+} #pi^{+}) [GeV]; Entries / "+ftoa1(binWidth3*1000)+" MeV").c_str());
  TCanvas *c_MassFit_DPlus = new TCanvas("c_MassFit_DPlus", "",125,111,700,600);
  TF1 *f_function_DPlus = new TF1("f_function_DPlus", "fitFunction2", 1.7, 2, 8);
  //f_function_DPlus->SetParameters(40, 11, 31, 1.865, 0.05, 40, 1.866, 0.01);
  //f_function_DPlus->SetParameters(41, 9, 45, 1.863, 0.05, 49, 1.866, 0.01);
  f_function_DPlus->SetParameters(41, 10, 65, 1.863, 0.05, 49, 1.866, 0.01);
  //f_function_DPlus->SetParameters(41, 11, 55, 1.865, 0.005, 70, 1.865, 0.001);
  h3_mass1->Fit("f_function_DPlus", "L", "ep");
  Double_t parameters_3[8];
  f_function_DPlus->GetParameters(parameters_3);

  TF1 *f_background_DPlus = new TF1("f_background_DPlus", exponential, 1.7, 2, 2);
  f_background_DPlus->SetParameters(parameters_3);
  TF1 *f_signal_DPlus = new TF1("f_signal_DPlus", doubleGaussianPeak, 1.7, 2, 6);
  f_signal_DPlus->SetParameters(parameters_3+2);

  // Beautify everything before putting on canvas
  beautifyHist(h3_mass1);
  h3_mass1->SetMarkerColor(kBlack);
  h3_mass2->SetFillStyle(3002);
  h3_mass2->SetFillColor(1);
  h3_mass2->SetLineColor(0);
  beautifyFunc(f_function_DPlus);
  beautifyFunc(f_background_DPlus);
  f_function_DPlus->SetLineColor(kBlue);
  f_background_DPlus->SetLineColor(kRed);

  h3_mass1->SetMaximum(120);
  h3_mass1->GetYaxis()->SetNdivisions(509);
  h3_mass1->Draw("HIST pe");
  h3_mass2->Scale(0.5);
  h3_mass2->Draw("HIST same");
  f_function_DPlus->Draw("same");
  f_background_DPlus->Draw("same");
  TLegend *legend_3=new TLegend(0.14,0.72,0.47,0.86);
  legend_3->SetLineColor(0);
  legend_3->SetTextSize(0.03);
  //legend_3->SetHeader("Reconstructed Mass Combinations","C");
  legend_3->AddEntry(h3_mass1,"RS","lpe");
  legend_3->AddEntry(h3_mass2,"WS","f");
  legend_3->Draw();

  cout<<"Integral = "<<f_signal_DPlus->Integral(1.7, 2)<<endl;

  c_MassFit_DPlus->SaveAs("Plots/Data/c_MassFit_DPlus.pdf");

  // *************************************************************** c_MassFit_DsPlus **************************************************************** //

  double binWidth4 = h4_mass1->GetBinWidth(1);
  std::cout<<"Bin Width = "<<binWidth4<<std::endl;

  h4_mass1->SetTitle(("; M(D_{s}^{+} #rightarrow K^{+} K^{-} #pi^{+}) [GeV]; Entries / "+ftoa1(binWidth4*1000)+" MeV").c_str());
  TCanvas *c_MassFit_DsPlus = new TCanvas("c_MassFit_DsPlus", "",125,111,700,600);
  TF1 *f_function_DsPlus = new TF1("f_function_DsPlus", "fitFunction2", 1.8, 2.1, 8);
  f_function_DsPlus->SetParameters(10, 10, 9, 1.87, 0.002, 16, 1.965, 0.005);
  //f_function_DsPlus->SetParameters(10, 10, 23, 1.965, 0.005, 25, 1.968, 0.004);
  h4_mass1->Fit("f_function_DsPlus", "L", "ep");
  Double_t parameters_4[8];
  f_function_DsPlus->GetParameters(parameters_4);

  TF1 *f_background_DsPlus = new TF1("f_background_DsPlus", exponential, 1.8, 2.1, 2);
  f_background_DsPlus->SetParameters(parameters_4);
  TF1 *f_signal_DsPlus = new TF1("f_signal_DsPlus", doubleGaussianPeak, 1.8, 2.1, 6);
  f_signal_DsPlus->SetParameters(parameters_4+2);

  // Beautify everything before putting on canvas
  beautifyHist(h4_mass1);
  h4_mass1->SetLineColor(kBlack);
  h4_mass2->SetFillStyle(3002);
  h4_mass2->SetFillColor(1);
  h4_mass2->SetLineColor(0);
  beautifyFunc(f_function_DsPlus);
  beautifyFunc(f_background_DsPlus);
  f_function_DsPlus->SetLineColor(kBlue);
  f_background_DsPlus->SetLineColor(kRed);

  h4_mass1->SetMaximum(50);
  h4_mass1->GetXaxis()->SetTitleOffset(0.88);
  h4_mass1->GetYaxis()->SetNdivisions(509);
  h4_mass1->Draw("HIST pe");
  h4_mass2->Scale(0.5);
  h4_mass2->Draw("HIST same");
  f_function_DsPlus->Draw("same");
  f_background_DsPlus->Draw("same");
  TLegend *legend_4=new TLegend(0.14,0.72,0.47,0.86);
  legend_4->SetLineColor(0);
  legend_4->SetTextSize(0.03);
  //legend_4->SetHeader("Reconstructed Mass Combinations","C");
  legend_4->AddEntry(h4_mass1,"RS","lpe");
  legend_4->AddEntry(h4_mass2,"WS","f");
  legend_4->Draw();

  cout<<"Integral = "<<f_signal_DsPlus->Integral(1.9, 2.1)<<endl;

  c_MassFit_DsPlus->SaveAs("Plots/Data/c_MassFit_DsPlus.pdf");

  // ********************************************************* c_MassFit_LambdacPlus ************************************************************** //

  double binWidth5 = h5_mass1->GetBinWidth(1);
  std::cout<<"Bin Width = "<<ftoa1(binWidth5)<<std::endl;

  h5_mass1->SetTitle(("; M(#Lambda_{c}^{+} #rightarrow K^{-} #pi^{+} p) [GeV]; Entries / "+ftoa1(binWidth5*1000)+" MeV").c_str());
  TCanvas *c_MassFit_LambdacPlus = new TCanvas("c_MassFit_LambdacPlus", "",125,111,700,600);
  TF1 *f_function_LambdacPlus = new TF1("f_function_LambdacPlus", "fitFunction2", 2.1, 2.4, 8);
  f_function_LambdacPlus->SetParameters(20, 1, 15, 2.25, 0.005, 10, 2.29, 0.005);
  h5_mass1->Fit("f_function_LambdacPlus", "L", "ep");
  Double_t parameters_5[8];
  f_function_LambdacPlus->GetParameters(parameters_5);

  TF1 *f_background_LambdacPlus = new TF1("f_background_LambdacPlus", background, 2.1, 2.4, 2);
  f_background_LambdacPlus->SetParameters(parameters_5);
  TF1 *f_signal_LambdacPlus = new TF1("f_signal_LambdacPlus", doubleGaussianPeak, 2.1, 2.4, 6);
  f_signal_LambdacPlus->SetParameters(parameters_5+2);

  // Beautify everything before putting on canvas
  beautifyHist(h5_mass1);
  h5_mass1->SetLineColor(kBlack);
  h5_mass2->SetFillStyle(3002);
  h5_mass2->SetFillColor(1);
  h5_mass2->SetLineColor(0);
  beautifyFunc(f_function_LambdacPlus);
  beautifyFunc(f_background_LambdacPlus);
  f_function_LambdacPlus->SetLineColor(kBlue);
  f_background_LambdacPlus->SetLineColor(kRed);

  h5_mass1->SetMaximum(40);
  h5_mass1->GetYaxis()->SetNdivisions(509);
  h5_mass1->Draw("HIST pe");
  h5_mass2->Scale(0.5);
  h5_mass2->Draw("HIST same");
  f_function_LambdacPlus->Draw("same");
  f_background_LambdacPlus->Draw("same");
  TLegend *legend_5=new TLegend(0.14,0.68,0.51,0.86);
  legend_5->SetLineColor(0);
  legend_5->SetTextSize(0.03);
  //legend_5->SetHeader("Reconstructed Mass Combinations","C");
  legend_5->AddEntry(h5_mass1,"RS","lpe");
  legend_5->AddEntry(h5_mass2,"WS","f");
  legend_5->Draw();

  cout<<"Integral = "<<f_signal_LambdacPlus->Integral(2.27, 2.3)<<endl;

  c_MassFit_LambdacPlus->SaveAs("Plots/Data/c_MassFit_LambdacPlus.pdf");
}

Double_t exponential(Double_t *x, Double_t *par) {
   return par[0] + par[1]*x[0];
}

Double_t background(Double_t *x, Double_t *par) {
   return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}

Double_t gaussianPeak(Double_t *x, Double_t *par) {
  return par[0]*TMath::Gaus(x[0],par[1],par[2]);
}

Double_t doubleGaussianPeak(Double_t *x, Double_t *par) {
  return par[0]*TMath::Gaus(x[0],par[1],par[2]) + par[3]*TMath::Gaus(x[0],par[4],par[5]);
}

// Sum of background and peak function
Double_t fitFunction1(Double_t *x, Double_t *par) {
  return exponential(x,par) + gaussianPeak(x,&par[2]);
}

Double_t fitFunction2(Double_t *x, Double_t *par) {
  return exponential(x,par) + doubleGaussianPeak(x,&par[2]);
}

Double_t fitFunction3(Double_t *x, Double_t *par) {
  return background(x,par) + doubleGaussianPeak(x,&par[3]);
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
  h->GetXaxis()->SetTitleOffset(0.9);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetLabelSize(0.045);
  h->GetYaxis()->SetTitleOffset(0.9);
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
