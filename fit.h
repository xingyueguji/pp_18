#ifndef FITVARIABLE_H
#define FITVARIABLE_H

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooExponential.h" // For RooExponential
#include "RooAddPdf.h"      // For RooAddPdf
#include "RooPlot.h"
#include "TFile.h"     // For TFile
#include "TCanvas.h"   // For TCanvas
#include "TPaveText.h" // For TPaveText

RooRealVar *x;

// This is BW

RooRealVar *bwmean;
RooRealVar *width;
RooBreitWigner *bw;

// This is CrystalBall
RooRealVar *cbmean;
RooRealVar *cbsigma;
RooRealVar *cbalpha;
RooRealVar *cbn;
RooCBShape *cb;

// This is EXP
RooRealVar *decayVar;
RooExponential *expo;

// This is combined
RooFFTConvPdf *newconvpdf;

// This is combined + exp
RooRealVar *fsig;
RooAddPdf *purepdf;

// This is for hist stat
RooRealVar *nEvents;
RooExtendPdf *expExtended;

RooPlot *frame;
RooPlot *framecheck;

RooPlot *pullFrame;
RooHist *residuals;
TPaveText *textBox;

TCanvas *c1;
TCanvas *c2;
TCanvas *c3;

RooDataSet *mass_array_eta[22];
RooDataSet *mass_array_raw[22];
TH1D *h_mass_array_raw[22];
TH1D *h_mass_array_eta[22];
TH1D *bksub_hist_pp_raw[22];
TH1D *bksub_hist_pp_eta[22];
RooDataHist *dataHistarray_raw[22];
RooDataHist *dataHistarray_eta[22];
TH1 *fitHist_pp_raw[22];
TH1 *fitHist_pp_eta[22];

TH1D *h_mass_array_raw_inclusive;
RooDataSet *mass_array_raw_inclusive;
TH1D *h_mass_array_eta_inclusive;
RooDataSet *mass_array_eta_inclusive;

Double_t dMass_raw[22];
Double_t dMass_Err_raw[22];
Double_t dWidth_raw[22];
Double_t dWidth_Err_raw[22];
Double_t Alpha_raw[22];
Double_t Alpha_Err_raw[22];
Double_t N_raw[22];
Double_t N_Err_raw[22];
Double_t Std_raw[22];
Double_t Std_Err_raw[22];

Double_t dMass_eta[22];
Double_t dMass_Err_eta[22];
Double_t dWidth_eta[22];
Double_t dWidth_Err_eta[22];
Double_t Alpha_eta[22];
Double_t Alpha_Err_eta[22];
Double_t N_eta[22];
Double_t N_Err_eta[22];
Double_t Std_eta[22];
Double_t Std_Err_eta[22];

Double_t xposition[22];
Double_t xposition_err[22];

// Here's HI area

RooDataSet *HI_mass_array_raw[6];
RooDataSet *HI_mass_array_eta[6];
TH1D *h_HI_mass_array_raw[6];
TH1D *h_HI_mass_array_eta[6];
RooDataHist *HI_dataHistarray[6];
TH1D *bksub_hist_HI_raw[6];
TH1D *bksub_hist_HI_eta[6];
RooDataHist *HI_dataHistarray_raw[6];
RooDataHist *HI_dataHistarray_eta[6];
TH1 *fitHist_HI_raw[6];
TH1 *fitHist_HI_eta[6];

Double_t HI_dMass_raw[6];
Double_t HI_dMass_Err_raw[6];
Double_t HI_dWidth_raw[6];
Double_t HI_dWidth_Err_raw[6];
Double_t HI_Alpha_raw[6];
Double_t HI_Alpha_Err_raw[6];
Double_t HI_N_raw[6];
Double_t HI_N_Err_raw[6];
Double_t HI_Std_raw[6];
Double_t HI_Std_Err_raw[6];

Double_t HI_dMass_eta[6];
Double_t HI_dMass_Err_eta[6];
Double_t HI_dWidth_eta[6];
Double_t HI_dWidth_Err_eta[6];
Double_t HI_Alpha_eta[6];
Double_t HI_Alpha_Err_eta[6];
Double_t HI_N_eta[6];
Double_t HI_N_Err_eta[6];
Double_t HI_Std_eta[6];
Double_t HI_Std_Err_eta[6];

Double_t HI_xposition[6];
Double_t HI_xposition_err[6];

Int_t fittype;

#endif