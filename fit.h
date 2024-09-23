#ifndef FITVARIABLE_H
#define FITVARIABLE_H

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooExponential.h"  // For RooExponential
#include "RooAddPdf.h"       // For RooAddPdf
#include "RooPlot.h"
#include "TFile.h"           // For TFile
#include "TCanvas.h"         // For TCanvas
#include "TPaveText.h"       // For TPaveText


RooRealVar *x;

// This is BW

RooRealVar *bwmean;
RooRealVar *width;
RooBreitWigner *bw;

//This is CrystalBall
RooRealVar *cbmean;
RooRealVar *cbsigma;
RooRealVar *cbalpha;
RooRealVar *cbn;
RooCBShape *cb;

//This is EXP
RooRealVar *decayVar;
RooExponential *expo;

//This is combined
RooFFTConvPdf* newconvpdf;

//This is combined + exp
RooRealVar *fsig;
RooAddPdf *purepdf;

RooPlot* frame;
RooPlot* framecheck;

RooPlot* pullFrame;
RooHist* residuals;
TPaveText *textBox;

TCanvas *c1;
TCanvas *c2;
TCanvas *c3;




#endif