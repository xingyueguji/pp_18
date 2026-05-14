#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"

#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooBreitWigner.h"
#include "RooCBShape.h"
#include "RooFFTConvPdf.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooArgList.h"
#include "RooHist.h"

#include <iostream>
#include <string>

using namespace RooFit;

void fit_antiiso_correct(const char *fname = "./new_pp_data_file_with_iso_0.25.root",
                         double fitLo = 85.0,
                         double fitHi = 120.0,
                         const char *outDir = "./antiiso_fit_plots")
{
    gSystem->mkdir(outDir, kTRUE);

    TFile *f = TFile::Open(fname, "UPDATE");
    if (!f || f->IsZombie())
    {
        std::cerr << "[ERROR] Cannot open " << fname << "\n";
        return;
    }

    TH1D *h_anti = dynamic_cast<TH1D *>(f->Get("FA_nominal_inclusive_antiiso"));
    TH1D *h_incl = dynamic_cast<TH1D *>(f->Get("FA_nominal_inclusive"));
    if (!h_anti || !h_incl)
    {
        std::cerr << "[ERROR] Missing histograms.\n";
        return;
    }
    if (h_anti->GetNbinsX() != h_incl->GetNbinsX())
    {
        std::cerr << "[ERROR] Bin count mismatch.\n";
        return;
    }

    // ---------------- RooFit model ----------------
    RooRealVar mass("mass", "m_{ll} (GeV)", 60.0, 120.0); // full variable range
    mass.setBins(10000, "fft");
    mass.setRange("peakOnly", fitLo, fitHi); // restricted fit window

    RooDataHist data("data", "antiiso data", mass, RooFit::Import(*h_anti));

    RooRealVar mZ("mZ", "mZ", 91.1876);
    // mZ.setConstant(true);
    RooRealVar gZ("gZ", "gZ", 2.4955);
    // gZ.setConstant(true);
    RooBreitWigner bw("bw", "bw", mass, mZ, gZ);

    RooRealVar cb_mean("cb_mean", "cb_mean", 0.0, -5.0, 5.0);
    RooRealVar cb_sigma("cb_sigma", "cb_sigma", 2.0, 0.5, 6.0);
    RooRealVar cb_alpha("cb_alpha", "cb_alpha", 1, 0.1, 5.0);
    RooRealVar cb_n("cb_n", "cb_n", 2.0, 0.5, 10.0);
    RooCBShape cb("cb", "cb", mass, cb_mean, cb_sigma, cb_alpha, cb_n);

    // Smooth background
    RooRealVar bkg_lambda("bkg_lambda", "bkg_lambda", -0.02, -0.5, 0.0);
    RooExponential bkg("bkg", "bkg", mass, bkg_lambda);

    RooFFTConvPdf sig("sig", "BW (x) CB", mass, bw, cb);
    sig.setBufferFraction(5.0); // default is 0.1, try 1.0 to 5.0

    const double Ntot = h_anti->Integral();
    RooRealVar nsig("nsig", "nsig", 0.90 * Ntot, 0.0, 2.0 * Ntot);
    RooRealVar nbkg("nbkg", "nbkg", 0.10 * Ntot, 0.0, 2.0 * Ntot);
    RooAddPdf model("model", "sig+bkg", RooArgList(sig, bkg), RooArgList(nsig, nbkg));

    // ---------------- Fit ----------------
    RooFitResult *fr = model.fitTo(data,
                                   Save(true),
                                   Extended(true),
                                   Range("peakOnly"),
                                   PrintLevel(-1), SumW2Error(true));
    std::cout << "\n[INFO] Fit status: " << fr->status()
              << "   covQual: " << fr->covQual() << "\n";
    fr->Print("v");

    // ---------------- Build extrapolated signal histogram (ONCE) ----------------
    const int nb = h_incl->GetNbinsX();
    const double xlo = h_incl->GetXaxis()->GetXmin();
    const double xhi = h_incl->GetXaxis()->GetXmax();

    // --- Build histograms from each PDF (unit-area shapes over [xlo, xhi]) ---
    TH1 *h_sig_pdf = sig.createHistogram("h_sig_extracted_tmp", mass,
                                         Binning(nb, xlo, xhi));
    TH1 *h_bkg_pdf = bkg.createHistogram("h_bkg_extracted_tmp", mass,
                                         Binning(nb, xlo, xhi));
    if (!h_sig_pdf || !h_bkg_pdf)
    {
        std::cerr << "[ERROR] createHistogram failed.\n";
        return;
    }
    h_sig_pdf->SetDirectory(nullptr);
    h_bkg_pdf->SetDirectory(nullptr);

    // above is new

    // --- PDF fractions inside the fit window (each PDF has its own) ---
    std::unique_ptr<RooAbsReal> I_sig_fit{
        sig.createIntegral(mass, NormSet(mass), Range("peakOnly"))};
    std::unique_ptr<RooAbsReal> I_bkg_fit{
        bkg.createIntegral(mass, NormSet(mass), Range("peakOnly"))};

    const double frac_sig_fit = I_sig_fit->getVal();
    cout << "frac sig is " << frac_sig_fit << endl;
    const double frac_bkg_fit = I_bkg_fit->getVal();
    cout << "frac bkg is " << frac_bkg_fit << endl;

    // --- Total events extrapolated to the full range ---
    const double N_sig_total = nsig.getVal() / frac_sig_fit;
    const double N_bkg_total = nbkg.getVal() / frac_bkg_fit;

    // --- Scale each histogram so its total integral matches N_*_total ---
    const double area_sig = h_sig_pdf->Integral();
    const double area_bkg = h_bkg_pdf->Integral();
    if (area_sig > 0)
        h_sig_pdf->Scale(N_sig_total / area_sig);
    if (area_bkg > 0)
        h_bkg_pdf->Scale(N_bkg_total / area_bkg);

    // --- Per-bin errors using each component's own relative uncertainty ---
    const double rel_sig = nsig.getError() / nsig.getVal();
    const double rel_bkg = nbkg.getError() / nbkg.getVal();
    for (int i = 1; i <= h_sig_pdf->GetNbinsX(); ++i)
    {
        h_sig_pdf->SetBinError(i, h_sig_pdf->GetBinContent(i) * rel_sig);
        h_bkg_pdf->SetBinError(i, h_bkg_pdf->GetBinContent(i) * rel_bkg);
    }

    const double n_total = h_sig_pdf->Integral();
    std::cout << "[INFO] N_sig in fit window  : " << nsig.getVal() << "\n";
    std::cout << "[INFO] N_sig extrapolated   : " << n_total
              << "  (factor " << n_total / nsig.getVal() << ")\n";

    //TH1 *h_diff = (TH1 *)h_sig_pdf->Clone("h_sig_minus_bkg");
    //h_diff->Add(h_bkg_pdf, -1.0);

    // ---------------- Plot ----------------
    RooPlot *frame = mass.frame(Title(""));
    data.plotOn(frame, Name("dat"), MarkerStyle(20), MarkerSize(0.9));

    // extrapolated curve over full range, normalized to fit window
    model.plotOn(frame, Name("ext"),
                 NormRange("peakOnly"), Range(xlo, xhi),
                 LineColor(kRed + 1), LineWidth(2));
    model.plotOn(frame, Components("sig"), Range(xlo, xhi), LineColor(kGreen + 2),
                 LineStyle(kDashed), LineWidth(2));

    // overlay the in-fit-window portion
    model.plotOn(frame, Name("fit"),
                 NormRange("peakOnly"), Range("peakOnly"),
                 LineColor(kBlue + 1), LineWidth(2));
    model.plotOn(frame, Components("bkg"), Range(xlo, xhi), LineColor(kGreen + 2),
                 LineWidth(2));

    RooHist *pullHist = frame->pullHist("dat", "ext"); // pulls of data wrt extrapolated curve

    TCanvas *c = new TCanvas("c_fit", "antiiso fit", 800, 800);
    TPad *p1 = new TPad("p1", "", 0, 0.30, 1, 1.00);
    TPad *p2 = new TPad("p2", "", 0, 0.00, 1, 0.30);
    p1->SetBottomMargin(0.02);
    p1->SetTicks(1, 1);
    p2->SetTopMargin(0.04);
    p2->SetBottomMargin(0.32);
    p2->SetTicks(1, 1);
    p1->Draw();
    p2->Draw();

    // top
    p1->cd();
    frame->GetXaxis()->SetLabelSize(0.0);
    frame->GetXaxis()->SetTitleSize(0.0);
    frame->GetYaxis()->SetTitle("Events");
    frame->GetYaxis()->SetTitleSize(0.05);
    frame->GetYaxis()->SetTitleOffset(1.1);
    frame->SetTitle("");
    frame->Draw();

    TLegend *leg = new TLegend(0.62, 0.62, 0.88, 0.85);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);
    leg->AddEntry(frame->findObject("dat"), "Data", "lep");
    leg->AddEntry(frame->findObject("fit"), "Fit (peakOnly)", "l");
    leg->AddEntry(frame->findObject("ext"), "Extrapolation", "l");
    leg->Draw();

    TLatex tl;
    tl.SetNDC();
    tl.SetTextFont(42);
    tl.SetTextSize(0.035);
    tl.DrawLatex(0.18, 0.86, Form("N_{sig}^{fit}  = %.0f #pm %.0f", nsig.getVal(), nsig.getError()));
    tl.DrawLatex(0.18, 0.81, Form("N_{sig}^{full} = %.0f", n_total));
    tl.DrawLatex(0.18, 0.66, Form("N_{bkg} = %.0f #pm %.0f", nbkg.getVal(), nbkg.getError()));
    tl.DrawLatex(0.18, 0.76, Form("CB #sigma = %.2f #pm %.2f", cb_sigma.getVal(), cb_sigma.getError()));
    tl.DrawLatex(0.18, 0.71, Form("CB #alpha = %.2f #pm %.2f", cb_alpha.getVal(), cb_alpha.getError()));

    // bottom
    p2->cd();
    RooPlot *frame2 = mass.frame(Title(""));
    frame2->addPlotable(pullHist, "P");
    frame2->SetTitle("");
    frame2->GetYaxis()->SetTitle("Pull");
    frame2->GetYaxis()->SetTitleSize(0.12);
    frame2->GetYaxis()->SetTitleOffset(0.45);
    frame2->GetYaxis()->SetLabelSize(0.10);
    frame2->GetYaxis()->SetNdivisions(505);
    frame2->GetYaxis()->SetRangeUser(-5, 5);
    frame2->GetXaxis()->SetTitle("m_{ll} (GeV)");
    frame2->GetXaxis()->SetTitleSize(0.12);
    frame2->GetXaxis()->SetTitleOffset(1.1);
    frame2->GetXaxis()->SetLabelSize(0.10);
    frame2->GetXaxis()->SetRangeUser(fitLo, fitHi); // clip pull display to fit window
    frame2->Draw();

    TLine *l0 = new TLine(fitLo, 0, fitHi, 0);
    l0->SetLineStyle(2);
    l0->Draw("SAME");
    TLine *lp = new TLine(fitLo, 2, fitHi, 2);
    lp->SetLineStyle(3);
    lp->SetLineColor(kRed + 1);
    lp->Draw("SAME");
    TLine *lm = new TLine(fitLo, -2, fitHi, -2);
    lm->SetLineStyle(3);
    lm->SetLineColor(kRed + 1);
    lm->Draw("SAME");

    std::string stem = std::string(outDir) + "/antiiso_fit";
    c->SaveAs((stem + ".png").c_str());
    c->SaveAs((stem + ".pdf").c_str());

    // ---------------- Add extrapolated signal to inclusive ----------------
    TH1D *h_corr = dynamic_cast<TH1D *>(h_incl->Clone("FA_nominal_inclusive_corrected"));
    h_corr->SetDirectory(nullptr);
    h_corr->Add(h_sig_pdf);

    std::cout << "[INFO] Inclusive integral (orig)      : " << h_incl->Integral() << "\n";
    std::cout << "[INFO] Inclusive integral (corrected) : " << h_corr->Integral() << "\n";

    // ---------------- Write back ----------------
    f->cd();
    h_corr->Write("FA_nominal_inclusive_corrected", TObject::kOverwrite);
    //h_diff->Write("FA_nominal_Fitted_sig", TObject::kOverwrite);
    h_sig_pdf->Write("signal",2);
    h_bkg_pdf->Write("bkg",2);
    f->Close();
    delete f;

    std::cout << "[INFO] Wrote 'FA_nominal_inclusive_corrected' into " << fname << "\n";
    std::cout << "[INFO] Plot at " << stem << ".png/.pdf\n";
}