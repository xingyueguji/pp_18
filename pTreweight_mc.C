#include "plotting_helper.h"

void drawandsave(TFile *f1, TH1D *h_PbPb, TH1D *h_pp, TH1D *h_ratio)
{
    // Create Canvas and Pads
    TCanvas *c = new TCanvas("c", "PbPb / pp MC", 800, 800);

    TPad *pad1 = new TPad("pad1", "Top Pad", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0.02); // Reduce bottom margin for top plot
    pad1->Draw();
    pad1->cd();

    // Upper Pad: Draw PbPb and pp histograms
    h_PbPb->SetLineColor(kBlue);
    h_PbPb->Draw("E");
    h_pp->SetLineColor(kRed);
    h_pp->Draw("E SAME");

    TLegend *leg1 = new TLegend(0.65, 0.75, 0.85, 0.85);
    leg1->AddEntry(h_PbPb, "Pb+Pb MC (5.02 TeV)", "l");
    leg1->AddEntry(h_pp, "p+p MC (13 TeV)", "l");
    leg1->Draw();

    // Lower Pad: Draw the ratio
    c->cd();
    TPad *pad2 = new TPad("pad2", "Bottom Pad", 0, 0.0, 1, 0.3);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.25);
    pad2->Draw();
    pad2->cd();

    h_ratio->SetTitle("");
    h_ratio->GetYaxis()->SetTitle("PbPb / pp p_{T} Ratio");
    h_ratio->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_ratio->GetYaxis()->SetTitleSize(0.08);
    h_ratio->GetYaxis()->SetTitleOffset(0.5);
    h_ratio->GetXaxis()->SetTitleSize(0.1);
    h_ratio->GetXaxis()->SetLabelSize(0.08);
    h_ratio->GetYaxis()->SetLabelSize(0.08);
    h_ratio->Draw("E");

    // Fit Function
    TF1 *fitFunc = new TF1("fitFunc", "[0]/(x + [1])^[2]", 0, 200);
    fitFunc->SetParLimits(1, 5.0, 6.0);
    fitFunc->SetParameters(1.4, 5.1, 0.1); // Initial parameters
    h_ratio->Fit(fitFunc, "R");

    fitFunc->SetLineColor(kOrange);
    fitFunc->Draw("SAME");

    // Fit parameters
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.05);
    latex.DrawLatex(0.2, 0.8, Form("A = %.1f #pm %.1f", fitFunc->GetParameter(0), fitFunc->GetParError(0)));
    latex.DrawLatex(0.2, 0.7, Form("B = %.1f #pm %.1f", fitFunc->GetParameter(1), fitFunc->GetParError(1)));
    latex.DrawLatex(0.2, 0.6, Form("n = %.1f #pm %.1f", fitFunc->GetParameter(2), fitFunc->GetParError(2)));

    f1->cd();
    h_ratio->Write("", 2);

    fitFunc->Write("FA_ratio_fit", 2);
    c->SaveAs("./pTratio/ratio_plot_FA_mc.pdf");
}

void pTreweight_mc()
{

    plotting_helper *ovo = new plotting_helper();

    TFile *f1 = new TFile("../ZBoson_18/rootfile/mc_pTratio.root", "UPDATE");

    TH1D *PbPb_FA = (TH1D *)f1->Get("mc_PbPb_pT");
    TH1D *pp_FA = (TH1D *)f1->Get("mc_pp_pT");

    PbPb_FA->Rebin(4);
    pp_FA->Rebin(4);

    ovo->areanormalize(PbPb_FA);
    ovo->areanormalize(pp_FA);

    // Clone PbPb histograms to store weights
    TH1D *weight_FA = new TH1D("weight_FA", "", 50, 0, 200);

    // Compute weights for FA
    for (int i = 1; i <= PbPb_FA->GetNbinsX(); i++)
    {
        double PbPb_val = PbPb_FA->GetBinContent(i);
        double pp_val = pp_FA->GetBinContent(i);

        double PbPb_err = PbPb_FA->GetBinError(i);
        double pp_err = pp_FA->GetBinError(i);

        if (pp_val > 0 && PbPb_val > 0) // Avoid division by zero
        {
            double ratio = PbPb_val / pp_val;
            double ratio_err = ratio * sqrt(pow(PbPb_err / PbPb_val, 2) + pow(pp_err / pp_val, 2));

            weight_FA->SetBinContent(i, ratio);
            weight_FA->SetBinError(i, ratio_err);
        }
        else
        {
            weight_FA->SetBinContent(i, 0); // Set weight to 0 if denominator is 0
            weight_FA->SetBinError(i, 0);
        }
    }

    drawandsave(f1, PbPb_FA, pp_FA, weight_FA);

    std::cout << "pT reweighting histograms saved in pT_reweight.root" << std::endl;
}