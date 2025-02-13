#include "tdrStyle.C"
#include "CMS_lumi.C"

TH1D *convertGraphToHist(TGraphErrors *graph, const char *histName)
{
    int nPoints = graph->GetN();
    if (nPoints == 0)
    {
        std::cerr << "Error: Graph " << histName << " has no points!" << std::endl;
        return nullptr;
    }

    // Extract bin edges (assuming uniform bin width)
    double x, y;

    // Create histogram with the same binning as the graph
    TH1D *hist = new TH1D(histName, histName, nPoints, 0, 5);

    for (int i = 0; i < nPoints; i++)
    {
        graph->GetPoint(i, x, y);
        double errY = graph->GetErrorY(i);

        int bin = i + 1;
        hist->SetBinContent(bin, y);
        hist->SetBinError(bin, errY);
    }

    return hist;
}

void plotSystematicUncertainty(TCanvas *c1, TGraphErrors *nominal, TGraphErrors *tnpU, TGraphErrors *tnpD,
                               TGraphErrors *acooff, TGraphErrors *no_bk, TGraphErrors *uniform_rebin,
                               TGraphErrors *mass_range, bool isdM, bool iseta)
{
    // Convert all TGraphErrors to TH1D
    TH1D *h_nominal = convertGraphToHist(nominal, "h_nominal");
    TH1D *h_tnpU = convertGraphToHist(tnpU, "h_tnpU");
    TH1D *h_tnpD = convertGraphToHist(tnpD, "h_tnpD");
    TH1D *h_acooff = convertGraphToHist(acooff, "h_acooff");
    TH1D *h_no_bk = convertGraphToHist(no_bk, "h_no_bk");
    TH1D *h_uniform_rebin = convertGraphToHist(uniform_rebin, "h_uniform_rebin");
    TH1D *h_mass_range = convertGraphToHist(mass_range, "h_mass_range");

    h_nominal->SetDirectory(0);
    h_tnpU->SetDirectory(0);
    h_tnpD->SetDirectory(0);
    h_acooff->SetDirectory(0);
    h_no_bk->SetDirectory(0);
    h_uniform_rebin->SetDirectory(0);
    h_mass_range->SetDirectory(0);

    // Create a histogram for total systematic uncertainty
    TH1D *h_nominalStat = (TH1D *)h_nominal->Clone("h_nominalStat");
    TH1D *h_totalSyst = (TH1D *)h_nominal->Clone("h_totalSyst");
    TH1D *h_tnpUSyst = (TH1D *)h_tnpU->Clone("h_tnpUSyst");
    TH1D *h_tnpDSyst = (TH1D *)h_tnpD->Clone("h_tnpDSyst");
    TH1D *h_acooffSyst = (TH1D *)h_acooff->Clone("h_acooffSyst");
    TH1D *h_no_bkSyst = (TH1D *)h_no_bk->Clone("h_no_bkSyst");
    TH1D *h_uniform_rebinSyst = (TH1D *)h_uniform_rebin->Clone("h_uniform_rebinSyst");
    TH1D *h_mass_rangeSyst = (TH1D *)h_mass_range->Clone("h_mass_rangeSyst");

    h_nominalStat->Reset();
    h_totalSyst->Reset(); // Clear bin contents
    h_tnpUSyst->Reset();
    h_tnpDSyst->Reset();
    h_acooffSyst->Reset();
    h_no_bkSyst->Reset();
    h_uniform_rebinSyst->Reset();
    h_mass_rangeSyst->Reset();

    // Compute total systematic uncertainty bin by bin using quadrature sum
    for (int i = 1; i <= h_nominal->GetNbinsX(); i++)
    {
        double nominalerror = h_nominal->GetBinError(i);
        double nominalvalue = h_nominal->GetBinContent(i);
        double tnpUvalue = abs(nominalvalue - h_tnpU->GetBinContent(i));
        double tnpDvalue = abs(nominalvalue - h_tnpD->GetBinContent(i));
        double acooffvalue = abs(nominalvalue - h_acooff->GetBinContent(i));
        double no_bkvalue = abs(nominalvalue - h_no_bk->GetBinContent(i));
        double uniformvalue = abs(nominalvalue - h_uniform_rebin->GetBinContent(i));
        double mass_rangevalue = abs(nominalvalue - h_mass_range->GetBinContent(i));

        cout << "Syst are " << tnpUvalue << " " << tnpDvalue << " " << acooffvalue << " " << no_bkvalue << " " << uniformvalue << " " << mass_rangevalue << endl;

        double totalSyst = tnpUvalue * tnpUvalue + tnpDvalue * tnpDvalue + acooffvalue * acooffvalue + no_bkvalue * no_bkvalue + mass_rangevalue * mass_rangevalue;
        totalSyst = sqrt(totalSyst);

        h_nominalStat->SetBinContent(i, nominalerror);
        h_totalSyst->SetBinContent(i, totalSyst);
        h_tnpUSyst->SetBinContent(i, tnpUvalue);
        h_tnpDSyst->SetBinContent(i, tnpDvalue);
        h_acooffSyst->SetBinContent(i, acooffvalue);
        h_no_bkSyst->SetBinContent(i, no_bkvalue);
        h_uniform_rebinSyst->SetBinContent(i, uniformvalue);
        h_mass_rangeSyst->SetBinContent(i, mass_rangevalue);
    }

    c1->cd();

    // Plot nominal as a histogram
    /*h_nominal->SetLineColor(kBlack);
    h_nominal->SetLineWidth(2);
    h_nominal->Draw("HIST");*/

    // Plot variations as lines
    // Plot total systematic uncertainty
    h_nominalStat->SetLineColor(kBlack); // Gray shaded region
    h_nominalStat->SetLineWidth(2);

    h_totalSyst->SetLineColor(kTeal); // Gray shaded region
    h_totalSyst->SetLineWidth(2);
    h_totalSyst->SetTitle("");
    h_totalSyst->GetYaxis()->SetTitle("Absolute Syst. Err (GeV)");
    h_totalSyst->GetXaxis()->SetTitle("Centrality");
    h_totalSyst->GetXaxis()->SetLabelSize(0);

    h_totalSyst->GetXaxis()->SetNdivisions(5, 0, 0, kFALSE);
    h_totalSyst->GetXaxis()->ChangeLabel(1, 0, 0.04, 11, -1, -1, "  0-100%");
    h_totalSyst->GetXaxis()->ChangeLabel(2, 0, 0.04, 11, -1, -1, "   0-10%");
    h_totalSyst->GetXaxis()->ChangeLabel(3, 0, 0.04, 11, -1, -1, "  10-20%");
    h_totalSyst->GetXaxis()->ChangeLabel(4, 0, 0.04, 11, -1, -1, "  20-30%");
    h_totalSyst->GetXaxis()->ChangeLabel(5, 0, 0.04, 11, -1, -1, " 30-100%");
    h_totalSyst->GetXaxis()->SetLabelOffset(0.03);
    if (iseta)
    {
        if (isdM)
        {
            h_totalSyst->GetYaxis()->SetRangeUser(0, 0.2);
        }
        else
        {
            h_totalSyst->GetYaxis()->SetRangeUser(0, 0.3);
        }
    }
    if (!iseta)
    {
        if (isdM)
        {
            h_totalSyst->GetYaxis()->SetRangeUser(0, 0.2);
        }
        else
        {
            h_totalSyst->GetYaxis()->SetRangeUser(0, 0.3);
        }
    }
    h_totalSyst->Draw("L");
    h_nominalStat->Draw("L SAME");

    h_tnpUSyst->SetLineWidth(2);
    h_tnpUSyst->SetLineColor(kRed);
    h_tnpUSyst->SetMarkerColor(kRed);

    h_tnpDSyst->SetLineWidth(2);
    h_tnpDSyst->SetLineColor(kBlue);
    h_tnpDSyst->SetMarkerColor(kBlue);

    h_acooffSyst->SetLineWidth(2);
    h_acooffSyst->SetLineColor(kGreen + 2);
    h_acooffSyst->SetMarkerColor(kGreen + 2);

    h_no_bkSyst->SetLineWidth(2);
    h_no_bkSyst->SetLineColor(kMagenta);
    h_no_bkSyst->SetMarkerColor(kMagenta);

    h_uniform_rebinSyst->SetLineWidth(2);
    h_uniform_rebinSyst->SetLineColor(kGray + 1);
    h_uniform_rebinSyst->SetMarkerColor(kGray + 1);

    h_mass_rangeSyst->SetLineWidth(2);
    h_mass_rangeSyst->SetLineColor(kOrange + 7);
    h_mass_rangeSyst->SetMarkerColor(kOrange + 7);

    h_tnpUSyst->Draw("L SAME");
    h_tnpDSyst->Draw("L SAME");
    h_acooffSyst->Draw("L SAME");
    h_no_bkSyst->Draw("L SAME");
    //h_uniform_rebinSyst->Draw("L SAME");
    h_mass_rangeSyst->Draw("L SAME");

    // Add legend
    TLegend *legend = new TLegend(0.6, 0.6, 0.88, 0.88);
    // legend->AddEntry(h_nominal, "Nominal", "L");
    legend->AddEntry(h_nominalStat, "Nominal stat. uncertainty", "L");
    legend->AddEntry(h_totalSyst, "Total syst. uncertainty", "L");
    legend->AddEntry(h_tnpUSyst, "tnpU", "L");
    legend->AddEntry(h_tnpDSyst, "tnpD", "L");
    legend->AddEntry(h_acooffSyst, "acooff", "L");
    legend->AddEntry(h_no_bkSyst, "bk_off", "L");
    //legend->AddEntry(h_uniform_rebinSyst, "Uniform rebin", "L");
    legend->AddEntry(h_mass_rangeSyst, "Mass range", "L");
    legend->SetTextSize(0.03);              // Set text size
    legend->SetTextFont(42);                // Use a modern, clean font
    legend->SetBorderSize(0);               // Set border size (0 for no border)
    legend->SetLineColor(kBlack);           // Border color (if any)
    legend->SetLineWidth(0);                // Border line width
    legend->SetFillColorAlpha(kWhite, 0.1); // Background color with transparency
    legend->Draw("SAME");

    // Update canvas
    gPad->RedrawAxis();
}

void get_systematic_stack()
{
    setTDRStyle();

    TFile *f1 = new TFile("./nominal_syst_tgraph.root", "UPDATE");

    TGraphErrors *HI_sub_pp_dM_chi2_raw_nominal = (TGraphErrors *)f1->Get("HI_sub_pp_dM_chi2_raw_nominal");
    TGraphErrors *HI_sub_pp_dM_chi2_raw_tnpU = (TGraphErrors *)f1->Get("HI_sub_pp_dM_chi2_raw_tnpU");
    TGraphErrors *HI_sub_pp_dM_chi2_raw_tnpD = (TGraphErrors *)f1->Get("HI_sub_pp_dM_chi2_raw_tnpD");
    TGraphErrors *HI_sub_pp_dM_chi2_raw_acooff = (TGraphErrors *)f1->Get("HI_sub_pp_dM_chi2_raw_acooff");
    TGraphErrors *HI_sub_pp_dM_chi2_raw_nominal_no_bk = (TGraphErrors *)f1->Get("HI_sub_pp_dM_chi2_raw_nominal_no_bk");
    TGraphErrors *HI_sub_pp_dM_chi2_raw_nominal_uniform_rebin = (TGraphErrors *)f1->Get("HI_sub_pp_dM_chi2_raw_nominal_uniform_rebin");
    TGraphErrors *HI_sub_pp_dM_chi2_raw_nominal_mass_range = (TGraphErrors *)f1->Get("HI_sub_pp_dM_chi2_raw_nominal_mass_range");

    TGraphErrors *HI_sub_pp_dWidth_chi2_raw_nominal = (TGraphErrors *)f1->Get("HI_sub_pp_dWidth_chi2_raw_nominal");
    TGraphErrors *HI_sub_pp_dWidth_chi2_raw_tnpU = (TGraphErrors *)f1->Get("HI_sub_pp_dWidth_chi2_raw_tnpU");
    TGraphErrors *HI_sub_pp_dWidth_chi2_raw_tnpD = (TGraphErrors *)f1->Get("HI_sub_pp_dWidth_chi2_raw_tnpD");
    TGraphErrors *HI_sub_pp_dWidth_chi2_raw_acooff = (TGraphErrors *)f1->Get("HI_sub_pp_dWidth_chi2_raw_acooff");
    TGraphErrors *HI_sub_pp_dWidth_chi2_raw_nominal_no_bk = (TGraphErrors *)f1->Get("HI_sub_pp_dWidth_chi2_raw_nominal_no_bk");
    TGraphErrors *HI_sub_pp_dWidth_chi2_raw_nominal_uniform_rebin = (TGraphErrors *)f1->Get("HI_sub_pp_dWidth_chi2_raw_nominal_uniform_rebin");
    TGraphErrors *HI_sub_pp_dWidth_chi2_raw_nominal_mass_range = (TGraphErrors *)f1->Get("HI_sub_pp_dWidth_chi2_raw_nominal_mass_range");

    TGraphErrors *HI_sub_pp_dM_chi2_eta_nominal = (TGraphErrors *)f1->Get("HI_sub_pp_dM_chi2_eta_nominal");
    TGraphErrors *HI_sub_pp_dM_chi2_eta_tnpU = (TGraphErrors *)f1->Get("HI_sub_pp_dM_chi2_eta_tnpU");
    TGraphErrors *HI_sub_pp_dM_chi2_eta_tnpD = (TGraphErrors *)f1->Get("HI_sub_pp_dM_chi2_eta_tnpD");
    TGraphErrors *HI_sub_pp_dM_chi2_eta_acooff = (TGraphErrors *)f1->Get("HI_sub_pp_dM_chi2_eta_acooff");
    TGraphErrors *HI_sub_pp_dM_chi2_eta_nominal_no_bk = (TGraphErrors *)f1->Get("HI_sub_pp_dM_chi2_eta_nominal_no_bk");
    TGraphErrors *HI_sub_pp_dM_chi2_eta_nominal_uniform_rebin = (TGraphErrors *)f1->Get("HI_sub_pp_dM_chi2_eta_nominal_uniform_rebin");
    TGraphErrors *HI_sub_pp_dM_chi2_eta_nominal_mass_range = (TGraphErrors *)f1->Get("HI_sub_pp_dM_chi2_eta_nominal_mass_range");

    TGraphErrors *HI_sub_pp_dWidth_chi2_eta_nominal = (TGraphErrors *)f1->Get("HI_sub_pp_dWidth_chi2_eta_nominal");
    TGraphErrors *HI_sub_pp_dWidth_chi2_eta_tnpU = (TGraphErrors *)f1->Get("HI_sub_pp_dWidth_chi2_eta_tnpU");
    TGraphErrors *HI_sub_pp_dWidth_chi2_eta_tnpD = (TGraphErrors *)f1->Get("HI_sub_pp_dWidth_chi2_eta_tnpD");
    TGraphErrors *HI_sub_pp_dWidth_chi2_eta_acooff = (TGraphErrors *)f1->Get("HI_sub_pp_dWidth_chi2_eta_acooff");
    TGraphErrors *HI_sub_pp_dWidth_chi2_eta_nominal_no_bk = (TGraphErrors *)f1->Get("HI_sub_pp_dWidth_chi2_eta_nominal_no_bk");
    TGraphErrors *HI_sub_pp_dWidth_chi2_eta_nominal_uniform_rebin = (TGraphErrors *)f1->Get("HI_sub_pp_dWidth_chi2_eta_nominal_uniform_rebin");
    TGraphErrors *HI_sub_pp_dWidth_chi2_eta_nominal_mass_range = (TGraphErrors *)f1->Get("HI_sub_pp_dWidth_chi2_eta_nominal_mass_range");

    TCanvas *c_HI_sub_pp_dM_chi2_raw_stack = new TCanvas("c_HI_sub_pp_dM_chi2_raw_stack", "", 800, 800);
    TCanvas *c_HI_sub_pp_dWidth_chi2_raw_stack = new TCanvas("c_HI_sub_pp_dWidth_chi2_raw_stack", "", 800, 800);
    TCanvas *c_HI_sub_pp_dM_chi2_eta_stack = new TCanvas("c_HI_sub_pp_dM_chi2_eta_stack", "", 800, 800);
    TCanvas *c_HI_sub_pp_dWidth_chi2_eta_stack = new TCanvas("c_HI_sub_pp_dWidth_chi2_eta_stack", "", 800, 800);

    plotSystematicUncertainty(c_HI_sub_pp_dM_chi2_raw_stack, HI_sub_pp_dM_chi2_raw_nominal, HI_sub_pp_dM_chi2_raw_tnpU, HI_sub_pp_dM_chi2_raw_tnpD, HI_sub_pp_dM_chi2_raw_acooff,
                              HI_sub_pp_dM_chi2_raw_nominal_no_bk, HI_sub_pp_dM_chi2_raw_nominal_uniform_rebin, HI_sub_pp_dM_chi2_raw_nominal_mass_range, 1, 0);
    CMS_lumi(c_HI_sub_pp_dM_chi2_raw_stack, 13, 10);

    plotSystematicUncertainty(c_HI_sub_pp_dWidth_chi2_raw_stack, HI_sub_pp_dWidth_chi2_raw_nominal, HI_sub_pp_dWidth_chi2_raw_tnpU, HI_sub_pp_dWidth_chi2_raw_tnpD, HI_sub_pp_dWidth_chi2_raw_acooff,
                              HI_sub_pp_dWidth_chi2_raw_nominal_no_bk, HI_sub_pp_dWidth_chi2_raw_nominal_uniform_rebin, HI_sub_pp_dWidth_chi2_raw_nominal_mass_range, 0, 0);
    CMS_lumi(c_HI_sub_pp_dWidth_chi2_raw_stack, 13, 10);

    plotSystematicUncertainty(c_HI_sub_pp_dM_chi2_eta_stack, HI_sub_pp_dM_chi2_eta_nominal, HI_sub_pp_dM_chi2_eta_tnpU, HI_sub_pp_dM_chi2_eta_tnpD, HI_sub_pp_dM_chi2_eta_acooff,
                              HI_sub_pp_dM_chi2_eta_nominal_no_bk, HI_sub_pp_dM_chi2_eta_nominal_uniform_rebin, HI_sub_pp_dM_chi2_eta_nominal_mass_range, 1, 1);
    CMS_lumi(c_HI_sub_pp_dM_chi2_eta_stack, 13, 10);

    plotSystematicUncertainty(c_HI_sub_pp_dWidth_chi2_eta_stack, HI_sub_pp_dWidth_chi2_eta_nominal, HI_sub_pp_dWidth_chi2_eta_tnpU, HI_sub_pp_dWidth_chi2_eta_tnpD, HI_sub_pp_dWidth_chi2_eta_acooff,
                              HI_sub_pp_dWidth_chi2_eta_nominal_no_bk, HI_sub_pp_dWidth_chi2_eta_nominal_uniform_rebin, HI_sub_pp_dWidth_chi2_eta_nominal_mass_range, 0, 1);
    CMS_lumi(c_HI_sub_pp_dWidth_chi2_eta_stack, 13, 10);

    c_HI_sub_pp_dM_chi2_raw_stack->SaveAs("./systematic/FA_dM_diff_stack.png");
    c_HI_sub_pp_dWidth_chi2_raw_stack->SaveAs("./systematic/FA_dW_diff_stack.png");
    c_HI_sub_pp_dM_chi2_eta_stack->SaveAs("./systematic/Eta_dM_diff_stack.png");
    c_HI_sub_pp_dWidth_chi2_eta_stack->SaveAs("./systematic/Eta_dW_diff_stack.png");
}