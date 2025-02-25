#include "tdrStyle.C"
#include "CMS_lumi.C"
void cosmetic(TGraphErrors *g1, TGraphErrors *g2, TGraphErrors *g3, TGraphErrors *g4, TGraphErrors *g5, TGraphErrors *g6, TGraphErrors *g7, TGraphErrors *g8, TGraphErrors *g9, bool isdM, bool isdiff = false)
{
    g1->SetMarkerStyle(20); // Circle marker
    g1->SetMarkerSize(1.2);
    g1->SetLineWidth(2);
    g1->SetLineColor(kBlack);
    g1->SetMarkerColor(kBlack);
    g1->SetTitle("");
    if (!isdiff)
    {
        if (isdM)
        {
            g1->GetYaxis()->SetTitle("#DeltaM = M - M_{PDG} (GeV)");
        }
        else
        {
            g1->GetYaxis()->SetTitle("#Delta#Gamma = #Gamma - #Gamma_{PDG} (GeV)");
        }
    }
    else
    {
        if (isdM)
        {
            g1->GetYaxis()->SetTitle("#DeltaM = M_{HI} - M_{pp} (GeV)");
        }
        else
        {
            g1->GetYaxis()->SetTitle("#Delta#Gamma = #Gamma_{HI} - #Gamma_{pp} (GeV)");
        }
    }

    g1->GetXaxis()->SetTitle("Centrality");

    g1->GetXaxis()->SetLimits(0.5, 6.5);
    g1->GetXaxis()->SetRangeUser(0.5, 6.5);
    g1->GetXaxis()->SetNdivisions(6, 0, 0, kFALSE);

    if (isdiff)
    {
        g1->GetXaxis()->SetLimits(0.5, 5.5);
        g1->GetXaxis()->SetRangeUser(0.5, 5.5);
        g1->GetXaxis()->SetNdivisions(5, 0, 0, kFALSE);
    }
    g1->GetXaxis()->SetLabelSize(0);

    g1->GetXaxis()->ChangeLabel(1, 0, 0.04, 11, -1, -1, " 0-100%");
    g1->GetXaxis()->ChangeLabel(2, 0, 0.04, 11, -1, -1, "  0-10%");
    g1->GetXaxis()->ChangeLabel(3, 0, 0.04, 11, -1, -1, " 10-20%");
    g1->GetXaxis()->ChangeLabel(4, 0, 0.04, 11, -1, -1, " 20-30%");
    g1->GetXaxis()->ChangeLabel(5, 0, 0.04, 11, -1, -1, " 30-100%");
    g1->GetXaxis()->ChangeLabel(6, 0, 0.04, 11, -1, -1, "    pp");

    if (isdiff)
    {
        g1->GetXaxis()->ChangeLabel(1, 0, 0.04, 11, -1, -1, "  0-100%");
        g1->GetXaxis()->ChangeLabel(2, 0, 0.04, 11, -1, -1, "   0-10%");
        g1->GetXaxis()->ChangeLabel(3, 0, 0.04, 11, -1, -1, "  10-20%");
        g1->GetXaxis()->ChangeLabel(4, 0, 0.04, 11, -1, -1, "  20-30%");
        g1->GetXaxis()->ChangeLabel(5, 0, 0.04, 11, -1, -1, " 30-100%");
    }
    g1->GetXaxis()->SetLabelOffset(0.03);

    g2->SetMarkerStyle(21); // Square marker
    g2->SetMarkerSize(1.2);
    g2->SetLineWidth(2);
    g2->SetLineColor(kRed);
    g2->SetMarkerColor(kRed);

    g3->SetMarkerStyle(24); // Triangle marker
    g3->SetMarkerSize(1.2);
    g3->SetLineWidth(2);
    g3->SetLineColor(kBlue);
    g3->SetMarkerColor(kBlue);

    g4->SetMarkerStyle(25); // Cross marker
    g4->SetMarkerSize(1.2);
    g4->SetLineWidth(2);
    g4->SetLineColor(kGreen + 2);
    g4->SetMarkerColor(kGreen + 2);

    g5->SetMarkerStyle(33); // Star marker
    g5->SetMarkerSize(1.4);
    g5->SetLineWidth(2);
    g5->SetLineColor(kMagenta);
    g5->SetMarkerColor(kMagenta);

    g6->SetMarkerStyle(28);
    g6->SetMarkerSize(1.4);
    g6->SetLineWidth(2);
    g6->SetLineColor(kGray + 2);
    g6->SetMarkerColor(kGray + 2);

    g7->SetMarkerStyle(29);
    g7->SetMarkerSize(1.4);
    g7->SetLineWidth(2);
    g7->SetLineColor(kOrange + 7);
    g7->SetMarkerColor(kOrange + 7);

    g8->SetMarkerStyle(32);
    g8->SetMarkerSize(1.4);
    g8->SetLineWidth(2);
    g8->SetLineColor(kCyan);
    g8->SetMarkerColor(kCyan);

    g9->SetMarkerStyle(32);
    g9->SetMarkerSize(1.4);
    g9->SetLineWidth(2);
    g9->SetLineColor(kYellow);
    g9->SetMarkerColor(kYellow);

    // Define small shifts for each graph
    double shift2 = 0.1;  // Small right shift for g2
    double shift3 = -0.1; // Slightly more shift for g3
    double shift4 = 0.2;  // More shift for g4
    double shift5 = -0.2; // Largest shift for g5
    double shift7 = -0.3; // 6
    double shift6 = 0.3;  // 7
    double shift8 = 0.4;
    double shift9 = -0.4;

    // Function to apply X shifts to a TGraphErrors
    auto shiftX = [](TGraphErrors *g, double shift)
    {
        for (int i = 0; i < g->GetN(); i++)
        {
            double x, y;
            g->GetPoint(i, x, y);
            g->SetPoint(i, x + shift, y); // Apply shift
        }
    };

    // Apply shifts to g2, g3, g4, g5
    shiftX(g2, shift2);
    shiftX(g3, shift3);
    shiftX(g4, shift4);
    shiftX(g5, shift5);
    shiftX(g6, shift6);
    shiftX(g7, shift7);
    shiftX(g8, shift8);
    shiftX(g9, shift9);
}

void cosmeticpp(TGraphErrors *g1, TGraphErrors *g2, TGraphErrors *g3, TGraphErrors *g4, TGraphErrors *g5, TGraphErrors *g6, TGraphErrors *g7, bool isdM)
{
    g1->SetMarkerStyle(20); // Circle marker
    g1->SetMarkerSize(1.2);
    g1->SetLineWidth(2);
    g1->SetLineColor(kBlack);
    g1->SetMarkerColor(kBlack);
    g1->SetTitle("");
    g1->GetYaxis()->SetTitleOffset(0.7);

    if (isdM)
    {
        g1->GetYaxis()->SetTitle("#DeltaM = M - M_{PDG} (GeV)");
    }
    else
    {
        g1->GetYaxis()->SetTitle("#Delta#Gamma = #Gamma - #Gamma_{PDG} (GeV)");
    }
    g1->GetXaxis()->SetTitle("Period");
    g1->GetXaxis()->SetLimits(0.5, 22.5);
    g1->GetXaxis()->SetRangeUser(0.5, 22.5);
    g1->GetXaxis()->SetNdivisions(22, 0, 0, kFALSE);
    g1->GetXaxis()->SetLabelSize(0);

    for (int i = 1; i <= 22; i++)
    {
        if (i == 22)
            g1->GetXaxis()->ChangeLabel(i, 0, 0.04, 11, -1, -1, "inclusive");
        else if (i >= 10)
        {
            g1->GetXaxis()->ChangeLabel(i, 0, 0.04, 11, -1, -1, Form("  %i", i));
        }
        else
            g1->GetXaxis()->ChangeLabel(i, 0, 0.04, 11, -1, -1, Form("   %i", i));
    }

    g1->GetXaxis()->SetLabelOffset(0.03);

    g2->SetMarkerStyle(21); // Square marker
    g2->SetMarkerSize(1.2);
    g2->SetLineWidth(2);
    g2->SetLineColor(kRed);
    g2->SetMarkerColor(kRed);

    g3->SetMarkerStyle(24); // Triangle marker
    g3->SetMarkerSize(1.2);
    g3->SetLineWidth(2);
    g3->SetLineColor(kBlue);
    g3->SetMarkerColor(kBlue);

    g4->SetMarkerStyle(25); // Cross marker
    g4->SetMarkerSize(1.2);
    g4->SetLineWidth(2);
    g4->SetLineColor(kGreen + 2);
    g4->SetMarkerColor(kGreen + 2);

    g5->SetMarkerStyle(33); // Star marker
    g5->SetMarkerSize(1.4);
    g5->SetLineWidth(2);
    g5->SetLineColor(kMagenta);
    g5->SetMarkerColor(kMagenta);

    g6->SetMarkerStyle(28);
    g6->SetMarkerSize(1.4);
    g6->SetLineWidth(2);
    g6->SetLineColor(kGray + 2);
    g6->SetMarkerColor(kGray + 2);

    g7->SetMarkerStyle(29);
    g7->SetMarkerSize(1.4);
    g7->SetLineWidth(2);
    g7->SetLineColor(kOrange + 7);
    g7->SetMarkerColor(kOrange + 7);

    /*g8->SetMarkerStyle(32);
    g8->SetMarkerSize(1.4);
    g8->SetLineWidth(2);
    g8->SetLineColor(kCyan);
    g8->SetMarkerColor(kCyan);*/

    // Define small shifts for each graph
    double shift2 = 0.1;  // Small right shift for g2
    double shift3 = -0.1; // Slightly more shift for g3
    double shift4 = 0.2;  // More shift for g4
    double shift5 = -0.2; // Largest shift for g5
    double shift7 = -0.3;
    double shift6 = 0.3;

    // Function to apply X shifts to a TGraphErrors
    auto shiftX = [](TGraphErrors *g, double shift)
    {
        for (int i = 0; i < g->GetN(); i++)
        {
            double x, y;
            g->GetPoint(i, x, y);
            g->SetPoint(i, x + shift, y); // Apply shift
        }
    };

    // Apply shifts to g2, g3, g4, g5
    shiftX(g2, shift2);
    shiftX(g3, shift3);
    shiftX(g4, shift4);
    shiftX(g5, shift5);
    shiftX(g6, shift6);
    shiftX(g7, shift7);
}

TGraphErrors *mergeGraphs(TGraphErrors *g1, TGraphErrors *g2)
{
    if (!g1 || g1->GetN() == 0 || !g2 || g2->GetN() == 0)
        return nullptr;

    int n1 = g1->GetN(); // Number of points in g1

    // Create a new TGraphErrors to store the merged result
    TGraphErrors *g_merged = new TGraphErrors();

    // Copy all points from g1
    for (int i = 0; i < n1; i++)
    {
        double x, y;
        g1->GetPoint(i, x, y);
        double ex = g1->GetErrorX(i);
        double ey = g1->GetErrorY(i);

        g_merged->SetPoint(i, x, y);
        g_merged->SetPointError(i, ex, ey);
    }

    // Get the last point from g2 and attach it to g1, setting x = 6
    int n2 = g2->GetN();
    double last_x, last_y;
    g2->GetPoint(n2 - 1, last_x, last_y);
    double ex_last = g2->GetErrorX(n2 - 1);
    double ey_last = g2->GetErrorY(n2 - 1);

    // Add the last point from g2 with x = 6
    g_merged->SetPoint(n1, 6, last_y);
    g_merged->SetPointError(n1, ex_last, ey_last);

    return g_merged;
}
void draw(TCanvas *c1, TGraphErrors *g1, TGraphErrors *g2, TGraphErrors *g3, TGraphErrors *g4, TGraphErrors *g5, TGraphErrors *g6, TGraphErrors *g7, TGraphErrors *g8, TGraphErrors *g9, bool isdM, bool iseta)
{
    c1->cd();
    /*if (isdM && !iseta)
        g1->GetYaxis()->SetRangeUser(-0.28, 0.1);
    if (isdM && iseta)
        g1->GetYaxis()->SetRangeUser(-0.4, 0.08);
    if (!isdM && !iseta)
        g1->GetYaxis()->SetRangeUser(-0.5, 0.8);
    if (!isdM && iseta)
        g1->GetYaxis()->SetRangeUser(-0.25, 0.6);*/

    if (isdM)
        g1->GetYaxis()->SetRangeUser(-0.4, 0.2);
    if (!isdM)
        g1->GetYaxis()->SetRangeUser(-0.5, 1.0);

    g1->Draw("AP");
    g2->Draw("P SAME");
    g3->Draw("P SAME");
    g4->Draw("P SAME");
    g5->Draw("P SAME");
    // g6->Draw("P SAME");
    g7->Draw("P SAME");
    if (g8 != nullptr && g9 != nullptr)
    {
        g8->Draw("P SAME");
        g9->Draw("P SAME");
    }

    TLegend *legend = new TLegend(0.2, 0.7, 0.4, 0.9);
    legend->AddEntry(g1, "Nominal", "PL");
    legend->AddEntry(g2, "tnpU", "PL");
    legend->AddEntry(g3, "tnpD", "PL");
    legend->AddEntry(g4, "acooff", "PL");
    legend->AddEntry(g5, "bk_off", "PL");
    // legend->AddEntry(g6, "nominal uniform", "PL");
    legend->AddEntry(g7, "nominal range", "PL");
    if (g8 != nullptr && g9 != nullptr)
    {
        legend->AddEntry(g8, "HF up", "PL");
        legend->AddEntry(g9, "HF down", "PL");
    }

    legend->SetTextSize(0.03);              // Set text size
    legend->SetTextFont(42);                // Use a modern, clean font
    legend->SetBorderSize(0);               // Set border size (0 for no border)
    legend->SetLineColor(kBlack);           // Border color (if any)
    legend->SetLineWidth(0);                // Border line width
    legend->SetFillColorAlpha(kWhite, 0.1); // Background color with transparency
    legend->Draw("SAME");

    TPaveText *pave = new TPaveText(0.65, 0.8, 0.88, 0.85, "NDC"); // Normalized coordinates
    pave->SetFillColor(0);                                         // Transparent background
    pave->SetBorderSize(0);                                        // No border
    pave->SetTextSize(0.04);                                       // Adjust text size
    pave->SetTextFont(42);                                         // CMS-like font
    if (iseta)
        pave->AddText("|#eta| < 1");
    if (!iseta)
        pave->AddText("|#eta| < 2.4");
    pave->Draw("same");
}
void drawdiff(TCanvas *c1, TGraphErrors *g1, TGraphErrors *g2, TGraphErrors *g3, TGraphErrors *g4, TGraphErrors *g5, TGraphErrors *g6, TGraphErrors *g7, TGraphErrors *g8, TGraphErrors *g9, bool isdM, bool iseta)
{
    c1->cd();

    /*if (isdM && !iseta)
        g1->GetYaxis()->SetRangeUser(-0.28, 0.1);
    if (isdM && iseta)
        g1->GetYaxis()->SetRangeUser(-0.4, 0.08);
    if (!isdM && !iseta)
        g1->GetYaxis()->SetRangeUser(-0.5, 0.8);
    if (!isdM && iseta)
        g1->GetYaxis()->SetRangeUser(-0.25, 0.6);*/

    if (isdM)
        g1->GetYaxis()->SetRangeUser(-0.4, 0.2);
    if (!isdM)
        g1->GetYaxis()->SetRangeUser(-0.5, 1.0);

    g1->Draw("AP");
    g2->Draw("P SAME");
    g3->Draw("P SAME");
    g4->Draw("P SAME");
    g5->Draw("P SAME");
    // g6->Draw("P SAME");
    g7->Draw("P SAME");
    if (g8 != nullptr && g9 != nullptr)
    {
        g8->Draw("P SAME");
        g9->Draw("P SAME");
    }

    TLegend *legend = new TLegend(0.2, 0.7, 0.4, 0.9);
    legend->AddEntry(g1, "Nominal", "PL");
    legend->AddEntry(g2, "tnpU", "PL");
    legend->AddEntry(g3, "tnpD", "PL");
    legend->AddEntry(g4, "acooff", "PL");
    legend->AddEntry(g5, "bk_off", "PL");
    // legend->AddEntry(g6, "nominal uniform", "PL");
    legend->AddEntry(g7, "nominal range", "PL");
    if (g8 != nullptr && g9 != nullptr)
    {
        legend->AddEntry(g8, "HF up", "PL");
        legend->AddEntry(g9, "HF down", "PL");
    }

    legend->SetTextSize(0.03);              // Set text size
    legend->SetTextFont(42);                // Use a modern, clean font
    legend->SetBorderSize(0);               // Set border size (0 for no border)
    legend->SetLineColor(kBlack);           // Border color (if any)
    legend->SetLineWidth(0);                // Border line width
    legend->SetFillColorAlpha(kWhite, 0.1); // Background color with transparency
    legend->Draw("SAME");

    TPaveText *pave = new TPaveText(0.65, 0.8, 0.88, 0.85, "NDC"); // Normalized coordinates
    pave->SetFillColor(0);                                         // Transparent background
    pave->SetBorderSize(0);                                        // No border
    pave->SetTextSize(0.04);                                       // Adjust text size
    pave->SetTextFont(42);                                         // CMS-like font
    if (iseta)
        pave->AddText("|#eta| < 1");
    if (!iseta)
        pave->AddText("|#eta| < 2.4");
    pave->Draw("same");
}

TGraphErrors *getdiffplot(TGraphErrors *HI_PP, int typeofsystematic = 0)
{
    // 1 = nominal, 2 = tnpU, 3 = tnpD, 4 = acooff, 5 = nominal_no_bk, 6 = uniform_rebin, 7 = massrange

    int nPoints = 5;
    double xVals[nPoints], yVals[nPoints], xErrs[nPoints], yErrs[nPoints];

    double x_last, y_last, err_x_last, err_y_last;
    HI_PP->GetPoint(5, x_last, y_last);
    err_x_last = HI_PP->GetErrorX(5);
    err_y_last = HI_PP->GetErrorY(5);

    for (int i = 0; i < nPoints; ++i)
    {
        double x, y, err_x, err_y;
        HI_PP->GetPoint(i, x, y);
        err_x = HI_PP->GetErrorX(i);
        err_y = HI_PP->GetErrorY(i);

        xVals[i] = x;
        yVals[i] = y - y_last;

        yErrs[i] = sqrt(err_y * err_y + err_y_last * err_y_last); // Quadratic sum
        xErrs[i] = err_x;
    }
    TGraphErrors *g_diff = new TGraphErrors(nPoints, xVals, yVals, xErrs, yErrs);
    return g_diff;
}

void get_systematic()
{
    setTDRStyle();

    TFile *f1 = new TFile("All_plots.root", "READ");

    TGraphErrors *HI_dM_chi2_raw_nominal = (TGraphErrors *)f1->Get("HI_dM_chi2_raw_nominal");
    TGraphErrors *HI_dM_chi2_raw_tnpU = (TGraphErrors *)f1->Get("HI_dM_chi2_raw_tnpU");
    TGraphErrors *HI_dM_chi2_raw_tnpD = (TGraphErrors *)f1->Get("HI_dM_chi2_raw_tnpD");
    TGraphErrors *HI_dM_chi2_raw_acooff = (TGraphErrors *)f1->Get("HI_dM_chi2_raw_acooff");
    TGraphErrors *HI_dM_chi2_raw_nominal_no_bk = (TGraphErrors *)f1->Get("HI_dM_chi2_raw_nominal_no_bk");
    TGraphErrors *HI_dM_chi2_raw_nominal_uniform_rebin = (TGraphErrors *)f1->Get("HI_dM_chi2_raw_nominal_uniform_rebin");
    TGraphErrors *HI_dM_chi2_raw_nominal_mass_range = (TGraphErrors *)f1->Get("HI_dM_chi2_raw_nominal_mass_range");
    TGraphErrors *HI_dM_chi2_raw_HF_up = (TGraphErrors *)f1->Get("HI_dM_chi2_raw_HF_up");
    TGraphErrors *HI_dM_chi2_raw_HF_down = (TGraphErrors *)f1->Get("HI_dM_chi2_raw_HF_down");

    TGraphErrors *HI_dWidth_chi2_raw_nominal = (TGraphErrors *)f1->Get("HI_dWidth_chi2_raw_nominal");
    TGraphErrors *HI_dWidth_chi2_raw_tnpU = (TGraphErrors *)f1->Get("HI_dWidth_chi2_raw_tnpU");
    TGraphErrors *HI_dWidth_chi2_raw_tnpD = (TGraphErrors *)f1->Get("HI_dWidth_chi2_raw_tnpD");
    TGraphErrors *HI_dWidth_chi2_raw_acooff = (TGraphErrors *)f1->Get("HI_dWidth_chi2_raw_acooff");
    TGraphErrors *HI_dWidth_chi2_raw_nominal_no_bk = (TGraphErrors *)f1->Get("HI_dWidth_chi2_raw_nominal_no_bk");
    TGraphErrors *HI_dWidth_chi2_raw_nominal_uniform_rebin = (TGraphErrors *)f1->Get("HI_dWidth_chi2_raw_nominal_uniform_rebin");
    TGraphErrors *HI_dWidth_chi2_raw_nominal_mass_range = (TGraphErrors *)f1->Get("HI_dWidth_chi2_raw_nominal_mass_range");
    TGraphErrors *HI_dWidth_chi2_raw_HF_up = (TGraphErrors *)f1->Get("HI_dWidth_chi2_raw_HF_up");
    TGraphErrors *HI_dWidth_chi2_raw_HF_down = (TGraphErrors *)f1->Get("HI_dWidth_chi2_raw_HF_down");

    TGraphErrors *HI_dM_chi2_eta_nominal = (TGraphErrors *)f1->Get("HI_dM_chi2_eta_nominal");
    TGraphErrors *HI_dM_chi2_eta_tnpU = (TGraphErrors *)f1->Get("HI_dM_chi2_eta_tnpU");
    TGraphErrors *HI_dM_chi2_eta_tnpD = (TGraphErrors *)f1->Get("HI_dM_chi2_eta_tnpD");
    TGraphErrors *HI_dM_chi2_eta_acooff = (TGraphErrors *)f1->Get("HI_dM_chi2_eta_acooff");
    TGraphErrors *HI_dM_chi2_eta_nominal_no_bk = (TGraphErrors *)f1->Get("HI_dM_chi2_eta_nominal_no_bk");
    TGraphErrors *HI_dM_chi2_eta_nominal_uniform_rebin = (TGraphErrors *)f1->Get("HI_dM_chi2_eta_nominal_uniform_rebin");
    TGraphErrors *HI_dM_chi2_eta_nominal_mass_range = (TGraphErrors *)f1->Get("HI_dM_chi2_eta_nominal_mass_range");
    TGraphErrors *HI_dM_chi2_eta_HF_up = (TGraphErrors *)f1->Get("HI_dM_chi2_eta_HF_up");
    TGraphErrors *HI_dM_chi2_eta_HF_down = (TGraphErrors *)f1->Get("HI_dM_chi2_eta_HF_down");

    TGraphErrors *HI_dWidth_chi2_eta_nominal = (TGraphErrors *)f1->Get("HI_dWidth_chi2_eta_nominal");
    TGraphErrors *HI_dWidth_chi2_eta_tnpU = (TGraphErrors *)f1->Get("HI_dWidth_chi2_eta_tnpU");
    TGraphErrors *HI_dWidth_chi2_eta_tnpD = (TGraphErrors *)f1->Get("HI_dWidth_chi2_eta_tnpD");
    TGraphErrors *HI_dWidth_chi2_eta_acooff = (TGraphErrors *)f1->Get("HI_dWidth_chi2_eta_acooff");
    TGraphErrors *HI_dWidth_chi2_eta_nominal_no_bk = (TGraphErrors *)f1->Get("HI_dWidth_chi2_eta_nominal_no_bk");
    TGraphErrors *HI_dWidth_chi2_eta_nominal_uniform_rebin = (TGraphErrors *)f1->Get("HI_dWidth_chi2_eta_nominal_uniform_rebin");
    TGraphErrors *HI_dWidth_chi2_eta_nominal_mass_range = (TGraphErrors *)f1->Get("HI_dWidth_chi2_eta_nominal_mass_range");
    TGraphErrors *HI_dWidth_chi2_eta_HF_up = (TGraphErrors *)f1->Get("HI_dWidth_chi2_eta_HF_up");
    TGraphErrors *HI_dWidth_chi2_eta_HF_down = (TGraphErrors *)f1->Get("HI_dWidth_chi2_eta_HF_down");

    TGraphErrors *pp_dM_chi2_raw_nominal = (TGraphErrors *)f1->Get("pp_dM_chi2_raw_nominal");
    TGraphErrors *pp_dM_chi2_raw_tnpU = (TGraphErrors *)f1->Get("pp_dM_chi2_raw_tnpU");
    TGraphErrors *pp_dM_chi2_raw_tnpD = (TGraphErrors *)f1->Get("pp_dM_chi2_raw_tnpD");
    TGraphErrors *pp_dM_chi2_raw_acooff = (TGraphErrors *)f1->Get("pp_dM_chi2_raw_acooff");
    TGraphErrors *pp_dM_chi2_raw_nominal_no_bk = (TGraphErrors *)f1->Get("pp_dM_chi2_raw_nominal_no_bk");
    TGraphErrors *pp_dM_chi2_raw_nominal_uniform_rebin = (TGraphErrors *)f1->Get("pp_dM_chi2_raw_nominal_uniform_rebin");
    TGraphErrors *pp_dM_chi2_raw_nominal_mass_range = (TGraphErrors *)f1->Get("pp_dM_chi2_raw_nominal_mass_range");

    TGraphErrors *pp_dWidth_chi2_raw_nominal = (TGraphErrors *)f1->Get("pp_dWidth_chi2_raw_nominal");
    TGraphErrors *pp_dWidth_chi2_raw_tnpU = (TGraphErrors *)f1->Get("pp_dWidth_chi2_raw_tnpU");
    TGraphErrors *pp_dWidth_chi2_raw_tnpD = (TGraphErrors *)f1->Get("pp_dWidth_chi2_raw_tnpD");
    TGraphErrors *pp_dWidth_chi2_raw_acooff = (TGraphErrors *)f1->Get("pp_dWidth_chi2_raw_acooff");
    TGraphErrors *pp_dWidth_chi2_raw_nominal_no_bk = (TGraphErrors *)f1->Get("pp_dWidth_chi2_raw_nominal_no_bk");
    TGraphErrors *pp_dWidth_chi2_raw_nominal_uniform_rebin = (TGraphErrors *)f1->Get("pp_dWidth_chi2_raw_nominal_uniform_rebin");
    TGraphErrors *pp_dWidth_chi2_raw_nominal_mass_range = (TGraphErrors *)f1->Get("pp_dWidth_chi2_raw_nominal_mass_range");

    TGraphErrors *pp_dM_chi2_eta_nominal = (TGraphErrors *)f1->Get("pp_dM_chi2_eta_nominal");
    TGraphErrors *pp_dM_chi2_eta_tnpU = (TGraphErrors *)f1->Get("pp_dM_chi2_eta_tnpU");
    TGraphErrors *pp_dM_chi2_eta_tnpD = (TGraphErrors *)f1->Get("pp_dM_chi2_eta_tnpD");
    TGraphErrors *pp_dM_chi2_eta_acooff = (TGraphErrors *)f1->Get("pp_dM_chi2_eta_acooff");
    TGraphErrors *pp_dM_chi2_eta_nominal_no_bk = (TGraphErrors *)f1->Get("pp_dM_chi2_eta_nominal_no_bk");
    TGraphErrors *pp_dM_chi2_eta_nominal_uniform_rebin = (TGraphErrors *)f1->Get("pp_dM_chi2_eta_nominal_uniform_rebin");
    TGraphErrors *pp_dM_chi2_eta_nominal_mass_range = (TGraphErrors *)f1->Get("pp_dM_chi2_eta_nominal_mass_range");

    TGraphErrors *pp_dWidth_chi2_eta_nominal = (TGraphErrors *)f1->Get("pp_dWidth_chi2_eta_nominal");
    TGraphErrors *pp_dWidth_chi2_eta_tnpU = (TGraphErrors *)f1->Get("pp_dWidth_chi2_eta_tnpU");
    TGraphErrors *pp_dWidth_chi2_eta_tnpD = (TGraphErrors *)f1->Get("pp_dWidth_chi2_eta_tnpD");
    TGraphErrors *pp_dWidth_chi2_eta_acooff = (TGraphErrors *)f1->Get("pp_dWidth_chi2_eta_acooff");
    TGraphErrors *pp_dWidth_chi2_eta_nominal_no_bk = (TGraphErrors *)f1->Get("pp_dWidth_chi2_eta_nominal_no_bk");
    TGraphErrors *pp_dWidth_chi2_eta_nominal_uniform_rebin = (TGraphErrors *)f1->Get("pp_dWidth_chi2_eta_nominal_uniform_rebin");
    TGraphErrors *pp_dWidth_chi2_eta_nominal_mass_range = (TGraphErrors *)f1->Get("pp_dWidth_chi2_eta_nominal_mass_range");

    TCanvas *c_FA_dM = new TCanvas("c_FA_dM", "", 800, 800);
    TCanvas *c_FA_dW = new TCanvas("c_FA_dW", "", 800, 800);

    TCanvas *c_Eta_dM = new TCanvas("c_Eta_dM", "", 800, 800);
    TCanvas *c_Eta_dW = new TCanvas("c_Eta_dW", "", 800, 800);

    TCanvas *c_FA_dM_diff = new TCanvas("c_FA_dM_diff", "", 800, 800);
    TCanvas *c_FA_dW_diff = new TCanvas("c_FA_dW_diff", "", 800, 800);

    TCanvas *c_Eta_dM_diff = new TCanvas("c_Eta_dM_diff", "", 800, 800);
    TCanvas *c_Eta_dW_diff = new TCanvas("c_Eta_dW_diff", "", 800, 800);

    TCanvas *c_pp_FA_dM = new TCanvas("c_pp_FA_dM", "", 1600, 800);
    TCanvas *c_pp_FA_dW = new TCanvas("c_pp_FA_dW", "", 1600, 800);

    TCanvas *c_pp_Eta_dM = new TCanvas("c_pp_Eta_dM", "", 1600, 800);
    TCanvas *c_pp_Eta_dW = new TCanvas("c_pp_Eta_dW", "", 1600, 800);

    c_pp_FA_dM->SetLeftMargin(0.12);
    c_pp_FA_dW->SetLeftMargin(0.12);
    c_pp_Eta_dM->SetLeftMargin(0.12);
    c_pp_Eta_dW->SetLeftMargin(0.12);

    // dM, FA, HI
    TGraphErrors *HI_pp_dM_chi2_raw_nominal = mergeGraphs(HI_dM_chi2_raw_nominal, pp_dM_chi2_raw_nominal);
    TGraphErrors *HI_pp_dM_chi2_raw_tnpU = mergeGraphs(HI_dM_chi2_raw_tnpU, pp_dM_chi2_raw_tnpU);
    TGraphErrors *HI_pp_dM_chi2_raw_tnpD = mergeGraphs(HI_dM_chi2_raw_tnpD, pp_dM_chi2_raw_tnpD);
    TGraphErrors *HI_pp_dM_chi2_raw_acooff = mergeGraphs(HI_dM_chi2_raw_acooff, pp_dM_chi2_raw_acooff);
    TGraphErrors *HI_pp_dM_chi2_raw_nominal_no_bk = mergeGraphs(HI_dM_chi2_raw_nominal_no_bk, pp_dM_chi2_raw_nominal_no_bk);
    TGraphErrors *HI_pp_dM_chi2_raw_nominal_uniform_rebin = mergeGraphs(HI_dM_chi2_raw_nominal_uniform_rebin, pp_dM_chi2_raw_nominal_uniform_rebin);
    TGraphErrors *HI_pp_dM_chi2_raw_nominal_mass_range = mergeGraphs(HI_dM_chi2_raw_nominal_mass_range, pp_dM_chi2_raw_nominal_mass_range);
    TGraphErrors *HI_pp_dM_chi2_raw_HF_up = mergeGraphs(HI_dM_chi2_raw_HF_up, pp_dM_chi2_raw_nominal);
    TGraphErrors *HI_pp_dM_chi2_raw_HF_down = mergeGraphs(HI_dM_chi2_raw_HF_down, pp_dM_chi2_raw_nominal);

    TGraphErrors *HI_sub_pp_dM_chi2_raw_nominal = getdiffplot(HI_pp_dM_chi2_raw_nominal);
    TGraphErrors *HI_sub_pp_dM_chi2_raw_tnpU = getdiffplot(HI_pp_dM_chi2_raw_tnpU);
    TGraphErrors *HI_sub_pp_dM_chi2_raw_tnpD = getdiffplot(HI_pp_dM_chi2_raw_tnpD);
    TGraphErrors *HI_sub_pp_dM_chi2_raw_acooff = getdiffplot(HI_pp_dM_chi2_raw_acooff);
    TGraphErrors *HI_sub_pp_dM_chi2_raw_nominal_no_bk = getdiffplot(HI_pp_dM_chi2_raw_nominal_no_bk);
    TGraphErrors *HI_sub_pp_dM_chi2_raw_nominal_uniform_rebin = getdiffplot(HI_pp_dM_chi2_raw_nominal_uniform_rebin);
    TGraphErrors *HI_sub_pp_dM_chi2_raw_nominal_mass_range = getdiffplot(HI_pp_dM_chi2_raw_nominal_mass_range);
    TGraphErrors *HI_sub_pp_dM_chi2_raw_HF_up = getdiffplot(HI_pp_dM_chi2_raw_HF_up);
    TGraphErrors *HI_sub_pp_dM_chi2_raw_HF_down = getdiffplot(HI_pp_dM_chi2_raw_HF_down);

    cosmetic(HI_pp_dM_chi2_raw_nominal, HI_pp_dM_chi2_raw_tnpU, HI_pp_dM_chi2_raw_tnpD, HI_pp_dM_chi2_raw_acooff, HI_pp_dM_chi2_raw_nominal_no_bk,
             HI_pp_dM_chi2_raw_nominal_uniform_rebin, HI_pp_dM_chi2_raw_nominal_mass_range, HI_pp_dM_chi2_raw_HF_up, HI_pp_dM_chi2_raw_HF_down, 1);
    draw(c_FA_dM, HI_pp_dM_chi2_raw_nominal, HI_pp_dM_chi2_raw_tnpU, HI_pp_dM_chi2_raw_tnpD, HI_pp_dM_chi2_raw_acooff, HI_pp_dM_chi2_raw_nominal_no_bk,
         HI_pp_dM_chi2_raw_nominal_uniform_rebin, HI_pp_dM_chi2_raw_nominal_mass_range, HI_pp_dM_chi2_raw_HF_up, HI_pp_dM_chi2_raw_HF_down, 1, 0);

    cosmetic(HI_sub_pp_dM_chi2_raw_nominal, HI_sub_pp_dM_chi2_raw_tnpU, HI_sub_pp_dM_chi2_raw_tnpD, HI_sub_pp_dM_chi2_raw_acooff, HI_sub_pp_dM_chi2_raw_nominal_no_bk,
             HI_sub_pp_dM_chi2_raw_nominal_uniform_rebin, HI_sub_pp_dM_chi2_raw_nominal_mass_range, HI_sub_pp_dM_chi2_raw_HF_up, HI_sub_pp_dM_chi2_raw_HF_down, 1, 1);
    drawdiff(c_FA_dM_diff, HI_sub_pp_dM_chi2_raw_nominal, HI_sub_pp_dM_chi2_raw_tnpU, HI_sub_pp_dM_chi2_raw_tnpD, HI_sub_pp_dM_chi2_raw_acooff, HI_sub_pp_dM_chi2_raw_nominal_no_bk,
             HI_sub_pp_dM_chi2_raw_nominal_uniform_rebin, HI_sub_pp_dM_chi2_raw_nominal_mass_range, HI_sub_pp_dM_chi2_raw_HF_up, HI_sub_pp_dM_chi2_raw_HF_down, 1, 0);

    // dW, FA, HI
    TGraphErrors *HI_pp_dWidth_chi2_raw_nominal = mergeGraphs(HI_dWidth_chi2_raw_nominal, pp_dWidth_chi2_raw_nominal);
    TGraphErrors *HI_pp_dWidth_chi2_raw_tnpU = mergeGraphs(HI_dWidth_chi2_raw_tnpU, pp_dWidth_chi2_raw_tnpU);
    TGraphErrors *HI_pp_dWidth_chi2_raw_tnpD = mergeGraphs(HI_dWidth_chi2_raw_tnpD, pp_dWidth_chi2_raw_tnpD);
    TGraphErrors *HI_pp_dWidth_chi2_raw_acooff = mergeGraphs(HI_dWidth_chi2_raw_acooff, pp_dWidth_chi2_raw_acooff);
    TGraphErrors *HI_pp_dWidth_chi2_raw_nominal_no_bk = mergeGraphs(HI_dWidth_chi2_raw_nominal_no_bk, pp_dWidth_chi2_raw_nominal_no_bk);
    TGraphErrors *HI_pp_dWidth_chi2_raw_nominal_uniform_rebin = mergeGraphs(HI_dWidth_chi2_raw_nominal_uniform_rebin, pp_dWidth_chi2_raw_nominal_uniform_rebin);
    TGraphErrors *HI_pp_dWidth_chi2_raw_nominal_mass_range = mergeGraphs(HI_dWidth_chi2_raw_nominal_mass_range, pp_dWidth_chi2_raw_nominal_mass_range);
    TGraphErrors *HI_pp_dWidth_chi2_raw_HF_up = mergeGraphs(HI_dWidth_chi2_raw_HF_up, pp_dWidth_chi2_raw_nominal);
    TGraphErrors *HI_pp_dWidth_chi2_raw_HF_down = mergeGraphs(HI_dWidth_chi2_raw_HF_down, pp_dWidth_chi2_raw_nominal);

    TGraphErrors *HI_sub_pp_dWidth_chi2_raw_nominal = getdiffplot(HI_pp_dWidth_chi2_raw_nominal);
    TGraphErrors *HI_sub_pp_dWidth_chi2_raw_tnpU = getdiffplot(HI_pp_dWidth_chi2_raw_tnpU);
    TGraphErrors *HI_sub_pp_dWidth_chi2_raw_tnpD = getdiffplot(HI_pp_dWidth_chi2_raw_tnpD);
    TGraphErrors *HI_sub_pp_dWidth_chi2_raw_acooff = getdiffplot(HI_pp_dWidth_chi2_raw_acooff);
    TGraphErrors *HI_sub_pp_dWidth_chi2_raw_nominal_no_bk = getdiffplot(HI_pp_dWidth_chi2_raw_nominal_no_bk);
    TGraphErrors *HI_sub_pp_dWidth_chi2_raw_nominal_uniform_rebin = getdiffplot(HI_pp_dWidth_chi2_raw_nominal_uniform_rebin);
    TGraphErrors *HI_sub_pp_dWidth_chi2_raw_nominal_mass_range = getdiffplot(HI_pp_dWidth_chi2_raw_nominal_mass_range);
    TGraphErrors *HI_sub_pp_dWidth_chi2_raw_HF_up = getdiffplot(HI_pp_dWidth_chi2_raw_HF_up);
    TGraphErrors *HI_sub_pp_dWidth_chi2_raw_HF_down = getdiffplot(HI_pp_dWidth_chi2_raw_HF_down);

    cosmetic(HI_pp_dWidth_chi2_raw_nominal, HI_pp_dWidth_chi2_raw_tnpU, HI_pp_dWidth_chi2_raw_tnpD, HI_pp_dWidth_chi2_raw_acooff, HI_pp_dWidth_chi2_raw_nominal_no_bk,
             HI_pp_dWidth_chi2_raw_nominal_uniform_rebin, HI_pp_dWidth_chi2_raw_nominal_mass_range, HI_pp_dWidth_chi2_raw_HF_up, HI_pp_dWidth_chi2_raw_HF_down, 0);
    draw(c_FA_dW, HI_pp_dWidth_chi2_raw_nominal, HI_pp_dWidth_chi2_raw_tnpU, HI_pp_dWidth_chi2_raw_tnpD, HI_pp_dWidth_chi2_raw_acooff, HI_pp_dWidth_chi2_raw_nominal_no_bk,
         HI_pp_dWidth_chi2_raw_nominal_uniform_rebin, HI_pp_dWidth_chi2_raw_nominal_mass_range, HI_pp_dWidth_chi2_raw_HF_up, HI_pp_dWidth_chi2_raw_HF_down, 0, 0);

    cosmetic(HI_sub_pp_dWidth_chi2_raw_nominal, HI_sub_pp_dWidth_chi2_raw_tnpU, HI_sub_pp_dWidth_chi2_raw_tnpD, HI_sub_pp_dWidth_chi2_raw_acooff, HI_sub_pp_dWidth_chi2_raw_nominal_no_bk,
             HI_sub_pp_dWidth_chi2_raw_nominal_uniform_rebin, HI_sub_pp_dWidth_chi2_raw_nominal_mass_range, HI_sub_pp_dM_chi2_raw_HF_up, HI_sub_pp_dM_chi2_raw_HF_down, 0, 1);
    drawdiff(c_FA_dW_diff, HI_sub_pp_dWidth_chi2_raw_nominal, HI_sub_pp_dWidth_chi2_raw_tnpU, HI_sub_pp_dWidth_chi2_raw_tnpD, HI_sub_pp_dWidth_chi2_raw_acooff, HI_sub_pp_dWidth_chi2_raw_nominal_no_bk,
             HI_sub_pp_dWidth_chi2_raw_nominal_uniform_rebin, HI_sub_pp_dWidth_chi2_raw_nominal_mass_range, HI_sub_pp_dM_chi2_raw_HF_up, HI_sub_pp_dM_chi2_raw_HF_down, 0, 0);

    // dM, eta, HI
    TGraphErrors *HI_pp_dM_chi2_eta_nominal = mergeGraphs(HI_dM_chi2_eta_nominal, pp_dM_chi2_eta_nominal);
    TGraphErrors *HI_pp_dM_chi2_eta_tnpU = mergeGraphs(HI_dM_chi2_eta_tnpU, pp_dM_chi2_eta_tnpU);
    TGraphErrors *HI_pp_dM_chi2_eta_tnpD = mergeGraphs(HI_dM_chi2_eta_tnpD, pp_dM_chi2_eta_tnpD);
    TGraphErrors *HI_pp_dM_chi2_eta_acooff = mergeGraphs(HI_dM_chi2_eta_acooff, pp_dM_chi2_eta_acooff);
    TGraphErrors *HI_pp_dM_chi2_eta_nominal_no_bk = mergeGraphs(HI_dM_chi2_eta_nominal_no_bk, pp_dM_chi2_eta_nominal_no_bk);
    TGraphErrors *HI_pp_dM_chi2_eta_nominal_uniform_rebin = mergeGraphs(HI_dM_chi2_eta_nominal_uniform_rebin, pp_dM_chi2_eta_nominal_uniform_rebin);
    TGraphErrors *HI_pp_dM_chi2_eta_nominal_mass_range = mergeGraphs(HI_dM_chi2_eta_nominal_mass_range, pp_dM_chi2_eta_nominal_mass_range);
    TGraphErrors *HI_pp_dM_chi2_eta_HF_up = mergeGraphs(HI_dM_chi2_eta_HF_up, pp_dM_chi2_eta_nominal);
    TGraphErrors *HI_pp_dM_chi2_eta_HF_down = mergeGraphs(HI_dM_chi2_eta_HF_down, pp_dM_chi2_eta_nominal);

    TGraphErrors *HI_sub_pp_dM_chi2_eta_nominal = getdiffplot(HI_pp_dM_chi2_eta_nominal);
    TGraphErrors *HI_sub_pp_dM_chi2_eta_tnpU = getdiffplot(HI_pp_dM_chi2_eta_tnpU);
    TGraphErrors *HI_sub_pp_dM_chi2_eta_tnpD = getdiffplot(HI_pp_dM_chi2_eta_tnpD);
    TGraphErrors *HI_sub_pp_dM_chi2_eta_acooff = getdiffplot(HI_pp_dM_chi2_eta_acooff);
    TGraphErrors *HI_sub_pp_dM_chi2_eta_nominal_no_bk = getdiffplot(HI_pp_dM_chi2_eta_nominal_no_bk);
    TGraphErrors *HI_sub_pp_dM_chi2_eta_nominal_uniform_rebin = getdiffplot(HI_pp_dM_chi2_eta_nominal_uniform_rebin);
    TGraphErrors *HI_sub_pp_dM_chi2_eta_nominal_mass_range = getdiffplot(HI_pp_dM_chi2_eta_nominal_mass_range);
    TGraphErrors *HI_sub_pp_dM_chi2_eta_HF_up = getdiffplot(HI_pp_dM_chi2_eta_HF_up);
    TGraphErrors *HI_sub_pp_dM_chi2_eta_HF_down = getdiffplot(HI_pp_dM_chi2_eta_HF_down);

    cosmetic(HI_pp_dM_chi2_eta_nominal, HI_pp_dM_chi2_eta_tnpU, HI_pp_dM_chi2_eta_tnpD, HI_pp_dM_chi2_eta_acooff, HI_pp_dM_chi2_eta_nominal_no_bk,
             HI_pp_dM_chi2_eta_nominal_uniform_rebin, HI_pp_dM_chi2_eta_nominal_mass_range, HI_pp_dM_chi2_eta_HF_up, HI_pp_dM_chi2_eta_HF_down, 1);
    draw(c_Eta_dM, HI_pp_dM_chi2_eta_nominal, HI_pp_dM_chi2_eta_tnpU, HI_pp_dM_chi2_eta_tnpD, HI_pp_dM_chi2_eta_acooff, HI_pp_dM_chi2_eta_nominal_no_bk,
         HI_pp_dM_chi2_eta_nominal_uniform_rebin, HI_pp_dM_chi2_eta_nominal_mass_range, HI_pp_dM_chi2_eta_HF_up, HI_pp_dM_chi2_eta_HF_down, 1, 1);

    cosmetic(HI_sub_pp_dM_chi2_eta_nominal, HI_sub_pp_dM_chi2_eta_tnpU, HI_sub_pp_dM_chi2_eta_tnpD, HI_sub_pp_dM_chi2_eta_acooff, HI_sub_pp_dM_chi2_eta_nominal_no_bk,
             HI_sub_pp_dM_chi2_eta_nominal_uniform_rebin, HI_sub_pp_dM_chi2_eta_nominal_mass_range, HI_sub_pp_dM_chi2_eta_HF_up, HI_sub_pp_dM_chi2_eta_HF_down, 1, 1);
    drawdiff(c_Eta_dM_diff, HI_sub_pp_dM_chi2_eta_nominal, HI_sub_pp_dM_chi2_eta_tnpU, HI_sub_pp_dM_chi2_eta_tnpD, HI_sub_pp_dM_chi2_eta_acooff, HI_sub_pp_dM_chi2_eta_nominal_no_bk,
             HI_sub_pp_dM_chi2_eta_nominal_uniform_rebin, HI_sub_pp_dM_chi2_eta_nominal_mass_range, HI_sub_pp_dM_chi2_eta_HF_up, HI_sub_pp_dM_chi2_eta_HF_down, 1, 1);

    // dW, eta, HI
    TGraphErrors *HI_pp_dWidth_chi2_eta_nominal = mergeGraphs(HI_dWidth_chi2_eta_nominal, pp_dWidth_chi2_eta_nominal);
    TGraphErrors *HI_pp_dWidth_chi2_eta_tnpU = mergeGraphs(HI_dWidth_chi2_eta_tnpU, pp_dWidth_chi2_eta_tnpU);
    TGraphErrors *HI_pp_dWidth_chi2_eta_tnpD = mergeGraphs(HI_dWidth_chi2_eta_tnpD, pp_dWidth_chi2_eta_tnpD);
    TGraphErrors *HI_pp_dWidth_chi2_eta_acooff = mergeGraphs(HI_dWidth_chi2_eta_acooff, pp_dWidth_chi2_eta_acooff);
    TGraphErrors *HI_pp_dWidth_chi2_eta_nominal_no_bk = mergeGraphs(HI_dWidth_chi2_eta_nominal_no_bk, pp_dWidth_chi2_eta_nominal_no_bk);
    TGraphErrors *HI_pp_dWidth_chi2_eta_nominal_uniform_rebin = mergeGraphs(HI_dWidth_chi2_eta_nominal_uniform_rebin, pp_dWidth_chi2_eta_nominal_uniform_rebin);
    TGraphErrors *HI_pp_dWidth_chi2_eta_nominal_mass_range = mergeGraphs(HI_dWidth_chi2_eta_nominal_mass_range, pp_dWidth_chi2_eta_nominal_mass_range);
    TGraphErrors *HI_pp_dWidth_chi2_eta_HF_up = mergeGraphs(HI_dWidth_chi2_eta_HF_up, pp_dWidth_chi2_eta_nominal);
    TGraphErrors *HI_pp_dWidth_chi2_eta_HF_down = mergeGraphs(HI_dWidth_chi2_eta_HF_down, pp_dWidth_chi2_eta_nominal);

    TGraphErrors *HI_sub_pp_dWidth_chi2_eta_nominal = getdiffplot(HI_pp_dWidth_chi2_eta_nominal);
    TGraphErrors *HI_sub_pp_dWidth_chi2_eta_tnpU = getdiffplot(HI_pp_dWidth_chi2_eta_tnpU);
    TGraphErrors *HI_sub_pp_dWidth_chi2_eta_tnpD = getdiffplot(HI_pp_dWidth_chi2_eta_tnpD);
    TGraphErrors *HI_sub_pp_dWidth_chi2_eta_acooff = getdiffplot(HI_pp_dWidth_chi2_eta_acooff);
    TGraphErrors *HI_sub_pp_dWidth_chi2_eta_nominal_no_bk = getdiffplot(HI_pp_dWidth_chi2_eta_nominal_no_bk);
    TGraphErrors *HI_sub_pp_dWidth_chi2_eta_nominal_uniform_rebin = getdiffplot(HI_pp_dWidth_chi2_eta_nominal_uniform_rebin);
    TGraphErrors *HI_sub_pp_dWidth_chi2_eta_nominal_mass_range = getdiffplot(HI_pp_dWidth_chi2_eta_nominal_mass_range);
    TGraphErrors *HI_sub_pp_dWidth_chi2_eta_HF_up = getdiffplot(HI_pp_dWidth_chi2_eta_HF_up);
    TGraphErrors *HI_sub_pp_dWidth_chi2_eta_HF_down = getdiffplot(HI_pp_dWidth_chi2_eta_HF_down);

    cosmetic(HI_pp_dWidth_chi2_eta_nominal, HI_pp_dWidth_chi2_eta_tnpU, HI_pp_dWidth_chi2_eta_tnpD, HI_pp_dWidth_chi2_eta_acooff, HI_pp_dWidth_chi2_eta_nominal_no_bk,
             HI_pp_dWidth_chi2_eta_nominal_uniform_rebin, HI_pp_dWidth_chi2_eta_nominal_mass_range, HI_pp_dWidth_chi2_eta_HF_up, HI_pp_dWidth_chi2_eta_HF_down, 0);
    draw(c_Eta_dW, HI_pp_dWidth_chi2_eta_nominal, HI_pp_dWidth_chi2_eta_tnpU, HI_pp_dWidth_chi2_eta_tnpD, HI_pp_dWidth_chi2_eta_acooff, HI_pp_dWidth_chi2_eta_nominal_no_bk,
         HI_pp_dWidth_chi2_eta_nominal_uniform_rebin, HI_pp_dWidth_chi2_eta_nominal_mass_range, HI_pp_dWidth_chi2_eta_HF_up, HI_pp_dWidth_chi2_eta_HF_down, 0, 1);

    cosmetic(HI_sub_pp_dWidth_chi2_eta_nominal, HI_sub_pp_dWidth_chi2_eta_tnpU, HI_sub_pp_dWidth_chi2_eta_tnpD, HI_sub_pp_dWidth_chi2_eta_acooff, HI_sub_pp_dWidth_chi2_eta_nominal_no_bk,
             HI_sub_pp_dWidth_chi2_eta_nominal_uniform_rebin, HI_sub_pp_dWidth_chi2_eta_nominal_mass_range, HI_sub_pp_dWidth_chi2_eta_HF_up, HI_sub_pp_dWidth_chi2_eta_HF_down, 0, 1);
    drawdiff(c_Eta_dW_diff, HI_sub_pp_dWidth_chi2_eta_nominal, HI_sub_pp_dWidth_chi2_eta_tnpU, HI_sub_pp_dWidth_chi2_eta_tnpD, HI_sub_pp_dWidth_chi2_eta_acooff, HI_sub_pp_dWidth_chi2_eta_nominal_no_bk,
             HI_sub_pp_dWidth_chi2_eta_nominal_uniform_rebin, HI_sub_pp_dWidth_chi2_eta_nominal_mass_range, HI_sub_pp_dWidth_chi2_eta_HF_up, HI_sub_pp_dWidth_chi2_eta_HF_down, 0, 1);

    // dM, FA, pp

    cosmeticpp(pp_dM_chi2_raw_nominal, pp_dM_chi2_raw_tnpU, pp_dM_chi2_raw_tnpD, pp_dM_chi2_raw_acooff, pp_dM_chi2_raw_nominal_no_bk, pp_dM_chi2_raw_nominal_uniform_rebin, pp_dM_chi2_raw_nominal_mass_range, 1);

    draw(c_pp_FA_dM, pp_dM_chi2_raw_nominal, pp_dM_chi2_raw_tnpU, pp_dM_chi2_raw_tnpD, pp_dM_chi2_raw_acooff, pp_dM_chi2_raw_nominal_no_bk, pp_dM_chi2_raw_nominal_uniform_rebin, pp_dM_chi2_raw_nominal_mass_range, nullptr, nullptr, 1, 0);

    // dW, FA, pp

    cosmeticpp(pp_dWidth_chi2_raw_nominal, pp_dWidth_chi2_raw_tnpU, pp_dWidth_chi2_raw_tnpD, pp_dWidth_chi2_raw_acooff, pp_dWidth_chi2_raw_nominal_no_bk, pp_dWidth_chi2_raw_nominal_uniform_rebin, pp_dWidth_chi2_raw_nominal_mass_range, 0);

    draw(c_pp_FA_dW, pp_dWidth_chi2_raw_nominal, pp_dWidth_chi2_raw_tnpU, pp_dWidth_chi2_raw_tnpD, pp_dWidth_chi2_raw_acooff, pp_dWidth_chi2_raw_nominal_no_bk, pp_dWidth_chi2_raw_nominal_uniform_rebin, pp_dWidth_chi2_raw_nominal_mass_range, nullptr, nullptr, 0, 0);

    // dM, eta, pp

    cosmeticpp(pp_dM_chi2_eta_nominal, pp_dM_chi2_eta_tnpU, pp_dM_chi2_eta_tnpD, pp_dM_chi2_eta_acooff, pp_dM_chi2_eta_nominal_no_bk, pp_dM_chi2_eta_nominal_uniform_rebin, pp_dM_chi2_eta_nominal_mass_range, 1);

    draw(c_pp_Eta_dM, pp_dM_chi2_eta_nominal, pp_dM_chi2_eta_tnpU, pp_dM_chi2_eta_tnpD, pp_dM_chi2_eta_acooff, pp_dM_chi2_eta_nominal_no_bk, pp_dM_chi2_eta_nominal_uniform_rebin, pp_dM_chi2_eta_nominal_mass_range, nullptr, nullptr, 1, 1);

    // dW, eta, pp

    cosmeticpp(pp_dWidth_chi2_eta_nominal, pp_dWidth_chi2_eta_tnpU, pp_dWidth_chi2_eta_tnpD, pp_dWidth_chi2_eta_acooff, pp_dWidth_chi2_eta_nominal_no_bk, pp_dWidth_chi2_eta_nominal_uniform_rebin, pp_dWidth_chi2_eta_nominal_mass_range, 0);

    draw(c_pp_Eta_dW, pp_dWidth_chi2_eta_nominal, pp_dWidth_chi2_eta_tnpU, pp_dWidth_chi2_eta_tnpD, pp_dWidth_chi2_eta_acooff, pp_dWidth_chi2_eta_nominal_no_bk, pp_dWidth_chi2_eta_nominal_uniform_rebin, pp_dWidth_chi2_eta_nominal_mass_range, nullptr, nullptr, 0, 1);

    c_FA_dM->SaveAs("./systematic/FA_dM.png");
    c_FA_dM_diff->SaveAs("./systematic/FA_dM_diff.png");
    c_FA_dW->SaveAs("./systematic/FA_dW.png");
    c_FA_dW_diff->SaveAs("./systematic/FA_dW_diff.png");
    c_Eta_dM->SaveAs("./systematic/Eta_dM.png");
    c_Eta_dM_diff->SaveAs("./systematic/Eta_dM_diff.png");
    c_Eta_dW->SaveAs("./systematic/Eta_dW.png");
    c_Eta_dW_diff->SaveAs("./systematic/Eta_dW_diff.png");

    c_pp_FA_dM->SaveAs("./systematic/pp_FA_dM.png");
    c_pp_FA_dW->SaveAs("./systematic/pp_FA_dW.png");
    c_pp_Eta_dM->SaveAs("./systematic/pp_Eta_dM.png");
    c_pp_Eta_dW->SaveAs("./systematic/pp_Eta_dW.png");

    TFile *nominal_syst = new TFile("./nominal_syst_tgraph.root", "RECREATE");
    nominal_syst->cd();

    HI_sub_pp_dM_chi2_raw_nominal->Write("HI_sub_pp_dM_chi2_raw_nominal", 2);
    HI_sub_pp_dM_chi2_raw_tnpU->Write("HI_sub_pp_dM_chi2_raw_tnpU", 2);
    HI_sub_pp_dM_chi2_raw_tnpD->Write("HI_sub_pp_dM_chi2_raw_tnpD", 2);
    HI_sub_pp_dM_chi2_raw_acooff->Write("HI_sub_pp_dM_chi2_raw_acooff", 2);
    HI_sub_pp_dM_chi2_raw_nominal_no_bk->Write("HI_sub_pp_dM_chi2_raw_nominal_no_bk", 2);
    HI_sub_pp_dM_chi2_raw_nominal_uniform_rebin->Write("HI_sub_pp_dM_chi2_raw_nominal_uniform_rebin", 2);
    HI_sub_pp_dM_chi2_raw_nominal_mass_range->Write("HI_sub_pp_dM_chi2_raw_nominal_mass_range", 2);
    HI_sub_pp_dM_chi2_raw_HF_up->Write("HI_sub_pp_dM_chi2_raw_HF_up", 2);
    HI_sub_pp_dM_chi2_raw_HF_down->Write("HI_sub_pp_dM_chi2_raw_HF_down", 2);

    HI_pp_dM_chi2_raw_nominal->Write("HI_pp_dM_chi2_raw_nominal", 2);
    HI_pp_dM_chi2_raw_tnpU->Write("HI_pp_dM_chi2_raw_tnpU", 2);
    HI_pp_dM_chi2_raw_tnpD->Write("HI_pp_dM_chi2_raw_tnpD", 2);
    HI_pp_dM_chi2_raw_acooff->Write("HI_pp_dM_chi2_raw_acooff", 2);
    HI_pp_dM_chi2_raw_nominal_no_bk->Write("HI_pp_dM_chi2_raw_nominal_no_bk", 2);
    HI_pp_dM_chi2_raw_nominal_uniform_rebin->Write("HI_pp_dM_chi2_raw_nominal_uniform_rebin", 2);
    HI_pp_dM_chi2_raw_nominal_mass_range->Write("HI_pp_dM_chi2_raw_nominal_mass_range", 2);
    HI_pp_dM_chi2_raw_HF_up->Write("HI_pp_dM_chi2_raw_HF_up", 2);
    HI_pp_dM_chi2_raw_HF_down->Write("HI_pp_dM_chi2_raw_HF_down", 2);

    HI_sub_pp_dWidth_chi2_raw_nominal->Write("HI_sub_pp_dWidth_chi2_raw_nominal", 2);
    HI_sub_pp_dWidth_chi2_raw_tnpU->Write("HI_sub_pp_dWidth_chi2_raw_tnpU", 2);
    HI_sub_pp_dWidth_chi2_raw_tnpD->Write("HI_sub_pp_dWidth_chi2_raw_tnpD", 2);
    HI_sub_pp_dWidth_chi2_raw_acooff->Write("HI_sub_pp_dWidth_chi2_raw_acooff", 2);
    HI_sub_pp_dWidth_chi2_raw_nominal_no_bk->Write("HI_sub_pp_dWidth_chi2_raw_nominal_no_bk", 2);
    HI_sub_pp_dWidth_chi2_raw_nominal_uniform_rebin->Write("HI_sub_pp_dWidth_chi2_raw_nominal_uniform_rebin", 2);
    HI_sub_pp_dWidth_chi2_raw_nominal_mass_range->Write("HI_sub_pp_dWidth_chi2_raw_nominal_mass_range", 2);
    HI_sub_pp_dWidth_chi2_raw_HF_up->Write("HI_sub_pp_dWidth_chi2_raw_HF_up", 2);
    HI_sub_pp_dWidth_chi2_raw_HF_down->Write("HI_sub_pp_dWidth_chi2_raw_HF_down", 2);

    HI_pp_dWidth_chi2_raw_nominal->Write("HI_pp_dWidth_chi2_raw_nominal", 2);
    HI_pp_dWidth_chi2_raw_tnpU->Write("HI_pp_dWidth_chi2_raw_tnpU", 2);
    HI_pp_dWidth_chi2_raw_tnpD->Write("HI_pp_dWidth_chi2_raw_tnpD", 2);
    HI_pp_dWidth_chi2_raw_acooff->Write("HI_pp_dWidth_chi2_raw_acooff", 2);
    HI_pp_dWidth_chi2_raw_nominal_no_bk->Write("HI_pp_dWidth_chi2_raw_nominal_no_bk", 2);
    HI_pp_dWidth_chi2_raw_nominal_uniform_rebin->Write("HI_pp_dWidth_chi2_raw_nominal_uniform_rebin", 2);
    HI_pp_dWidth_chi2_raw_nominal_mass_range->Write("HI_pp_dWidth_chi2_raw_nominal_mass_range", 2);
    HI_pp_dWidth_chi2_raw_HF_up->Write("HI_pp_dWidth_chi2_raw_HF_up", 2);
    HI_pp_dWidth_chi2_raw_HF_down->Write("HI_pp_dWidth_chi2_raw_HF_down", 2);

    HI_sub_pp_dM_chi2_eta_nominal->Write("HI_sub_pp_dM_chi2_eta_nominal", 2);
    HI_sub_pp_dM_chi2_eta_tnpU->Write("HI_sub_pp_dM_chi2_eta_tnpU", 2);
    HI_sub_pp_dM_chi2_eta_tnpD->Write("HI_sub_pp_dM_chi2_eta_tnpD", 2);
    HI_sub_pp_dM_chi2_eta_acooff->Write("HI_sub_pp_dM_chi2_eta_acooff", 2);
    HI_sub_pp_dM_chi2_eta_nominal_no_bk->Write("HI_sub_pp_dM_chi2_eta_nominal_no_bk", 2);
    HI_sub_pp_dM_chi2_eta_nominal_uniform_rebin->Write("HI_sub_pp_dM_chi2_eta_nominal_uniform_rebin", 2);
    HI_sub_pp_dM_chi2_eta_nominal_mass_range->Write("HI_sub_pp_dM_chi2_eta_nominal_mass_range", 2);
    HI_sub_pp_dM_chi2_eta_HF_up->Write("HI_sub_pp_dM_chi2_eta_HF_up", 2);
    HI_sub_pp_dM_chi2_eta_HF_down->Write("HI_sub_pp_dM_chi2_eta_HF_down", 2);

    HI_pp_dM_chi2_eta_nominal->Write("HI_pp_dM_chi2_eta_nominal", 2);
    HI_pp_dM_chi2_eta_tnpU->Write("HI_pp_dM_chi2_eta_tnpU", 2);
    HI_pp_dM_chi2_eta_tnpD->Write("HI_pp_dM_chi2_eta_tnpD", 2);
    HI_pp_dM_chi2_eta_acooff->Write("HI_pp_dM_chi2_eta_acooff", 2);
    HI_pp_dM_chi2_eta_nominal_no_bk->Write("HI_pp_dM_chi2_eta_nominal_no_bk", 2);
    HI_pp_dM_chi2_eta_nominal_uniform_rebin->Write("HI_pp_dM_chi2_eta_nominal_uniform_rebin", 2);
    HI_pp_dM_chi2_eta_nominal_mass_range->Write("HI_pp_dM_chi2_eta_nominal_mass_range", 2);
    HI_pp_dM_chi2_eta_HF_up->Write("HI_pp_dM_chi2_eta_HF_up", 2);
    HI_pp_dM_chi2_eta_HF_down->Write("HI_pp_dM_chi2_eta_HF_down", 2);

    HI_sub_pp_dWidth_chi2_eta_nominal->Write("HI_sub_pp_dWidth_chi2_eta_nominal", 2);
    HI_sub_pp_dWidth_chi2_eta_tnpU->Write("HI_sub_pp_dWidth_chi2_eta_tnpU", 2);
    HI_sub_pp_dWidth_chi2_eta_tnpD->Write("HI_sub_pp_dWidth_chi2_eta_tnpD", 2);
    HI_sub_pp_dWidth_chi2_eta_acooff->Write("HI_sub_pp_dWidth_chi2_eta_acooff", 2);
    HI_sub_pp_dWidth_chi2_eta_nominal_no_bk->Write("HI_sub_pp_dWidth_chi2_eta_nominal_no_bk", 2);
    HI_sub_pp_dWidth_chi2_eta_nominal_uniform_rebin->Write("HI_sub_pp_dWidth_chi2_eta_nominal_uniform_rebin", 2);
    HI_sub_pp_dWidth_chi2_eta_nominal_mass_range->Write("HI_sub_pp_dWidth_chi2_eta_nominal_mass_range", 2);
    HI_sub_pp_dWidth_chi2_eta_HF_up->Write("HI_sub_pp_dWidth_chi2_eta_HF_up", 2);
    HI_sub_pp_dWidth_chi2_eta_HF_down->Write("HI_sub_pp_dWidth_chi2_eta_HF_down", 2);

    HI_pp_dWidth_chi2_eta_nominal->Write("HI_pp_dWidth_chi2_eta_nominal", 2);
    HI_pp_dWidth_chi2_eta_tnpU->Write("HI_pp_dWidth_chi2_eta_tnpU", 2);
    HI_pp_dWidth_chi2_eta_tnpD->Write("HI_pp_dWidth_chi2_eta_tnpD", 2);
    HI_pp_dWidth_chi2_eta_acooff->Write("HI_pp_dWidth_chi2_eta_acooff", 2);
    HI_pp_dWidth_chi2_eta_nominal_no_bk->Write("HI_pp_dWidth_chi2_eta_nominal_no_bk", 2);
    HI_pp_dWidth_chi2_eta_nominal_uniform_rebin->Write("HI_pp_dWidth_chi2_eta_nominal_uniform_rebin", 2);
    HI_pp_dWidth_chi2_eta_nominal_mass_range->Write("HI_pp_dWidth_chi2_eta_nominal_mass_range", 2);
    HI_pp_dWidth_chi2_eta_HF_up->Write("HI_pp_dWidth_chi2_eta_HF_up", 2);
    HI_pp_dWidth_chi2_eta_HF_down->Write("HI_pp_dWidth_chi2_eta_HF_down", 2);
}