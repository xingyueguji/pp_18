#include "tdrStyle.C"
#include "CMS_lumi.C"
void cosmetic(TGraphErrors *g1, TGraphErrors *g2, TGraphErrors *g3, TGraphErrors *g4, TGraphErrors *g5, bool isdM)
{
    g1->SetMarkerStyle(20); // Circle marker
    g1->SetMarkerSize(1.2);
    g1->SetLineWidth(2);
    g1->SetLineColor(kBlack);
    g1->SetMarkerColor(kBlack);
    g1->SetTitle("");
    if (isdM)
    {
        g1->GetYaxis()->SetTitle("#DeltaM (GeV)");
    }
    else
    {
        g1->GetYaxis()->SetTitle("#DeltaWidth (GeV)");
    }

    g1->GetXaxis()->SetLimits(0.5, 5.5);
    g1->GetXaxis()->SetRangeUser(0.5, 5.5);
    g1->GetXaxis()->SetNdivisions(5, 0, 0, kFALSE);
    g1->GetXaxis()->SetLabelSize(0);

    g1->GetXaxis()->ChangeLabel(1, 0, 0.04, 11, -1, -1, "     0-100%");
    g1->GetXaxis()->ChangeLabel(2, 0, 0.04, 11, -1, -1, "      0-10%");
    g1->GetXaxis()->ChangeLabel(3, 0, 0.04, 11, -1, -1, "     10-20%");
    g1->GetXaxis()->ChangeLabel(4, 0, 0.04, 11, -1, -1, "     20-30%");
    g1->GetXaxis()->ChangeLabel(5, 0, 0.04, 11, -1, -1, "    30-100%");
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

    // Define small shifts for each graph
    double shift2 = 0.1;  // Small right shift for g2
    double shift3 = -0.1; // Slightly more shift for g3
    double shift4 = 0.2;  // More shift for g4
    double shift5 = -0.2; // Largest shift for g5

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

    TGraphErrors *HI_dWidth_chi2_raw_nominal = (TGraphErrors *)f1->Get("HI_dWidth_chi2_raw_nominal");
    TGraphErrors *HI_dWidth_chi2_raw_tnpU = (TGraphErrors *)f1->Get("HI_dWidth_chi2_raw_tnpU");
    TGraphErrors *HI_dWidth_chi2_raw_tnpD = (TGraphErrors *)f1->Get("HI_dWidth_chi2_raw_tnpD");
    TGraphErrors *HI_dWidth_chi2_raw_acooff = (TGraphErrors *)f1->Get("HI_dWidth_chi2_raw_acooff");
    TGraphErrors *HI_dWidth_chi2_raw_nominal_no_bk = (TGraphErrors *)f1->Get("HI_dWidth_chi2_raw_nominal_no_bk");

    TGraphErrors *HI_dM_chi2_eta_nominal = (TGraphErrors *)f1->Get("HI_dM_chi2_eta_nominal");
    TGraphErrors *HI_dM_chi2_eta_tnpU = (TGraphErrors *)f1->Get("HI_dM_chi2_eta_tnpU");
    TGraphErrors *HI_dM_chi2_eta_tnpD = (TGraphErrors *)f1->Get("HI_dM_chi2_eta_tnpD");
    TGraphErrors *HI_dM_chi2_eta_acooff = (TGraphErrors *)f1->Get("HI_dM_chi2_eta_acooff");
    TGraphErrors *HI_dM_chi2_eta_nominal_no_bk = (TGraphErrors *)f1->Get("HI_dM_chi2_eta_nominal_no_bk");

    TGraphErrors *HI_dWidth_chi2_eta_nominal = (TGraphErrors *)f1->Get("HI_dWidth_chi2_eta_nominal");
    TGraphErrors *HI_dWidth_chi2_eta_tnpU = (TGraphErrors *)f1->Get("HI_dWidth_chi2_eta_tnpU");
    TGraphErrors *HI_dWidth_chi2_eta_tnpD = (TGraphErrors *)f1->Get("HI_dWidth_chi2_eta_tnpD");
    TGraphErrors *HI_dWidth_chi2_eta_acooff = (TGraphErrors *)f1->Get("HI_dWidth_chi2_eta_acooff");
    TGraphErrors *HI_dWidth_chi2_eta_nominal_no_bk = (TGraphErrors *)f1->Get("HI_dWidth_chi2_eta_nominal_no_bk");

    TGraphErrors *pp_dM_chi2_raw_nominal = (TGraphErrors *)f1->Get("pp_dM_chi2_raw_nominal");
    TGraphErrors *pp_dM_chi2_raw_tnpU = (TGraphErrors *)f1->Get("pp_dM_chi2_raw_tnpU");
    TGraphErrors *pp_dM_chi2_raw_tnpD = (TGraphErrors *)f1->Get("pp_dM_chi2_raw_tnpD");
    TGraphErrors *pp_dM_chi2_raw_acooff = (TGraphErrors *)f1->Get("pp_dM_chi2_raw_acooff");
    TGraphErrors *pp_dM_chi2_raw_nominal_no_bk = (TGraphErrors *)f1->Get("pp_dM_chi2_raw_nominal_no_bk");

    TGraphErrors *pp_dWidth_chi2_raw_nominal = (TGraphErrors *)f1->Get("pp_dWidth_chi2_raw_nominal");
    TGraphErrors *pp_dWidth_chi2_raw_tnpU = (TGraphErrors *)f1->Get("pp_dWidth_chi2_raw_tnpU");
    TGraphErrors *pp_dWidth_chi2_raw_tnpD = (TGraphErrors *)f1->Get("pp_dWidth_chi2_raw_tnpD");
    TGraphErrors *pp_dWidth_chi2_raw_acooff = (TGraphErrors *)f1->Get("pp_dWidth_chi2_raw_acooff");
    TGraphErrors *pp_dWidth_chi2_raw_nominal_no_bk = (TGraphErrors *)f1->Get("pp_dWidth_chi2_raw_nominal_no_bk");

    TGraphErrors *pp_dM_chi2_eta_nominal = (TGraphErrors *)f1->Get("pp_dM_chi2_eta_nominal");
    TGraphErrors *pp_dM_chi2_eta_tnpU = (TGraphErrors *)f1->Get("pp_dM_chi2_eta_tnpU");
    TGraphErrors *pp_dM_chi2_eta_tnpD = (TGraphErrors *)f1->Get("pp_dM_chi2_eta_tnpD");
    TGraphErrors *pp_dM_chi2_eta_acooff = (TGraphErrors *)f1->Get("pp_dM_chi2_eta_acooff");
    TGraphErrors *pp_dM_chi2_eta_nominal_no_bk = (TGraphErrors *)f1->Get("pp_dM_chi2_eta_nominal_no_bk");

    TGraphErrors *pp_dWidth_chi2_eta_nominal = (TGraphErrors *)f1->Get("pp_dWidth_chi2_eta_nominal");
    TGraphErrors *pp_dWidth_chi2_eta_tnpU = (TGraphErrors *)f1->Get("pp_dWidth_chi2_eta_tnpU");
    TGraphErrors *pp_dWidth_chi2_eta_tnpD = (TGraphErrors *)f1->Get("pp_dWidth_chi2_eta_tnpD");
    TGraphErrors *pp_dWidth_chi2_eta_acooff = (TGraphErrors *)f1->Get("pp_dWidth_chi2_eta_acooff");
    TGraphErrors *pp_dWidth_chi2_eta_nominal_no_bk = (TGraphErrors *)f1->Get("pp_dWidth_chi2_eta_nominal_no_bk");

    TCanvas *c_FA_dM = new TCanvas("c_FA_dM", "", 800, 600);
    TCanvas *c_FA_dW = new TCanvas("c_FA_dW", "", 800, 600);

    TCanvas *c_Eta_dM = new TCanvas("c_Eta_dM", "", 800, 600);
    TCanvas *c_Eta_dW = new TCanvas("c_Eta_dW", "", 800, 600);

    TCanvas *c_pp_FA_dM = new TCanvas("c_pp_FA_dM", "", 800, 600);
    TCanvas *c_pp_FA_dW = new TCanvas("c_pp_FA_dW", "", 800, 600);

    TCanvas *c_pp_Eta_dM = new TCanvas("c_pp_Eta_dM", "", 800, 600);
    TCanvas *c_pp_Eta_dW = new TCanvas("c_pp_Eta_dW", "", 800, 600);

    c_FA_dM->cd();

    cosmetic(HI_dM_chi2_raw_nominal, HI_dM_chi2_raw_tnpU, HI_dM_chi2_raw_tnpD, HI_dM_chi2_raw_acooff, HI_dM_chi2_raw_nominal_no_bk,1);
    HI_dM_chi2_raw_nominal->Draw("AP");
    HI_dM_chi2_raw_tnpU->Draw("P SAME");
    HI_dM_chi2_raw_tnpD->Draw("P SAME");
    HI_dM_chi2_raw_acooff->Draw("P SAME");
    HI_dM_chi2_raw_nominal_no_bk->Draw("P SAME");

    TLegend *legend = new TLegend(0.2, 0.7, 0.4, 0.9);
    legend->AddEntry(HI_dM_chi2_raw_nominal, "Nominal", "PL");
    legend->AddEntry(HI_dM_chi2_raw_tnpU, "tnpU", "PL");
    legend->AddEntry(HI_dM_chi2_raw_tnpD, "tnpD", "PL");
    legend->AddEntry(HI_dM_chi2_raw_acooff, "acooff", "PL");
    legend->AddEntry(HI_dM_chi2_raw_nominal_no_bk, "bk_off", "PL");
    legend->SetTextSize(0.03);              // Set text size
    legend->SetTextFont(42);                // Use a modern, clean font
    legend->SetBorderSize(0);               // Set border size (0 for no border)
    legend->SetLineColor(kBlack);           // Border color (if any)
    legend->SetLineWidth(0);                // Border line width
    legend->SetFillColorAlpha(kWhite, 0.1); // Background color with transparency
    legend->Draw("SAME");

    c_FA_dW->cd();
    cosmetic(HI_dWidth_chi2_raw_nominal, HI_dWidth_chi2_raw_tnpU, HI_dWidth_chi2_raw_tnpD, HI_dWidth_chi2_raw_acooff, HI_dWidth_chi2_raw_nominal_no_bk,0);
    HI_dWidth_chi2_raw_nominal->Draw("AP");
    HI_dWidth_chi2_raw_tnpU->Draw("P SAME");
    HI_dWidth_chi2_raw_tnpD->Draw("P SAME");
    HI_dWidth_chi2_raw_acooff->Draw("P SAME");
    HI_dWidth_chi2_raw_nominal_no_bk->Draw("P SAME");

    TLegend *legend1 = new TLegend(0.2, 0.7, 0.4, 0.9);
    legend1->AddEntry(HI_dWidth_chi2_raw_nominal, "Nominal", "PL");
    legend1->AddEntry(HI_dWidth_chi2_raw_tnpU, "tnpU", "PL");
    legend1->AddEntry(HI_dWidth_chi2_raw_tnpD, "tnpD", "PL");
    legend1->AddEntry(HI_dWidth_chi2_raw_acooff, "acooff", "PL");
    legend1->AddEntry(HI_dWidth_chi2_raw_nominal_no_bk, "bk_off", "PL");
    legend1->SetTextSize(0.03);              // Set text size
    legend1->SetTextFont(42);                // Use a modern, clean font
    legend1->SetBorderSize(0);               // Set border size (0 for no border)
    legend1->SetLineColor(kBlack);           // Border color (if any)
    legend1->SetLineWidth(0);                // Border line width
    legend1->SetFillColorAlpha(kWhite, 0.1); // Background color with transparency
    legend1->Draw("SAME");

    c_Eta_dM->cd();
    cosmetic(HI_dM_chi2_eta_nominal, HI_dM_chi2_eta_tnpU, HI_dM_chi2_eta_tnpD, HI_dM_chi2_eta_acooff, HI_dM_chi2_eta_nominal_no_bk,1);
    HI_dM_chi2_eta_nominal->Draw("AP");
    HI_dM_chi2_eta_tnpU->Draw("P SAME");
    HI_dM_chi2_eta_tnpD->Draw("P SAME");
    HI_dM_chi2_eta_acooff->Draw("P SAME");
    HI_dM_chi2_eta_nominal_no_bk->Draw("P SAME");

    TLegend *legend2 = new TLegend(0.2, 0.7, 0.4, 0.9);
    legend2->AddEntry(HI_dM_chi2_eta_nominal, "Nominal", "PL");
    legend2->AddEntry(HI_dM_chi2_eta_tnpU, "tnpU", "PL");
    legend2->AddEntry(HI_dM_chi2_eta_tnpD, "tnpD", "PL");
    legend2->AddEntry(HI_dM_chi2_eta_acooff, "acooff", "PL");
    legend2->AddEntry(HI_dM_chi2_eta_nominal_no_bk, "bk_off", "PL");
    legend2->SetTextSize(0.03);              // Set text size
    legend2->SetTextFont(42);                // Use a modern, clean font
    legend2->SetBorderSize(0);               // Set border size (0 for no border)
    legend2->SetLineColor(kBlack);           // Border color (if any)
    legend2->SetLineWidth(0);                // Border line width
    legend2->SetFillColorAlpha(kWhite, 0.1); // Background color with transparency
    legend2->Draw("SAME");

    c_Eta_dW->cd();
    cosmetic(HI_dWidth_chi2_eta_nominal, HI_dWidth_chi2_eta_tnpU, HI_dWidth_chi2_eta_tnpD, HI_dWidth_chi2_eta_acooff, HI_dWidth_chi2_eta_nominal_no_bk,0);
    HI_dWidth_chi2_eta_nominal->Draw("AP");
    HI_dWidth_chi2_eta_tnpU->Draw("P SAME");
    HI_dWidth_chi2_eta_tnpD->Draw("P SAME");
    HI_dWidth_chi2_eta_acooff->Draw("P SAME");
    HI_dWidth_chi2_eta_nominal_no_bk->Draw("P SAME");
    TLegend *legend3 = new TLegend(0.2, 0.7, 0.4, 0.9);
    legend3->AddEntry(HI_dWidth_chi2_eta_nominal, "Nominal", "PL");
    legend3->AddEntry(HI_dWidth_chi2_eta_tnpU, "tnpU", "PL");
    legend3->AddEntry(HI_dWidth_chi2_eta_tnpD, "tnpD", "PL");
    legend3->AddEntry(HI_dWidth_chi2_eta_acooff, "acooff", "PL");
    legend3->AddEntry(HI_dWidth_chi2_eta_nominal_no_bk, "bk_off", "PL");
    legend3->SetTextSize(0.03);              // Set text size
    legend3->SetTextFont(42);                // Use a modern, clean font
    legend3->SetBorderSize(0);               // Set border size (0 for no border)
    legend3->SetLineColor(kBlack);           // Border color (if any)
    legend3->SetLineWidth(0);                // Border line width
    legend3->SetFillColorAlpha(kWhite, 0.1); // Background color with transparency
    legend3->Draw("SAME");

    c_pp_FA_dM->cd();

    cosmetic(pp_dM_chi2_raw_nominal, pp_dM_chi2_raw_tnpU, pp_dM_chi2_raw_tnpD, pp_dM_chi2_raw_acooff, pp_dM_chi2_raw_nominal_no_bk,1);
    pp_dM_chi2_raw_nominal->Draw("AP");
    pp_dM_chi2_raw_tnpU->Draw("P SAME");
    pp_dM_chi2_raw_tnpD->Draw("P SAME");
    pp_dM_chi2_raw_acooff->Draw("P SAME");
    pp_dM_chi2_raw_nominal_no_bk->Draw("P SAME");

    TLegend *legend4 = new TLegend(0.2, 0.7, 0.4, 0.9);
    legend4->AddEntry(pp_dM_chi2_raw_nominal, "Nominal", "PL");
    legend4->AddEntry(pp_dM_chi2_raw_tnpU, "tnpU", "PL");
    legend4->AddEntry(pp_dM_chi2_raw_tnpD, "tnpD", "PL");
    legend4->AddEntry(pp_dM_chi2_raw_acooff, "acooff", "PL");
    legend4->AddEntry(pp_dM_chi2_raw_nominal_no_bk, "bk_off", "PL");
    legend4->SetTextSize(0.03);              // Set text size
    legend4->SetTextFont(42);                // Use a modern, clean font
    legend4->SetBorderSize(0);               // Set border size (0 for no border)
    legend4->SetLineColor(kBlack);           // Border color (if any)
    legend4->SetLineWidth(0);                // Border line width
    legend4->SetFillColorAlpha(kWhite, 0.1); // Background color with transparency
    legend4->Draw("SAME");

    c_pp_FA_dW->cd();
    cosmetic(pp_dWidth_chi2_raw_nominal, pp_dWidth_chi2_raw_tnpU, pp_dWidth_chi2_raw_tnpD, pp_dWidth_chi2_raw_acooff, pp_dWidth_chi2_raw_nominal_no_bk,0);
    pp_dWidth_chi2_raw_nominal->GetYaxis()->SetRangeUser(0.2,0.34);
    pp_dWidth_chi2_raw_nominal->Draw("AP");
    pp_dWidth_chi2_raw_tnpU->Draw("P SAME");
    pp_dWidth_chi2_raw_tnpD->Draw("P SAME");
    pp_dWidth_chi2_raw_acooff->Draw("P SAME");
    pp_dWidth_chi2_raw_nominal_no_bk->Draw("P SAME");

    TLegend *legend5 = new TLegend(0.2, 0.7, 0.4, 0.9);
    legend5->AddEntry(pp_dWidth_chi2_raw_nominal, "Nominal", "PL");
    legend5->AddEntry(pp_dWidth_chi2_raw_tnpU, "tnpU", "PL");
    legend5->AddEntry(pp_dWidth_chi2_raw_tnpD, "tnpD", "PL");
    legend5->AddEntry(pp_dWidth_chi2_raw_acooff, "acooff", "PL");
    legend5->AddEntry(pp_dWidth_chi2_raw_nominal_no_bk, "bk_off", "PL");
    legend5->SetTextSize(0.03);              // Set text size
    legend5->SetTextFont(42);                // Use a modern, clean font
    legend5->SetBorderSize(0);               // Set border size (0 for no border)
    legend5->SetLineColor(kBlack);           // Border color (if any)
    legend5->SetLineWidth(0);                // Border line width
    legend5->SetFillColorAlpha(kWhite, 0.1); // Background color with transparency
    legend5->Draw("SAME");

    c_pp_Eta_dM->cd();
    cosmetic(pp_dM_chi2_eta_nominal, pp_dM_chi2_eta_tnpU, pp_dM_chi2_eta_tnpD, pp_dM_chi2_eta_acooff, pp_dM_chi2_eta_nominal_no_bk,1);
    pp_dM_chi2_eta_nominal->Draw("AP");
    pp_dM_chi2_eta_tnpU->Draw("P SAME");
    pp_dM_chi2_eta_tnpD->Draw("P SAME");
    pp_dM_chi2_eta_acooff->Draw("P SAME");
    pp_dM_chi2_eta_nominal_no_bk->Draw("P SAME");

    TLegend *legend6 = new TLegend(0.2, 0.7, 0.4, 0.9);
    legend6->AddEntry(pp_dM_chi2_eta_nominal, "Nominal", "PL");
    legend6->AddEntry(pp_dM_chi2_eta_tnpU, "tnpU", "PL");
    legend6->AddEntry(pp_dM_chi2_eta_tnpD, "tnpD", "PL");
    legend6->AddEntry(pp_dM_chi2_eta_acooff, "acooff", "PL");
    legend6->AddEntry(pp_dM_chi2_eta_nominal_no_bk, "bk_off", "PL");
    legend6->SetTextSize(0.03);              // Set text size
    legend6->SetTextFont(42);                // Use a modern, clean font
    legend6->SetBorderSize(0);               // Set border size (0 for no border)
    legend6->SetLineColor(kBlack);           // Border color (if any)
    legend6->SetLineWidth(0);                // Border line width
    legend6->SetFillColorAlpha(kWhite, 0.1); // Background color with transparency
    legend6->Draw("SAME");

    c_pp_Eta_dW->cd();
    cosmetic(pp_dWidth_chi2_eta_nominal, pp_dWidth_chi2_eta_tnpU, pp_dWidth_chi2_eta_tnpD, pp_dWidth_chi2_eta_acooff, pp_dWidth_chi2_eta_nominal_no_bk,0);
    pp_dWidth_chi2_eta_nominal->GetYaxis()->SetRangeUser(0.2,0.34);
    pp_dWidth_chi2_eta_nominal->Draw("AP");
    pp_dWidth_chi2_eta_tnpU->Draw("P SAME");
    pp_dWidth_chi2_eta_tnpD->Draw("P SAME");
    pp_dWidth_chi2_eta_acooff->Draw("P SAME");
    pp_dWidth_chi2_eta_nominal_no_bk->Draw("P SAME");
    TLegend *legend7 = new TLegend(0.2, 0.7, 0.4, 0.9);
    legend7->AddEntry(pp_dWidth_chi2_eta_nominal, "Nominal", "PL");
    legend7->AddEntry(pp_dWidth_chi2_eta_tnpU, "tnpU", "PL");
    legend7->AddEntry(pp_dWidth_chi2_eta_tnpD, "tnpD", "PL");
    legend7->AddEntry(pp_dWidth_chi2_eta_acooff, "acooff", "PL");
    legend7->AddEntry(pp_dWidth_chi2_eta_nominal_no_bk, "bk_off", "PL");
    legend7->SetTextSize(0.03);              // Set text size
    legend7->SetTextFont(42);                // Use a modern, clean font
    legend7->SetBorderSize(0);               // Set border size (0 for no border)
    legend7->SetLineColor(kBlack);           // Border color (if any)
    legend7->SetLineWidth(0);                // Border line width
    legend7->SetFillColorAlpha(kWhite, 0.1); // Background color with transparency
    legend7->Draw("SAME");

    c_FA_dM->SaveAs("./systematic/FA_dM.png");
    c_FA_dW->SaveAs("./systematic/FA_dW.png");
    c_Eta_dM->SaveAs("./systematic/Eta_dM.png");
    c_Eta_dW->SaveAs("./systematic/Eta_dW.png");

    c_pp_FA_dM->SaveAs("./systematic/pp_FA_dM.png");
    c_pp_FA_dW->SaveAs("./systematic/pp_FA_dW.png");
    c_pp_Eta_dM->SaveAs("./systematic/pp_Eta_dM.png");
    c_pp_Eta_dW->SaveAs("./systematic/pp_Eta_dW.png");
}