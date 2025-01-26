#include "tdrStyle.C"
#include "CMS_lumi.C"

TGraphErrors *SubtractTGraphErrors(const TGraphErrors *g1, const TGraphErrors *g2)
{
    if (g1->GetN() != g2->GetN())
    {
        throw std::runtime_error("Graphs must have the same number of points");
    }

    int nPoints = g1->GetN();
    TGraphErrors *result = new TGraphErrors(nPoints);

    for (int i = 0; i < nPoints; ++i)
    {
        double x1, y1, x2, y2;
        double ex1, ey1, ex2, ey2;

        g1->GetPoint(i, x1, y1);
        g2->GetPoint(i, x2, y2);

        // Ensure the x-coordinates match
        if (TMath::Abs(x1 - x2) > 1e-6)
        {
            throw std::runtime_error("Graphs have mismatched x-coordinates");
        }

        ex1 = g1->GetErrorX(i);
        ey1 = g1->GetErrorY(i);
        ex2 = g2->GetErrorX(i);
        ey2 = g2->GetErrorY(i);

        cout << "Error 1 is " << ey1 << " Error 2 is " << ey2 << endl;

        double yDiff = y1 - y2;
        double eYDiff = TMath::Sqrt(ey1 * ey1 + ey2 * ey2);

        result->SetPoint(i, x1, yDiff);
        result->SetPointError(i, ex1, eYDiff); // Assume x-errors are the same
    }

    return result;
}

void drawandsave(bool israw, bool isdm, TGraphErrors *g1, TGraphErrors *g2, TGraphErrors *g3, TGraphErrors *g4)
{
    // Create an array of TGraphErrors pointers for easy looping
    TGraphErrors *graphs[] = {g1, g2, g3, g4};
    TString Name[] = {"raw-eff", "raw-bk", "eff-bk+eff", "bk-bk+eff"};
    TString Title[] = {"no correction - eff correction", "no correction - bk subtraction", "eff correction - bk with eff correction", "bk subtraction - bk with eff correction"};
    const int nGraphs = 4; // Number of graphs

    // Create a canvas
    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 800);

    // Loop through the graphs
    for (int i = 0; i < nGraphs; i++)
    {
        // Clear the canvas for the next graph
        c1->Clear();

        // Access the current graph
        TGraphErrors *currentGraph = graphs[i];
        TString currentName = Name[i];
        TString currentTitle = Title[i];

        // Customize the graph's appearance (optional)
        currentGraph->SetMarkerStyle(20);   // Different marker style for each graph
        currentGraph->SetLineColor(kBlack); // Different color for each graph

        currentGraph->SetTitle("");

        currentGraph->GetXaxis()->SetLimits(0.5, 5.5);
        currentGraph->GetXaxis()->SetRangeUser(0.5, 5.5);
        currentGraph->GetXaxis()->SetNdivisions(5, 0, 0, kFALSE);
        currentGraph->GetXaxis()->SetLabelSize(0);
        if (isdm)
            currentGraph->GetYaxis()->SetRangeUser(-0.2, 0.2);
        if (!isdm)
            currentGraph->GetYaxis()->SetRangeUser(-0.5, 0.5);
        if (isdm)
            currentGraph->GetYaxis()->SetTitle("dM");
        if (!isdm)
            currentGraph->GetYaxis()->SetTitle("dWidth");

        currentGraph->GetXaxis()->ChangeLabel(1, 0, 0.04, 11, -1, -1, "   0-100%");
        currentGraph->GetXaxis()->ChangeLabel(2, 0, 0.04, 11, -1, -1, "    0-10%");
        currentGraph->GetXaxis()->ChangeLabel(3, 0, 0.04, 11, -1, -1, "  10-20%");
        currentGraph->GetXaxis()->ChangeLabel(4, 0, 0.04, 11, -1, -1, "  20-30%");
        currentGraph->GetXaxis()->ChangeLabel(5, 0, 0.04, 11, -1, -1, "  30-100%");
        currentGraph->GetXaxis()->SetLabelOffset(0.03);

        // Draw the graph
        currentGraph->Draw("AP");

        TPaveText *pt1 = new TPaveText(0.15, 0.85, 0.7, 0.9, "NDC");
		pt1->AddText(currentTitle);
		pt1->SetTextSize(0.03);
        pt1->SetTextFont(42);
		pt1->SetTextAlign(22); // Center alignment
		pt1->SetFillStyle(0);  // Make the background transparent
		pt1->SetBorderSize(0); // Remove the border
		pt1->SetLineColor(0);  // Remove the border line (optional)
		pt1->SetTextColor(1);  // Set text color (default: black)
		pt1->Draw();

        // Save the canvas as an image file
        TString filename;
        if (israw)
        {
            if (isdm)
            {
                filename = Form("./backgroundeff/raw/%s_dm.png", currentName.Data());
            }
            else
            {
                filename = Form("./backgroundeff/raw/%s_dw.png", currentName.Data());
            }
        }
        if (!israw)
        {
            if (isdm)
            {
                filename = Form("./backgroundeff/eta/%s_dm.png", currentName.Data());
            }
            else
            {
                filename = Form("./backgroundeff/eta/%s_dw.png", currentName.Data());
            }
        }
        c1->SaveAs(filename);
    }

    // Delete the canvas to free memory
    delete c1;
}

void backgroundcheck()
{
    setTDRStyle();
    TFile *f_1 = new TFile("./All_plots.root", "READ");

    // This is very very tedious.... I am going to do two version
    // 1. raw
    // raw vs raw_eff, raw vs raw_bk, raw_eff vs raw_eff_bk, raw_bk vs raw_bk_eff
    TGraphErrors *raw_dm = (TGraphErrors *)f_1->Get("HI_dM_chi2_raw");
    TGraphErrors *raw_dw = (TGraphErrors *)f_1->Get("HI_dWidth_chi2_raw");
    TGraphErrors *raw_eff_dm = (TGraphErrors *)f_1->Get("HI_dM_chi2_raw_eff");
    TGraphErrors *raw_eff_dw = (TGraphErrors *)f_1->Get("HI_dWidth_chi2_raw_eff");
    TGraphErrors *raw_bk_dm = (TGraphErrors *)f_1->Get("HI_dM_chi2_raw_bksub");
    TGraphErrors *raw_bk_dw = (TGraphErrors *)f_1->Get("HI_dWidth_chi2_raw_bksub");
    TGraphErrors *raw_eff_bk_dm = (TGraphErrors *)f_1->Get("HI_dM_chi2_raw_bksub_eff");
    TGraphErrors *raw_eff_bk_dw = (TGraphErrors *)f_1->Get("HI_dWidth_chi2_raw_bksub_eff");

    // 2. eta
    // eta vs eta_eff, eta vs eta_bk, eta_eff vs eta_eff_bk, eta_bk vs eta_bk_eff
    TGraphErrors *eta_dm = (TGraphErrors *)f_1->Get("HI_dM_chi2_eta");
    TGraphErrors *eta_dw = (TGraphErrors *)f_1->Get("HI_dWidth_chi2_eta");
    TGraphErrors *eta_eff_dm = (TGraphErrors *)f_1->Get("HI_dM_chi2_eta_eff");
    TGraphErrors *eta_eff_dw = (TGraphErrors *)f_1->Get("HI_dWidth_chi2_eta_eff");
    TGraphErrors *eta_bk_dm = (TGraphErrors *)f_1->Get("HI_dM_chi2_eta_bksub");
    TGraphErrors *eta_bk_dw = (TGraphErrors *)f_1->Get("HI_dWidth_chi2_eta_bksub");
    TGraphErrors *eta_eff_bk_dm = (TGraphErrors *)f_1->Get("HI_dM_chi2_eta_bksub_eff");
    TGraphErrors *eta_eff_bk_dw = (TGraphErrors *)f_1->Get("HI_dWidth_chi2_eta_bksub_eff");

    TGraphErrors *raw_raw_eff_dm = SubtractTGraphErrors(raw_dm, raw_eff_dm);
    TGraphErrors *raw_raw_bk_dm = SubtractTGraphErrors(raw_dm, raw_bk_dm);
    TGraphErrors *raw_eff_raw_eff_bk_dm = SubtractTGraphErrors(raw_eff_dm, raw_eff_bk_dm);
    TGraphErrors *raw_bk_raw_eff_bk_dm = SubtractTGraphErrors(raw_bk_dm, raw_eff_bk_dm);

    TGraphErrors *raw_raw_eff_dw = SubtractTGraphErrors(raw_dw, raw_eff_dw);
    TGraphErrors *raw_raw_bk_dw = SubtractTGraphErrors(raw_dw, raw_bk_dw);
    TGraphErrors *raw_eff_raw_eff_bk_dw = SubtractTGraphErrors(raw_eff_dw, raw_eff_bk_dw);
    TGraphErrors *raw_bk_raw_eff_bk_dw = SubtractTGraphErrors(raw_bk_dw, raw_eff_bk_dw);

    TGraphErrors *eta_eta_eff_dm = SubtractTGraphErrors(eta_dm, eta_eff_dm);
    TGraphErrors *eta_eta_bk_dm = SubtractTGraphErrors(eta_dm, eta_bk_dm);
    TGraphErrors *eta_eff_eta_eff_bk_dm = SubtractTGraphErrors(eta_eff_dm, eta_eff_bk_dm);
    TGraphErrors *eta_bk_eta_eff_bk_dm = SubtractTGraphErrors(eta_bk_dm, eta_eff_bk_dm);

    TGraphErrors *eta_eta_eff_dw = SubtractTGraphErrors(eta_dw, eta_eff_dw);
    TGraphErrors *eta_eta_bk_dw = SubtractTGraphErrors(eta_dw, eta_bk_dw);
    TGraphErrors *eta_eff_eta_eff_bk_dw = SubtractTGraphErrors(eta_eff_dw, eta_eff_bk_dw);
    TGraphErrors *eta_bk_eta_eff_bk_dw = SubtractTGraphErrors(eta_bk_dw, eta_eff_bk_dw);

    drawandsave(1, 1, raw_raw_eff_dm, raw_raw_bk_dm, raw_eff_raw_eff_bk_dm, raw_bk_raw_eff_bk_dm);
    drawandsave(1, 0, raw_raw_eff_dw, raw_raw_bk_dw, raw_eff_raw_eff_bk_dw, raw_bk_raw_eff_bk_dw);
    drawandsave(0, 1, eta_eta_eff_dm, eta_eta_bk_dm, eta_eff_eta_eff_bk_dm, eta_bk_eta_eff_bk_dm);
    drawandsave(0, 0, eta_eta_eff_dw, eta_eta_bk_dw, eta_eff_eta_eff_bk_dw, eta_bk_eta_eff_bk_dw);
}