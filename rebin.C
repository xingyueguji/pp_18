std::vector<double> createCustomBinning(
    TH1D *hist,
    double range1_min, double range1_max, int rebin1,
    double range2_min, double range2_max, int rebin2,
    double range3_min, double range3_max, int rebin3,
    double range4_min, double range4_max, int rebin4)
{
    // Get the original histogram range and binning
    double xMin = hist->GetXaxis()->GetXmin();
    double xMax = hist->GetXaxis()->GetXmax();
    int nBins = hist->GetNbinsX();
    double binWidth = (xMax - xMin) / nBins; // Original bin width

    std::vector<double> newBins;

    // Function to check if the rebinning factor is proper
    auto is_proper_rebin = [](double range_min, double range_max, double binWidth, int rebin) -> bool
    {
        double nOriginalBins = (range_max - range_min) / binWidth;
        // cout << std::fmod(nOriginalBins, rebin) << endl;
        return std::fmod(nOriginalBins, rebin) == 0; // Proper if nOriginalBins is divisible by rebin
    };

    // Check if the rebinning factors are proper for each range
    if (!is_proper_rebin(range1_min, range1_max, binWidth, rebin1))
    {
        std::cerr << "Warning: Rebinning factor for range [" << range1_min << ", " << range1_max << "] is improper." << std::endl;
    }
    if (!is_proper_rebin(range2_min, range2_max, binWidth, rebin2))
    {
        std::cerr << "Warning: Rebinning factor for range [" << range2_min << ", " << range2_max << "] is improper." << std::endl;
    }
    if (!is_proper_rebin(range3_min, range3_max, binWidth, rebin3))
    {
        std::cerr << "Warning: Rebinning factor for range [" << range3_min << ", " << range3_max << "] is improper." << std::endl;
    }
    if (!is_proper_rebin(range4_min, range4_max, binWidth, rebin4))
    {
        std::cerr << "Warning: Rebinning factor for range [" << range4_min << ", " << range4_max << "] is improper." << std::endl;
    }

    // Step 1: Loop over the histogram bins and create bin edges
    double currentEdge = xMin;

    while (currentEdge < xMax)
    {
        // If we are in the first rebinning range
        if (currentEdge >= range1_min && currentEdge < range1_max)
        {
            newBins.push_back(currentEdge);
            currentEdge += rebin1 * binWidth;
        }
        // If we are in the second rebinning range
        else if (currentEdge >= range2_min && currentEdge < range2_max)
        {
            newBins.push_back(currentEdge);
            currentEdge += rebin2 * binWidth;
        }
        // If we are in the third rebinning range
        else if (currentEdge >= range3_min && currentEdge < range3_max)
        {
            newBins.push_back(currentEdge);
            currentEdge += rebin3 * binWidth;
        }
        // If we are in the fourth rebinning range
        else if (currentEdge >= range4_min && currentEdge < range4_max)
        {
            newBins.push_back(currentEdge);
            currentEdge += rebin4 * binWidth;
        }
        // If we are outside the custom rebinning ranges
        else
        {
            newBins.push_back(currentEdge);
            currentEdge += binWidth; // Regular bin width
        }
    }

    newBins.push_back(xMax); // Make sure to include the upper edge of the last bin

    return newBins;
}

void ConvertToDNdx(TH1D* h2) {
    // Check if the histogram exists
    if (!h2) {
        std::cerr << "Error: Histogram does not exist!" << std::endl;
        return;
    }

    // Loop over all bins in the histograms
    int nBins = h2->GetNbinsX();
    for (int bin = 1; bin <= nBins; ++bin) { // Loop over bins (1 to nBins)
        double binContent = h2->GetBinContent(bin);
        double binWidth = h2->GetBinWidth(bin);

        if (binWidth > 0) {
            // Normalize bin content by bin width
            double normalizedContent = binContent / binWidth;
            h2->SetBinContent(bin, normalizedContent);

            // Scale the error as well (if applicable)
            double binError = h2->GetBinError(bin);
            double normalizedError = binError / binWidth;
            h2->SetBinError(bin, normalizedError);
        }
    }

    // Optionally, update the histogram y-axis title
    h2->GetYaxis()->SetTitle("dN/dx");
    std::cout << "Histogram converted to dN/dx." << std::endl;
}

void rebin()
{
    TFile *f1 = new TFile("./data_file.root");

    Int_t cenlowlimit[11] = {0, 10, 20, 30, 30, 0, 15, 50, 0, 14, 0};
    Int_t cenhighlimit[11] = {10, 20, 30, 100, 50, 15, 100, 100, 14, 100, 100};

    TH1D *mass_array_data_raw[11];
    TH1D *mass_array_data_eta[11];

    TH1D *modified_mc_signal_raw[11];
    TH1D *modified_mc_signal_eta[11];

    TH1D *mc_bk_raw[11];
    TH1D *mc_bk_eta[11];

    // We have 120 bins before, now want to combine two tails first to check the result?
    // Used eta cut as ref. Cent based.

    // 0-10% 60 - 84 rebin(3), 96-120 rebin(8)

    for (int i = 0; i < 11; i++)
    {
        mass_array_data_raw[i] = (TH1D *)f1->Get(Form("mass_array_data_%i", i));
        mass_array_data_eta[i] = (TH1D *)f1->Get(Form("mass_array_data_witheta_%i", i));
    }

    std::vector<double> newBinEdges_0_10 = createCustomBinning(mass_array_data_eta[0], 60, 76, 4, 76, 83, 2, 96, 102, 3, 102, 120, 18);
    std::vector<double> newBinEdges_10_20 = createCustomBinning(mass_array_data_eta[0], 60, 76, 8, 76, 84, 2, 96, 100, 2, 100, 120, 20);
    std::vector<double> newBinEdges_20_30 = createCustomBinning(mass_array_data_eta[0], 60, 76, 16, 76, 85, 3, 94, 98, 2, 98, 120, 22);
    std::vector<double> newBinEdges_30_50 = createCustomBinning(mass_array_data_eta[0], 60, 78, 12, 78, 86, 4, 95, 100, 5, 100, 120, 40);
    std::vector<double> newBinEdges_50_100 = createCustomBinning(mass_array_data_eta[0], 60, 85, 25, 85, 89, 4, 92.5, 96.5, 4, 96.5, 120, 47);
    std::vector<double> newBinEdges_0_100 = createCustomBinning(mass_array_data_eta[0], 60, 80, 10, 80, 86, 4, 95, 100, 5, 100, 120, 40);
    //60, 80, 20, 80, 86, 2, 95.5, 100.5, 10, 100.5, 120, 39);

    int nNewBins[6];

    nNewBins[0] = newBinEdges_0_100.size() - 1;
    nNewBins[1] = newBinEdges_0_100.size() - 1;
    nNewBins[2] = newBinEdges_0_100.size() - 1;
    nNewBins[3] = newBinEdges_0_100.size() - 1;
    nNewBins[4] = newBinEdges_0_100.size() - 1;
    nNewBins[5] = newBinEdges_0_100.size() - 1;

    Double_t *binEdgesArray[6];

    binEdgesArray[0] = &newBinEdges_0_100[0];
    binEdgesArray[1] = &newBinEdges_0_100[0];
    binEdgesArray[2] = &newBinEdges_0_100[0];
    binEdgesArray[3] = &newBinEdges_0_100[0];
    binEdgesArray[4] = &newBinEdges_0_100[0];
    binEdgesArray[5] = &newBinEdges_0_100[0];

    TH1D *h2[11];

    for (int i = 0; i < 11; i++)
    {
        TCanvas *c1 = new TCanvas("", "", 1600, 1200);
        TCanvas *c2 = new TCanvas("", "", 1600, 1200);
        TCanvas *c3 = new TCanvas("", "", 1600, 1200);

        c1->cd();

        mass_array_data_raw[i]->Draw("HIST");
        mass_array_data_raw[i]->Draw("textsame");

        c2->cd();

        mass_array_data_eta[i]->Draw("HIST");
        mass_array_data_eta[i]->Draw("textsame");

        c3->cd();

        if (i < 4)
        {
            h2[i] = (TH1D *)mass_array_data_eta[i]->Rebin(nNewBins[i], Form("mass_array_data_witheta_rebin_%i", i), binEdgesArray[i]);
            ConvertToDNdx(h2[i]);
            h2[i]->Draw("HIST");
            h2[i]->Draw("textsame");

            cout << "# of bins - 2 - 1 for " << cenlowlimit[i] << " " << cenhighlimit[i] << " is " << h2[i]->GetNbinsX() - 3 << endl;
        }
        else if (i == 10)
        {
            h2[i] = (TH1D *)mass_array_data_eta[i]->Rebin(nNewBins[5], Form("mass_array_data_witheta_rebin_%i", i), binEdgesArray[5]);
            ConvertToDNdx(h2[i]);
            h2[i]->Draw("HIST");
            h2[i]->Draw("textsame");

            cout << "# of bins - 2 - 1 for " << cenlowlimit[i] << " " << cenhighlimit[i] << " is " << h2[i]->GetNbinsX() - 3 << endl;
        }
        else
        {
            continue;
        }

        c1->SaveAs(Form("./rebintest/raw/cent_%i_%i.png", cenlowlimit[i], cenhighlimit[i]));
        c2->SaveAs(Form("./rebintest/eta/cent_%i_%i.png", cenlowlimit[i], cenhighlimit[i]));
        c3->SaveAs(Form("./rebintest/eta/cent_rebin_%i_%i.png", cenlowlimit[i], cenhighlimit[i]));

        delete c1;
        delete c2;
        delete c3;
    }
}