void areanormalize(TH1D *h_1)
{
    if (h_1 == nullptr)
    {
        return;
    }
    double normalization_factor = h_1->Integral("width");
    h_1->Scale(1 / normalization_factor);
}

Double_t myownfunctionchi2(TH1D *h1, TH1D *h2)
{
    int nbinsh1 = h1->GetNbinsX();
    int nbinsh2 = h2->GetNbinsX();

    if (nbinsh1 != nbinsh2)
    {
        std::cerr << "Warning: h1 and h2 have different bin numbers!" << nbinsh1 << " " << nbinsh2 << std::endl;
        return -99; // Optionally, return an error value
    }

    Double_t chi2 = 0;
    for (int i = 1; i <= nbinsh1; i++)
    {
        double num = TMath::Power((h1->GetBinContent(i) - h2->GetBinContent(i)), 2);
        double den = TMath::Power((h1->GetBinError(i)), 2) + TMath::Power((h2->GetBinError(i)), 2);
        // double den = TMath::Power((h1->GetBinError(i)), 2) + TMath::Power((h2->GetBinError(i)), 2);

        if (TMath::Power((h1->GetBinContent(i)), 1) <= 0)
        {
            std::cerr << "Warning: Bin " << i << " has zero errors/missing bins/ missing bins after bksub for data. Not Skipped" << std::endl;
            // continue; // Skip this bin to avoid division by zero
        }

        chi2 += num / den;
    }

    Double_t ndf = nbinsh1 - 2;

    return chi2;
}

std::vector<std::pair<double, double>> getContour(TH2D *h1)
{
    const int nbinsx = h1->GetNbinsX();
    const int nbinsy = h1->GetNbinsY();

    // 1. Find the minimum bin
    int minBinX = 0, minBinY = 0, minBinZ = 0;
    h1->GetMinimumBin(minBinX, minBinY, minBinZ);

    double minBinContent = h1->GetBinContent(minBinX, minBinY);

    cout << "minBinContent is " << minBinContent << endl;

    // 2. Compute threshold for 1σ contour (Δχ² = 2.30)
    const double threshold = minBinContent + 2.30;

    std::vector<std::pair<double, double>> results;
    results.reserve(nbinsy);

    // 3. Scan along x for each y-bin
    for (int yindex = 1; yindex <= nbinsy; ++yindex)
    {
        bool isLeft = true;
        double left = NAN, right = NAN;

        for (int xindex = 1; xindex <= nbinsx; ++xindex)
        {
            double binContent = h1->GetBinContent(xindex, yindex);

            // find left crossing
            if (binContent < threshold && isLeft && xindex > 1)
            {
                double c1 = h1->GetBinContent(xindex - 1, yindex);
                double c2 = binContent;
                double x1 = h1->GetXaxis()->GetBinCenter(xindex - 1);
                double x2 = h1->GetXaxis()->GetBinCenter(xindex);
                left = (x1 * c1 + x2 * c2) / (c1 + c2);
                isLeft = false;
            }

            // find right crossing
            if (binContent > threshold && !isLeft && xindex > 1)
            {
                double c1 = h1->GetBinContent(xindex - 1, yindex);
                double c2 = binContent;
                double x1 = h1->GetXaxis()->GetBinCenter(xindex - 1);
                double x2 = h1->GetXaxis()->GetBinCenter(xindex);
                right = (x1 * c1 + x2 * c2) / (c1 + c2);
                break;
            }
        }

        results.emplace_back(left, right);
    }

    // optional: print where the minimum is found
    double minX = h1->GetXaxis()->GetBinCenter(minBinX);
    double minY = h1->GetYaxis()->GetBinCenter(minBinY);
    // printf("Minimum bin found at (x=%.4f, y=%.4f) with content=%.4f\n",
    //        minX, minY, minBinContent);

    return results;
}

Int_t GetPhiBin(double phi)
{
    while (phi <= -TMath::Pi())
        phi += 2 * TMath::Pi();
    while (phi > TMath::Pi())
        phi -= 2 * TMath::Pi();

    double xmin = -TMath::Pi();
    double xmax = TMath::Pi();
    int nbins = 8;
    double width = (xmax - xmin) / nbins;

    int bin = static_cast<int>((phi - xmin) / width);
    if (bin < 0)
        bin = 0;
    if (bin >= nbins)
        bin = nbins - 1;

    return bin;
}

void cosmetic(TH2D *h_1)
{
    gStyle->SetPalette(kRainBow);
    h_1->SetTitle("");
    h_1->SetTitleFont(42);
    h_1->Draw("COLZ");
    // h_1->Draw("TEXTSAME");
    h_1->GetXaxis()->SetNdivisions(21, 0, 0);
    h_1->GetYaxis()->SetNdivisions(21, 0, 0);
    h_1->GetXaxis()->SetLabelSize(0.02); // Change this value to make the labels smaller
    h_1->GetYaxis()->SetLabelSize(0.02);
    h_1->GetXaxis()->SetTitle("Mass Shifted Amount (GeV)");
    h_1->GetYaxis()->SetTitle("Width Smeared Amount (GeV)");

    for (int j = 1; j <= 42; j++)
    {
        for (int k = 1; k <= 42; k++)
        {
            double xlow = h_1->GetXaxis()->GetBinLowEdge(k);
            double xup = h_1->GetXaxis()->GetBinUpEdge(k);
            double ylow = h_1->GetYaxis()->GetBinLowEdge(j);
            double yup = h_1->GetYaxis()->GetBinUpEdge(j);

            TBox *box = new TBox(xlow, ylow, xup, yup);
            box->SetFillStyle(0);      // No fill
            box->SetLineColor(kBlack); // Black border
            box->SetLineWidth(1);      // Border width

            box->Draw("same");
        }
    }

    Int_t minBinX = -1, minBinY = -1;
    Double_t minContent = h_1->GetMaximum();
    for (Int_t binX = 1; binX <= h_1->GetNbinsX(); ++binX)
    {
        for (Int_t binY = 1; binY <= h_1->GetNbinsY(); ++binY)
        {
            Double_t content = h_1->GetBinContent(binX, binY);
            if (content < minContent)
            {
                minContent = content;
                minBinX = binX;
                minBinY = binY;
            }
        }
    }
    Double_t xMin = h_1->GetXaxis()->GetBinLowEdge(minBinX);
    Double_t xMax = h_1->GetXaxis()->GetBinUpEdge(minBinX);
    Double_t yMin = h_1->GetYaxis()->GetBinLowEdge(minBinY);
    Double_t yMax = h_1->GetYaxis()->GetBinUpEdge(minBinY);

    Double_t xCenter = h_1->GetXaxis()->GetBinCenter(minBinX);
    Double_t yCenter = h_1->GetYaxis()->GetBinCenter(minBinY);

    TBox *box1 = new TBox(xMin, yMin, xMax, yMax);
    box1->SetLineColor(kRed);
    box1->SetLineWidth(4);
    box1->SetFillStyle(0);
    box1->Draw("same");

    auto contourPoints = getContour(h_1);

    for (size_t i = 0; i < contourPoints.size(); ++i)
    {
        double y = h_1->GetYaxis()->GetBinCenter(i + 1);
        double left = contourPoints[i].first;
        double right = contourPoints[i].second;

        if (!std::isnan(left))
        {
            TMarker *mL = new TMarker(left, y, 20);
            mL->SetMarkerColor(kRed);
            mL->SetMarkerStyle(20);
            mL->SetMarkerSize(2);
            mL->Draw("same");
        }

        if (!std::isnan(right))
        {
            TMarker *mR = new TMarker(right, y, 20);
            mR->SetMarkerColor(kRed);
            mR->SetMarkerStyle(20);
            mR->SetMarkerSize(2);
            mR->Draw("same");
        }
    }
}

std::pair<TGraphErrors *, TGraphErrors *>
BuildMassAndWidthShiftGraphs(std::vector<TH2D *> &h_chi2_maps)
{
    const int nPhi = h_chi2_maps.size();

    std::vector<double> x_phi(nPhi);
    std::vector<double> y_massShift(nPhi);
    std::vector<double> y_widthShift(nPhi);
    std::vector<double> ex_phi(nPhi, 0.0);
    std::vector<double> ey_massShift(nPhi);
    std::vector<double> ey_widthShift(nPhi);

    for (int i = 0; i < nPhi; ++i)
    {
        TH2D *h2 = h_chi2_maps[i];
        if (!h2)
            continue;

        // --- Locate minimum bin ---
        int minBinX = 0, minBinY = 0, dummyZ = 0;
        h2->GetMinimumBin(minBinX, minBinY, dummyZ);
        double minX = h2->GetXaxis()->GetBinCenter(minBinX);
        double minY = h2->GetYaxis()->GetBinCenter(minBinY);

        // --- Get horizontal contour (mass shift direction) ---
        auto contourX = getContour(h2);
        double left = contourX[minBinY - 1].first;
        double right = contourX[minBinY - 1].second;
        double errLeftX = fabs(minX - left);
        double errRightX = fabs(right - minX);
        double errMaxX = std::max(errLeftX, errRightX);

        // --- Get vertical contour (width shift direction) ---
        // We just rotate the logic: scan along y for fixed x slice
        const int nbinsy = h2->GetNbinsY();
        const double minContent = h2->GetBinContent(minBinX, minBinY);
        const double threshold = minContent + 2.30;
        bool isDown = true;
        double down = NAN, up = NAN;
        for (int y = 1; y <= nbinsy; ++y)
        {
            double binContent = h2->GetBinContent(minBinX, y);
            if (binContent < threshold && isDown && y > 1)
            {
                double c1 = h2->GetBinContent(minBinX, y - 1);
                double c2 = binContent;
                double y1 = h2->GetYaxis()->GetBinCenter(y - 1);
                double y2 = h2->GetYaxis()->GetBinCenter(y);
                down = (y1 * c1 + y2 * c2) / (c1 + c2);
                isDown = false;
            }
            if (binContent > threshold && !isDown && y > 1)
            {
                double c1 = h2->GetBinContent(minBinX, y - 1);
                double c2 = binContent;
                double y1 = h2->GetYaxis()->GetBinCenter(y - 1);
                double y2 = h2->GetYaxis()->GetBinCenter(y);
                up = (y1 * c1 + y2 * c2) / (c1 + c2);
                break;
            }
        }

        double errDownY = fabs(minY - down);
        double errUpY = fabs(up - minY);
        double errMaxY = std::max(errDownY, errUpY);

        // --- Phi coordinate ---
        double phiWidth = (2 * TMath::Pi()) / nPhi;
        double phiCenter = -TMath::Pi() + (i + 0.5) * phiWidth;

        x_phi[i] = phiCenter;
        y_massShift[i] = minX;
        y_widthShift[i] = minY;
        ey_massShift[i] = errMaxX;
        ey_widthShift[i] = errMaxY;

        std::cout << Form("phi=%+.3f | mass=%.5f ± %.5f | width=%.5f ± %.5f\n",
                          phiCenter, minX, errMaxX, minY, errMaxY);
    }

    // --- Build TGraphErrors ---
    TGraphErrors *gMass = new TGraphErrors(nPhi, x_phi.data(), y_massShift.data(), ex_phi.data(), ey_massShift.data());
    gMass->SetTitle("Best-fit mass shift vs #phi;#phi (rad);Mass shift (GeV)");
    gMass->SetMarkerStyle(20);
    gMass->SetMarkerColor(kRed + 1);
    gMass->SetLineColor(kRed + 1);
    gMass->SetLineWidth(2);

    TGraphErrors *gWidth = new TGraphErrors(nPhi, x_phi.data(), y_widthShift.data(), ex_phi.data(), ey_widthShift.data());
    gWidth->SetTitle("Best-fit width smear vs #phi;#phi (rad);Width smear (GeV)");
    gWidth->SetMarkerStyle(21);
    gWidth->SetMarkerColor(kBlue + 1);
    gWidth->SetLineColor(kBlue + 1);
    gWidth->SetLineWidth(2);

    return {gMass, gWidth};
}

void get_scan_on_PbPb_and_pp_data_old()
{
    gStyle->SetOptStat(0);
    TFile *PbPb_data = new TFile("../ZBoson_18/rootfile/data_file.root", "READ");
    TFile *pp_data = new TFile("./new_pp_data_file_stability_readonly.root", "READ");
    TFile *PbPb_mc_template_PbPb_use = new TFile("../ZBoson_18/rootfile/version_2/template_PbPb_pp_data_phi_PbPb.root", "READ");
    TFile *PbPb_mc_template_pp_use = new TFile("../ZBoson_18/rootfile/version_2/template_PbPb_pp_data_phi_pp.root", "READ");

    TH1D *PbPb_data_in_Phi[8];
    TH1D *pp_data_in_Phi[8];
    TH1D *pp_data_in_Phi_without_pT_reweight[8];
    TH1D *PbPb_mc_nominal_template_PbPb[42][42];
    TH1D *PbPb_mc_nominal_template_pp[42][42];

    TH2D *h_chisquare_pp_PbPb[8];
    TH2D *h_chisquare_PbPb_PbPb[8];
    TH2D *h_chisquare_pp_PbPb_no_pT[8];

    double highbin_mass_shift = 0.0;
    double highbin_smear = 0.4;
    double lowbin_mass_shift = -0.5;
    double lowbin_smear = -0.3;

    double h_low_mass_shift = lowbin_mass_shift - ((highbin_mass_shift - lowbin_mass_shift) / (42 - 1)) / 2;
    double h_high_mass_shift = highbin_mass_shift + ((highbin_mass_shift - lowbin_mass_shift) / (42 - 1)) / 2;

    double h_low_smear = lowbin_smear - ((highbin_smear - lowbin_smear) / (42 - 1)) / 2;
    double h_high_smear = highbin_smear + ((highbin_smear - lowbin_smear) / (42 - 1)) / 2;

    double highbin_mass_shift_pp = 0.1;
    double highbin_smear_pp = 0.2;
    double lowbin_mass_shift_pp = -0.3;
    double lowbin_smear_pp = 0.0;

    double h_low_mass_shift_pp = lowbin_mass_shift_pp - ((highbin_mass_shift_pp - lowbin_mass_shift_pp) / (42 - 1)) / 2;
    double h_high_mass_shift_pp = highbin_mass_shift_pp + ((highbin_mass_shift_pp - lowbin_mass_shift_pp) / (42 - 1)) / 2;

    double h_low_smear_pp = lowbin_smear_pp - ((highbin_smear_pp - lowbin_smear_pp) / (42 - 1)) / 2;
    double h_high_smear_pp = highbin_smear_pp + ((highbin_smear_pp - lowbin_smear_pp) / (42 - 1)) / 2;

    for (int i = 0; i < 8; i++)
    {
        PbPb_data_in_Phi[i] = (TH1D *)PbPb_data->Get(Form("FA_nominal_phi_plus_%i", i));
        pp_data_in_Phi[i] = (TH1D *)pp_data->Get(Form("pp_FA_nominal_phi_plus_%i", i));
        pp_data_in_Phi_without_pT_reweight[i] = (TH1D *)pp_data->Get(Form("pp_FA_nominal_phi_plus_without_pT_reweight_%i", i));
        h_chisquare_pp_PbPb[i] = new TH2D(Form("h_chisquare_pp_PbPb_%i", i), "", 42, h_low_mass_shift_pp, h_high_mass_shift_pp, 42, h_low_smear_pp, h_high_smear_pp);
        h_chisquare_PbPb_PbPb[i] = new TH2D(Form("h_chisquare_PbPb_PbPb_%i", i), "", 42, h_low_mass_shift, h_high_mass_shift, 42, h_low_smear, h_high_smear);
        h_chisquare_pp_PbPb_no_pT[i] = new TH2D(Form("h_chisquare_pp_PbPb_no_pT_%i", i), "", 42, h_low_mass_shift_pp, h_high_mass_shift_pp, 42, h_low_smear_pp, h_high_smear_pp);

        PbPb_data_in_Phi[i]->Rebin(4);
        pp_data_in_Phi[i]->Rebin(4);
        pp_data_in_Phi_without_pT_reweight[i]->Rebin(4);

        areanormalize(PbPb_data_in_Phi[i]);
        areanormalize(pp_data_in_Phi[i]);
        areanormalize(pp_data_in_Phi_without_pT_reweight[i]);
    }

    for (int i = 0; i < 42; i++)
    {
        for (int j = 0; j < 42; j++)
        {
            PbPb_mc_nominal_template_PbPb[i][j] = (TH1D *)PbPb_mc_template_PbPb_use->Get(Form("template_FA_nominal_%i_%i_%i", i, j, 10));
            PbPb_mc_nominal_template_PbPb[i][j]->Rebin(4);
            areanormalize(PbPb_mc_nominal_template_PbPb[i][j]);

            PbPb_mc_nominal_template_pp[i][j] = (TH1D *)PbPb_mc_template_pp_use->Get(Form("template_FA_nominal_%i_%i_%i", i, j, 10));
            PbPb_mc_nominal_template_pp[i][j]->Rebin(4);
            areanormalize(PbPb_mc_nominal_template_pp[i][j]);
        }
    }

    for (int phi = 0; phi < 8; phi++)
    {
        for (int i = 0; i < 42; i++)
        {
            for (int j = 0; j < 42; j++)
            {
                double chisquarevalue_PbPb_PbPb = 0;
                double chisquarevalue_PbPb_pp = 0;
                double chisquarevalue_PbPb_pp_without_pT_reweight = 0;

                chisquarevalue_PbPb_PbPb = myownfunctionchi2(PbPb_data_in_Phi[phi], PbPb_mc_nominal_template_PbPb[i][j]);
                chisquarevalue_PbPb_pp = myownfunctionchi2(pp_data_in_Phi[phi], PbPb_mc_nominal_template_pp[i][j]);
                chisquarevalue_PbPb_pp_without_pT_reweight = myownfunctionchi2(pp_data_in_Phi_without_pT_reweight[phi], PbPb_mc_nominal_template_pp[i][j]);

                std::ostringstream stream;
                stream << std::fixed << std::setprecision(4) << chisquarevalue_PbPb_PbPb;
                double formatted_PbPb_PbPb = std::stod(stream.str());

                std::ostringstream stream_1;
                stream_1 << std::fixed << std::setprecision(4) << chisquarevalue_PbPb_pp;
                double formatted_PbPb_pp = std::stod(stream_1.str());

                std::ostringstream stream_2;
                stream_2 << std::fixed << std::setprecision(4) << chisquarevalue_PbPb_pp_without_pT_reweight;
                double formatted_PbPb_pp_wihout_pT_reweight = std::stod(stream_2.str());

                h_chisquare_pp_PbPb[phi]->SetBinContent(i + 1, j + 1, formatted_PbPb_pp);
                h_chisquare_PbPb_PbPb[phi]->SetBinContent(i + 1, j + 1, formatted_PbPb_PbPb);
                h_chisquare_pp_PbPb_no_pT[phi]->SetBinContent(i + 1, j + 1, formatted_PbPb_pp_wihout_pT_reweight);
            }
        }
    }

    for (int i = 0; i < 8; i++)
    {
        TCanvas *c_PbPb_PbPb = new TCanvas("c_PbPb_PbPb", "", 3200, 2400);
        TCanvas *c_PbPb_pp = new TCanvas("c_PbPb_pp", "", 3200, 2400);
        TCanvas *c_PbPb_pp_no_pT = new TCanvas("c_PbPb_pp_no_pT", "", 3200, 2400);

        c_PbPb_PbPb->SetLeftMargin(0.10);   // smaller left margin
        c_PbPb_PbPb->SetRightMargin(0.12);  // smaller right margin
        c_PbPb_PbPb->SetTopMargin(0.08);    // smaller top margin
        c_PbPb_PbPb->SetBottomMargin(0.10); // smaller bottom margin

        c_PbPb_pp->SetLeftMargin(0.10);   // smaller left margin
        c_PbPb_pp->SetRightMargin(0.12);  // smaller right margin
        c_PbPb_pp->SetTopMargin(0.08);    // smaller top margin
        c_PbPb_pp->SetBottomMargin(0.10); // smaller bottom margin

        c_PbPb_pp_no_pT->SetLeftMargin(0.10);   // smaller left margin
        c_PbPb_pp_no_pT->SetRightMargin(0.12);  // smaller right margin
        c_PbPb_pp_no_pT->SetTopMargin(0.08);    // smaller top margin
        c_PbPb_pp_no_pT->SetBottomMargin(0.10); // smaller bottom margin

        c_PbPb_pp->cd();
        h_chisquare_pp_PbPb[i]->Draw("COLZ");
        cosmetic(h_chisquare_pp_PbPb[i]);
        c_PbPb_pp->Update();
        c_PbPb_pp->SaveAs(Form("./ScaninPhi/data/pp_%i.png", i));

        c_PbPb_PbPb->cd();
        h_chisquare_PbPb_PbPb[i]->Draw("COLZ");
        cosmetic(h_chisquare_PbPb_PbPb[i]);
        c_PbPb_PbPb->Update();
        c_PbPb_PbPb->SaveAs(Form("./ScaninPhi/data/PbPb_%i.png", i));

        c_PbPb_pp_no_pT->cd();
        h_chisquare_pp_PbPb_no_pT[i]->Draw("COLZ");
        cosmetic(h_chisquare_pp_PbPb_no_pT[i]);
        c_PbPb_pp_no_pT->Update();
        c_PbPb_pp_no_pT->SaveAs(Form("./ScaninPhi/data/pp_no_pT_%i.png", i));

        delete c_PbPb_PbPb;
        delete c_PbPb_pp;
        delete c_PbPb_pp_no_pT;
    }
    // === Build TGraphErrors for pp_PbPb, PbPb_PbPb, and pp_PbPb_no_pT ===
    std::vector<TH2D *> vec_pp_PbPb, vec_PbPb_PbPb, vec_pp_no_pT;
    for (int i = 0; i < 8; ++i)
    {
        vec_pp_PbPb.push_back(h_chisquare_pp_PbPb[i]);
        vec_PbPb_PbPb.push_back(h_chisquare_PbPb_PbPb[i]);
        vec_pp_no_pT.push_back(h_chisquare_pp_PbPb_no_pT[i]);
    }

    // Build mass- and width-shift graphs for each case
    auto [gMass_ppPbPb, gWidth_ppPbPb] = BuildMassAndWidthShiftGraphs(vec_pp_PbPb);
    auto [gMass_PbPbPbPb, gWidth_PbPbPbPb] = BuildMassAndWidthShiftGraphs(vec_PbPb_PbPb);
    auto [gMass_pp_no_pT, gWidth_pp_no_pT] = BuildMassAndWidthShiftGraphs(vec_pp_no_pT);

    //---------------------------------------------------------------
    // --- Compute Y-axis ranges from PbPb→PbPb (blue) ---
    double yMinMass = gMass_PbPbPbPb->GetY()[0];
    double yMaxMass = gMass_PbPbPbPb->GetY()[0];
    double yMinWidth = gWidth_PbPbPbPb->GetY()[0];
    double yMaxWidth = gWidth_PbPbPbPb->GetY()[0];

    for (int i = 0; i < gMass_PbPbPbPb->GetN(); ++i)
    {
        double ym = gMass_PbPbPbPb->GetY()[i];
        double ye = gMass_PbPbPbPb->GetEY()[i];
        yMinMass = std::min(yMinMass, ym - ye);
        yMaxMass = std::max(yMaxMass, ym + ye);
    }
    for (int i = 0; i < gWidth_PbPbPbPb->GetN(); ++i)
    {
        double yw = gWidth_PbPbPbPb->GetY()[i];
        double ye = gWidth_PbPbPbPb->GetEY()[i];
        yMinWidth = std::min(yMinWidth, yw - ye);
        yMaxWidth = std::max(yMaxWidth, yw + ye);
    }
    // Add 10% padding
    yMinMass -= 0.2 * fabs(yMaxMass - yMinMass);
    yMaxMass += 0.2 * fabs(yMaxMass - yMinMass);
    yMinWidth -= 0.4 * fabs(yMaxWidth - yMinWidth);
    yMaxWidth += 0.4 * fabs(yMaxWidth - yMinWidth);

    // --- CMS style setup ---
    gStyle->SetOptStat(0);
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTitleSize(0.045, "XYZ");
    gStyle->SetLabelSize(0.04, "XYZ");
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetFrameLineWidth(2);

    // --- Separate canvas: Mass shift ---

    gStyle->SetPadTopMargin(0.05);
    gStyle->SetPadBottomMargin(0.13);
    gStyle->SetPadLeftMargin(0.16);
    gStyle->SetPadRightMargin(0.04);

    TCanvas *cMass = new TCanvas("cMass", "Mass shift vs phi", 800, 800);
    gMass_ppPbPb->SetLineColor(kRed + 1);
    gMass_PbPbPbPb->SetLineColor(kBlue + 1);
    gMass_pp_no_pT->SetLineColor(kGreen + 2);
    gMass_ppPbPb->SetMarkerColor(kRed + 1);
    gMass_PbPbPbPb->SetMarkerColor(kBlue + 1);
    gMass_pp_no_pT->SetMarkerColor(kGreen + 2);
    gMass_PbPbPbPb->GetYaxis()->SetRangeUser(-0.5, 0.4);
    gMass_PbPbPbPb->SetTitle(";#phi (rad);dMass (GeV)");
    gMass_PbPbPbPb->Draw("APL");
    gMass_ppPbPb->Draw("PLsame");
    gMass_pp_no_pT->Draw("PLsame");
    auto legM = new TLegend(0.2, 0.75, 0.45, 0.9);
    legM->AddEntry(gMass_ppPbPb, "pp#rightarrowPbPb", "lp");
    legM->AddEntry(gMass_PbPbPbPb, "PbPb#rightarrowPbPb", "lp");
    legM->AddEntry(gMass_pp_no_pT, "pp no p_{T} reweight", "lp");
    legM->SetBorderSize(0);
    legM->SetTextFont(42);
    legM->Draw();
    TLatex latexM;
    latexM.SetNDC();
    latexM.SetTextFont(62);
    latexM.SetTextSize(0.045);
    // latexM.DrawLatex(0.16, 0.93, "#bf{CMS}  #it{Preliminary}");
    cMass->SaveAs("./ScaninPhi/data/massShift_vs_phi_CMS.png");

    // --- Separate canvas: Width smear ---
    TCanvas *cWidth = new TCanvas("cWidth", "Width smear vs phi", 800, 800);
    gWidth_ppPbPb->SetLineColor(kRed + 1);
    gWidth_PbPbPbPb->SetLineColor(kBlue + 1);
    gWidth_pp_no_pT->SetLineColor(kGreen + 2);
    gWidth_ppPbPb->SetMarkerColor(kRed + 1);
    gWidth_PbPbPbPb->SetMarkerColor(kBlue + 1);
    gWidth_pp_no_pT->SetMarkerColor(kGreen + 2);
    gWidth_PbPbPbPb->GetYaxis()->SetRangeUser(-0.5, 0.8);
    gWidth_PbPbPbPb->SetTitle(";#phi (rad); dWidth (GeV)");
    gWidth_PbPbPbPb->Draw("APL");
    gWidth_ppPbPb->Draw("PLsame");
    gWidth_pp_no_pT->Draw("PLsame");
    auto legW = new TLegend(0.2, 0.75, 0.45, 0.9);
    legW->AddEntry(gWidth_ppPbPb, "pp#rightarrowPbPb", "lp");
    legW->AddEntry(gWidth_PbPbPbPb, "PbPb#rightarrowPbPb", "lp");
    legW->AddEntry(gWidth_pp_no_pT, "pp no p_{T} reweight", "lp");
    legW->SetBorderSize(0);
    legW->SetTextFont(42);
    legW->Draw();
    TLatex latexW;
    latexW.SetNDC();
    latexW.SetTextFont(62);
    latexW.SetTextSize(0.045);
    // latexW.DrawLatex(0.16, 0.93, "#bf{CMS}  #it{Preliminary}");
    cWidth->SaveAs("./ScaninPhi/data/widthShift_vs_phi_CMS.png");
}