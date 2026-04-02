void areanormalize(TH1D *h_1)
{
    if (h_1 == nullptr)
    {
        return;
    }
    double normalization_factor = h_1->Integral("width");
    h_1->Scale(1 / normalization_factor);
}

void shiftpp(TGraphErrors *g1, double add)
{
    double y = g1->GetY()[0];
    double ex = g1->GetEX()[0];
    double ey = g1->GetEY()[0];

    g1->SetPoint(0, 5.0 + add, y);
    g1->SetPointError(0, ex, ey);

    cout << "X is " << g1->GetX()[0] << endl;
}

void outputresult(TH2D *h)
{
    // Get minimum bin content
    double minVal = h->GetMinimum();

    // Get global bin index of minimum
    int gbin = h->GetMinimumBin();

    // Convert global bin index → (ix, iy, iz)
    int ix, iy, iz;
    h->GetBinXYZ(gbin, ix, iy, iz);

    // Get axis coordinates of that bin
    double x_center = h->GetXaxis()->GetBinCenter(ix);
    double y_center = h->GetYaxis()->GetBinCenter(iy);
}

std::vector<std::pair<double, double>> getContour(TH2D *h1)
{
    const int nbinsx = h1->GetNbinsX();
    const int nbinsy = h1->GetNbinsY();

    // 1. Find the minimum bin
    int minBinX = 0, minBinY = 0, minBinZ = 0;
    h1->GetMinimumBin(minBinX, minBinY, minBinZ);

    double minBinContent = h1->GetBinContent(minBinX, minBinY);

    // cout << "minBinContent is " << minBinContent << endl;

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

std::pair<TGraphErrors *, TGraphErrors *>
BuildMassAndWidthShiftGraphs(std::vector<TH2D *> &h_chi2_maps)
{
    const int nAll = h_chi2_maps.size();
    if (nAll == 0)
        return {nullptr, nullptr};

    const int nPhi = nAll - 1; // last one is inclusive

    std::vector<double> x_phi(nAll);
    std::vector<double> y_massShift(nAll);
    std::vector<double> y_widthShift(nAll);
    std::vector<double> ex_phi(nAll, 0.0);
    std::vector<double> ey_massShift(nAll);
    std::vector<double> ey_widthShift(nAll);

    const double phiWidth = (2 * TMath::Pi()) / nPhi;
    const double halfBin = phiWidth / 2.0;

    for (int i = 0; i < nAll; ++i)
    {
        TH2D *h2 = h_chi2_maps[i];
        if (!h2)
            continue;

        // locate minimum bin
        int minBinX = 0, minBinY = 0, dummyZ = 0;
        h2->GetMinimumBin(minBinX, minBinY, dummyZ);
        double minX = h2->GetXaxis()->GetBinCenter(minBinX);
        double minY = h2->GetYaxis()->GetBinCenter(minBinY);

        // 1σ error estimates
        auto contourX = getContour(h2);
        double left = contourX[minBinY - 1].first;
        double right = contourX[minBinY - 1].second;
        double errMaxX = std::max(fabs(minX - left), fabs(right - minX));

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
        double errMaxY = std::max(fabs(minY - down), fabs(up - minY));

        x_phi[i] = i;
        ex_phi[i] = 0;
        y_massShift[i] = minX;
        y_widthShift[i] = minY;
        ey_massShift[i] = errMaxX;
        ey_widthShift[i] = errMaxY;

        cout << "X is " << x_phi[i] << "Y is " << y_massShift[i] << "Y Error is " << ey_massShift[i] << endl;

        const char *tag = (i == nPhi ? "inclusive" : Form("%2d", i));
        // std::cout << Form("phi=%s (%+.3f) | mass=%.5f ± %.5f | width=%.5f ± %.5f\n",
        //                   tag, phiCenter, minX, errMaxX, minY, errMaxY);
    }

    // build graphs
    auto gMass = new TGraphErrors(nAll, x_phi.data(), y_massShift.data(),
                                  ex_phi.data(), ey_massShift.data());
    gMass->SetTitle("Best-fit mass shift vs #phi;#phi (rad);Mass shift (GeV)");
    gMass->SetMarkerStyle(20);
    gMass->SetMarkerColor(kRed + 1);
    gMass->SetLineColor(kRed + 1);
    gMass->SetLineWidth(2);

    auto gWidth = new TGraphErrors(nAll, x_phi.data(), y_widthShift.data(),
                                   ex_phi.data(), ey_widthShift.data());
    gWidth->SetTitle("Best-fit width smear vs #phi;#phi (rad);Width smear (GeV)");
    gWidth->SetMarkerStyle(21);
    gWidth->SetMarkerColor(kBlue + 1);
    gWidth->SetLineColor(kBlue + 1);
    gWidth->SetLineWidth(2);

    return {gMass, gWidth};
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
}

void get_scan_on_PbPb_and_pp_mc()
{
    gStyle->SetOptStat(0);
    TFile *PbPb_mc = new TFile("../ZBoson_18/rootfile/mc_signal.root", "READ");
    TFile *pp_mc = new TFile("./new_pp_data_file_stability_readonly.root", "READ");
    TFile *PbPb_mc_template = new TFile("../ZBoson_18/rootfile/version_2/template_PbPb_pp_mc_special.root", "READ");

    TH1D *h_PbPb_mc[11];
    for (int i = 0; i < 11; i++)
    {
        if (!((i < 4) || (i == 10)))
            continue;
        if (i == 10)
        {
            h_PbPb_mc[i] = (TH1D *)PbPb_mc->Get(Form("FA_nominal_%i", 4));
        }
        else
        {
            h_PbPb_mc[i] = (TH1D *)PbPb_mc->Get(Form("FA_nominal_%i", i));
        }

        areanormalize(h_PbPb_mc[i]);
    }
    TH1D *pp_mc_nominal_inclusive = (TH1D *)pp_mc->Get("pp_mc_FA_nominal_phi_plus_inclusive");
    TH1D *pp_mc_nominal_inclusive_1D_pT = (TH1D *)pp_mc->Get("pp_mc_FA_nominal_phi_plus_inclusive_1D_pT");
    TH1D *pp_mc_nominal_inclusive_no_pT = (TH1D *)pp_mc->Get("pp_mc_FA_nominal_phi_plus_inclusive_without_pT_reweight");

    areanormalize(pp_mc_nominal_inclusive);
    areanormalize(pp_mc_nominal_inclusive_1D_pT);
    areanormalize(pp_mc_nominal_inclusive_no_pT);

    TH1D *PbPb_mc_nominal_template[42][42][11];

    double highbin_mass_shift = 0.01;
    double highbin_smear = 0.01;
    double lowbin_mass_shift = -0.01;
    double lowbin_smear = -0.01;

    double h_low_mass_shift = lowbin_mass_shift - ((highbin_mass_shift - lowbin_mass_shift) / (42 - 1)) / 2;
    double h_high_mass_shift = highbin_mass_shift + ((highbin_mass_shift - lowbin_mass_shift) / (42 - 1)) / 2;

    double h_low_smear = lowbin_smear - ((highbin_smear - lowbin_smear) / (42 - 1)) / 2;
    double h_high_smear = highbin_smear + ((highbin_smear - lowbin_smear) / (42 - 1)) / 2;

    TH2D *h_chisquare_pp_PbPb = new TH2D("h_chisquare_pp_PbPb", "", 42, h_low_mass_shift, h_high_mass_shift, 42, h_low_smear, h_high_smear);
    TH2D *h_chisquare_pp_PbPb_1D_pT = new TH2D("h_chisquare_pp_PbPb_1D_pT", "", 42, h_low_mass_shift, h_high_mass_shift, 42, h_low_smear, h_high_smear);
    TH2D *h_chisquare_pp_PbPb_no_pT = new TH2D("h_chisquare_pp_PbPb_no_pT", "", 42, h_low_mass_shift, h_high_mass_shift, 42, h_low_smear, h_high_smear);

    TH2D *h_chisquare_PbPb_PbPb[11];

    for (int i = 0; i < 11; i++)
    {
        h_chisquare_PbPb_PbPb[i] = new TH2D(Form("h_chisquare_PbPb_PbPb_%i", i), "", 42, h_low_mass_shift, h_high_mass_shift, 42, h_low_smear, h_high_smear);
    }

    for (int cent = 0; cent < 11; cent++)
    {
        if (!((cent < 4) || (cent == 10)))
            continue;

        // Here's still old cent binning choice since special, I never re-run it.

        for (int i = 0; i < 42; i++)
        {
            for (int j = 0; j < 42; j++)
            {
                PbPb_mc_nominal_template[i][j][cent] = (TH1D *)PbPb_mc_template->Get(Form("template_FA_nominal_%i_%i_%i", i, j, cent));
                areanormalize(PbPb_mc_nominal_template[i][j][cent]);
            }
        }
    }

    for (int i = 0; i < 42; i++)
    {
        for (int j = 0; j < 42; j++)
        {
            double chisquarevalue_PbPb_pp = 0;
            double chisquarevalue_PbPb_pp_1D = 0;
            double chisquarevalue_PbPb_pp_no_pT = 0;

            chisquarevalue_PbPb_pp = myownfunctionchi2(pp_mc_nominal_inclusive, PbPb_mc_nominal_template[i][j][10]);
            chisquarevalue_PbPb_pp_1D = myownfunctionchi2(pp_mc_nominal_inclusive_1D_pT, PbPb_mc_nominal_template[i][j][10]);
            chisquarevalue_PbPb_pp_no_pT = myownfunctionchi2(pp_mc_nominal_inclusive_no_pT, PbPb_mc_nominal_template[i][j][10]);

            std::ostringstream stream_1;
            stream_1 << std::fixed << std::setprecision(6) << chisquarevalue_PbPb_pp;
            double formatted_PbPb_pp = std::stod(stream_1.str());

            std::ostringstream stream_2;
            stream_2 << std::fixed << std::setprecision(6) << chisquarevalue_PbPb_pp_1D;
            double formatted_PbPb_pp_1D = std::stod(stream_2.str());

            std::ostringstream stream_3;
            stream_3 << std::fixed << std::setprecision(6) << chisquarevalue_PbPb_pp_no_pT;
            double formatted_PbPb_pp_no_pT = std::stod(stream_3.str());

            h_chisquare_pp_PbPb->SetBinContent(i + 1, j + 1, formatted_PbPb_pp);
            h_chisquare_pp_PbPb_1D_pT->SetBinContent(i + 1, j + 1, formatted_PbPb_pp_1D);
            h_chisquare_pp_PbPb_no_pT->SetBinContent(i + 1, j + 1, formatted_PbPb_pp_no_pT);
        }
    }

    for (int cent = 0; cent < 11; cent++)
    {
        if (!((cent < 4) || (cent == 10)))
            continue;

        for (int i = 0; i < 42; i++)
        {
            for (int j = 0; j < 42; j++)
            {
                double chisquarevalue_PbPb_PbPb = 0;

                chisquarevalue_PbPb_PbPb = myownfunctionchi2(h_PbPb_mc[cent], PbPb_mc_nominal_template[i][j][cent]);

                std::ostringstream stream;
                stream << std::fixed << std::setprecision(4) << chisquarevalue_PbPb_PbPb;
                double formatted_PbPb_PbPb = std::stod(stream.str());

                h_chisquare_PbPb_PbPb[cent]->SetBinContent(i + 1, j + 1, formatted_PbPb_PbPb);
            }
        }
    }

    TCanvas *c_PbPb_PbPb = new TCanvas("c_PbPb_PbPb", "", 3200, 2400);
    TCanvas *c_PbPb_pp = new TCanvas("c_PbPb_pp", "", 3200, 2400);

    c_PbPb_PbPb->SetLeftMargin(0.10);   // smaller left margin
    c_PbPb_PbPb->SetRightMargin(0.12);  // smaller right margin
    c_PbPb_PbPb->SetTopMargin(0.08);    // smaller top margin
    c_PbPb_PbPb->SetBottomMargin(0.10); // smaller bottom margin

    c_PbPb_pp->SetLeftMargin(0.10);   // smaller left margin
    c_PbPb_pp->SetRightMargin(0.12);  // smaller right margin
    c_PbPb_pp->SetTopMargin(0.08);    // smaller top margin
    c_PbPb_pp->SetBottomMargin(0.10); // smaller bottom margin

    for (int i = 0; i < 11; i++)
    {
        if (!((i < 4) || (i == 10)))
            continue;
        c_PbPb_PbPb->cd();
        cosmetic(h_chisquare_PbPb_PbPb[i]);
        c_PbPb_PbPb->SaveAs(Form("./PbPb_pp_mc_compare/PbPb_fit_to_PbPb_mc_zoomin_%i.pdf", i));
        c_PbPb_PbPb->Clear();
    }

    c_PbPb_pp->cd();
    cosmetic(h_chisquare_pp_PbPb);
    outputresult(h_chisquare_pp_PbPb);
    c_PbPb_pp->SaveAs("./PbPb_pp_mc_compare/pp_fit_to_PbPb_mc_zoomin.pdf");
    c_PbPb_pp->Clear();

    c_PbPb_pp->cd();
    cosmetic(h_chisquare_pp_PbPb_1D_pT);
    outputresult(h_chisquare_pp_PbPb_1D_pT);

    c_PbPb_pp->SaveAs("./PbPb_pp_mc_compare/pp_fit_to_PbPb_mc_zoomin_1D_pT.pdf");
    c_PbPb_pp->Clear();

    c_PbPb_pp->cd();
    cosmetic(h_chisquare_pp_PbPb_no_pT);
    outputresult(h_chisquare_pp_PbPb_no_pT);

    c_PbPb_pp->SaveAs("./PbPb_pp_mc_compare/pp_fit_to_PbPb_mc_zoomin_no_pT.pdf");
    c_PbPb_pp->Clear();

    std::vector<TH2D *> vec_PbPb, vec_pp, vec_pp_1D_pT, vec_pp_no_pT;

    for (int i = 0; i < 11; i++)
    {
        if (!((i < 4) || (i == 10)))
            continue;
        vec_PbPb.push_back(h_chisquare_PbPb_PbPb[i]);
    }

    vec_pp.push_back(h_chisquare_pp_PbPb);
    vec_pp_1D_pT.push_back(h_chisquare_pp_PbPb_1D_pT);
    vec_pp_no_pT.push_back(h_chisquare_pp_PbPb_no_pT);

    auto [gMass_PbPb, gWidth_PbPb] = BuildMassAndWidthShiftGraphs(vec_PbPb);
    auto [gMass_pp, gWidth_pp] = BuildMassAndWidthShiftGraphs(vec_pp);
    auto [gMass_pp_1D_pT, gWidth_pp_1D_pT] = BuildMassAndWidthShiftGraphs(vec_pp_1D_pT);
    auto [gMass_pp_no_pT, gWidth_pp_no_pT] = BuildMassAndWidthShiftGraphs(vec_pp_no_pT);

    shiftpp(gMass_pp, 0);
    shiftpp(gMass_pp_1D_pT, -0.2);
    shiftpp(gMass_pp_no_pT, +0.2);
    shiftpp(gWidth_pp, 0);
    shiftpp(gWidth_pp_1D_pT, -0.2);
    shiftpp(gWidth_pp_no_pT, +0.2);

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
    gMass_pp->SetLineColor(kRed + 1);
    gMass_PbPb->SetLineColor(kBlue + 1);
    gMass_pp_no_pT->SetLineColor(kGreen + 2);
    gMass_pp_1D_pT->SetLineColor(kBlack);
    gMass_pp->SetMarkerColor(kRed + 1);
    gMass_PbPb->SetMarkerColor(kBlue + 1);
    gMass_pp_no_pT->SetMarkerColor(kGreen + 2);
    gMass_pp_1D_pT->SetMarkerColor(kBlack);
    gMass_PbPb->GetYaxis()->SetRangeUser(-0.05, 0.05);
    gMass_PbPb->GetXaxis()->SetLimits(-0.5, 5.5);
    gMass_PbPb->SetTitle(";Centrality;dMass (GeV)");
    gMass_PbPb->Draw("AP");
    // LabelLastPoint(gMass_PbPbPbPb, gMass_PbPbPbPb->GetMarkerColor());
    gMass_pp->Draw("Psame");
    // LabelLastPoint(gMass_ppPbPb, gMass_ppPbPb->GetMarkerColor());
    gMass_pp_no_pT->Draw("Psame");
    gMass_pp_1D_pT->Draw("Psame");
    auto legM = new TLegend(0.2, 0.75, 0.45, 0.9);
    legM->AddEntry(gMass_pp, "pp#rightarrowPbPb", "lp");
    legM->AddEntry(gMass_PbPb, "PbPb#rightarrowPbPb", "lp");
    legM->AddEntry(gMass_pp_no_pT, "pp no p_{T} reweight", "lp");
    legM->AddEntry(gMass_pp_1D_pT, "pp 1D p_{T} reweight", "lp");
    legM->SetBorderSize(0);
    legM->SetTextFont(42);
    legM->Draw();
    TLatex latexM;
    latexM.SetNDC();
    latexM.SetTextFont(62);
    latexM.SetTextSize(0.045);
    // FitAndAnnotateThreeGraphs(cMass, isflattened);
    //  latexM.DrawLatex(0.16, 0.93, "#bf{CMS}  #it{Preliminary}");
    cMass->SaveAs("./PbPb_pp_mc_compare/PbPb_and_pp_result_dM.pdf");

    // --- Separate canvas: Width smear ---
    TCanvas *cWidth = new TCanvas("cWidth", "Width smear vs phi", 800, 800);
    gWidth_pp->SetLineColor(kRed + 1);
    gWidth_PbPb->SetLineColor(kBlue + 1);
    gWidth_pp_no_pT->SetLineColor(kGreen + 2);
    gWidth_pp_1D_pT->SetLineColor(kBlack);
    gWidth_pp->SetMarkerColor(kRed + 1);
    gWidth_PbPb->SetMarkerColor(kBlue + 1);
    gWidth_pp_no_pT->SetMarkerColor(kGreen + 2);
    gWidth_pp_1D_pT->SetMarkerColor(kBlack);
    gWidth_PbPb->GetYaxis()->SetRangeUser(-0.05, 0.05);
    gWidth_PbPb->GetXaxis()->SetLimits(-0.5, 5.5);
    gWidth_PbPb->SetTitle(";Centrality; dWidth (GeV)");
    gWidth_PbPb->Draw("AP");
    // LabelLastPoint(gWidth_PbPbPbPb, gWidth_PbPbPbPb->GetMarkerColor());
    gWidth_pp->Draw("Psame");
    // LabelLastPoint(gWidth_ppPbPb, gWidth_ppPbPb->GetMarkerColor());
    gWidth_pp_no_pT->Draw("Psame");
    gWidth_pp_1D_pT->Draw("Psame");
    auto legW = new TLegend(0.2, 0.75, 0.45, 0.9);
    legW->AddEntry(gWidth_pp, "pp#rightarrowPbPb", "lp");
    legW->AddEntry(gWidth_PbPb, "PbPb#rightarrowPbPb", "lp");
    legW->AddEntry(gWidth_pp_no_pT, "pp no p_{T} reweight", "lp");
    legW->AddEntry(gWidth_pp_1D_pT, "pp 1D p_{T} reweight", "lp");
    legW->SetBorderSize(0);
    legW->SetTextFont(42);
    legW->Draw();
    TLatex latexW;
    latexW.SetNDC();
    latexW.SetTextFont(62);
    latexW.SetTextSize(0.045);
    // latexW.DrawLatex(0.16, 0.93, "#bf{CMS}  #it{Preliminary}");
    // FitAndAnnotateThreeGraphs(cWidth, isflattened, "width");
    cWidth->SaveAs("./PbPb_pp_mc_compare/PbPb_and_pp_result_dW.pdf");

    auto gMass_diff_pp = new TGraphErrors();
    auto gMass_diff_ppNoPT = new TGraphErrors();
    auto gMass_diff_pp1DPT = new TGraphErrors();
    auto gWidth_diff_pp = new TGraphErrors();
    auto gWidth_diff_ppNoPT = new TGraphErrors();
    auto gWidth_diff_pp1DPT = new TGraphErrors();

    int nPts = gMass_PbPb->GetN();
    gMass_diff_pp->Set(nPts);
    gMass_diff_ppNoPT->Set(nPts);
    gMass_diff_pp1DPT->Set(nPts);
    gWidth_diff_pp->Set(nPts);
    gWidth_diff_ppNoPT->Set(nPts);
    gWidth_diff_pp1DPT->Set(nPts);

    for (int i = 0; i < nPts; ++i)
    {
        double x = gMass_PbPb->GetX()[i];
        double ex = gMass_PbPb->GetEX()[i];

        double yPbPb = gMass_PbPb->GetY()[i];
        double yPP = gMass_pp->GetY()[0];
        double yPPnoPT = gMass_pp_no_pT->GetY()[0];
        double yPP1DPT = gMass_pp_1D_pT->GetY()[0];

        double ePbPb = gMass_PbPb->GetEY()[i];
        double ePP = gMass_pp->GetEY()[0];
        double ePPnoPT = gMass_pp_no_pT->GetEY()[0];
        double ePP1DPT = gMass_pp_1D_pT->GetEY()[0];

        // === MASS DIFF ===
        gMass_diff_pp->SetPoint(i, x, yPbPb - yPP);
        gMass_diff_pp->SetPointError(i, ex, std::hypot(ePbPb, ePP));

        gMass_diff_ppNoPT->SetPoint(i, x, yPbPb - yPPnoPT);
        gMass_diff_ppNoPT->SetPointError(i, ex, std::hypot(ePbPb, ePPnoPT));

        gMass_diff_pp1DPT->SetPoint(i, x, yPbPb - yPP1DPT);
        gMass_diff_pp1DPT->SetPointError(i, ex, std::hypot(ePbPb, ePP1DPT));

        // === WIDTH DIFF ===
        double wPbPb = gWidth_PbPb->GetY()[i];
        double wPP = gWidth_pp->GetY()[0];
        double wPPnoPT = gWidth_pp_no_pT->GetY()[0];
        double wPP1DPT = gWidth_pp_1D_pT->GetY()[0];

        double ewPbPb = gWidth_PbPb->GetEY()[i];
        double ewPP = gWidth_pp->GetEY()[0];
        double ewPPnoPT = gWidth_pp_no_pT->GetEY()[0];
        double ewPP1DPT = gWidth_pp_1D_pT->GetEY()[0];

        gWidth_diff_pp->SetPoint(i, x, wPbPb - wPP);
        gWidth_diff_pp->SetPointError(i, ex, std::hypot(ewPbPb, ewPP));

        gWidth_diff_ppNoPT->SetPoint(i, x, wPbPb - wPPnoPT);
        gWidth_diff_ppNoPT->SetPointError(i, ex, std::hypot(ewPbPb, ewPPnoPT));

        gWidth_diff_pp1DPT->SetPoint(i, x, wPbPb - wPP1DPT);
        gWidth_diff_pp1DPT->SetPointError(i, ex, std::hypot(ewPbPb, ewPP1DPT));
    }

    // style
    gMass_diff_pp->SetMarkerStyle(20);
    gMass_diff_pp->SetMarkerColor(kRed + 1);
    gMass_diff_pp->SetLineColor(kRed + 1);
    gMass_diff_pp->SetLineWidth(2);

    gMass_diff_ppNoPT->SetMarkerStyle(21);
    gMass_diff_ppNoPT->SetMarkerColor(kBlue + 1);
    gMass_diff_ppNoPT->SetLineColor(kBlue + 1);
    gMass_diff_ppNoPT->SetLineWidth(2);

    gMass_diff_pp1DPT->SetMarkerStyle(21);
    gMass_diff_pp1DPT->SetMarkerColor(kBlack);
    gMass_diff_pp1DPT->SetLineColor(kBlack);
    gMass_diff_pp1DPT->SetLineWidth(2);

    gWidth_diff_pp->SetMarkerStyle(20);
    gWidth_diff_pp->SetMarkerColor(kRed + 1);
    gWidth_diff_pp->SetLineColor(kRed + 1);
    gWidth_diff_pp->SetLineWidth(2);

    gWidth_diff_ppNoPT->SetMarkerStyle(21);
    gWidth_diff_ppNoPT->SetMarkerColor(kBlue + 1);
    gWidth_diff_ppNoPT->SetLineColor(kBlue + 1);
    gWidth_diff_ppNoPT->SetLineWidth(2);

    gWidth_diff_pp1DPT->SetMarkerStyle(21);
    gWidth_diff_pp1DPT->SetMarkerColor(kBlack);
    gWidth_diff_pp1DPT->SetLineColor(kBlack);
    gWidth_diff_pp1DPT->SetLineWidth(2);

    TCanvas *cMassDiff = new TCanvas("cMassDiff", "Mass difference vs phi", 800, 800);
    gMass_diff_pp->SetTitle(";Centrality; #Delta Mass (PbPb - pp) (GeV)");
    gMass_diff_pp->GetYaxis()->SetRangeUser(-0.05, 0.05);
    gMass_diff_pp->GetXaxis()->SetLimits(-0.5, 4.5);
    gMass_diff_pp->Draw("AP");
    // LabelLastPoint(gMass_diff_pp, gMass_diff_pp->GetMarkerColor());
    gMass_diff_ppNoPT->Draw("Psame");
    gMass_diff_pp1DPT->Draw("Psame");

    auto legMD = new TLegend(0.2, 0.75, 0.45, 0.9);
    legMD->AddEntry(gMass_diff_pp, "PbPb - pp", "lp");
    legMD->AddEntry(gMass_diff_ppNoPT, "PbPb - pp (no p_{T})", "lp");
    legMD->AddEntry(gMass_diff_pp1DPT, "PbPb - pp (1D p_{T})", "lp");
    legMD->SetBorderSize(0);
    legMD->SetTextFont(42);
    legMD->Draw();

    // FitAndAnnotateThreeGraphs(cMassDiff, isflattened, "I dont want to save");
    cMassDiff->SaveAs("./PbPb_pp_mc_compare/PbPb_-_pp_dM.pdf");

    TCanvas *cWidthDiff = new TCanvas("cWidthDiff", "Width difference vs phi", 800, 800);
    gWidth_diff_pp->SetTitle(";Centrality; #Delta Width (PbPb - pp) (GeV)");
    gWidth_diff_pp->GetYaxis()->SetRangeUser(-0.05, 0.05);
    gWidth_diff_pp->GetXaxis()->SetLimits(-0.5, 4.5);
    gWidth_diff_pp->Draw("AP");
    // LabelLastPoint(gWidth_diff_pp, gWidth_diff_pp->GetMarkerColor());
    gWidth_diff_ppNoPT->Draw("Psame");
    gWidth_diff_pp1DPT->Draw("Psame");

    auto legWD = new TLegend(0.2, 0.75, 0.45, 0.9);
    legWD->AddEntry(gWidth_diff_pp, "PbPb - pp", "lp");
    legWD->AddEntry(gWidth_diff_ppNoPT, "PbPb - pp (no p_{T})", "lp");
    legWD->AddEntry(gWidth_diff_pp1DPT, "PbPb - pp (1D p_{T})", "lp");
    legWD->SetBorderSize(0);
    legWD->SetTextFont(42);
    legWD->Draw();

    // FitAndAnnotateThreeGraphs(cWidthDiff, isflattened, "I dont want to save");
    cWidthDiff->SaveAs("./PbPb_pp_mc_compare/PbPb_-_pp_dW.pdf");

    TFile *send_to_frank = new TFile("./cos_fit_save/slides_13.root", "UPDATE");
    send_to_frank->cd();

    gMass_diff_pp->Write("2D_dM", 2);
    gMass_diff_pp1DPT->Write("1D_dM", 2);
    gMass_diff_ppNoPT->Write("nopT_dM", 2);
    gWidth_diff_pp->Write("2D_dW", 2);
    gWidth_diff_pp1DPT->Write("1D_dW", 2);
    gWidth_diff_ppNoPT->Write("nopT_dW", 2);

    /*TFile *send_to_frank = new TFile("./cos_fit_save/send_to_frank.root", "UPDATE");
    send_to_frank->cd();

    cMassDiff->Write(Form("Data_Mass_Diff%s", isflattened ? "_flatten" : ""), 2);
    cWidthDiff->Write(Form("Data_Width_Diff%s", isflattened ? "_flatten" : ""), 2);
    cWidth->Write(Form("Data_Width%s", isflattened ? "_flatten" : ""), 2);
    cMass->Write(Form("Data_Mass%s", isflattened ? "_flatten" : ""), 2);*/
}