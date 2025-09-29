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
    TFile *PbPb_mc_template = new TFile("../ZBoson_18/rootfile/version_2/template_PbPb_pp_mc_special_zoomin.root", "READ");

    TH1D *PbPb_mc_nominal_inclusive = (TH1D *)PbPb_mc->Get("FA_nominal_10");
    TH1D *pp_mc_nominal_inclusive = (TH1D *)pp_mc->Get("pp_mc_inclusive_test_with_pt_reweight");
    TH1D *PbPb_mc_nominal_template[42][42];

    double highbin_mass_shift = 0.03;
    double highbin_smear = 0.03;
    double lowbin_mass_shift = -0.03;
    double lowbin_smear = -0.03;

    double h_low_mass_shift = lowbin_mass_shift - ((highbin_mass_shift - lowbin_mass_shift) / (42 - 1)) / 2;
    double h_high_mass_shift = highbin_mass_shift + ((highbin_mass_shift - lowbin_mass_shift) / (42 - 1)) / 2;

    double h_low_smear = lowbin_smear - ((highbin_smear - lowbin_smear) / (42 - 1)) / 2;
    double h_high_smear = highbin_smear + ((highbin_smear - lowbin_smear) / (42 - 1)) / 2;

    TH2D *h_chisquare_pp_PbPb = new TH2D("h_chisquare_pp_PbPb", "", 42, h_low_mass_shift, h_high_mass_shift, 42, h_low_smear, h_high_smear);
    TH2D *h_chisquare_PbPb_PbPb = new TH2D("h_chisquare_PbPb_PbPb", "", 42, h_low_mass_shift, h_high_mass_shift, 42, h_low_smear, h_high_smear);

    for (int i = 0; i < 42; i++)
    {
        for (int j = 0; j < 42; j++)
        {
            PbPb_mc_nominal_template[i][j] = (TH1D *)PbPb_mc_template->Get(Form("template_FA_nominal_%i_%i_%i", i, j, 10));
            areanormalize(PbPb_mc_nominal_template[i][j]);
        }
    }
    areanormalize(PbPb_mc_nominal_inclusive);
    areanormalize(pp_mc_nominal_inclusive);

    for (int i = 0; i < 42; i++)
    {
        for (int j = 0; j < 42; j++)
        {
            double chisquarevalue_PbPb_PbPb = 0;
            double chisquarevalue_PbPb_pp = 0;

            chisquarevalue_PbPb_PbPb = myownfunctionchi2(PbPb_mc_nominal_inclusive, PbPb_mc_nominal_template[i][j]);
            chisquarevalue_PbPb_pp = myownfunctionchi2(pp_mc_nominal_inclusive, PbPb_mc_nominal_template[i][j]);

            std::ostringstream stream;
            stream << std::fixed << std::setprecision(4) << chisquarevalue_PbPb_PbPb;
            double formatted_PbPb_PbPb = std::stod(stream.str());

            std::ostringstream stream_1;
            stream_1 << std::fixed << std::setprecision(4) << chisquarevalue_PbPb_pp;
            double formatted_PbPb_pp = std::stod(stream_1.str());

            h_chisquare_pp_PbPb->SetBinContent(i + 1, j + 1, formatted_PbPb_pp);
            h_chisquare_PbPb_PbPb->SetBinContent(i + 1, j + 1, formatted_PbPb_PbPb);
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

    cosmetic(h_chisquare_PbPb_PbPb);

    c_PbPb_pp->cd();
    h_chisquare_pp_PbPb->Draw("COLZ");
    cosmetic(h_chisquare_pp_PbPb);
    c_PbPb_PbPb->Update();

    c_PbPb_PbPb->cd();
    h_chisquare_PbPb_PbPb->Draw("COLZ");
    cosmetic(h_chisquare_PbPb_PbPb);
    c_PbPb_PbPb->Update();

    c_PbPb_pp->SaveAs("./PbPb_pp_mc_compare/pp_fit_to_PbPb_mc_zoomin_ptReweight.pdf");
    c_PbPb_PbPb->SaveAs("./PbPb_pp_mc_compare/PbPb_fit_to_PbPb_mc_zoomin_ptReweight.pdf");
}