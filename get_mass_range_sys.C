void areanormalize(TH1D *h_1)
{
    if (h_1 == nullptr)
    {
        return;
    }
    double normalization_factor = h_1->Integral("width");
    h_1->Scale(1 / normalization_factor);
}

TH1D *CloneWithNewRangeAndBins(const TH1D *h_orig, const char *newname = "h_new")
{
    if (!h_orig)
        return nullptr;

    // Define new histogram parameters
    const int nBins_new = 20;
    const double xMin_new = 70.0;
    const double xMax_new = 110.0;

    // Create the new histogram
    TH1D *h_new = new TH1D(newname, h_orig->GetTitle(), nBins_new, xMin_new, xMax_new);
    h_new->SetDirectory(0); // Detach from file/directory

    // Fill the new histogram by sampling from the original
    for (int i = 1; i <= h_orig->GetNbinsX(); ++i)
    {
        double x = h_orig->GetBinCenter(i);
        double y = h_orig->GetBinContent(i);
        double err = h_orig->GetBinError(i);

        if (x >= xMin_new && x < xMax_new)
        {
            int new_bin = h_new->FindBin(x);
            h_new->SetBinContent(new_bin, y);
            h_new->SetBinError(new_bin, err);
        }
    }

    return h_new;
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
            //std::cerr << "Warning: Bin " << i << " has zero errors/missing bins/ missing bins after bksub for data. Not Skipped" << std::endl;
            // continue; // Skip this bin to avoid division by zero
        }

        chi2 += num / den;
    }

    Double_t ndf = nbinsh1 - 2;

    return chi2;
}

void get_mass_range_sys()
{
    const int nbins_cent = 11;
    const int ntrials = 10000;
    const int nbins_mass_shift = 42;
    const int nbins_smear = 42;

    TFile *f_pseudo_experiment = new TFile("./bestfittemplaterootfile/pseudo_experiment.root", "READ");
    TFile *f_PbPb_template = new TFile("../ZBoson_18/rootfile/version_2/template_PbPb.root", "READ");
    TFile *f_pp_template = new TFile("../ZBoson_18/rootfile/version_2/template_pp_bk.root", "READ");
    //TFile *f_background = new TFile("../ZBoson_18/rootfile/normalized/FA_nominal.root", "READ");

    TH1D *h_PbPb_pseudo_data[nbins_cent][ntrials];
    TH1D *h_PbPb_pseudo_data_reduced[nbins_cent][ntrials];
    TH1D *h_pp_pseudo_data[ntrials];
    TH1D *h_PbPb_template[nbins_mass_shift][nbins_smear][nbins_cent];
    TH1D *h_PbPb_template_reduced[nbins_mass_shift][nbins_smear][nbins_cent];
    TH1D *h_pp_template[nbins_mass_shift][nbins_smear];
    TH1D *h_mc_bk[nbins_cent];

    // Nasty loop for reading TH1Ds.
    for (int cent = 0; cent < nbins_cent; cent++)
    {
        if (!((cent < 4) || (cent == 10)))
            continue;

        // I don't think I need to do bk sub again
        // h_mc_bk[nbins_cent] = (TH1D *)f_background->Get(Form("Normalized_mc_bk_%i", cent));

        for (int mass = 0; mass < nbins_mass_shift; mass++)
        {
            for (int smear = 0; smear < nbins_smear; smear++)
            {
                h_PbPb_template[mass][smear][cent] = (TH1D *)f_PbPb_template->Get(Form("template_FA_nominal_%i_%i_%i", mass, smear, cent));
                areanormalize(h_PbPb_template[mass][smear][cent]);
                h_PbPb_template[mass][smear][cent]->Rebin(4);
                h_PbPb_template_reduced[mass][smear][cent] = (TH1D *)CloneWithNewRangeAndBins(h_PbPb_template[mass][smear][cent], Form("h_PbPb_template_reduced_%i_%i_%i", mass, smear, cent));

                if (cent == 0)
                {
                    //h_pp_template[mass][smear] = (TH1D *)f_pp_template->Get(Form("template_FA_nominal_%i_%i_%i", mass, smear, 10));
                    //areanormalize(h_pp_template[mass][smear]);
                    //h_pp_template[mass][smear]->Rebin(4);
                }
                for (int trials = 0; trials < ntrials; trials++)
                {
                    if (mass == 0 && smear == 0)
                    {
                        h_PbPb_pseudo_data[cent][trials] = (TH1D *)f_pseudo_experiment->Get(Form("pseudo_cent_%i_trial_%i", cent, trials));
                        areanormalize(h_PbPb_pseudo_data[cent][trials]);
                        h_PbPb_pseudo_data[cent][trials]->Rebin(4);
                        h_PbPb_pseudo_data_reduced[cent][trials] = (TH1D *)CloneWithNewRangeAndBins(h_PbPb_pseudo_data[cent][trials], Form("h_PbPb_pseudo_data_%i_%i", cent, trials));

                        if (cent == 0)
                        {
                            //h_pp_pseudo_data[trials] = (TH1D *)f_pseudo_experiment->Get(Form("pseudo_pp_trial_%i", trials));
                            //areanormalize(h_pp_pseudo_data[trials]);
                            //h_pp_pseudo_data[trials]->Rebin(4);
                        }
                    }
                }
            }
        }
    }

    // Here's getting the reduced version

    // Now start calculating chi2, find the minimum, get the corresponding mass shift and width shift.
    // PbPb

    Double_t minimumvalue[nbins_cent][ntrials];
    Double_t minimumvalue_mass[nbins_cent][ntrials];
    Double_t minimumvalue_width[nbins_cent][ntrials];

    Double_t minimumvalue_reduced[nbins_cent][ntrials];
    Double_t minimumvalue_mass_reduced[nbins_cent][ntrials];
    Double_t minimumvalue_width_reduced[nbins_cent][ntrials];

    double raw_mass_shift_array_low[nbins_cent] = {-0.24, -0.32, -0.32, -0.25, 0, 0, 0, 0, 0, 0, -0.22};
    double raw_mass_shift_array_high[nbins_cent] = {-0.02, -0.06, -0.04, 0.02, 0, 0, 0, 0, 0, 0, -0.08};
    double raw_mass_smear_array_low[nbins_cent] = {-0.05, 0.075, -0.2, 0.15, 0, 0, 0, 0, 0, 0, 0.13};
    double raw_mass_smear_array_high[nbins_cent] = {0.3, 0.5, 0.3, 0.65, 0, 0, 0, 0, 0, 0, 0.36};

    for (int i = 0; i < ntrials; i++)
    {
        for (int cent = 0; cent < nbins_cent; cent++)
        {
            minimumvalue[cent][i] = 9999;
            minimumvalue_mass[cent][i] = 9999;
            minimumvalue_width[cent][i] = 9999;

            minimumvalue_reduced[cent][i] = 9999;
            minimumvalue_mass_reduced[cent][i] = 9999;
            minimumvalue_width_reduced[cent][i] = 9999;
        }
    }
    for (int trials = 0; trials < ntrials; trials++)
    {
        for (int cent = 0; cent < nbins_cent; cent++)
        {
            if (!((cent < 4) || (cent == 10)))
                continue;

            for (int mass = 0; mass < nbins_mass_shift; mass++)
            {
                for (int smear = 0; smear < nbins_smear; smear++)
                {
                    double chi2value = myownfunctionchi2(h_PbPb_pseudo_data[cent][trials], h_PbPb_template[mass][smear][cent]);
                    double chi2value_reduced = myownfunctionchi2(h_PbPb_pseudo_data_reduced[cent][trials], h_PbPb_template_reduced[mass][smear][cent]);
                    if (chi2value < minimumvalue[cent][trials])
                    {
                        minimumvalue[cent][trials] = chi2value;
                        minimumvalue_mass[cent][trials] = mass;
                        minimumvalue_width[cent][trials] = smear;
                    }
                    if (chi2value_reduced < minimumvalue_reduced[cent][trials])
                    {
                        minimumvalue_reduced[cent][trials] = chi2value_reduced;
                        minimumvalue_mass_reduced[cent][trials] = mass;
                        minimumvalue_width_reduced[cent][trials] = smear;
                    }
                }
            }
        }
    }

    TH1D *h_PbPb_mass_distribution[nbins_cent];
    TH1D *h_PbPb_width_distribution[nbins_cent];
    TH1D *h_PbPb_mass_distribution_reduced[nbins_cent];
    TH1D *h_PbPb_width_distribution_reduced[nbins_cent];
    TH1D *h_PbPb_mass_diff[nbins_cent];
    TH1D *h_PbPb_width_diff[nbins_cent];

    TFile *f1 = new TFile("./bestfittemplaterootfile/pseudoresult.root", "RECREATE");
    f1->cd();

    for (int cent = 0; cent < nbins_cent; cent++)
    {
        if (!((cent < 4) || (cent == 10)))
            continue;

        h_PbPb_mass_distribution[cent] = new TH1D(Form("h_PbPb_mass_distribution_%i", cent), "", 42, 0, 42);
        h_PbPb_width_distribution[cent] = new TH1D(Form("h_PbPb_width_distribution_%i", cent), "", 42, 0, 42);
        h_PbPb_mass_distribution_reduced[cent] = new TH1D(Form("h_PbPb_mass_distribution_reduced_%i", cent), "", 42, 0, 42);
        h_PbPb_width_distribution_reduced[cent] = new TH1D(Form("h_PbPb_width_distribution_reduced_%i", cent), "", 42, 0, 42);

        for (int trials = 0; trials < ntrials; trials++)
        {
            h_PbPb_mass_distribution[cent]->Fill(minimumvalue_mass[cent][trials]);
            h_PbPb_width_distribution[cent]->Fill(minimumvalue_width[cent][trials]);
            h_PbPb_mass_distribution_reduced[cent]->Fill(minimumvalue_mass_reduced[cent][trials]);
            h_PbPb_width_distribution_reduced[cent]->Fill(minimumvalue_width_reduced[cent][trials]);
        }

        h_PbPb_mass_diff[cent] = (TH1D *)h_PbPb_mass_distribution[cent]->Clone(Form("h_PbPb_mass_diff_%i", cent));
        h_PbPb_width_diff[cent] = (TH1D *)h_PbPb_width_distribution[cent]->Clone(Form("h_PbPb_width_diff_%i", cent));

        h_PbPb_mass_diff[cent]->Add(h_PbPb_mass_distribution_reduced[cent], -1);
        h_PbPb_width_diff[cent]->Add(h_PbPb_width_distribution_reduced[cent], -1);

        h_PbPb_mass_distribution[cent]->Write("", 2);
        h_PbPb_width_distribution[cent]->Write("", 2);
        h_PbPb_mass_distribution_reduced[cent]->Write("", 2);
        h_PbPb_width_distribution_reduced[cent]->Write("", 2);
        h_PbPb_mass_diff[cent]->Write("",2);
        h_PbPb_width_diff[cent]->Write("",2);

    }
}