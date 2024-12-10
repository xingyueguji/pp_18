void areanormalize(TH1D *h_1)
{
    double normalization_factor = h_1->Integral("width");
    h_1->Scale(1 / normalization_factor);
}
void mcvsdata(int type = 1, bool isbk = 0)
{
    // type 0 == raw , 1 == eta
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);

    TH1::SetDefaultSumw2();
    TString bkfilepath;
    TString datafilepath;
    TString mcfilepath;

    TH1D *h_mc_signal[21][21][11];
    TH1D *h_data[11];
    TH1D *h_mc_bk[11];
    TH1D *h_data_bksub[11];

    TFile *datafile;
    TFile *bkfile;
    TFile *mcfile;
    Int_t cenlowlimit[11] = {0, 10, 20, 30, 30, 0, 15, 50, 0, 14, 0};
    Int_t cenhighlimit[11] = {10, 20, 30, 100, 50, 15, 100, 100, 14, 100, 100};

    TCanvas *c_2d_chisquare_ndf[11];
    // What I need: invariant mass comparison data vs mc raw, data vs data bk sub.

    TCanvas *c_data_mc_raw[21][21][11];

    if (type == 0)
        bkfilepath = "rawfile.root";
    if (type == 1)
        bkfilepath = "etacut_file.root";

    datafilepath = "data_file.root";

    datafile = new TFile(datafilepath, "READ");
    bkfile = new TFile(bkfilepath, "READ");

    for (int cent = 0; cent < 11; cent++)
    {
        if (!((cent < 4) || (cent == 10)))
            continue;

        cout << "cent is " << cent << endl;

        if (type == 1)
        {
            if (cent == 0)
            {
                // mcfilepath = "../ZBoson_18/rootfile/eta_0_10_shift_-0.44_-0.33_smear_0_0.0085_modified_signal_21_21_100.root";
                mcfilepath = "../ZBoson_18/rootfile/eta_0_10_shift_-0.45_-0.05_smear_0_0.01_modified_signal_21_21_1000.root";
            }
            if (cent == 1)
            {
                mcfilepath = "../ZBoson_18/rootfile/eta_10_20_shift_-0.3_0.2_smear_0_0.01_modified_signal_21_21_1000.root";
            }
            if (cent == 2)
            {
                mcfilepath = "../ZBoson_18/rootfile/eta_20_30_shift_-0.45_0.1_smear_0_0.01_modified_signal_21_21_1000.root";
            }
            if (cent == 3)
            {
                mcfilepath = "../ZBoson_18/rootfile/eta_30_100_shift_-0.5_0.0_smear_0_0.01_modified_signal_21_21_1000.root";
            }
            if (cent == 10)
            {
                mcfilepath = "../ZBoson_18/rootfile/eta_0_100_shift_-0.34_-0.08_smear_0_0.007_modified_signal_21_21_1000.root";
            }
        }
        if (type == 0)
        {
            if (cent == 0)
            {
                mcfilepath = "../ZBoson_18/rootfile/raw_0_10_shift_-0.3_0.0_smear_0_0.008_modified_signal_21_21_1000.root";
            }
            if (cent == 1)
            {
                mcfilepath = "../ZBoson_18/rootfile/raw_10_20_shift_-0.4_-0.05_smear_0.001_0.012_modified_signal_21_21_1000.root";
            }
            if (cent == 2)
            {
                mcfilepath = "../ZBoson_18/rootfile/raw_20_30_shift_-0.4_0.0_smear_0.0_0.008_modified_signal_21_21_1000.root";
            }
            if (cent == 3)
            {
                mcfilepath = "../ZBoson_18/rootfile/raw_30_100_shift_-0.4_0.05_smear_0.002_0.013_modified_signal_21_21_1000.root";
            }
            if (cent == 10)
            {
                // This is the same as 0-10
                mcfilepath = "../ZBoson_18/rootfile/raw_0_10_shift_-0.3_0.0_smear_0_0.008_modified_signal_21_21_1000.root";
            }
        }

        mcfile = new TFile(mcfilepath, "READ");

        if (type == 1)
        {
        }
        if (type == 0)
        {
        }

        if (type == 0)
            h_data[cent] = (TH1D *)datafile->Get(Form("mass_array_data_%i", cent));
        if (type == 1)
            h_data[cent] = (TH1D *)datafile->Get(Form("mass_array_data_witheta_%i", cent));

        h_mc_bk[cent] = (TH1D *)bkfile->Get(Form("Normalized_mc_bk_%i", cent));

        areanormalize(h_data[cent]);

        h_data_bksub[cent] = (TH1D *)h_data[cent]->Clone();
        h_data_bksub[cent]->Add(h_mc_bk[cent], -1);

        // cout << "nbins is " << h_data_bksub[cent]->GetNbinsX() << endl;

        for (int shift = 0; shift < 21; shift++)
        {
            for (int smear = 0; smear < 21; smear++)
            {
                c_data_mc_raw[shift][smear][cent] = new TCanvas(Form("shift_%i_smear_%i_cent_%i", shift, smear, cent), "", 800, 800);
                if (type == 0)
                {
                    if (true)
                    {
                        h_mc_signal[shift][smear][cent] = (TH1D *)mcfile->Get(Form("modifiedmass_raw_without_eff_%i_%i_%i_new", shift, smear, cent));
                    }
                    else
                    {
                        h_mc_signal[shift][smear][cent] = (TH1D *)mcfile->Get(Form("modifiedmass_raw_without_eff_%i_%i_%i", shift, smear, cent));
                    }
                }
                if (type == 1)
                {
                    if (true)
                    {
                        h_mc_signal[shift][smear][cent] = (TH1D *)mcfile->Get(Form("modifiedmass_eta_without_eff_%i_%i_%i_new", shift, smear, cent));
                    }
                    else
                    {
                        h_mc_signal[shift][smear][cent] = (TH1D *)mcfile->Get(Form("modifiedmass_eta_without_eff_%i_%i_%i", shift, smear, cent));
                    }
                }
                areanormalize(h_mc_signal[shift][smear][cent]);
            }
        }
    }

    for (int cent = 0; cent < 11; cent++)
    {
        if (!((cent < 4) || (cent == 10)))
            continue;

        for (int shift = 0; shift < 21; shift++)
        {
            for (int smear = 0; smear < 21; smear++)
            {
                c_data_mc_raw[shift][smear][cent]->cd();
                c_data_mc_raw[shift][smear][cent]->SetLogy();

                if (isbk)
                {
                    // h_data_bksub[cent]->SetTitle(Form(data_mc_title, this->cenlowlimit[cent], this->cenhighlimit[cent]));
                    h_data_bksub[cent]->SetMarkerColor(kRed);
                    h_data_bksub[cent]->SetMarkerStyle(kFullCircle);
                    h_data_bksub[cent]->Draw("A P");
                    h_mc_signal[shift][smear][cent]->SetMarkerColor(kGreen);
                    h_mc_signal[shift][smear][cent]->SetMarkerStyle(kFullDotLarge);
                    h_mc_signal[shift][smear][cent]->Draw("P SAME");
                }
                if (!isbk)
                {
                    // h_data[cent]->SetTitle(Form(data_mc_title, this->cenlowlimit[cent], this->cenhighlimit[cent]));
                    h_data[cent]->SetMarkerColor(kRed);
                    h_data[cent]->SetMarkerStyle(kFullCircle);
                    h_data[cent]->Draw("P");
                    h_mc_signal[shift][smear][cent]->SetMarkerColor(kGreen);
                    h_mc_signal[shift][smear][cent]->SetMarkerStyle(kFullDotLarge);
                    h_mc_signal[shift][smear][cent]->Draw("P SAME");
                }

                TPaveText *pt1 = new TPaveText(0.1, 0.8, 0.5, 0.9, "NDC");
                pt1->AddText("Red is data, green is template");
                pt1->SetTextSize(0.04);
                pt1->SetTextAlign(22); // Center alignment
                pt1->Draw();

                c_data_mc_raw[shift][smear][cent]->SaveAs(Form("./syst/cent_%i_%i/data_mc_%i_%i.png", cenlowlimit[cent], cenhighlimit[cent], shift, smear));
            }
        }
    }
}