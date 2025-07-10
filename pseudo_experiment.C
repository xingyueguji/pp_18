void pseudo_experiment()
{
    const int ntrials = 10000;
    TFile *f1 = new TFile("./bestfittemplaterootfile/template.root", "READ");
    f1->cd();

    TFile *PbPb_data = new TFile("../ZBoson_18/rootfile/data_file.root", "READ");
    TFile *pp_data = new TFile("./new_pp_data_file_stability_readonly.root", "READ");

    TH1D *h_PbPb_data[11];
    TH1D *h_PbPb_besttemplate[11];
    TH1D *h_pp_data;
    TH1D *h_pp_besttemplate;

    Int_t nentries_PbPb[11] = {};
    Int_t nentries_pp = 0;

    // Get nEntries first

    for (int cent = 0; cent < 11; cent++)
    {
        if (!((cent < 4) || (cent == 10)))
            continue;

        h_PbPb_data[cent] = (TH1D *)PbPb_data->Get(Form("FA_nominal_%i", cent));
        nentries_PbPb[cent] = h_PbPb_data[cent]->GetEntries();
        cout << "PbPb cent " << cent << " has " << nentries_PbPb[cent] << " entries" << endl;
    }

    h_pp_data = (TH1D *)pp_data->Get("FA_nominal_inclusive");
    nentries_pp = h_pp_data->GetEntries();
    cout << "pp inclusive" << " has " << nentries_pp << " entries" << endl;

    // READ best template

    h_PbPb_besttemplate[0] = (TH1D *)f1->Get("template_FA_nominal_zoomin_clone_26_17_0");
    h_PbPb_besttemplate[1] = (TH1D *)f1->Get("template_FA_nominal_zoomin_clone_20_16_1");
    h_PbPb_besttemplate[2] = (TH1D *)f1->Get("template_FA_nominal_zoomin_clone_16_16_2");
    h_PbPb_besttemplate[3] = (TH1D *)f1->Get("template_FA_nominal_zoomin_clone_20_16_3");
    h_PbPb_besttemplate[10] = (TH1D *)f1->Get("template_FA_nominal_zoomin_clone_18_26_10");
    h_pp_besttemplate = (TH1D *)f1->Get("h_mc_signal_pp_zoomin_not_rebinned_clone_14_16");

    TH1D *h_pseudo_PbPb[11][ntrials];
    TH1D *h_pseudo_pp[ntrials];

    for (int cent = 0; cent < 11; cent++)
    {
        if (!((cent < 4) || (cent == 10)))
            continue;

        if (!h_PbPb_besttemplate[cent])
            continue; // Safety check

        // Normalize a clone of the original to get a PDF
        TH1D *pdf = (TH1D *)h_PbPb_besttemplate[cent]->Clone(Form("pdf_cent_%i", cent));
        pdf->SetDirectory(0);                     // Detach from TFile
        pdf->Scale(1.0 / pdf->Integral("width")); // Normalize shape

        for (int i = 0; i < ntrials; i++)
        {
            cout << "we are at PbPb trial " << i << " for cent " << cent << endl;
            h_pseudo_PbPb[cent][i] = (TH1D *)h_PbPb_besttemplate[cent]->Clone(Form("pseudo_cent_%i_trial_%i", cent, i));
            h_pseudo_PbPb[cent][i]->Reset();         // Clear contents
            h_pseudo_PbPb[cent][i]->SetDirectory(0); // Optional: avoid attaching to gDirectory

            for (int j = 0; j < nentries_PbPb[cent]; j++)
            {
                // cout << nentries_PbPb[cent] << endl;
                double val = pdf->GetRandom(); // Randomly sample from shape
                h_pseudo_PbPb[cent][i]->Fill(val);
            }
        }

        delete pdf; // Clean up normalized shape
    }

    if (!h_pp_besttemplate)
        cout << "pp template read error" << endl; // Safety check

    // Normalize a clone of the original to get a PDF
    TH1D *pdf = (TH1D *)h_pp_besttemplate->Clone("pdf_pp");
    pdf->SetDirectory(0);                     // Detach from TFile
    pdf->Scale(1.0 / pdf->Integral("width")); // Normalize shape

    /*for (int i = 0; i < ntrials; i++)
    {
        cout << "we are at pp trial " << i << endl;
        h_pseudo_pp[i] = (TH1D *)h_pp_besttemplate->Clone(Form("pseudo_pp_trial_%i", i));
        h_pseudo_pp[i]->Reset();         // Clear contents
        h_pseudo_pp[i]->SetDirectory(0); // Optional: avoid attaching to gDirectory

        for (int j = 0; j < nentries_pp; j++)
        {
            double val = pdf->GetRandom(); // Randomly sample from shape
            h_pseudo_pp[i]->Fill(val);
        }
    }*/
    //delete pdf;

    TFile *f2 = new TFile("./bestfittemplaterootfile/pseudo_experiment.root", "RECREATE");
    f2->cd();

    for (int cent = 0; cent < 11; cent++)
    {
        if (!((cent < 4) || (cent == 10)))
            continue;

        if (cent == 0)
        {
            for (int i = 0; i < ntrials; i++)
            {
                h_pseudo_PbPb[cent][i]->Write("", 2);
                //h_pseudo_pp[i]->Write("", 2);
            }
        }
        else
        {
            for (int i = 0; i < ntrials; i++)
            {
                h_pseudo_PbPb[cent][i]->Write("", 2);
            }
        }
    }
}