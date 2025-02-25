#include "plotting_helper.h"

void pTreweight()
{

    plotting_helper *ovo = new plotting_helper();

    TFile *f1 = new TFile("../ZBoson_18/rootfile/pt_File.root", "UPDATE");

    TH1D *PbPb_FA = (TH1D *)f1->Get("pT_spec_FA");
    TH1D *pp_FA = (TH1D *)f1->Get("pT_spec_pp_FA");

    TH1D *PbPb_Eta = (TH1D *)f1->Get("pT_spec_Eta");
    TH1D *pp_Eta = (TH1D *)f1->Get("pT_spec_pp_Eta");

    ovo->areanormalize(PbPb_FA);
    ovo->areanormalize(pp_FA);
    ovo->areanormalize(PbPb_Eta);
    ovo->areanormalize(pp_Eta);

    // Clone PbPb histograms to store weights
    TH1D *weight_FA = new TH1D("weight_FA", "", 200, 0, 200);
    TH1D *weight_Eta = new TH1D("weight_Eta", "", 200, 0, 200);

    // Compute weights for FA
    for (int i = 1; i <= PbPb_FA->GetNbinsX(); i++)
    {
        double PbPb_val = PbPb_FA->GetBinContent(i);
        double pp_val = pp_FA->GetBinContent(i);

        if (pp_val > 0) // Avoid division by zero
        {
            double ratio = PbPb_val / pp_val;
            weight_FA->SetBinContent(i, ratio);
        }
        else
        {
            weight_FA->SetBinContent(i, 0); // Set weight to 0 if denominator is 0
        }
    }

    // Compute weights for Eta
    for (int i = 1; i <= PbPb_Eta->GetNbinsX(); i++)
    {
        double PbPb_val = PbPb_Eta->GetBinContent(i);
        double pp_val = pp_Eta->GetBinContent(i);

        if (pp_val > 0)
        {
            double ratio = PbPb_val / pp_val;

            weight_Eta->SetBinContent(i, ratio);
        }
        else
        {
            weight_Eta->SetBinContent(i, 0);
        }
    }

    // Save to output file
    f1->cd();
    weight_FA->Write("", 2);
    weight_Eta->Write("", 2);

    std::cout << "pT reweighting histograms saved in pT_reweight.root" << std::endl;
}