#include "plotting_helper.h"
void pTreweight_2D()
{

    plotting_helper *ovo = new plotting_helper();

    TFile *f1 = new TFile("../ZBoson_18/rootfile/pt_File.root", "UPDATE");

    TH2D *PbPb_FA = (TH2D *)f1->Get("pT_y_spec_FA");
    TH2D *pp_FA = (TH2D *)f1->Get("pT_y_spec_pp_FA");

    ovo->areanormalize(PbPb_FA);
    ovo->areanormalize(pp_FA);

    TH2D *ratio_FA = (TH2D *)PbPb_FA->Clone("ratio_FA");
    ratio_FA->SetTitle("PbPb / pp Ratio");
    ratio_FA->GetXaxis()->SetTitle("y");
    ratio_FA->GetYaxis()->SetTitle("pT (GeV)");

    ratio_FA->Divide(pp_FA);

    ratio_FA->SetStats(false);
    ratio_FA->GetXaxis()->SetTitle("y");
    ratio_FA->GetYaxis()->SetTitle("p_{T} [GeV]");

    TCanvas *c1 = new TCanvas("c1", "Ratio Plot", 800, 700);
    ratio_FA->Draw("COLZ");

    c1->SaveAs("./pTratio/2D_ratio_plot_FA_data.pdf");
    f1->cd();
    ratio_FA->Write("FA_2D_ratio_data", 2);
    f1->Close();
}