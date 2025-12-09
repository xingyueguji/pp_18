#include "plotting_helper.h"
void pTreweight_2D_mc()
{

    plotting_helper *ovo = new plotting_helper();

    TFile *f1 = new TFile("../ZBoson_18/rootfile/mc_pTratio.root", "UPDATE");

    TH2D *PbPb_FA = (TH2D *)f1->Get("mc_PbPb_pT_y");
    TH2D *pp_FA = (TH2D *)f1->Get("mc_pp_pT_y");

    ovo->areanormalize(PbPb_FA);
    ovo->areanormalize(pp_FA);

    TH2D *ratio_FA = (TH2D *)PbPb_FA->Clone("ratio_FA");
    ratio_FA->SetTitle("PbPb / pp Ratio (mc)");
    ratio_FA->Divide(pp_FA);

    ratio_FA->SetStats(false);
    ratio_FA->GetXaxis()->SetTitle("y");
    ratio_FA->GetYaxis()->SetTitle("p_{T} [GeV]");

    TCanvas *c1 = new TCanvas("c1", "Ratio Plot", 800, 700);
    ratio_FA->Draw("COLZ");

    c1->SaveAs("./pTratio/2D_ratio_plot_FA_mc.pdf");

    f1->cd();
    ratio_FA->Write("FA_2D_ratio_mc", 2);
    f1->Close();
}