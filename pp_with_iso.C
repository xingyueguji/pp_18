#define pp_with_iso_cxx
#include "pp_with_iso.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TMath.h>
#include <iostream>
#include <bitset>

using std::cout;
using std::endl;

double getWeightFromHist(TF1 *weightfucntion, double pt)
{
   double ratio = weightfucntion->Eval(pt);
   return ratio;
}

// Relative PF isolation for muon at index i.
// Adapted from the std::vector<float>* version to work with the
// fixed-size arrays in this tree (Reco_mu_PF*Iso[Reco_mu_size]).
// pt is passed in directly since this tree stores muon kinematics
// inside a TClonesArray (Reco_mu_4mom), not a Float_t pT array.
// Returns (chargedHad + neutralHad + photon) / pt.
static double RelIsoPF(int i,
                       double pt,
                       const Float_t *muPFChIso,
                       const Float_t *muPFNeuIso,
                       const Float_t *muPFPhoIso,
                       const Float_t *muPFPUIso)
{
   if (!muPFChIso || !muPFNeuIso || !muPFPhoIso || !muPFPUIso)
      return 999.0;
   if (pt <= 0)
      return 999.0;
   const double neutralSum = muPFNeuIso[i] + muPFPhoIso[i] - 0.5 * muPFPUIso[i];
   const double isoAbs = muPFChIso[i] + std::max(0.0, neutralSum);
   return isoAbs / pt;
}

Int_t GetPhiBin(double phi, const int numberofphibins)
{
   while (phi <= -TMath::Pi())
      phi += 2 * TMath::Pi();
   while (phi > TMath::Pi())
      phi -= 2 * TMath::Pi();

   double xmin = -TMath::Pi();
   double xmax = TMath::Pi();
   double width = (xmax - xmin) / numberofphibins;

   int bin = static_cast<int>((phi - xmin) / width);

   if (bin < 0)
      bin = 0;
   if (bin >= numberofphibins)
   {
      bin = numberofphibins - 1;
      cout << "Max range exceeded? you should not able to see this" << endl;
   }
   return bin;
}

// Helper used to be in SkimNew.h; here we put it in the .C since
// pp_with_iso.h was auto-generated and doesn't declare it.
static double getEfficiency(TEfficiency *e, double y, double pt)
{
   double originalPt = pt;
   if (pt >= 200)
      pt = 199.9;

   int bin = e->FindFixBin(y, pt);
   float efficiency = e->GetEfficiency(bin);

   if (efficiency > 0 && efficiency <= 1)
      return efficiency;
   std::cout << "efficiency not in the range [0,1], returning 1!" << std::endl;
   std::cout << "Rapidity: " << y << " Pt: " << originalPt << std::endl;
   return 1;
}

static Bool_t CheckTrigBit(ULong64_t num, int bitPosition)
{
   return ((num >> bitPosition) & 1ULL) == 1ULL;
}

void pp_with_iso::Loop(bool useIsolation)
{
   // Standard CMS tight-muon PF relative-isolation working point.
   // Tweak this if you want a different threshold.
   const double kRelIsoCut = 0.2;

   if (fChain == 0)
      return;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   const int numberofphibins = 8;

   Double_t runupperlimit[22] = {320916, 321012, 321140, 321232, 321393, 321436, 321735, 321820, 321909, 322014, 322118, 322319, 322381, 322617, 323525, 323778, 324021, 324245, 324772, 324897, 325001, 325175};
   Double_t runlowerlimit[22] = {320500, 320917, 321051, 321149, 321233, 321396, 321457, 321755, 321831, 321917, 322022, 322179, 322322, 322407, 322625, 323526, 323790, 324022, 324293, 324785, 324970, 325022};

   TH1D *FA_nominal[22];
   TH1D *FA_AcoUp[22];
   TH1D *FA_AcoDown[22];
   TH1D *FA_tnpU[22];
   TH1D *FA_tnpD[22];
   TH1D *FA_mass_range[22];
   TH1D *FA_1D_pT[22];
   TH1D *FA_nominal_inclusive = new TH1D("FA_nominal_inclusive", "", 120, 60, 120);
   TH1D *FA_nominal_inclusive_antiiso = new TH1D("FA_nominal_inclusive_antiiso", "", 120, 60, 120);
   TH1D *FA_nominal_inclusive_without_eff = new TH1D("FA_nominal_inclusive_without_eff", "", 120, 60, 120);
   TH1D *FA_nominal_inclusive_antiiso_without_eff = new TH1D("FA_nominal_inclusive_antiiso_without_eff", "", 120, 60, 120);

   TH1D *FA_nominal_inclusive_1D_pT = new TH1D("FA_nominal_inclusive_1D_pT", "", 120, 60, 120);
   TH1D *FA_nominal_inclusive_no_pT = new TH1D("FA_nominal_inclusive_no_pT", "", 120, 60, 120);
   TH1D *FA_AcoUp_inclusive = new TH1D("FA_AcoUp_inclusive", "", 120, 60, 120);
   TH1D *FA_AcoDown_inclusive = new TH1D("FA_AcoDown_inclusive", "", 120, 60, 120);
   TH1D *FA_tnpU_inclusive = new TH1D("FA_tnpU_inclusive", "", 120, 60, 120);
   TH1D *FA_tnpD_inclusive = new TH1D("FA_tnpD_inclusive", "", 120, 60, 120);
   TH1D *FA_mass_range_inclusive = new TH1D("FA_mass_range_inclusive", "", 80, 70, 110);

   TH1D *pT_spec_pp_FA = new TH1D("pT_spec_pp_FA", "", 200, 0, 200);
   double x_edges[] = {-2.4, -2.1, -1.8, -1.5, -1.2, -0.9, -0.6, -0.3,
                       0.0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4};
   double y_edges[] = {0.0, 1.0, 3.0, 5.0, 10.0, 20.0, 40.0, 70.0, 200.0};
   int nx = sizeof(x_edges) / sizeof(double) - 1;
   int ny = sizeof(y_edges) / sizeof(double) - 1;

   TH2D *pT_y_spec_pp_FA = new TH2D("pT_y_spec_pp_FA",
                                    "pT vs y (variable bins)",
                                    nx, x_edges,
                                    ny, y_edges);

   TH1D *FA_nominal_phi_plus[numberofphibins];
   TH1D *FA_nominal_phi_plus_1D_pT[numberofphibins];
   TH1D *FA_nominal_phi_plus_without_pT_reweight[numberofphibins];

   TH1D *FA_nominal_phi_plus_flattened[numberofphibins];
   TH1D *FA_nominal_phi_plus_1D_pT_flattened[numberofphibins];
   TH1D *FA_nominal_phi_plus_without_pT_reweight_flattened[numberofphibins];

   TH1D *FA_nominal_inclusive_flattened = new TH1D("FA_nominal_inclusive_flattened", "", 120, 60, 120);
   TH1D *FA_nominal_inclusive_1D_pT_flattened = new TH1D("FA_nominal_inclusive_1D_pT_flattened", "", 120, 60, 120);
   TH1D *FA_nominal_inclusive_no_pT_flattened = new TH1D("FA_nominal_inclusive_no_pT_flattened", "", 120, 60, 120);

   TFile *f_fit_function = new TFile("./cos_fit_save/data.root", "READ");
   TF1 *t_pp_PbPb_data = (TF1 *)f_fit_function->Get("pp_PbPb_data");
   TF1 *t_pp_PbPb_no_pT_data = (TF1 *)f_fit_function->Get("pp_PbPb_no_pT_data");
   TF1 *t_pp_PbPb_1D_pT_data = (TF1 *)f_fit_function->Get("pp_PbPb_1D_pT_data");

   for (int i = 0; i < 22; i++)
   {
      FA_nominal[i] = new TH1D(Form("FA_nominal_%i", i), "", 120, 60, 120);
      FA_AcoUp[i] = new TH1D(Form("FA_AcoUp_%i", i), "", 120, 60, 120);
      FA_AcoDown[i] = new TH1D(Form("FA_AcoDown_%i", i), "", 120, 60, 120);
      FA_tnpU[i] = new TH1D(Form("FA_tnpU_%i", i), "", 120, 60, 120);
      FA_tnpD[i] = new TH1D(Form("FA_tnpD_%i", i), "", 120, 60, 120);
      FA_mass_range[i] = new TH1D(Form("FA_mass_range_%i", i), "", 80, 70, 110);
      FA_1D_pT[i] = new TH1D(Form("FA_1D_pT_%i", i), "", 120, 60, 120);
   }

   for (int i = 0; i < numberofphibins; i++)
   {
      FA_nominal_phi_plus[i] = new TH1D(Form("pp_FA_nominal_phi_plus_%i", i), "", 120, 60, 120);
      FA_nominal_phi_plus_1D_pT[i] = new TH1D(Form("pp_FA_nominal_phi_plus_1D_pT_%i", i), "", 120, 60, 120);
      FA_nominal_phi_plus_without_pT_reweight[i] = new TH1D(Form("pp_FA_nominal_phi_plus_without_pT_reweight_%i", i), "", 120, 60, 120);
      FA_nominal_phi_plus_flattened[i] = new TH1D(Form("pp_FA_nominal_phi_plus_flattened_%i", i), "", 120, 60, 120);
      FA_nominal_phi_plus_1D_pT_flattened[i] = new TH1D(Form("pp_FA_nominal_phi_plus_1D_pT_flattened_%i", i), "", 120, 60, 120);
      FA_nominal_phi_plus_without_pT_reweight_flattened[i] = new TH1D(Form("pp_FA_nominal_phi_plus_without_pT_reweight_flattened_%i", i), "", 120, 60, 120);
   }

   TEfficiency *e, *e_up, *e_down, *e_acoup, *e_acodown;

   TFile *eff_f1 = new TFile("../ZBoson_18/rootfile/mc_eff.root", "READ");
   e = (TEfficiency *)eff_f1->Get("eff_0.0_100.0");
   e_up = (TEfficiency *)eff_f1->Get("eff_U_0.0_100.0");
   e_down = (TEfficiency *)eff_f1->Get("eff_D_0.0_100.0");
   e_acoup = (TEfficiency *)eff_f1->Get("eff_Acoup_0.0_100.0");
   e_acodown = (TEfficiency *)eff_f1->Get("eff_Acodown_0.0_100.0");

   TFile *pT_PbPb_weight = new TFile("../ZBoson_18/rootfile/pT_file.root", "READ");
   TF1 *pTweight_FA = (TF1 *)pT_PbPb_weight->Get("FA_ratio_fit");
   TH2D *pT_y_weight_FA = (TH2D *)pT_PbPb_weight->Get("FA_2D_ratio_data");

   for (Long64_t jentry = 0; jentry < nentries; jentry++)
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0)
         break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;

      double percentage = 100.0 * jentry / nentries;
      if (jentry % 100000 == 0)
         std::cout << "Progress: " << percentage << "% completed\r" << std::flush;

      if (abs(zVtx) > 15)
         continue;

      for (int znum = 0; znum < Reco_QQ_size; znum++)
      {
         Int_t muonindexplus = Reco_QQ_mupl_idx[znum];
         Int_t muonindexminus = Reco_QQ_mumi_idx[znum];

         Int_t nZs = Reco_QQ_4mom->GetEntriesFast();
         if (nZs != Reco_QQ_size)
            cout << "Error: Not Sure why Reco_QQ_Size != Size of TCloneArray" << endl;

         TLorentzVector *Z_momentum = (TLorentzVector *)Reco_QQ_4mom->At(znum);
         TLorentzVector *muonplus_momentum = (TLorentzVector *)Reco_mu_4mom->At(muonindexplus);
         TLorentzVector *muonminus_momentum = (TLorentzVector *)Reco_mu_4mom->At(muonindexminus);

         // Cut on Z
         if (Z_momentum->M() < 60 || Z_momentum->M() > 120)
            continue;
         if (Reco_QQ_VtxProb[znum] < 0.001)
            continue;
         if (abs(Z_momentum->Rapidity()) > 2.4)
            continue;

         // Cut on Muons
         if (muonplus_momentum->Pt() < 20 || muonminus_momentum->Pt() < 20)
            continue;
         if (abs(muonplus_momentum->Eta()) > 2.4 || abs(muonminus_momentum->Eta()) > 2.4)
            continue;

         // Cut on Charge
         Bool_t isOppositeSign = Reco_mu_charge[muonindexplus] != Reco_mu_charge[muonindexminus];
         if (!isOppositeSign)
         {
            if (Reco_QQ_sign[znum] == 0)
               cout << "You should never see this, charge matching failed!!!!" << endl;
            continue;
         }

         if (Reco_mu_isTightCutBased[muonindexplus] != 1 || Reco_mu_isTightCutBased[muonindexminus] != 1)
            cout << "we have a case that does not pass tightID selection" << endl;

         // ------------------------------------------------------------

         Bool_t isEtacutPassed = (abs(muonplus_momentum->Eta()) < 1) && (abs(muonminus_momentum->Eta()) < 1);
         Double_t ZMass = Z_momentum->M();

         double efficiency = getEfficiency(e, Z_momentum->Rapidity(), Z_momentum->Pt());
         double efficiency_U = getEfficiency(e_up, Z_momentum->Rapidity(), Z_momentum->Pt());
         double efficiency_D = getEfficiency(e_down, Z_momentum->Rapidity(), Z_momentum->Pt());
         double efficiency_acoup = getEfficiency(e_acoup, Z_momentum->Rapidity(), Z_momentum->Pt());
         double efficiency_acodown = getEfficiency(e_acodown, Z_momentum->Rapidity(), Z_momentum->Pt());

         float acoplanarity = 1 - TMath::Abs(TMath::ACos(TMath::Cos(muonplus_momentum->Phi() - muonminus_momentum->Phi()))) / TMath::Pi();
         bool passesAco[3] = {1, 1, 1};
         if (Z_momentum->Pt() < 1.25 && acoplanarity < 0.001)
            passesAco[0] = false;
         if (Z_momentum->Pt() < 1.25 && acoplanarity < 0.0015)
            passesAco[1] = false;
         if (Z_momentum->Pt() < 1.25 && acoplanarity < 0.0005)
            passesAco[2] = false;

         double FA_pTweight_1D = getWeightFromHist(pTweight_FA, Z_momentum->Pt());
         double FA_pTweight_2D = pT_y_weight_FA->GetBinContent(pT_y_weight_FA->FindBin(Z_momentum->Rapidity(), Z_momentum->Pt()));

         bool passisolation = true;

         // ----- Optional PF relative-isolation cut on both muons -----
         if (useIsolation)
         {
            double relIsoPlus = RelIsoPF(muonindexplus,
                                         muonplus_momentum->Pt(),
                                         Reco_mu_PFChIso,
                                         Reco_mu_PFNeuIso,
                                         Reco_mu_PFPhoIso,
                                         Reco_mu_PFPUIso);
            double relIsoMinus = RelIsoPF(muonindexminus,
                                          muonminus_momentum->Pt(),
                                          Reco_mu_PFChIso,
                                          Reco_mu_PFNeuIso,
                                          Reco_mu_PFPhoIso,
                                          Reco_mu_PFPUIso);
            if (relIsoPlus >= kRelIsoCut || relIsoMinus >= kRelIsoCut)
               passisolation = false;
         }

         if (passesAco[0])
         {
            Int_t phibin = GetPhiBin(muonplus_momentum->Phi(), numberofphibins);

            pT_spec_pp_FA->Fill(Z_momentum->Pt(), 1.0 / efficiency);
            pT_y_spec_pp_FA->Fill(Z_momentum->Rapidity(), Z_momentum->Pt(), 1.0 / efficiency);

            if (passisolation)
               FA_nominal_inclusive->Fill(ZMass, 1.0 * FA_pTweight_2D / efficiency);
            else
            {
               FA_nominal_inclusive_antiiso->Fill(ZMass, 1.0 * FA_pTweight_2D / efficiency);
            }

            if (passisolation)
               FA_nominal_inclusive_without_eff->Fill(ZMass, 1.0 * FA_pTweight_2D);
            else
            {
               FA_nominal_inclusive_antiiso_without_eff->Fill(ZMass, 1.0 * FA_pTweight_2D);
            }
         } // end QQ loop
      }
   } // end tree loop

   // ----- output filename depends on whether iso cut was applied -----
   const char *outputFileName = useIsolation
                                    ? Form("./new_pp_data_file_with_iso_%.2f.root",kRelIsoCut)
                                    : "./new_pp_data_file_no_iso.root";
   TFile *writeout = new TFile(outputFileName, "UPDATE");
   writeout->cd();

   FA_nominal_inclusive->Write("", 2);
   FA_nominal_inclusive_without_eff->Write("", 2);
   FA_nominal_inclusive_antiiso_without_eff->Write("", 2);
   FA_nominal_inclusive_antiiso->Write("", 2);

   writeout->Close();

   // pT spectra still go to the shared pT_file.root, as in skimnew.C.
   // If you want them gated on the iso flag too, change the path here.
   /*TFile *pt_File = new TFile("../ZBoson_18/rootfile/pT_file.root", "UPDATE");
   pt_File->cd();
   pT_spec_pp_FA->Write("", 2);
   pT_y_spec_pp_FA->Write("", 2);
   pt_File->Close();*/
}