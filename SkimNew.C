#define SkimNew_cxx
#include "SkimNew.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

double getWeightFromHist(TF1 *weightfucntion, double pt)
{
   double ratio = weightfucntion->Eval(pt); // Find the corresponding bin for pT

   return ratio;
}

Int_t GetPhiBin(double phi, const int numberofphibins)
{
   // wrap phi into [-pi, pi]
   while (phi <= -TMath::Pi())
      phi += 2 * TMath::Pi();
   while (phi > TMath::Pi())
      phi -= 2 * TMath::Pi();

   double xmin = -TMath::Pi();
   double xmax = TMath::Pi();
   double width = (xmax - xmin) / numberofphibins;

   int bin = static_cast<int>((phi - xmin) / width);

   // safety clamp
   if (bin < 0)
      bin = 0;
   if (bin >= numberofphibins)
   {
      bin = numberofphibins - 1;
      cout << "Max range exceeded? you should not able to see this" << endl;
   }

   return bin;
}

void SkimNew::Loop()
{
   //   In a ROOT session, you can do:
   //      root> .L SkimNew.C
   //      root> SkimNew t
   //      root> t.GetEntry(12); // Fill t data members with entry number 12
   //      root> t.Show();       // Show values of entry 12
   //      root> t.Show(16);     // Read and show values of entry 16
   //      root> t.Loop();       // Loop on all entries
   //

   //     This is the loop skeleton where:
   //    jentry is the global entry number in the chain
   //    ientry is the entry number in the current Tree
   //  Note that the argument to GetEntry must be:
   //    jentry for TChain::GetEntry
   //    ientry for TTree::GetEntry and TBranch::GetEntry
   //
   //       To read only selected branches, Insert statements like:
   // METHOD1:
   //    fChain->SetBranchStatus("*",0);  // disable all branches
   //    fChain->SetBranchStatus("branchname",1);  // activate branchname
   // METHOD2: replace line
   //    fChain->GetEntry(jentry);       //read all branches
   // by  b_branchname->GetEntry(ientry); //read only this branch
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
   TH1D *FA_nominal_inclusive = new TH1D("FA_nominal_inclusive", "", 120, 60, 120);
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

   // n_bins = (#edges - 1)
   int nx = sizeof(x_edges) / sizeof(double) - 1;
   int ny = sizeof(y_edges) / sizeof(double) - 1;

   TH2D *pT_y_spec_pp_FA = new TH2D("pT_y_spec_pp_FA",
                                    "pT vs y (variable bins)",
                                    nx, x_edges,
                                    ny, y_edges);

   TH1D *FA_nominal_phi_plus[numberofphibins];
   TH1D *FA_nominal_phi_plus_without_pT_reweight[numberofphibins];
   TH1D *FA_nominal_phi_plus_flattened[numberofphibins];
   TH1D *FA_nominal_phi_plus_without_pT_reweight_flattened[numberofphibins];
   TH1D *FA_nominal_inclusive_flattened = new TH1D("FA_nominal_inclusive_flattened", "", 120, 60, 120);
   TH1D *FA_nominal_inclusive_no_pT_flattened = new TH1D("FA_nominal_inclusive_no_pT_flattened", "", 120, 60, 120);

   TFile *f_fit_function = new TFile("../pp_18/cos_fit_save/data.root", "READ");
   TF1 *t_pp_PbPb_data = (TF1 *)f_fit_function->Get("pp_PbPb_data");
   TF1 *t_pp_PbPb_no_pT_data = (TF1 *)f_fit_function->Get("pp_PbPb_no_pT_data");

   for (int i = 0; i < 22; i++)
   {
      FA_nominal[i] = new TH1D(Form("FA_nominal_%i", i), "", 120, 60, 120);
      FA_AcoUp[i] = new TH1D(Form("FA_AcoUp_%i", i), "", 120, 60, 120);
      FA_AcoDown[i] = new TH1D(Form("FA_AcoDown_%i", i), "", 120, 60, 120);
      FA_tnpU[i] = new TH1D(Form("FA_tnpU_%i", i), "", 120, 60, 120);
      FA_tnpD[i] = new TH1D(Form("FA_tnpD_%i", i), "", 120, 60, 120);
      FA_mass_range[i] = new TH1D(Form("FA_mass_range_%i", i), "", 80, 70, 110);
   }

   for (int i = 0; i < numberofphibins; i++)
   {
      FA_nominal_phi_plus[i] = new TH1D(Form("pp_FA_nominal_phi_plus_%i", i), "", 120, 60, 120);
      FA_nominal_phi_plus_without_pT_reweight[i] = new TH1D(Form("pp_FA_nominal_phi_plus_without_pT_reweight_%i", i), "", 120, 60, 120);
      FA_nominal_phi_plus_flattened[i] = new TH1D(Form("pp_FA_nominal_phi_plus_flattened_%i", i), "", 120, 60, 120);
      FA_nominal_phi_plus_without_pT_reweight_flattened[i] = new TH1D(Form("pp_FA_nominal_phi_plus_without_pT_reweight_flattened_%i", i), "", 120, 60, 120);
   }

   TEfficiency *e;
   TEfficiency *e_up;
   TEfficiency *e_down;
   TEfficiency *e_acoup;
   TEfficiency *e_acodown;

   TFile *eff_f1 = new TFile("../ZBoson_18/rootfile/mc_eff.root", "READ");

   e = (TEfficiency *)eff_f1->Get("eff_0_100");
   e_up = (TEfficiency *)eff_f1->Get("eff_U_0_100");
   e_down = (TEfficiency *)eff_f1->Get("eff_D_0_100");
   e_acoup = (TEfficiency *)eff_f1->Get("eff_Acoup_0_100");
   e_acodown = (TEfficiency *)eff_f1->Get("eff_Acodown_0_100");

   TFile *pT_PbPb_weight = new TFile("../ZBoson_18/rootfile/pT_file.root", "READ");

   TF1 *pTweight_FA = (TF1 *)pT_PbPb_weight->Get("FA_ratio_fit");
   TF1 *pTweight_Eta = (TF1 *)pT_PbPb_weight->Get("Eta_ratio_fit");

   TH2D *pT_y_weight_FA = (TH2D*)pT_PbPb_weight->Get("FA_2D_ratio_data");

   for (Long64_t jentry = 0; jentry < nentries; jentry++)
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0)
         break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      double percentage = 100.0 * jentry / nentries;
      if (jentry % 100000 == 0)
         std::cout << "Progress: " << percentage << "% completed\r" << std::flush;

      // Event selection:

      if (abs(zVtx) > 15)
         continue;
      // if (!(CheckTrigBit(HLTriggers, 5)))
      //    continue;
      //  if ((CheckTrigBit(HLTriggers,8))) cout << "For this event we have HLT 6 == 1" << endl;
      //  cout << " HLT for event is " << std::bitset<18>(HLTriggers) << endl;

      // Now looping through all reco dimuon pairs
      for (int znum = 0; znum < Reco_QQ_size; znum++)
      {

         // Get Two candidate muons index first
         Int_t muonindexplus = Reco_QQ_mupl_idx[znum];
         Int_t muonindexminus = Reco_QQ_mumi_idx[znum];

         // Get TLorentzvectors
         Int_t nZs = Reco_QQ_4mom->GetEntriesFast();
         if (nZs != Reco_QQ_size)
            cout << "Error: " << " Not Sure why Reco_QQ_Size != Size of TCloneArray" << endl;
         TLorentzVector *Z_momentum = (TLorentzVector *)Reco_QQ_4mom->At(znum);
         TLorentzVector *muonplus_momentum = (TLorentzVector *)Reco_mu_4mom->At(muonindexplus);
         TLorentzVector *muonminus_momentum = (TLorentzVector *)Reco_mu_4mom->At(muonindexminus);

         // Apply cuts
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

         // if (Reco_mu_isTightCutBased[muonindexplus])
         //  Cut on Trigger of two candidate muons
         //  Bool_t isDaughter1Trigger = CheckTrigBit(Reco_mu_trig[muonindexplus], 5);
         //  Bool_t isDaughter2Trigger = CheckTrigBit(Reco_mu_trig[muonindexminus], 5);
         //  if (!(isDaughter1Trigger || isDaughter2Trigger))
         //    continue;
         /*if ((isDaughter1Trigger||isDaughter2Trigger)){
            cout << " We have one dimuon pair with HLT 6 or == 1 " << endl;
            cout << " Trigger Bit is * and * " << Reco_mu_trig[muonindexplus] << " " << Reco_mu_trig[muonindexminus] << endl;
         }*/

         // Cut on Charge
         Bool_t isOppositeSign = Reco_mu_charge[muonindexplus] != Reco_mu_charge[muonindexminus];
         if (!isOppositeSign)
         {
            if (Reco_QQ_sign[znum] == 0)
               cout << "You should never see this, charge matching failed!!!!" << endl;
            continue;
         }

         if (Reco_mu_isTightCutBased[muonindexplus] != 1 || Reco_mu_isTightCutBased[muonindexminus] != 1)
         {
            cout << "we have a case that does not pass tightID selection" << endl;
         }

         // Aco cut
         // Float_t acoplanarity = 1 - TMath::Abs(TMath::ACos(TMath::Cos(muonplus_momentum->Phi() - muonminus_momentum->Phi()))) / TMath::Pi();
         // if (acoplanarity < 0.001)
         // continue;

         // Eta < 1 cut I will apply this later when filling
         Bool_t isEtacutPassed = (abs(muonplus_momentum->Eta()) < 1) && (abs(muonminus_momentum->Eta()) < 1);

         // Get Mass for easier fill

         Double_t ZMass = Z_momentum->M();

         // Fill the incluive one first
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
         double Eta_pTweight_1D = getWeightFromHist(pTweight_Eta, Z_momentum->Pt());

         double FA_pTweight_2D = pT_y_weight_FA->GetBinContent(pT_y_weight_FA->FindBin(Z_momentum->Rapidity(), Z_momentum->Pt()));

         // Here for inclusive
         if (passesAco[0])
         {
            Int_t phibin = GetPhiBin(muonplus_momentum->Phi(), numberofphibins);

            pT_spec_pp_FA->Fill(Z_momentum->Pt(), 1.0 / efficiency);
            pT_y_spec_pp_FA->Fill(Z_momentum->Rapidity(), Z_momentum->Pt(), 1.0 / efficiency);
            FA_nominal_inclusive->Fill(ZMass, 1.0 * FA_pTweight_2D / efficiency);
            FA_tnpU_inclusive->Fill(ZMass, 1.0 * FA_pTweight_2D / efficiency_U);
            FA_tnpD_inclusive->Fill(ZMass, 1.0 * FA_pTweight_2D / efficiency_D);
            FA_mass_range_inclusive->Fill(ZMass, 1.0 * FA_pTweight_2D / efficiency);
            FA_nominal_inclusive_no_pT->Fill(ZMass, 1.0 / efficiency);

            FA_nominal_phi_plus[phibin]->Fill(ZMass, 1.0 * FA_pTweight_2D / efficiency);
            FA_nominal_phi_plus_without_pT_reweight[phibin]->Fill(ZMass, 1.0 / efficiency);

            double phi = TVector2::Phi_mpi_pi(muonplus_momentum->Phi());

            // read parameters explicitly
            double a = t_pp_PbPb_data->GetParameter(0);
            double b = t_pp_PbPb_data->GetParameter(1);
            double phi0 = t_pp_PbPb_data->GetParameter(2);

            // evaluate the fitted modulation
            double fphi = a + b * cos(phi - phi0);

            // subtract only the oscillatory component (keep the mean)
            double mass_corr = ZMass - (fphi - a);

            FA_nominal_phi_plus_flattened[phibin]->Fill(mass_corr, FA_pTweight_2D / efficiency);
            FA_nominal_inclusive_flattened->Fill(mass_corr, FA_pTweight_2D / efficiency);

            double a_no = t_pp_PbPb_no_pT_data->GetParameter(0);
            double b_no = t_pp_PbPb_no_pT_data->GetParameter(1);
            double phi0_no = t_pp_PbPb_no_pT_data->GetParameter(2);

            double fphi_no = a_no + b_no * cos(phi - phi0_no);
            double mass_corr_no_pT = ZMass - (fphi_no - a_no);
            FA_nominal_phi_plus_without_pT_reweight_flattened[phibin]->Fill(mass_corr_no_pT, 1.0 / efficiency);
            FA_nominal_inclusive_no_pT_flattened->Fill(mass_corr_no_pT, 1.0 / efficiency);
         }
         if (passesAco[1])
         {
            FA_AcoUp_inclusive->Fill(ZMass, 1.0 * FA_pTweight_2D / efficiency_acoup);
         }
         if (passesAco[2])
         {
            FA_AcoDown_inclusive->Fill(ZMass, 1.0 * FA_pTweight_2D / efficiency_acodown);
         }

         // Fill the run number based then

         for (int runindex = 0; runindex < 22; ++runindex)
         {
            if (runNb >= runlowerlimit[runindex] && runNb <= runupperlimit[runindex])
            {
               if (passesAco[0])
               {
                  FA_nominal[runindex]->Fill(ZMass, 1.0 * FA_pTweight_2D / efficiency);
                  FA_tnpU[runindex]->Fill(ZMass, 1.0 * FA_pTweight_2D / efficiency_U);
                  FA_tnpD[runindex]->Fill(ZMass, 1.0 * FA_pTweight_2D / efficiency_D);
                  FA_mass_range[runindex]->Fill(ZMass, 1.0 * FA_pTweight_2D / efficiency);
               }
               if (passesAco[1])
               {
                  FA_AcoUp[runindex]->Fill(ZMass, 1.0 * FA_pTweight_2D / efficiency_acoup);
               }
               if (passesAco[2])
               {
                  FA_AcoDown[runindex]->Fill(ZMass, 1.0 * FA_pTweight_2D / efficiency_acodown);
               }
            }
         }
      } // End of QQ loop
   } // End of tree loop

   TFile *writeout = new TFile("./new_pp_data_file_stability_readonly.root", "UPDATE");
   writeout->cd();

   FA_nominal_inclusive->Write("", 2);
   FA_AcoUp_inclusive->Write("", 2);
   FA_AcoDown_inclusive->Write("", 2);
   FA_tnpU_inclusive->Write("", 2);
   FA_tnpD_inclusive->Write("", 2);
   FA_mass_range_inclusive->Write("", 2);
   FA_nominal_inclusive_no_pT->Write("", 2);

   for (int Z = 0; Z < numberofphibins; Z++)
   {
      FA_nominal_phi_plus[Z]->Write("", 2);
      FA_nominal_phi_plus_without_pT_reweight[Z]->Write("", 2);
      FA_nominal_phi_plus_flattened[Z]->Write("", 2);
      FA_nominal_phi_plus_without_pT_reweight_flattened[Z]->Write("", 2);
   }

   FA_nominal_inclusive_flattened->Write("", 2);
   FA_nominal_inclusive_no_pT_flattened->Write("", 2);

   for (int j = 0; j < 22; j++)
   {
      FA_nominal[j]->Write("", 2);
      FA_AcoUp[j]->Write("", 2);
      FA_AcoDown[j]->Write("", 2);
      FA_tnpU[j]->Write("", 2);
      FA_tnpD[j]->Write("", 2);
      FA_mass_range[j]->Write("", 2);
   }

   writeout->Close();

   TFile *pt_File = new TFile("../ZBoson_18/rootfile/pT_file.root", "UPDATE");
   pt_File->cd();

   pT_spec_pp_FA->Write("", 2);
   pT_y_spec_pp_FA->Write("", 2);

   pt_File->Close();
}
