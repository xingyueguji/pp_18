#define ppmc_cxx
#include "ppmc.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

double getWeightFromHist(TF1 *weightfucntion, double pt)
{
   double ratio = weightfucntion->Eval(pt); // Find the corresponding bin for pT

   return ratio;
}

void ppmc::Loop()
{
   //   In a ROOT session, you can do:
   //      root> .L ppmc.C
   //      root> ppmc t
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

   TH1D *FA_nominal_inclusive = new TH1D("FA_nominal_inclusive", "", 120, 60, 120);

   TFile *pT_PbPb_weight = new TFile("../ZBoson_18/rootfile/pT_file.root", "READ");

   TF1 *pTweight_FA = (TF1 *)pT_PbPb_weight->Get("FA_ratio_fit");

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
         double FA_pTweight = getWeightFromHist(pTweight_FA, Z_momentum->Pt());

         float acoplanarity = 1 - TMath::Abs(TMath::ACos(TMath::Cos(muonplus_momentum->Phi() - muonminus_momentum->Phi()))) / TMath::Pi();
         bool passesAco[3] = {1, 1, 1};
         if (Z_momentum->Pt() < 1.25 && acoplanarity < 0.001)
            passesAco[0] = false;
         if (Z_momentum->Pt() < 1.25 && acoplanarity < 0.0015)
            passesAco[1] = false;
         if (Z_momentum->Pt() < 1.25 && acoplanarity < 0.0005)
            passesAco[2] = false;

         if (passesAco[0])
         {
            FA_nominal_inclusive->Fill(ZMass, 1.0 * FA_pTweight);
         }
      }
   }

   TFile *writeout = new TFile("./new_pp_data_file_stability_readonly.root", "UPDATE");
   writeout->cd();
   FA_nominal_inclusive->Write("pp_mc_inclusive_test_with_pt_reweight", 2);
   writeout->Close();
}
