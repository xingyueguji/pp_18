#define SkimNew_cxx
#include "SkimNew.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

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

   Double_t runupperlimit[22] = {320916, 321012, 321140, 321232, 321393, 321436, 321735, 321820, 321909, 322014, 322118, 322319, 322381, 322617, 323525, 323778, 324021, 324245, 324772, 324897, 325001, 325175};
   Double_t runlowerlimit[22] = {320500, 320917, 321051, 321149, 321233, 321396, 321457, 321755, 321831, 321917, 322022, 322179, 322322, 322407, 322625, 323526, 323790, 324022, 324293, 324785, 324970, 325022};

   // Preparing the needed roodataset
   RooRealVar *x = new RooRealVar("x", "x", 60, 120);

   //"Raw"
   TH1D *h_mass_array_raw_inclusive = new TH1D("h_mass_array_raw_inclusive", "h_mass_array_raw_inclusive", 120, 60, 120);
   TH1D *h_mass_array_raw[22];
   RooDataSet *mass_array_raw[22];
   RooDataSet *mass_array_raw_inclusive = new RooDataSet("mass_array_raw_inclusive", "mass_array_raw_inclusive", RooArgSet(*x));

   //"Etacut"
   TH1D *h_mass_array_eta_inclusive = new TH1D("h_mass_array_eta_inclusive", "h_mass_array_eta_inclusive", 120, 60, 120);
   TH1D *h_mass_array_eta[22];
   RooDataSet *mass_array_eta[22];
   RooDataSet *mass_array_eta_inclusive = new RooDataSet("mass_array_eta_inclusive", "mass_array_eta_inclusive", RooArgSet(*x));

   for (int i = 0; i < 22; i++)
   {

      h_mass_array_raw[i] = new TH1D(Form("h_mass_array_raw_%i", i), Form("h_mass_array_raw_%i", i), 120, 60, 120);
      h_mass_array_eta[i] = new TH1D(Form("h_mass_array_eta_%i", i), Form("h_mass_array_eta_%i", i), 120, 60, 120);

      mass_array_raw[i] = new RooDataSet(Form("mass_array_raw_%i", i), Form("mass_array_raw_%i", i), RooArgSet(*x));
      mass_array_eta[i] = new RooDataSet(Form("mass_array_eta_%i", i), Form("mass_array_eta_%i", i), RooArgSet(*x));
   }

   for (Long64_t jentry = 0; jentry < nentries; jentry++)
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0)
         break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      double progress = (static_cast<double>(jentry) / nentries) * 100;

      if (jentry % 100000 == 0)
         cout << "We are at event " << jentry << " Percentage done: " << progress << "%" << endl;

      // Event selection:

      if (abs(zVtx) > 15)
         continue;
      // if (!(CheckTrigBit(HLTriggers,6))) continue;

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

         // Cut on Trigger of two candidate muons
         /*Bool_t isDoughter1Trigger = CheckTrigBit(Reco_mu_trig(muonindexplus),6);
         Bool_t isDoughter2Trigger = CheckTrigBit(Reco_mu_trig(muonindexminus),6);
         if (!(isDoughter1Trigger||isDoughter2Trigger)) continue;*/

         // Cut on Charge
         Bool_t isOppositeSign = Reco_mu_charge[muonindexplus] != Reco_mu_charge[muonindexminus];
         if (!isOppositeSign)
         {
            if (Reco_QQ_sign[znum] == 0)
               cout << "You should never see this, charge matching failed!!!!" << endl;
            continue;
         }

         // Aco cut
         Float_t acoplanarity = 1 - TMath::Abs(TMath::ACos(TMath::Cos(muonplus_momentum->Phi() - muonminus_momentum->Phi()))) / TMath::Pi();
         if (acoplanarity < 0.001)
            continue;

         // Eta < 1 cut I will apply this later when filling
         Bool_t isEtacutPassed = abs(muonplus_momentum->Eta()) < 1 && abs(muonminus_momentum->Eta()) < 1;

         // Get Mass for easier fill

         Double_t ZMass = Z_momentum->M();
         x->setVal(ZMass);

         // Finally can fill the Roodataset

         // Fill the incluive one first
         h_mass_array_raw_inclusive->Fill(ZMass);
         mass_array_raw_inclusive->add(RooArgSet(*x));
         if (isEtacutPassed)
            h_mass_array_eta_inclusive->Fill(ZMass);
         if (isEtacutPassed)
            mass_array_eta_inclusive->add(RooArgSet(*x));

         // Fill the run number based then

         for (int runindex = 0; runindex < 22; ++runindex)
         {
            if (runNb >= runlowerlimit[runindex] && runNb <= runupperlimit[runindex])
            {
               h_mass_array_raw[runindex]->Fill(ZMass);
               mass_array_raw[runindex]->add(RooArgSet(*x));
               if (isEtacutPassed)
               {
                  h_mass_array_eta[runindex]->Fill(ZMass);
                  mass_array_eta[runindex]->add(RooArgSet(*x));
               }
            }
         }
      } // End of QQ loop
   } // End of tree loop

   TFile *writeout = new TFile("./new_pp_data_file_stability_readonly.root", "UPDATE");
   writeout->cd();
   h_mass_array_raw_inclusive->Write("", 2);
   mass_array_raw_inclusive->Write("", 2);
   h_mass_array_eta_inclusive->Write("", 2);
   mass_array_eta_inclusive->Write("", 2);

   for (int j = 0; j < 22; j++)
   {
      h_mass_array_raw[j]->Write("", 2);
      mass_array_raw[j]->Write("", 2);
      h_mass_array_eta[j]->Write("", 2);
      mass_array_eta[j]->Write("", 2);
   }

   writeout->Close();
}
