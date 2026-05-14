//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun May  3 21:00:04 2026 by ROOT version 6.32.06
// from TTree myTree/My TTree of dimuons
// found on file: combined_pp_with_iso_v2.root
//////////////////////////////////////////////////////////

#ifndef pp_with_iso_h
#define pp_with_iso_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TClonesArray.h"

class pp_with_iso
{
public:
   TTree *fChain;  //! pointer to the analyzed TTree or TChain
   Int_t fCurrent; //! current Tree number in a TChain

   // Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t eventNb;
   UInt_t runNb;
   UInt_t LS;
   Float_t zVtx;
   Short_t nPV;
   Short_t Ntracks;
   Int_t trigPrescale[18];
   ULong64_t HLTriggers;
   Short_t Reco_QQ_size;
   Short_t Reco_QQ_type[5]; //[Reco_QQ_size]
   Short_t Reco_QQ_sign[5]; //[Reco_QQ_size]
   TClonesArray *Reco_QQ_4mom;
   Short_t Reco_QQ_mupl_idx[5];   //[Reco_QQ_size]
   Short_t Reco_QQ_mumi_idx[5];   //[Reco_QQ_size]
   ULong64_t Reco_QQ_trig[5];     //[Reco_QQ_size]
   Bool_t Reco_QQ_isCowboy[5];    //[Reco_QQ_size]
   Float_t Reco_QQ_ctau[5];       //[Reco_QQ_size]
   Float_t Reco_QQ_ctauErr[5];    //[Reco_QQ_size]
   Float_t Reco_QQ_cosAlpha[5];   //[Reco_QQ_size]
   Float_t Reco_QQ_ctau3D[5];     //[Reco_QQ_size]
   Float_t Reco_QQ_ctauErr3D[5];  //[Reco_QQ_size]
   Float_t Reco_QQ_cosAlpha3D[5]; //[Reco_QQ_size]
   Float_t Reco_QQ_VtxProb[5];    //[Reco_QQ_size]
   Float_t Reco_QQ_dca[5];        //[Reco_QQ_size]
   Float_t Reco_QQ_MassErr[5];    //[Reco_QQ_size]
   TClonesArray *Reco_QQ_vtx;
   Short_t Reco_mu_size;
   Short_t Reco_mu_type[5];        //[Reco_mu_size]
   Int_t Reco_mu_SelectionType[5]; //[Reco_mu_size]
   Short_t Reco_mu_charge[5];      //[Reco_mu_size]
   TClonesArray *Reco_mu_4mom;
   ULong64_t Reco_mu_trig[5];          //[Reco_mu_size]
   Bool_t Reco_mu_InTightAcc[5];       //[Reco_mu_size]
   Bool_t Reco_mu_InLooseAcc[5];       //[Reco_mu_size]
   Bool_t Reco_mu_highPurity[5];       //[Reco_mu_size]
   Bool_t Reco_mu_isPF[5];             //[Reco_mu_size]
   Bool_t Reco_mu_isTracker[5];        //[Reco_mu_size]
   Bool_t Reco_mu_isGlobal[5];         //[Reco_mu_size]
   Bool_t Reco_mu_isSoftCutBased[5];   //[Reco_mu_size]
   Bool_t Reco_mu_isHybridSoft[5];     //[Reco_mu_size]
   Bool_t Reco_mu_isMedium[5];         //[Reco_mu_size]
   Bool_t Reco_mu_isTightCutBased[5];  //[Reco_mu_size]
   Short_t Reco_mu_candType[5];        //[Reco_mu_size]
   Int_t Reco_mu_nPixValHits[5];       //[Reco_mu_size]
   Int_t Reco_mu_nMuValHits[5];        //[Reco_mu_size]
   Int_t Reco_mu_nTrkHits[5];          //[Reco_mu_size]
   Float_t Reco_mu_normChi2_inner[5];  //[Reco_mu_size]
   Float_t Reco_mu_normChi2_global[5]; //[Reco_mu_size]
   Int_t Reco_mu_nPixWMea[5];          //[Reco_mu_size]
   Int_t Reco_mu_nTrkWMea[5];          //[Reco_mu_size]
   Int_t Reco_mu_StationsMatched[5];   //[Reco_mu_size]
   Float_t Reco_mu_dxy[5];             //[Reco_mu_size]
   Float_t Reco_mu_dxyErr[5];          //[Reco_mu_size]
   Float_t Reco_mu_dz[5];              //[Reco_mu_size]
   Float_t Reco_mu_dzErr[5];           //[Reco_mu_size]
   Float_t Reco_mu_ptErr_inner[5];     //[Reco_mu_size]
   Float_t Reco_mu_PFChIso[5];         //[Reco_mu_size]
   Float_t Reco_mu_PFNeuIso[5];        //[Reco_mu_size]
   Float_t Reco_mu_PFPhoIso[5];        //[Reco_mu_size]
   Float_t Reco_mu_PFPUIso[5];         //[Reco_mu_size]

   // List of branches
   TBranch *b_eventNb;                 //!
   TBranch *b_runNb;                   //!
   TBranch *b_LS;                      //!
   TBranch *b_zVtx;                    //!
   TBranch *b_nPV;                     //!
   TBranch *b_Ntracks;                 //!
   TBranch *b_trigPrescale;            //!
   TBranch *b_HLTriggers;              //!
   TBranch *b_Reco_QQ_size;            //!
   TBranch *b_Reco_QQ_type;            //!
   TBranch *b_Reco_QQ_sign;            //!
   TBranch *b_Reco_QQ_4mom;            //!
   TBranch *b_Reco_QQ_mupl_idx;        //!
   TBranch *b_Reco_QQ_mumi_idx;        //!
   TBranch *b_Reco_QQ_trig;            //!
   TBranch *b_Reco_QQ_isCowboy;        //!
   TBranch *b_Reco_QQ_ctau;            //!
   TBranch *b_Reco_QQ_ctauErr;         //!
   TBranch *b_Reco_QQ_cosAlpha;        //!
   TBranch *b_Reco_QQ_ctau3D;          //!
   TBranch *b_Reco_QQ_ctauErr3D;       //!
   TBranch *b_Reco_QQ_cosAlpha3D;      //!
   TBranch *b_Reco_QQ_VtxProb;         //!
   TBranch *b_Reco_QQ_dca;             //!
   TBranch *b_Reco_QQ_MassErr;         //!
   TBranch *b_Reco_QQ_vtx;             //!
   TBranch *b_Reco_mu_size;            //!
   TBranch *b_Reco_mu_type;            //!
   TBranch *b_Reco_mu_SelectionType;   //!
   TBranch *b_Reco_mu_charge;          //!
   TBranch *b_Reco_mu_4mom;            //!
   TBranch *b_Reco_mu_trig;            //!
   TBranch *b_Reco_mu_InTightAcc;      //!
   TBranch *b_Reco_mu_InLooseAcc;      //!
   TBranch *b_Reco_mu_highPurity;      //!
   TBranch *b_Reco_mu_isPF;            //!
   TBranch *b_Reco_mu_isTracker;       //!
   TBranch *b_Reco_mu_isGlobal;        //!
   TBranch *b_Reco_mu_isSoftCutBased;  //!
   TBranch *b_Reco_mu_isHybridSoft;    //!
   TBranch *b_Reco_mu_isMedium;        //!
   TBranch *b_Reco_mu_isTightCutBased; //!
   TBranch *b_Reco_mu_candType;        //!
   TBranch *b_Reco_mu_nPixValHits;     //!
   TBranch *b_Reco_mu_nMuValHits;      //!
   TBranch *b_Reco_mu_nTrkHits;        //!
   TBranch *b_Reco_mu_normChi2_inner;  //!
   TBranch *b_Reco_mu_normChi2_global; //!
   TBranch *b_Reco_mu_nPixWMea;        //!
   TBranch *b_Reco_mu_nTrkWMea;        //!
   TBranch *b_Reco_mu_StationsMatched; //!
   TBranch *b_Reco_mu_dxy;             //!
   TBranch *b_Reco_mu_dxyErr;          //!
   TBranch *b_Reco_mu_dz;              //!
   TBranch *b_Reco_mu_dzErr;           //!
   TBranch *b_Reco_mu_ptErr_inner;     //!
   TBranch *b_Reco_mu_PFChIso;         //!
   TBranch *b_Reco_mu_PFNeuIso;        //!
   TBranch *b_Reco_mu_PFPhoIso;        //!
   TBranch *b_Reco_mu_PFPUIso;         //!

   pp_with_iso(TTree *tree = 0);
   virtual ~pp_with_iso();
   virtual Int_t Cut(Long64_t entry);
   virtual Int_t GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void Init(TTree *tree);
   virtual void Loop(bool useIsolation = false);
   virtual bool Notify();
   virtual void Show(Long64_t entry = -1);
};

#endif

#ifdef pp_with_iso_cxx
pp_with_iso::pp_with_iso(TTree *tree) : fChain(0)
{
   // if parameter tree is not specified (or zero), connect the file
   // used to generate this class and read the Tree.
   TFile *f = TFile::Open("~/ROOTFILES/combined_pp_with_iso_v2.root");
   TDirectory *dir = (TDirectory *)f->Get("hionia");
   dir->GetObject("myTree", tree);
   Init(tree);
}

pp_with_iso::~pp_with_iso()
{
   if (!fChain)
      return;
   delete fChain->GetCurrentFile();
}

Int_t pp_with_iso::GetEntry(Long64_t entry)
{
   // Read contents of entry.
   if (!fChain)
      return 0;
   return fChain->GetEntry(entry);
}
Long64_t pp_with_iso::LoadTree(Long64_t entry)
{
   // Set the environment to read one entry
   if (!fChain)
      return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0)
      return centry;
   if (fChain->GetTreeNumber() != fCurrent)
   {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void pp_with_iso::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Reco_QQ_4mom = 0;
   Reco_QQ_vtx = 0;
   Reco_mu_4mom = 0;
   // Set branch addresses and branch pointers
   if (!tree)
      return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("eventNb", &eventNb, &b_eventNb);
   fChain->SetBranchAddress("runNb", &runNb, &b_runNb);
   fChain->SetBranchAddress("LS", &LS, &b_LS);
   fChain->SetBranchAddress("zVtx", &zVtx, &b_zVtx);
   fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
   fChain->SetBranchAddress("Ntracks", &Ntracks, &b_Ntracks);
   fChain->SetBranchAddress("trigPrescale", trigPrescale, &b_trigPrescale);
   fChain->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
   fChain->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
   fChain->SetBranchAddress("Reco_QQ_type", Reco_QQ_type, &b_Reco_QQ_type);
   fChain->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);
   fChain->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
   fChain->SetBranchAddress("Reco_QQ_mupl_idx", Reco_QQ_mupl_idx, &b_Reco_QQ_mupl_idx);
   fChain->SetBranchAddress("Reco_QQ_mumi_idx", Reco_QQ_mumi_idx, &b_Reco_QQ_mumi_idx);
   fChain->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
   fChain->SetBranchAddress("Reco_QQ_isCowboy", Reco_QQ_isCowboy, &b_Reco_QQ_isCowboy);
   fChain->SetBranchAddress("Reco_QQ_ctau", Reco_QQ_ctau, &b_Reco_QQ_ctau);
   fChain->SetBranchAddress("Reco_QQ_ctauErr", Reco_QQ_ctauErr, &b_Reco_QQ_ctauErr);
   fChain->SetBranchAddress("Reco_QQ_cosAlpha", Reco_QQ_cosAlpha, &b_Reco_QQ_cosAlpha);
   fChain->SetBranchAddress("Reco_QQ_ctau3D", Reco_QQ_ctau3D, &b_Reco_QQ_ctau3D);
   fChain->SetBranchAddress("Reco_QQ_ctauErr3D", Reco_QQ_ctauErr3D, &b_Reco_QQ_ctauErr3D);
   fChain->SetBranchAddress("Reco_QQ_cosAlpha3D", Reco_QQ_cosAlpha3D, &b_Reco_QQ_cosAlpha3D);
   fChain->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);
   fChain->SetBranchAddress("Reco_QQ_dca", Reco_QQ_dca, &b_Reco_QQ_dca);
   fChain->SetBranchAddress("Reco_QQ_MassErr", Reco_QQ_MassErr, &b_Reco_QQ_MassErr);
   fChain->SetBranchAddress("Reco_QQ_vtx", &Reco_QQ_vtx, &b_Reco_QQ_vtx);
   fChain->SetBranchAddress("Reco_mu_size", &Reco_mu_size, &b_Reco_mu_size);
   fChain->SetBranchAddress("Reco_mu_type", Reco_mu_type, &b_Reco_mu_type);
   fChain->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType, &b_Reco_mu_SelectionType);
   fChain->SetBranchAddress("Reco_mu_charge", Reco_mu_charge, &b_Reco_mu_charge);
   fChain->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom, &b_Reco_mu_4mom);
   fChain->SetBranchAddress("Reco_mu_trig", Reco_mu_trig, &b_Reco_mu_trig);
   fChain->SetBranchAddress("Reco_mu_InTightAcc", Reco_mu_InTightAcc, &b_Reco_mu_InTightAcc);
   fChain->SetBranchAddress("Reco_mu_InLooseAcc", Reco_mu_InLooseAcc, &b_Reco_mu_InLooseAcc);
   fChain->SetBranchAddress("Reco_mu_highPurity", Reco_mu_highPurity, &b_Reco_mu_highPurity);
   fChain->SetBranchAddress("Reco_mu_isPF", Reco_mu_isPF, &b_Reco_mu_isPF);
   fChain->SetBranchAddress("Reco_mu_isTracker", Reco_mu_isTracker, &b_Reco_mu_isTracker);
   fChain->SetBranchAddress("Reco_mu_isGlobal", Reco_mu_isGlobal, &b_Reco_mu_isGlobal);
   fChain->SetBranchAddress("Reco_mu_isSoftCutBased", Reco_mu_isSoftCutBased, &b_Reco_mu_isSoftCutBased);
   fChain->SetBranchAddress("Reco_mu_isHybridSoft", Reco_mu_isHybridSoft, &b_Reco_mu_isHybridSoft);
   fChain->SetBranchAddress("Reco_mu_isMedium", Reco_mu_isMedium, &b_Reco_mu_isMedium);
   fChain->SetBranchAddress("Reco_mu_isTightCutBased", Reco_mu_isTightCutBased, &b_Reco_mu_isTightCutBased);
   fChain->SetBranchAddress("Reco_mu_candType", Reco_mu_candType, &b_Reco_mu_candType);
   fChain->SetBranchAddress("Reco_mu_nPixValHits", Reco_mu_nPixValHits, &b_Reco_mu_nPixValHits);
   fChain->SetBranchAddress("Reco_mu_nMuValHits", Reco_mu_nMuValHits, &b_Reco_mu_nMuValHits);
   fChain->SetBranchAddress("Reco_mu_nTrkHits", Reco_mu_nTrkHits, &b_Reco_mu_nTrkHits);
   fChain->SetBranchAddress("Reco_mu_normChi2_inner", Reco_mu_normChi2_inner, &b_Reco_mu_normChi2_inner);
   fChain->SetBranchAddress("Reco_mu_normChi2_global", Reco_mu_normChi2_global, &b_Reco_mu_normChi2_global);
   fChain->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea, &b_Reco_mu_nPixWMea);
   fChain->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea, &b_Reco_mu_nTrkWMea);
   fChain->SetBranchAddress("Reco_mu_StationsMatched", Reco_mu_StationsMatched, &b_Reco_mu_StationsMatched);
   fChain->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy, &b_Reco_mu_dxy);
   fChain->SetBranchAddress("Reco_mu_dxyErr", Reco_mu_dxyErr, &b_Reco_mu_dxyErr);
   fChain->SetBranchAddress("Reco_mu_dz", Reco_mu_dz, &b_Reco_mu_dz);
   fChain->SetBranchAddress("Reco_mu_dzErr", Reco_mu_dzErr, &b_Reco_mu_dzErr);
   fChain->SetBranchAddress("Reco_mu_ptErr_inner", Reco_mu_ptErr_inner, &b_Reco_mu_ptErr_inner);
   fChain->SetBranchAddress("Reco_mu_PFChIso", Reco_mu_PFChIso, &b_Reco_mu_PFChIso);
   fChain->SetBranchAddress("Reco_mu_PFNeuIso", Reco_mu_PFNeuIso, &b_Reco_mu_PFNeuIso);
   fChain->SetBranchAddress("Reco_mu_PFPhoIso", Reco_mu_PFPhoIso, &b_Reco_mu_PFPhoIso);
   fChain->SetBranchAddress("Reco_mu_PFPUIso", Reco_mu_PFPUIso, &b_Reco_mu_PFPUIso);
   Notify();
}

bool pp_with_iso::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

void pp_with_iso::Show(Long64_t entry)
{
   // Print contents of entry.
   // If entry is not specified, print current entry
   if (!fChain)
      return;
   fChain->Show(entry);
}
Int_t pp_with_iso::Cut(Long64_t entry)
{
   // This function may be called from Loop.
   // returns  1 if entry is accepted.
   // returns -1 otherwise.
   return 1;
}
#endif // #ifdef pp_with_iso_cxx
