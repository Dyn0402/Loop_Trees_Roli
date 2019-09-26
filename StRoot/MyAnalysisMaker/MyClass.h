//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Aug  2 09:54:26 2019 by ROOT version 5.34/30
// from TTree nsmTree/nsmTree
// found on file: ../output/7GeV/auau200_C697A3EB48DE1A1E2BD155655A98189F_96.root
//////////////////////////////////////////////////////////

#ifndef MyClass_h
#define MyClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <TClonesArray.h>
#include <TObject.h>

// Fixed size dimensions of array or collections stored in the TTree if any.
   const Int_t kMaxProton = 100;

class MyClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           Proton_;
   UInt_t          Proton_fUniqueID[kMaxProton];   //[Proton_]
   UInt_t          Proton_fBits[kMaxProton];   //[Proton_]
   Float_t         Proton_pt[kMaxProton];   //[Proton_]
   Float_t         Proton_p[kMaxProton];   //[Proton_]
   Float_t         Proton_phi[kMaxProton];   //[Proton_]
   Float_t         Proton_eta[kMaxProton];   //[Proton_]
   Float_t         Proton_dca[kMaxProton];   //[Proton_]
   Float_t         Proton_nsigma[kMaxProton];   //[Proton_]
   Float_t         Proton_beta[kMaxProton];   //[Proton_]
   Float_t         Proton_charge[kMaxProton];   //[Proton_]
 //nsmEvent        *Event;
   UInt_t          fUniqueID;
   UInt_t          fBits;
   Float_t         vtx_x;
   Float_t         vtx_y;
   Float_t         vtx_z;
   UInt_t          Nprim;
   UInt_t          run;
   UInt_t          ref2;
   Float_t         psi;
   UInt_t          btof;

   // List of branches
   TBranch        *b_Proton_;   //!
   TBranch        *b_Proton_fUniqueID;   //!
   TBranch        *b_Proton_fBits;   //!
   TBranch        *b_Proton_pt;   //!
   TBranch        *b_Proton_p;   //!
   TBranch        *b_Proton_phi;   //!
   TBranch        *b_Proton_eta;   //!
   TBranch        *b_Proton_dca;   //!
   TBranch        *b_Proton_nsigma;   //!
   TBranch        *b_Proton_beta;   //!
   TBranch        *b_Proton_charge;   //!
   TBranch        *b_Event_fUniqueID;   //!
   TBranch        *b_Event_fBits;   //!
   TBranch        *b_Event_vtx_x;   //!
   TBranch        *b_Event_vtx_y;   //!
   TBranch        *b_Event_vtx_z;   //!
   TBranch        *b_Event_Nprim;   //!
   TBranch        *b_Event_run;   //!
   TBranch        *b_Event_ref2;   //!
   TBranch        *b_Event_psi;   //!
   TBranch        *b_Event_btof;   //!

   MyClass(TTree *tree=0);
   virtual ~MyClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MyClass_cxx
MyClass::MyClass(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../output/7GeV/auau200_C697A3EB48DE1A1E2BD155655A98189F_96.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../output/7GeV/auau200_C697A3EB48DE1A1E2BD155655A98189F_96.root");
      }
      f->GetObject("nsmTree",tree);

   }
   Init(tree);
}

MyClass::~MyClass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MyClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MyClass::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MyClass::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Proton", &Proton_, &b_Proton_);
   fChain->SetBranchAddress("Proton.fUniqueID", &Proton_fUniqueID, &b_Proton_fUniqueID);
   fChain->SetBranchAddress("Proton.fBits", &Proton_fBits, &b_Proton_fBits);
   fChain->SetBranchAddress("Proton.pt", &Proton_pt, &b_Proton_pt);
   fChain->SetBranchAddress("Proton.p", &Proton_p, &b_Proton_p);
   fChain->SetBranchAddress("Proton.phi", &Proton_phi, &b_Proton_phi);
   fChain->SetBranchAddress("Proton.eta", &Proton_eta, &b_Proton_eta);
   fChain->SetBranchAddress("Proton.dca", &Proton_dca, &b_Proton_dca);
   fChain->SetBranchAddress("Proton.nsigma", &Proton_nsigma, &b_Proton_nsigma);
   fChain->SetBranchAddress("Proton.beta", &Proton_beta, &b_Proton_beta);
   fChain->SetBranchAddress("Proton.charge", &Proton_charge, &b_Proton_charge);
   fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_Event_fUniqueID);
   fChain->SetBranchAddress("fBits", &fBits, &b_Event_fBits);
   fChain->SetBranchAddress("vtx_x", &vtx_x, &b_Event_vtx_x);
   fChain->SetBranchAddress("vtx_y", &vtx_y, &b_Event_vtx_y);
   fChain->SetBranchAddress("vtx_z", &vtx_z, &b_Event_vtx_z);
   fChain->SetBranchAddress("Nprim", &Nprim, &b_Event_Nprim);
   fChain->SetBranchAddress("run", &run, &b_Event_run);
   fChain->SetBranchAddress("ref2", &ref2, &b_Event_ref2);
   fChain->SetBranchAddress("psi", &psi, &b_Event_psi);
   fChain->SetBranchAddress("btof", &btof, &b_Event_btof);
   Notify();
}

Bool_t MyClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MyClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MyClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MyClass_cxx
