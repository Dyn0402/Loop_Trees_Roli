#ifndef MyAnalysisMaker_def
#define MyAnalysisMaker_def

#include <fstream>
#include <math.h>
#include <iostream>

#include "TNamed.h"
#include "TClonesArray.h"
#include "StMaker.h"
#include "TObject.h"
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TChain.h"
#include "TBranch.h"
#include "TMath.h"
#include "TString.h"
#include "TDirectory.h"
#include "TF1.h"

#include "MyClass.h"

class TObject;
class TTree;
class TBranch;
class TNtuple;
class TFile;
class TChain;
class TClass;
class StFile;
class TH1D;
class TH1F;
class TH2F;
class TH2D;
class TH3F;
class TH3D;
class TProfile;
class MyClass;


class MyAnalysisMaker : public StMaker
{
public:
    MyAnalysisMaker(const char *name="Loop") ;
     ~MyAnalysisMaker() ;
    
    Double_t  doLoop(char *inputfile, char *outname, int energy = 14);
    Int_t getCentrality(int refmult2, double energy);
    Int_t getPhiBin(double phi);
    void declareHistograms();
    void FillHistograms(float ratio, int refMult, int i, int cent);
    void writeHistograms();
    void makeMixedPairs(int bufferPointer);
    void copyCurrentToBuffer(int bufferPointer);
    
protected:
    
    TChain* chain;                //!
    TFile* n;             	 //!
    
    
private:
   
    TFile *histogram_output;
    
    UInt_t  mEventsProcessed;
    char name[60];
    
    //----- general -----
    float Pi;
    
    //----- event variables -------
    int refMult2,  nProton, nParticle, num_particle_1, num_particle_2, num_particle_3, num_particle_4, num_particle_5, num_particle_6, cent;
    float ratio, Psi, refMult, tofmult;
    
    //----- track variables ------
    int check, phibin;
    float pt, p, phi, beta, m, dca, nsigma, charge, eta;
    
    //----- mixing variables ------
    int CurrentEvent_centrality;
    
    float CurrentEvent_Psi;
    
    int CurrentEvent_nProton;
    float CurrentEvent_ProtonPhi[10000];
    
    int BufferPointer;
    
    int BufferEvent_NEvents[100]; 
    
      int BufferEvent_Full[100];
    
    float BufferEvent_Psi[10000];
    
    int BufferEvent_nProton[10000];
    float BufferEvent_ProtonPhi[10000][10000];
    
    int iran ;
    
    
    //------histograms pointer------------
    TH1F *EventCount;
    TH2F *htr;

	TH1I event_cut_hist;
	TH1I track_cut_hist;
	TH1I cent16_events;
	TH1I cent9_events;


    TH1F *hnParticle[5][16];
    TH1F *hratio[5][16];
    TH1F *hnParticle1[5][16];
    
    TProfile *h_p1[5][16];
    TProfile *h_p2[5][16];
    TProfile *h_p3[5][16];
    TProfile *h_p4[5][16];
    TProfile *h_p5[5][16];
    TProfile *h_p6[5][16];
    TProfile *h_p7[5][16];
    TProfile *h_p8[5][16];

    
    ClassDef(MyAnalysisMaker,1)                       //Macro for CINT compatability
    
};

#endif















