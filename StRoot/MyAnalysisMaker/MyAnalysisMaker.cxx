#include "MyAnalysisMaker.h"
#include "MyClass.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <fstream>
#include <cstdlib>

#include "math.h"
#include "string.h"

#include "StMaker.h"
#include "TStyle.h"
#include "TObject.h"
#include <TSystem.h>
#include "TROOT.h"
#include "TNamed.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TChain.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TRandom.h"
#include "TProfile.h"
#include "TObjArray.h"
#include "TClonesArray.h"

using namespace std;

ClassImp(MyAnalysisMaker)                       // Macro for CINT compatibility

MyAnalysisMaker::MyAnalysisMaker(const char *name) : StMaker(name){
    
}

Double_t MyAnalysisMaker::doLoop(char *inputfile, char* outname, int energy){

  for(int i = 0; i < 100; i++){
    BufferEvent_NEvents[i] = 0;
    BufferEvent_Full[i] = 0;
  }
    
    cout << "inputfile " << inputfile << endl;
    chain = new TChain("nsmTree");
    MyClass *s = new MyClass(chain);
    
    TString inputfilename = inputfile;
    string dirFile = inputfilename.Data();
    
    if (dirFile.find(".list")!=string::npos) {
        ifstream inputStream(dirFile.c_str());
        
        if(!(inputStream.good())) {
            LOG_ERROR << "ERROR: Cannot open list file " << dirFile << endm;
        }
        
        char line[512];
        int nFile=0;
        string ltest;
        
        while (inputStream.good()) {
            inputStream.getline(line,512);
            string aFile = line;
            
            if (inputStream.good() && aFile.find(".root")!=string::npos) {
                TFile *ftmp = new TFile(line);
                
                if(ftmp && ftmp->IsOpen() && ftmp->GetNkeys()) {
                    LOG_INFO << " Read in ntuple file " << line << endm;
                    chain->Add(line);
                    nFile++;
                }
            }
        }
        
        LOG_INFO << " Total " << nFile << " files have been read in. " << endm;
        
    } else if (dirFile.find("root")!=string::npos) {
        
        chain->Add(dirFile.c_str());
        
    } else {
        
        LOG_WARN << " No good input file to read ... " << endm;
    }
    
    //========== Define global variables =========
    
    mEventsProcessed = 0;
    Pi = 6.28318530718/2.;
    iran = 0;
    
    histogram_output = new TFile(outname,"RECREATE") ;
    declareHistograms();
    
    //============= Event Loop starts ============
    
    Int_t nentries = chain->GetEntries();
    
    int nbytes=0;
    for(int ievt = 0; ievt < nentries; ievt++){
        nbytes += chain->GetEntry(ievt);
        
        EventCount->Fill(0.2); //record # of events
        event_cut_hist.Fill(0);
        
	if((ievt+1)%1000000==1) cout<<"Processing entry == "<< ievt+1 <<" == out of "<<nentries<<" last event centrality " << cent <<".\n" ;
        
        // --- read event variables ---
		// Dylan Edit 10/9/19 Start
		int current_run = s->run;
		bool bad_run_bool = false;
		int bad_runs[1000] = {11116006, 11116017, 11116055, 11117094, 11120084, 11122060, 11122080, 11123050, 11123056, 11125099, 11125102, 11126022, 11126068, 11127007, 11127078, 11129047, 11130043, 11130064, 11131027, 11132025, 11132058, 11132067, 11133085, 11134019, 11134053, 11134066, 11135025, 11136059, 11136077, 11137111, 11137114, 11138005, 11138047, 11138061, 11138066, 11138082, 11141019, 11142008, 11142018, 11142036, 11143027, 11143047, 11145101, 11148001, 11148006, 11148014, 11148045, 11148055, 11149001, 11149010, 11149015, 11149017, 11149047, 11150029, 11151050, 11151051, 11152036, 11153032, 11153042, 11154026, 11156003, 11156009, 11156036, 11156043, 11156044, 11156045, 12114007, 12114092, 12114103, 12115014, 12115015, 12115016, 12115018, 12115019, 12115020, 12115022, 12115023, 12115073, 12120026, 12121017, 12121022, 12121034, 12122019, 12172049, 12172050, 12172051, 12173031, 12173034, 12174077, 12174085, 12175030, 12175040, 12175062, 12176044, 12176047, 12176067, 12176069, 12176104, 12177015, 12177101, 12178051, 12178120, 12179068, 11002145, 11003001, 11003002, 11003004, 11003005, 11003006, 11003012, 11003013, 11023046, 11023048, 11025034, 11025054, 11026010, 11026114, 11028004, 11028005, 11028006, 11028018, 11028019, 11028020, 11028021, 11028022, 11030041, 11030042, 11057012, 11057035, 11059060, 11060059, 11064003, 11070051, 11070052, 11070053, 11070054, 11070055, 11070056, 11070057, 11071001, 11071002, 11071003, 11071004, 11071006, 11071007, 11071013, 11071014, 11071015, 11071016, 11071017, 11071018, 11071019, 11071020, 11072030, 11072031, 11072032, 11072044, 11072045, 11073001, 11073003, 11073030, 11073031, 11080072, 11084009, 11084010, 11087001, 11087002, 11087003, 11087004, 11087057, 11087058, 11087059, 11089026, 11089028, 11089029, 11089048, 11089055, 11092012, 11092031, 11092032, 11092033, 11092034, 11095040, 11095048, 11095049, 11095050, 11095051, 11095061, 11095062, 11095063, 11095064, 15047023, 15047065, 15047071, 15047087, 15049027, 15051131, 15053034, 15053035, 15053052, 15053053, 15053054, 15053055, 15053056, 15053057, 15061001, 15066064, 15067026, 11099124, 11100045, 11101030, 11101046, 11101104, 11102054, 11102098, 11103008, 11103009, 11103058, 11103062, 11105011, 11105018, 11105029, 11105045, 11106026, 11106027, 11106028, 11106029, 11106030, 11106039, 11106040, 11106041, 11106049, 11109013, 11109088, 11109090, 11109092, 11109097, 11109098, 11109102, 11109103, 11109104, 11110034, 11110041, 11110073, 11110076, 0, 0};
        for(int b_run=0; b_run<1000; b_run++) {
        	if(bad_runs[b_run]==0) break;
        	if(current_run == bad_runs[b_run]) { bad_run_bool = true; }
        }
        if(bad_run_bool) { continue; }
        EventCount->Fill(1.5);
        event_cut_hist.Fill(1);
        // Dylan Edit 10/9/19 End
        
        nProton = s->Proton_;
        
        if (nProton < 2) continue;
        EventCount -> Fill(2.5);
        event_cut_hist.Fill(2);
        
        tofmult = s->btof;
        refMult = s->Nprim;
        htr->Fill(tofmult, refMult);
        // Dylan Edit 10/9/19 Start
        if(energy == 27) {
        	if((double)tofmult / refMult > (double)(1878-207)/(400+18) || (double)tofmult / refMult < (double)(1211+39)/(460-6)) { continue; }
        } else if(energy == 11) {
        	if((double)tofmult / refMult > (double)(1274-172)/(250+6) || (double)tofmult / refMult < (double)(660+20)/(307-7)) { continue; }
        } else if(energy == 7) {
			if((double)tofmult / refMult > (double)(1565-187)/(310+6) || (double)tofmult / refMult < (double)(924+34)/(418-7)) { continue; }
		} else if(energy == 19) {
			if((double)tofmult / refMult > (double)(1596-195)/(364+11) || (double)tofmult / refMult < (double)(1049+35)/(426-5)) { continue; }
		}
        // Dylan Edit 10/9/19 End
        EventCount->Fill(3.5);
        event_cut_hist.Fill(3);

        refMult2 = s->ref2;
		cent = getCentrality(refMult2,energy);
		if(cent < 0 || cent > 9) continue;
		CurrentEvent_centrality = cent;
		EventCount -> Fill(4.5);
        
        Psi = s->psi;
        if(Psi < 0) Psi = Psi + Pi;
        CurrentEvent_Psi = Psi;
        
        CurrentEvent_nProton = 0;
        nParticle = 0; num_particle_1 = 0; num_particle_2 = 0; num_particle_3 = 0; num_particle_4 = 0; num_particle_5 = 0; num_particle_6 = 0;
        
        for(int trk = 0; trk < nProton; trk++){
            
            check = 0;
            track_cut_hist.Fill(0);
            
            pt = s->Proton_pt[trk];
            p = s->Proton_p[trk];
            phi = s->Proton_phi[trk];
            beta = s->Proton_beta[trk];
            dca = s->Proton_dca[trk];
            nsigma = s->Proton_nsigma[trk];
            charge = s->Proton_charge[trk];
            
            if(p < 0.15) continue;
            track_cut_hist.Fill(1);
            if(pt < 0.4 || pt > 2.0) continue;
            track_cut_hist.Fill(2);
            if(charge != 1) continue;
            track_cut_hist.Fill(3);
            if(fabs(eta) > 0.5) continue;
            track_cut_hist.Fill(4);
            if(fabs(nsigma) > 2.0) continue;
            track_cut_hist.Fill(5);
            if(dca < 0 || dca > 2.0) continue;
            track_cut_hist.Fill(6);
            
            if (pt < 0.8){
                check = 1;
            }
            
            if (pt >= 0.8){
                m = -999.;
                if(beta>0)  m = p*p*(1./(beta*beta)-1.);
                if(m!=-999 && m>0.800 && m <1.0) check = 1;
            }
            
            if(check == 1){
            	track_cut_hist.Fill(7);
                
                CurrentEvent_ProtonPhi[CurrentEvent_nProton] = phi;
                CurrentEvent_nProton++;
                
                nParticle++ ;
                phibin = getPhiBin(phi);
                

                if(phibin==1)num_particle_1++;
                if(phibin==2)num_particle_2++;
                //if(phibin==3)num_particle_3++;
                //if(phibin==4)num_particle_4++;
                //if(phibin==5)num_particle_5++;
                //if(phibin==6)num_particle_6++;
            }
            
        }//==================Proton loop ends=========================

        if(CurrentEvent_nProton < 2) continue;
        if(nParticle < 2) continue;
        
        event_cut_hist.Fill(4);

        cent9_events.Fill(cent);

        if(cent > 0 && cent < 10 &&  nParticle > 1){
            
            hnParticle[1][cent]->Fill(nParticle);
            
            ratio = (double) num_particle_1/nParticle; FillHistograms(ratio, refMult2, 1, cent); hratio[1][cent]->Fill(ratio);
            if(ratio == 1) hnParticle1[1][cent]->Fill(nParticle);
            
            ratio = (double) num_particle_2/nParticle; FillHistograms(ratio, refMult2, 1, cent); hratio[1][cent]->Fill(ratio);
            if(ratio == 1) hnParticle1[1][cent]->Fill(nParticle);
	    /*            
            ratio = (double) num_particle_3/nParticle; FillHistograms(ratio, refMult2, 1, cent); hratio[1][cent]->Fill(ratio);
            if(ratio == 1) hnParticle1[1][cent]->Fill(nParticle);
            
            ratio = (double) num_particle_4/nParticle; FillHistograms(ratio, refMult2, 1, cent); hratio[1][cent]->Fill(ratio);
            if(ratio == 1) hnParticle1[1][cent]->Fill(nParticle);
            
            ratio = (double) num_particle_5/nParticle; FillHistograms(ratio, refMult2, 1, cent); hratio[1][cent]->Fill(ratio);
            if(ratio == 1) hnParticle1[1][cent]->Fill(nParticle);
            
            ratio = (double) num_particle_6/nParticle; FillHistograms(ratio, refMult2, 1, cent); hratio[1][cent]->Fill(ratio);
            if(ratio == 1) hnParticle1[1][cent]->Fill(nParticle);
            */
        }

        
        BufferPointer = -1;
        BufferPointer = (int)(CurrentEvent_centrality)*10+(int)((10*CurrentEvent_Psi)/Pi);
        if(BufferPointer>=100||BufferPointer<0) continue;

        copyCurrentToBuffer(BufferPointer);
        if (mEventsProcessed > 100)makeMixedPairs(BufferPointer);

        mEventsProcessed++ ;
        EventCount->Fill(7.2);
        
    }//==================Event loop ends=========================
    
    
    histogram_output->cd();
    writeHistograms();
    histogram_output->Close();
    
    cout <<"\n ======> All done <======"<<endl;
    cout<<" Acutal #Events Processed = " <<mEventsProcessed<<"\n###### Thank You ######\n"<< endl ;
    
    delete chain;
    return 1;
    
}

MyAnalysisMaker::~MyAnalysisMaker() {
}

Int_t MyAnalysisMaker::getCentrality(int mult, double energy){
    
    int central = -1;
    
    if(energy == 7){
        
        float centFull[9] = {32,41,51,64,78,95,114,137,165};
        if      (mult>=centFull[8]) central=9;
        else if (mult>=centFull[7]) central=8;
        else if (mult>=centFull[6]) central=7;
        else if (mult>=centFull[5]) central=6;
        else if (mult>=centFull[4]) central=5;
        else if (mult>=centFull[3]) central=4;
        else if (mult>=centFull[2]) central=3;
        else if (mult>=centFull[1]) central=2;
        else if (mult>=centFull[0]) central=1;
        
    }
    
    else if(energy == 11){
        
        float centFull[9] = {41,52,65,80,98,118,143,172,206};
        if      (mult>=centFull[8]) central=9;
        else if (mult>=centFull[7]) central=8;
        else if (mult>=centFull[6]) central=7;
        else if (mult>=centFull[5]) central=6;
        else if (mult>=centFull[4]) central=5;
        else if (mult>=centFull[3]) central=4;
        else if (mult>=centFull[2]) central=3;
        else if (mult>=centFull[1]) central=2;
        else if (mult>=centFull[0]) central=1;
        
    }
    
    else if(energy == 19){
        
        float centFull[9] = {51,65,81,100,123,149,180,215,258};
        if      (mult>=centFull[8]) central=9;
        else if (mult>=centFull[7]) central=8;
        else if (mult>=centFull[6]) central=7;
        else if (mult>=centFull[5]) central=6;
        else if (mult>=centFull[4]) central=5;
        else if (mult>=centFull[3]) central=4;
        else if (mult>=centFull[2]) central=3;
        else if (mult>=centFull[1]) central=2;
        else if (mult>=centFull[0]) central=1;
        
    }
    
    else if(energy == 27){
        
        float centFull[9] = {56,71,90,111,135,164,198,237,284};
        if      (mult>=centFull[8]) central=9;
        else if (mult>=centFull[7]) central=8;
        else if (mult>=centFull[6]) central=7;
        else if (mult>=centFull[5]) central=6;
        else if (mult>=centFull[4]) central=5;
        else if (mult>=centFull[3]) central=4;
        else if (mult>=centFull[2]) central=3;
        else if (mult>=centFull[1]) central=2;
        else if (mult>=centFull[0]) central=1;
        
    }
    
    else if(energy == 39){
        
        float centFull[9] = {61,78,97,121,147,179,215,257,307};
        if      (mult>=centFull[8]) central=9;
        else if (mult>=centFull[7]) central=8;
        else if (mult>=centFull[6]) central=7;
        else if (mult>=centFull[5]) central=6;
        else if (mult>=centFull[4]) central=5;
        else if (mult>=centFull[3]) central=4;
        else if (mult>=centFull[2]) central=3;
        else if (mult>=centFull[1]) central=2;
        else if (mult>=centFull[0]) central=1;
        
    }
    
    else if(energy == 62){
        
        float centFull[9] = {66,84,106,131,160,194,233,279,334};
        if      (mult>=centFull[8]) central=9;
        else if (mult>=centFull[7]) central=8;
        else if (mult>=centFull[6]) central=7;
        else if (mult>=centFull[5]) central=6;
        else if (mult>=centFull[4]) central=5;
        else if (mult>=centFull[3]) central=4;
        else if (mult>=centFull[2]) central=3;
        else if (mult>=centFull[1]) central=2;
        else if (mult>=centFull[0]) central=1;
        
    }
    
    else if(energy == 200){
        
        float centFull[9] = {85,108,135,167,204,247,297,355,421};
        if      (mult>=centFull[8]) central=9;
        else if (mult>=centFull[7]) central=8;
        else if (mult>=centFull[6]) central=7;
        else if (mult>=centFull[5]) central=6;
        else if (mult>=centFull[4]) central=5;
        else if (mult>=centFull[3]) central=4;
        else if (mult>=centFull[2]) central=3;
        else if (mult>=centFull[1]) central=2;
        else if (mult>=centFull[0]) central=1;
        
    }
    
    return central;
    
}


Int_t MyAnalysisMaker::getPhiBin(double phi){
    
    int N = 2;
    int bin=-1;
    double twoPi = 6.28318530718;
    
    if(phi > 0 && phi <= twoPi/N) bin=1;
    if(phi > twoPi/N && phi <= 2*twoPi/N) bin=2;
    if(phi > 2*twoPi/N && phi <= 3*twoPi/N) bin=3;
    if(phi > 3*twoPi/N && phi <= 4*twoPi/N) bin=4;
    if(phi > 4*twoPi/N && phi <= 5*twoPi/N) bin=5;
    if(phi > 5*twoPi/N && phi <= 6*twoPi/N) bin=6;
    
    return bin;
}

void MyAnalysisMaker::declareHistograms(){
    
    EventCount    = new TH1F("EventCount","EventCount",10,0,10);
    htr           = new TH2F("htr","tofmult vs refmult",3001,-0.5,3000.5,601,-0.5,600.5);
    
	event_cut_hist = TH1I("event_cut", "Event Cuts", 5, -0.5, 4.5);
	track_cut_hist = TH1I("track_cut", "Track Cuts", 8, -0.5, 7.5);
	cent16_events = TH1I("cent16_events", "Cent16 Events", 18, -1.5, 16.5);
	cent9_events = TH1I("cent9_events", "Cent9 Events", 11, -1.5, 9.5);

    const Double_t binSize=700.0, minBin=-0.5, maxBin=699.5;
    
    for(int i = 1; i < 3; i++){
        
        for(int j = 1; j < 10; j++){
            
            sprintf(name,"hnParticle_%d_%d",i,j);
            hnParticle[i][j] = new TH1F(name,"Count",100,0,100);
            
            sprintf(name,"hratio_%d_%d",i,j);
            hratio[i][j] = new TH1F(name,"Ratio",100,-0.05,1.05);
            
            sprintf(name,"hnParticle1_%d_%d",i,j);
            hnParticle1[i][j] = new TH1F(name,"Count_ratio1",20,0,20);
            
            sprintf(name,"p1_%d_%d",i,j);
            h_p1[i][j]          = new TProfile(name,"",binSize,minBin,maxBin);
            sprintf(name,"p2_%d_%d",i,j);
            h_p2[i][j]          = new TProfile(name,"",binSize,minBin,maxBin);
            sprintf(name,"p3_%d_%d",i,j);
            h_p3[i][j]          = new TProfile(name,"",binSize,minBin,maxBin);
            sprintf(name,"p4_%d_%d",i,j);
            h_p4[i][j]          = new TProfile(name,"",binSize,minBin,maxBin);
            sprintf(name,"p5_%d_%d",i,j);
            h_p5[i][j]          = new TProfile(name,"",binSize,minBin,maxBin);
            sprintf(name,"p6_%d_%d",i,j);
            h_p6[i][j]          = new TProfile(name,"",binSize,minBin,maxBin);
            sprintf(name,"p7_%d_%d",i,j);
            h_p7[i][j]          = new TProfile(name,"",binSize,minBin,maxBin);
            sprintf(name,"p8_%d_%d",i,j);
            h_p8[i][j]          = new TProfile(name,"",binSize,minBin,maxBin);
            
        }
        
    }
    
}


void MyAnalysisMaker::FillHistograms(float ratio, int refMult, int i, int cent){
    
    h_p1[i][cent]->Fill(refMult,ratio);
    h_p2[i][cent]->Fill(refMult,pow(ratio,2));
    h_p3[i][cent]->Fill(refMult,pow(ratio,3));
    h_p4[i][cent]->Fill(refMult,pow(ratio,4));
    h_p5[i][cent]->Fill(refMult,pow(ratio,5));
    h_p6[i][cent]->Fill(refMult,pow(ratio,6));
    h_p7[i][cent]->Fill(refMult,pow(ratio,7));
    h_p8[i][cent]->Fill(refMult,pow(ratio,8));
    
}

void MyAnalysisMaker::writeHistograms(){
    
    EventCount->Write();
    htr->Write();
    
    event_cut_hist.Write();
	track_cut_hist.Write();
	cent16_events.Write();
	cent9_events.Write();

    for(int i = 1; i < 3; i++){
        
        for(int j = 1; j < 10; j++){
            
            hnParticle[i][j]->Write();
            hratio[i][j]->Write();
            hnParticle1[i][j]->Write();
            
            h_p1[i][j]->Write();
            h_p2[i][j]->Write();
            h_p3[i][j]->Write();
            h_p4[i][j]->Write();
            h_p5[i][j]->Write();
            h_p6[i][j]->Write();
            h_p7[i][j]->Write();
            h_p8[i][j]->Write();
        }
    }
}

void MyAnalysisMaker::makeMixedPairs(int bufferPointer){
    
    TRandom ran(0);
    
    nParticle = 0; num_particle_1 = 0; num_particle_2 = 0; num_particle_3 = 0; num_particle_4 = 0; num_particle_5 = 0; num_particle_6 = 0;

   int nmax = CurrentEvent_nProton;

    while(nParticle != nmax){

        phi = -999;
        
        while(phi < -99){

	  int k = ran.Rndm() * BufferEvent_NEvents[bufferPointer];
	  int j = ran.Rndm() * BufferEvent_nProton[k*100+bufferPointer];
	  phi = BufferEvent_ProtonPhi[k*100+bufferPointer][j]; 
	    
        }
        
        nParticle++ ;
        phibin = getPhiBin(phi);
        
        if(phibin==1)num_particle_1++;
        if(phibin==2)num_particle_2++;
        //if(phibin==3)num_particle_3++;
        //if(phibin==4)num_particle_4++;
        //if(phibin==5)num_particle_5++;
        //if(phibin==6)num_particle_6++;
        
    }
    
    if(CurrentEvent_centrality > 0 && CurrentEvent_centrality < 10){
        
        if(nParticle != nmax) cout << "Tragedy!!" << endl;
        
        cent = CurrentEvent_centrality;
        hnParticle[2][cent]->Fill(nParticle);
        
        ratio = (double) num_particle_1/nParticle; FillHistograms(ratio, refMult2, 2, cent); hratio[2][cent]->Fill(ratio);
        if(ratio == 1) hnParticle1[2][cent]->Fill(nParticle);
        
        ratio = (double) num_particle_2/nParticle; FillHistograms(ratio, refMult2, 2, cent); hratio[2][cent]->Fill(ratio);
        if(ratio == 1) hnParticle1[2][cent]->Fill(nParticle);
        /*
        ratio = (double) num_particle_3/nParticle; FillHistograms(ratio, refMult2, 2, cent); hratio[2][cent]->Fill(ratio);
        if(ratio == 1) hnParticle1[2][cent]->Fill(nParticle);
        
        ratio = (double) num_particle_4/nParticle; FillHistograms(ratio, refMult2, 2, cent); hratio[2][cent]->Fill(ratio);
        if(ratio == 1) hnParticle1[2][cent]->Fill(nParticle);
        
        ratio = (double) num_particle_5/nParticle; FillHistograms(ratio, refMult2, 2, cent); hratio[2][cent]->Fill(ratio);
        if(ratio == 1) hnParticle1[2][cent]->Fill(nParticle);
        
        ratio = (double) num_particle_6/nParticle; FillHistograms(ratio, refMult2, 2, cent); hratio[2][cent]->Fill(ratio);
        if(ratio == 1) hnParticle1[2][cent]->Fill(nParticle);
        */
    }
}

void MyAnalysisMaker::copyCurrentToBuffer(int bufferPointer){
    
    if(BufferEvent_NEvents[bufferPointer] >= 100) BufferEvent_Full[bufferPointer] = 1;
    
    TRandom3 *gRandom = new TRandom3(iran++);
    int eventPointer = -1;
    
    if(BufferEvent_Full[bufferPointer]){ // full - random rewrite one
        
        do{ double rrr = gRandom->Rndm();
            eventPointer = (int)(100.*(1.0 - rrr));
        } while(eventPointer<0||eventPointer>=100);
    }
    else { // not full
        eventPointer = BufferEvent_NEvents[bufferPointer];
    }
    delete gRandom;
    
    BufferEvent_Psi[bufferPointer+100*eventPointer] = CurrentEvent_Psi;
    
    BufferEvent_nProton[bufferPointer+100*eventPointer] = CurrentEvent_nProton;
    
    for(int i = 0; i < CurrentEvent_nProton; i++){
        
        BufferEvent_ProtonPhi[bufferPointer+100*eventPointer][i]=CurrentEvent_ProtonPhi[i];
        
    }
    
    if(BufferEvent_NEvents[bufferPointer]<100) {
        BufferEvent_NEvents[bufferPointer]+=1;
    }

}



