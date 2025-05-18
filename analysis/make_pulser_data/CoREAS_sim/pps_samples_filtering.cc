#include <stdio.h>
#include <dirent.h>
#include <sys/types.h>
#include "utils.cc"

using namespace std;

void GetTChain(string directory, TChain *chaintree); 

void pps_samples_filtering(){


  Short_t T1H[1500];
  Short_t T2H[1500];
  Short_t T3H[1500];
  Short_t T4H[1500];
  Short_t T1V[1500];
  Short_t T2V[1500];
  Short_t T3V[1500];
  Short_t T4V[1500];
  TBits* triggerBits = 0;
  double eventTime = 0;
  Int_t gDate = 0;
  Int_t gEventNum = 0;
  Double_t gEventRate = 0.;
  Int_t gThresholdStep = 0;


	TChain *pps_samples_Tree = new TChain("t");
	pps_samples_Tree->AddFile("pps_samples_2109.root");
	pps_samples_Tree->AddFile("pps_samples_2110.root");
	pps_samples_Tree->AddFile("pps_samples_2111.root");
	pps_samples_Tree->AddFile("pps_samples_2112.root");
	pps_samples_Tree->AddFile("pps_samples_2201.root");
	pps_samples_Tree->AddFile("pps_samples_2202.root");
	pps_samples_Tree->AddFile("pps_samples_2203.root");
	pps_samples_Tree->AddFile("pps_samples_2204.root");
	pps_samples_Tree->AddFile("pps_samples_2205.root");
	
	
 	
 	TFile *new_pps_samples = new TFile("pps_samples_filtered.root", "RECREATE");
 	TTree *new_pps_samples_Tree = new TTree("t","pps samples");

  pps_samples_Tree->SetBranchAddress("date", &gDate);
  pps_samples_Tree->SetBranchAddress("EventNum", &gEventNum);
  pps_samples_Tree->SetBranchAddress("eventTime", &eventTime); //PC timestamp
  pps_samples_Tree->SetBranchAddress("EventRate", &gEventRate);
  pps_samples_Tree->SetBranchAddress("ThresholdStep", &gThresholdStep);
  pps_samples_Tree->SetBranchAddress("T1H", T1H);
  pps_samples_Tree->SetBranchAddress("T2H", T2H);
  pps_samples_Tree->SetBranchAddress("T3H", T3H);
  pps_samples_Tree->SetBranchAddress("T4H", T4H);
  pps_samples_Tree->SetBranchAddress("T1V", T1V);
  pps_samples_Tree->SetBranchAddress("T2V", T2V);
  pps_samples_Tree->SetBranchAddress("T3V", T3V);
  pps_samples_Tree->SetBranchAddress("T4V", T4V);
  
  
  new_pps_samples_Tree->Branch("date", &gDate, "date/I");
  new_pps_samples_Tree->Branch("EventNum", &gEventNum, "EventNum/I");
  new_pps_samples_Tree->Branch("eventTime", &eventTime, "eventTime/D"); //PC timestamp
  new_pps_samples_Tree->Branch("EventRate", &gEventRate, "EventRate/D");
  new_pps_samples_Tree->Branch("ThresholdStep", &gThresholdStep, "ThresholdStep/I");
  new_pps_samples_Tree->Branch("T1H", T1H, "T1H[1500]/S");
  new_pps_samples_Tree->Branch("T2H", T2H, "T2H[1500]/S");
  new_pps_samples_Tree->Branch("T3H", T3H, "T3H[1500]/S");
  new_pps_samples_Tree->Branch("T4H", T4H, "T4H[1500]/S");
  new_pps_samples_Tree->Branch("T1V", T1V, "T1V[1500]/S");
  new_pps_samples_Tree->Branch("T2V", T2V, "T2V[1500]/S");
  new_pps_samples_Tree->Branch("T3V", T3V, "T3V[1500]/S");
  new_pps_samples_Tree->Branch("T4V", T4V, "T4V[1500]/S");
	


	for(int entry=0; entry<pps_samples_Tree->GetEntries(); entry++){
	
		pps_samples_Tree->GetEntry(entry);
	
		double x[8][4096] = {0.};
		for(int i=0; i<1500; i++){
			x[0][i] = T1H[i] *500./32512.;
			x[2][i] = T2H[i] *500./32512.;
			x[4][i] = T3H[i] *500./32512.;
			x[6][i] = T4H[i] *500./32512.;
			x[1][i] = T1V[i] *500./32512.;
			x[3][i] = T2V[i] *500./32512.;
			x[5][i] = T3V[i] *500./32512.;
			x[7][i] = T4V[i] *500./32512.;
		}
		
		for(int ch=0; ch<8; ch++){
			filter_fre(x[ch],0);
		}

		bool bad_LNA = false;
		int count_over=0;
		for(int j=0; j<8; j++){
			for(int i=0; i<1500; i++){
				// rms=10, 3 sigma
				if(fabs(x[j][i])>30){count_over++;}
				if(fabs(x[j][i])>100){bad_LNA = true;}
			}
		}
		if(count_over>80){continue;}
		if(bad_LNA){continue;}
		if((gThresholdStep<480)||(gThresholdStep>1280)){continue;}
		//cout<<entry<<" cyc "<<count_over<<"  "<<gEventRate<<"   "<<gThresholdStep<<endl;
		//cout<<"ccc "<<gEventRate<<"   "<<gThresholdStep<<endl;
				
		for(int i=0; i<1500; i++){
			T1H[i] = (Short_t) round(x[0][i] * 32512./500.);
			T2H[i] = (Short_t) round(x[2][i] * 32512./500.);
			T3H[i] = (Short_t) round(x[4][i] * 32512./500.);
			T4H[i] = (Short_t) round(x[6][i] * 32512./500.);
			
			T1V[i] = (Short_t) round(x[1][i] * 32512./500.);
			T2V[i] = (Short_t) round(x[3][i] * 32512./500.);
			T3V[i] = (Short_t) round(x[5][i] * 32512./500.);
			T4V[i] = (Short_t) round(x[7][i] * 32512./500.);
			
		}
		
		new_pps_samples_Tree->Fill();
	}
	    

	new_pps_samples->cd();
    new_pps_samples_Tree->Write();
    new_pps_samples->Close();
    
	
}





