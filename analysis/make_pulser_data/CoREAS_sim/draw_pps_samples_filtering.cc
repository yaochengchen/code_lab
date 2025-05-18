#include <stdio.h>
#include <dirent.h>
#include <sys/types.h>
#include "utils.cc"

using namespace std;

void GetTChain(string directory, TChain *chaintree); 

void draw_pps_samples_filtering(){


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
	
	
 	
 	TFile *new_pps_samples = new TFile("pps_samples_filtered.root");
 	TTree *new_pps_samples_Tree = (TTree *) new_pps_samples->Get("t");

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
  
  
  new_pps_samples_Tree->SetBranchAddress("date", &gDate);
  new_pps_samples_Tree->SetBranchAddress("EventNum", &gEventNum);
  new_pps_samples_Tree->SetBranchAddress("eventTime", &eventTime); //PC timestamp
  new_pps_samples_Tree->SetBranchAddress("EventRate", &gEventRate);
  new_pps_samples_Tree->SetBranchAddress("ThresholdStep", &gThresholdStep);
  new_pps_samples_Tree->SetBranchAddress("T1H", T1H);
  new_pps_samples_Tree->SetBranchAddress("T2H", T2H);
  new_pps_samples_Tree->SetBranchAddress("T3H", T3H);
  new_pps_samples_Tree->SetBranchAddress("T4H", T4H);
  new_pps_samples_Tree->SetBranchAddress("T1V", T1V);
  new_pps_samples_Tree->SetBranchAddress("T2V", T2V);
  new_pps_samples_Tree->SetBranchAddress("T3V", T3V);
  new_pps_samples_Tree->SetBranchAddress("T4V", T4V);
	


	//for(int entry=0; entry<new_pps_samples_Tree->GetEntries(); entry++){
		int entry = 42776;
		new_pps_samples_Tree->GetEntry(entry);
	
		double x[8][4096] = {0.};
		double t[1500] = {0.};
		for(int i=0; i<1500; i++){
			t[i] = i*0.8;
			x[0][i] = T1H[i] *500./32512.;
			x[2][i] = T2H[i] *500./32512.;
			x[4][i] = T3H[i] *500./32512.;
			x[6][i] = T4H[i] *500./32512.;
			x[1][i] = T1V[i] *500./32512.;
			x[3][i] = T2V[i] *500./32512.;
			x[5][i] = T3V[i] *500./32512.;
			x[7][i] = T4V[i] *500./32512.;
		}
		
		TGraph *cyc = new TGraph(1500, t, x[1]);
		cyc->Draw();


   int N = 1500;
   TVirtualFFT *fft_forward = TVirtualFFT::FFT(1, &N,"R2C");
      
   TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &N, "C2RBACKWARD M K");


	double Re_ant[751];
	double Im_ant[751];

   fft_forward->SetPoints(x[1]);
   fft_forward->Transform();
   fft_forward->GetPointsComplex(Re_ant,Im_ant);
   
   double power[751] = {0.};
   for(int i=0; i<751; i++){
   		power[i] = 10*log10((pow(Re_ant[i], 2) + pow(Im_ant[i], 2)));
   	}
   	//TGraph *wyn = new TGraph(751, t, power);
	//wyn->Draw();
   	
   
   
   		
		
	//}

    
	
}





