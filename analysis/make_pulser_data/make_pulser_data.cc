#include "TTree.h"
#include <iostream>
#include "TGraph.h" 
#include "TVirtualFFT.h" 
#include "TFile.h"
#include <fstream>
#include <math.h>
#include "TMath.h"
#include "TVirtualFFT.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TGraph.h"
#include <iostream>
#include "TTree.h"
#include "TNtuple.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TChain.h"
#include <TBits.h>
#include <ctime>
using namespace std;
#define PI TMath::Pi()
#define Speed_Of_Light 2.99792458e8

#define HALF_WINDOW_SPACING 75
 
 

int filter_time(double x[1500]);
void filter_fre(double x[1500], bool save_100MHz);
bool CreateDataTreeBranches(TTree* tr);

 
 Double_t new_x[8][1500]={0};//8 antennas
 Double_t pps_x[8][1500]={0};//8 antennas
 Double_t pulser_x[8][1500]={0};//8 antennas
 
 

  Short_t pps_T1H[1500];
  Short_t pps_T2H[1500];
  Short_t pps_T3H[1500];
  Short_t pps_T4H[1500];
  Short_t pps_T1V[1500];
  Short_t pps_T2V[1500];
  Short_t pps_T3V[1500];
  Short_t pps_T4V[1500];

  Short_t pulser_T1H[1500];
  Short_t pulser_T2H[1500];
  Short_t pulser_T3H[1500];
  Short_t pulser_T4H[1500];
  Short_t pulser_T1V[1500];
  Short_t pulser_T2V[1500];
  Short_t pulser_T3V[1500];
  Short_t pulser_T4V[1500];
  
  Short_t new_T1H[1500];
  Short_t new_T2H[1500];
  Short_t new_T3H[1500];
  Short_t new_T4H[1500];
  Short_t new_T1V[1500];
  Short_t new_T2V[1500];
  Short_t new_T3V[1500];
  Short_t new_T4V[1500];
    
  TBits * gTrigBits = NULL;
  TBits * gOverVoltBits = NULL;
  //TBits gTrigBits(64);
  //TBits gOverVoltBits(8);
  
  int gEventNum;
  

  bool bGood;
  double eventTime = 0;
  Int_t timeStamp_FPGA;


//.x make_pulser_data.cc(0, "pulser_data_H_H_V.root", 20210408, 18969, 18970)
void make_pulser_data(int pol, string pulser_data_name, int date, int start_run, int end_run){

//cout<<"working here"<<endl;
//cout<<argv[0]<<endl;
//cout<<argv[1]<<endl;
//cout<<argv[2]<<endl;
//cout<<argv[3]<<endl;
//cout<<argv[4]<<endl;
//cout<<argv[5]<<endl;

//int pol = std::stoi(argv[1]); 
//string pulser_data_name = argv[2];
//int date = std::stoi(argv[3]);
//int start_run = std::stoi(argv[4]);
//int end_run = std::stoi(argv[5]);




    TFile * pulser_file = new TFile(pulser_data_name.c_str());
    TTree *Tree_pulser = (TTree*) pulser_file->Get("t");

	Tree_pulser->SetBranchAddress("T1H",pulser_T1H);
	Tree_pulser->SetBranchAddress("T2H",pulser_T2H);
	Tree_pulser->SetBranchAddress("T3H",pulser_T3H);
	Tree_pulser->SetBranchAddress("T4H",pulser_T4H);
	Tree_pulser->SetBranchAddress("T1V",pulser_T1V);
	Tree_pulser->SetBranchAddress("T2V",pulser_T2V);
	Tree_pulser->SetBranchAddress("T3V",pulser_T3V);
	Tree_pulser->SetBranchAddress("T4V",pulser_T4V);
	
	int pulser_number = Tree_pulser->GetEntries();
	int pulser_count = 0;
		  
 
 
 string out_put_root_name = "run20240102.root";//"make_pulser_data_" + std::to_string(pol) + ".root";
 TFile *make_pulser_data = new TFile(out_put_root_name.c_str(), "RECREATE");

 TTree *make_Tree = new TTree("t","TAROGE data");
 
 CreateDataTreeBranches(make_Tree);



  Int_t event_number;
  const char * fname_evt;
  TFile *file;
  TTree *Tree_Muon;


 


  for(int number = start_run; number<end_run+1; number++){


    //fname_evt =  Form( "/media/cyc/For_Linux/TAROGE4_DATA/data/20211209/run%08d.root",number);
    fname_evt =  Form( "/media/cyc/For_Linux/TAROGE4_DATA/data/%d/run%08d.root", date, number);
    file = new TFile(fname_evt);
    Tree_Muon = (TTree*) file->Get("t");
          Tree_Muon->SetBranchAddress("bGood",&bGood);
		  Tree_Muon->SetBranchAddress("timeStamp_FPGA",&timeStamp_FPGA);
		  Tree_Muon->SetBranchAddress("eventTime", &eventTime);
  		  Tree_Muon->SetBranchAddress("triggerBits", &gTrigBits);
  		  Tree_Muon->SetBranchAddress("overVoltBits", &gOverVoltBits);
		  Tree_Muon->SetBranchAddress("T1H",pps_T1H);
		  Tree_Muon->SetBranchAddress("T2H",pps_T2H);
		  Tree_Muon->SetBranchAddress("T3H",pps_T3H);
		  Tree_Muon->SetBranchAddress("T4H",pps_T4H);
		  Tree_Muon->SetBranchAddress("T1V",pps_T1V);
		  Tree_Muon->SetBranchAddress("T2V",pps_T2V);
		  Tree_Muon->SetBranchAddress("T3V",pps_T3V);
		  Tree_Muon->SetBranchAddress("T4V",pps_T4V);
		  
 event_number = Tree_Muon->GetEntries();

 for (int j=0;j<event_number;j++){

      Tree_Muon->GetEntry(j);
      if(gTrigBits->CountBits() != 0) {continue;}
      cout<<j<<endl;
      
   //also store the pps original data
	for (int i=0;i<1500;i++)
	{
		new_T1H[i] = pps_T1H[i];
		new_T2H[i] = pps_T2H[i];
		new_T3H[i] = pps_T3H[i];
		new_T4H[i] = pps_T4H[i];
		new_T1V[i] = pps_T1V[i];
		new_T2V[i] = pps_T2V[i];
		new_T3V[i] = pps_T3V[i];
		new_T4V[i] = pps_T4V[i];
	  	  
	}
	make_Tree->Fill();
	
	
    // store made pulser data 
    
      cout<<gTrigBits->CountBits()<<endl;
      //H-pol 0; V-pol 32
      gTrigBits->SetBitNumber(32*pol, kTRUE);//let it become not all zeros
      cout<<gTrigBits->CountBits()<<endl;
      
      Tree_pulser->GetEntry(pulser_count);
      pulser_count++;
      
    for (int i=0;i<1500;i++)
	{
	  pulser_x[0][i]= (double)pulser_T1H[i]*500/32512;
	  pulser_x[2][i]= (double)pulser_T2H[i]*500/32512;
	  pulser_x[4][i]= (double)pulser_T3H[i]*500/32512;
	  pulser_x[6][i]= (double)pulser_T4H[i]*500/32512;
	  pulser_x[1][i]= (double)pulser_T1V[i]*500/32512;
	  pulser_x[3][i]= (double)pulser_T2V[i]*500/32512;
	  pulser_x[5][i]= (double)pulser_T3V[i]*500/32512;
	  pulser_x[7][i]= (double)pulser_T4V[i]*500/32512;	  	  
	}

    for (int i=0;i<1500;i++)
	{
	  pps_x[0][i]= (double)pps_T1H[i]*500/32512;
	  pps_x[2][i]= (double)pps_T2H[i]*500/32512;
	  pps_x[4][i]= (double)pps_T3H[i]*500/32512;
	  pps_x[6][i]= (double)pps_T4H[i]*500/32512;
	  pps_x[1][i]= (double)pps_T1V[i]*500/32512;
	  pps_x[3][i]= (double)pps_T2V[i]*500/32512;
	  pps_x[5][i]= (double)pps_T3V[i]*500/32512;
	  pps_x[7][i]= (double)pps_T4V[i]*500/32512;
	  //cout<<pps_x[0][i]<<endl;	  
	}

	for(int ch=0; ch<8; ch++){
		filter_fre(pps_x[ch], 1);
		filter_fre(pulser_x[ch], 0);
		//cout<<max_i<<endl;
		}
			
	for(int ant=0; ant<4; ant++){
		int ch = 2*ant + pol;
		int o_ch = 2*ant + (pol+1)%2;
		
		int max_i = filter_time(pulser_x[ch]);
		

		for(int k=0; k<1500; k++){
			//1?3:4+3 is 7
			new_x[ch][k] = ((k<max_i-HALF_WINDOW_SPACING)||(k>max_i+HALF_WINDOW_SPACING))?pps_x[ch][k]:(pulser_x[ch][k]/5.6+pps_x[ch][k]);
			new_x[o_ch][k] = ((k<max_i-HALF_WINDOW_SPACING)||(k>max_i+HALF_WINDOW_SPACING))?pps_x[o_ch][k]:(pulser_x[o_ch][k]/5.6+pps_x[o_ch][k]);
			// / 5.6
			//if(fabs(new_x[ch][k])>100){cout<<"k: "<<k<<endl;}
		}
		
		//new_x[ch][500] += -100;
		//new_x[ch][503] += 200;
		//new_x[ch][506] += -300;
		//new_x[ch][509] += 200;
		//new_x[ch][512] += -100;
		//new_x[o_ch][500] += -100;
		//new_x[o_ch][503] += 200;
		//new_x[o_ch][506] += -300;
		//new_x[o_ch][509] += 200;
		//new_x[o_ch][512] += -100;
		
	}
	
	for (int i=0;i<1500;i++)
	{
		new_T1H[i] = (Short_t) round(new_x[0][i]*32512/500);
		new_T2H[i] = (Short_t) round(new_x[2][i]*32512/500);
		new_T3H[i] = (Short_t) round(new_x[4][i]*32512/500);
		new_T4H[i] = (Short_t) round(new_x[6][i]*32512/500);
		new_T1V[i] = (Short_t) round(new_x[1][i]*32512/500);
		new_T2V[i] = (Short_t) round(new_x[3][i]*32512/500);
		new_T3V[i] = (Short_t) round(new_x[5][i]*32512/500);
		new_T4V[i] = (Short_t) round(new_x[7][i]*32512/500);
	  	  
	}
	make_Tree->Fill();
  }//event loop
  }//run loop

   make_pulser_data->cd();
   make_Tree->Write();
   make_pulser_data->Close();
   delete make_pulser_data;

   
}//end of main function



int filter_time(double x[1500]){

int Max_i = 0;
double Max = 0.;

           for (int i=0;i<1500;i++)
	{

	  if(fabs(x[i])>Max) {Max = fabs(x[i]);  Max_i = i;} 	  

	}
	//cout<<Max_i<<endl;
	
           for (int i=0;i<1500;i++)
	{

   	       x[i] = ((i<Max_i-100)||(i>Max_i+100))?0:x[i];
 	                               
	}

return Max_i;

}




void filter_fre(double x[1500], bool save_100MHz){

int N = 1500;

   TVirtualFFT *fft_forward = TVirtualFFT::FFT(1, &N,"R2C");
      
   TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &N, "C2RBACKWARD M K");
   
   	
double Re_ant[751];
double Im_ant[751];



   fft_forward->SetPoints(x);
   fft_forward->Transform();
   fft_forward->GetPointsComplex(Re_ant,Im_ant);
   
   // remove DC part
   Re_ant[0] = 0.; Im_ant[0] = 0;

	if(!save_100MHz){
 		for(Int_t k=119; k<122; k++){//remove 100MHz
       		Re_ant[k] = 0.;
       		Im_ant[k] = 0;
                             }
                    }
  
          

  


Double_t * filtered_x = new Double_t [1500];

  fft_back->SetPointsComplex(Re_ant,Im_ant);
  fft_back->Transform();
  filtered_x = fft_back->GetPointsReal();

for(int i=0;i<1500;i++){
x[i] = filtered_x[i]/1500.;
  						}

delete fft_forward;
delete fft_back;
}




bool CreateDataTreeBranches(TTree* tr){

  printf("CreateDataTreeBranches(): %s\t%s\n", tr->GetName(), tr->GetTitle() );
  tr->Branch("event",&gEventNum,"event/I");
  tr->Branch("bGood",&bGood,"bGood/O");     
  tr->Branch("timeStamp_FPGA", &timeStamp_FPGA,"timeStamp_FPGA/I"); //TAROGE-4, FPGA 24-bit counter    
  tr->Branch("triggerBits",&gTrigBits ); //TAROGE-4  TBits, align with ChIDs    Hpol and Vpol
  tr->Branch("overVoltBits",&gOverVoltBits ); //TAROGE-4  TBits(kNCh)   overvoltage or not
  tr->Branch("eventTime",&eventTime,"eventTime/D"); //PC timestamp

  tr->Branch("T1H", &new_T1H, "T1H[1500]/S");
  tr->Branch("T2H", &new_T2H, "T2H[1500]/S");
  tr->Branch("T3H", &new_T3H, "T3H[1500]/S");
  tr->Branch("T4H", &new_T4H, "T4H[1500]/S");
  tr->Branch("T1V", &new_T1V, "T1V[1500]/S");
  tr->Branch("T2V", &new_T2V, "T2V[1500]/S");
  tr->Branch("T3V", &new_T3V, "T3V[1500]/S");
  tr->Branch("T4V", &new_T4V, "T4V[1500]/S");
  // /S: 2-byte signed short  /I: 4 byte signed int  /F double
 

 
  return true;
}

