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
#include "utils.cc"
#include "trigger_sim.cc"


#include <stdio.h>
#include <dirent.h>
#include <sys/types.h>

using namespace std;
#define PI TMath::Pi()
#define Speed_Of_Light 2.99792458e8

#define HALF_WINDOW_SPACING 75


void GetTChain(string directory, TChain *chaintree); 
 

int filter_time(double x[1500]);
void filter_fre(double x[1500], bool save_100MHz);
bool CreateDataTreeBranches(TTree* tr);
double loop_make_data(int run);
 
 Double_t new_x[8][1500]={0};//8 antennas
 Double_t pps_x[8][1500]={0};//8 antennas
 Double_t wCoreas_x[8][1500]={0};//8 antennas
 
 

  Short_t pps_T1H[1500];
  Short_t pps_T2H[1500];
  Short_t pps_T3H[1500];
  Short_t pps_T4H[1500];
  Short_t pps_T1V[1500];
  Short_t pps_T2V[1500];
  Short_t pps_T3V[1500];
  Short_t pps_T4V[1500];

  Short_t coreas_T1H[1500];
  Short_t coreas_T2H[1500];
  Short_t coreas_T3H[1500];
  Short_t coreas_T4H[1500];
  Short_t coreas_T1V[1500];
  Short_t coreas_T2V[1500];
  Short_t coreas_T3V[1500];
  Short_t coreas_T4V[1500];
  TVector3 *vec_pos = NULL;
  double radius = 0;
  double pattern_phi = 0;
  double ShowerTheta;
  double ShowerPhi;
  struct coreas_sample {
    double radius=-1.;
    double phi=-1.;
    int entry=-1;
  };


  Short_t new_T1H[1500];
  Short_t new_T2H[1500];
  Short_t new_T3H[1500];
  Short_t new_T4H[1500];
  Short_t new_T1V[1500];
  Short_t new_T2V[1500];
  Short_t new_T3V[1500];
  Short_t new_T4V[1500];
    
  TBits gTrigBits(64);
  TBits * gOverVoltBits = NULL;
  //TBits gTrigBits(64);
  //TBits gOverVoltBits(8);
  
  int gEventNum;
  

  bool bGood;
  double eventTime = 0;
  Int_t timeStamp_FPGA;

//TRandom *r3;

//.x make_coreas_data.cc(1, "", 20210408, 18969, 18970)
void make_coreas_data(){

	trigger_sim_initial();

	initial_coreas_coordinate();
	r3 = new TRandom3();
	r3->SetSeed(851);
	
	TFile *f_ntuple = new TFile("ntuple1.root","RECREATE");
	TNtuple *ntuple = new TNtuple("ntuple","t","Energy:ShowerTheta:ShowerPhi:EffectedArea");
	
	/*
	for(int run = 60000; run<=60019; run++){
		if(run%1000>57){continue;}
		if(run%1000<1){continue;}
		double EffectedArea = loop_make_data(run);
		double Energy = run/10000-2+17.5;
		cout<<run<<"  "<<EffectedArea<<endl;
		ntuple->Fill(Energy, ShowerTheta, ShowerPhi, EffectedArea);
	}
	*/
	
	
	for(int run = 56006; run<=56012; run+=2){
		if(run%1000>19){continue;}
		if(run%1000<1){continue;}
		int true_run = run;
		if(run/1000==30){true_run-=20000;}
		double EffectedArea = loop_make_data(true_run);
		double Energy = run/10000-2+17.5;
		cout<<run<<"  "<<EffectedArea<<endl;
		ntuple->Fill(Energy, ShowerTheta, ShowerPhi, EffectedArea);
	}
	
	//double EffectedArea = loop_make_data(31009);
	
	
	
	//ntuple->Close();
	f_ntuple->Write(0,TObject::kOverwrite);
	f_ntuple->Close();

}

double loop_make_data(int run){	
	//int run = 17015;
	const string coreas_root = Form("/media/cyc/For_Linux/CoREAS_Sim/proton/cyc/Coreas-t4-r%06d.root", run);
	TFile *coreas_file = new TFile(coreas_root.c_str());
	
	TTree *EfieldTree = (TTree*) coreas_file->Get("EfieldTree");
	EfieldTree->SetBranchAddress("radius", &radius);
	EfieldTree->SetBranchAddress("phi", &pattern_phi);
	EfieldTree->SetBranchAddress("ShowerTheta", &ShowerTheta);
	EfieldTree->SetBranchAddress("ShowerPhi", &ShowerPhi);
	
	EfieldTree->SetBranchAddress("vecPos", &vec_pos);
	const int coreas_entries = EfieldTree->GetEntries();
	
	// maximum 120 sample points, 8 angle * 15 radius
	coreas_sample E_sample[200];
	double max_radius=0.;
	
	for(int i=0; i<coreas_entries; i++){
		EfieldTree->GetEntry(i);
		E_sample[i].entry = i;
		E_sample[i].radius = radius;
		E_sample[i].phi = pattern_phi;
		if(radius>max_radius){max_radius=radius;}
	}
	

	EfieldTree->SetBranchAddress("T1V", coreas_T1V);
	EfieldTree->SetBranchAddress("T2V", coreas_T2V);
	EfieldTree->SetBranchAddress("T3V", coreas_T3V);
	EfieldTree->SetBranchAddress("T4V", coreas_T4V);
	EfieldTree->SetBranchAddress("T1H", coreas_T1H);
	EfieldTree->SetBranchAddress("T2H", coreas_T2H);
	EfieldTree->SetBranchAddress("T3H", coreas_T3H);
	EfieldTree->SetBranchAddress("T4H", coreas_T4H);
	
	double impact_point_x[2][12000] = {0.};
	double impact_point_y[2][12000] = {0.};
	
	int trigger_count = 0;
	for(int sp=0; sp<1200; sp++){
		double impact_radius = sqrt(r3->Rndm()) * max_radius;
		double impact_phi = r3->Rndm() * 360;
		
		double d_square = 1.0e13;
		int s_entry = 0;
		for(int i=0; i<coreas_entries; i++){
			double ds = impact_radius*impact_radius + E_sample[i].radius*E_sample[i].radius - 2*impact_radius*E_sample[i].radius*TMath::Cos((impact_phi-E_sample[i].phi)*TMath::DegToRad());
			
			if(ds<d_square){
				d_square = ds;
				s_entry = E_sample[i].entry;
			}
		}
			//cout<<s_entry<<endl;
		EfieldTree->GetEntry(s_entry);
		
		
		
		double x_vpp[8][4096] = {0.};
		for(int i=0; i<1500; i++){
			x_vpp[0][i] = coreas_T1H[i] *500./32512;
			x_vpp[2][i] = coreas_T2H[i] *500./32512;
			x_vpp[4][i] = coreas_T3H[i] *500./32512;
			x_vpp[6][i] = coreas_T4H[i] *500./32512;
			x_vpp[1][i] = coreas_T1V[i] *500./32512;
			x_vpp[3][i] = coreas_T2V[i] *500./32512;
			x_vpp[5][i] = coreas_T3V[i] *500./32512;
			x_vpp[7][i] = coreas_T4V[i] *500./32512;
		}
		
		bool need_trigger = false;
		for(int i=0; i<1500; i++){
			if((x_vpp[0][i]>10)||(x_vpp[1][i]>10)){
				need_trigger = true;
				break;
			}
		}
		
		//double max_T1H = *std::max_element(&x_vpp[0][0], &x_vpp[0][0]+1500);
		//double max_T1V = *std::max_element(&x_vpp[1][0], &x_vpp[1][0]+1500);
		
		bool is_trigger = false;
		if(need_trigger){
			is_trigger = IsTrigger(x_vpp, 650, gTrigBits);
		}
		//cout<<s_entry<<endl;
		double detected_x = impact_radius*TMath::Cos(impact_phi*TMath::DegToRad());
		double detected_y = impact_radius*TMath::Sin(impact_phi*TMath::DegToRad());
		//cout<<detected_x<<endl;
		detected_x /= TMath::Cos((90-ShowerTheta)*TMath::DegToRad());
		//cout<<ShowerTheta<<endl;
		//cout<<detected_x<<endl;
		//double detected_x = vec_pos->X();
		//double detected_y = vec_pos->Y();
		impact_point_x[0][sp] = detected_x;
		impact_point_y[0][sp] = detected_y;
		if(is_trigger){
			//cout<<x_vpp[0]<<"  "<<x_vpp[1]<<endl;
			impact_point_x[1][trigger_count] = detected_x;
			impact_point_y[1][trigger_count] = detected_y;
			trigger_count++;
			
			for(int i=0; i<8; i++){
			//cout<<impact_radius<<"   "<<i<<": "<<*std::max_element(&x_vpp[i][0], &x_vpp[i][0]+1500)<<endl;
			}
		}
		
		
	}
	
	
	
	
	TGraph * cyc = new TGraph(1200, impact_point_x[0], impact_point_y[0]);
	cyc->Draw("AP");
	cyc->SetMarkerColor(1);
	cyc->SetMarkerStyle(7);
	cyc->SetMarkerSize(2);
   
	TGraph * wyn = new TGraph(trigger_count, impact_point_x[1], impact_point_y[1]);
	wyn->Draw("P same");
	wyn->SetMarkerColor(2);
	wyn->SetMarkerStyle(7);
	wyn->SetMarkerSize(2);
	
	
	
	delete EfieldTree;
	coreas_file->Close();
	delete coreas_file;
	cout<<trigger_count<<"   "<<max_radius<<endl;
	
	return 3.14159*max_radius*max_radius*trigger_count/1200.;


}

/*
	TBits* triggerBits = 0;
	double eventTime = 0;

	TChain *pps_chain = new TChain("t");
	string pps_directory = "/media/cyc/For_Linux/TAROGE4_DATA/data/20201218/";
	GetTChain(pps_directory, pps_chain);
	pps_chain->SetBranchAddress("eventTime", &eventTime);
	pps_chain->SetBranchAddress("triggerBits",&triggerBits);
    const int pps_entries = pps_chain->GetEntries();


	int count_pps=0;
	double pps_eventTime_list[pps_entries];
	double entry_list[pps_entries];
	for(int i=0; i<pps_entries; i++){
		pps_chain->GetEntry(i);
		if(triggerBits->CountBits()==0){
			entry_list[count_pps] = i;
			pps_eventTime_list[count_pps] = eventTime;
			count_pps++;
		}
	}


	double eventRate = 0;
	double threshold = 0;
	
	TChain *trigger_chain = new TChain("t");
	string trigger_directory = "/media/cyc/For_Linux/TAROGE4_DATA/trigger/20201218/";
	GetTChain(trigger_directory, trigger_chain);
	trigger_chain->SetBranchAddress("eventTime", &eventTime);
	trigger_chain->SetBranchAddress("threshold",&threshold);
	trigger_chain->SetBranchAddress("eventRate",&eventRate);
    const int trigger_entries = trigger_chain->GetEntries();

	int count_trigger=0;
	double trigger_eventTime_list[trigger_entries];
	double threshold_list[trigger_entries];
	double eventrate_list[trigger_entries];
	for(int i=0; i<trigger_entries; i++){
		trigger_chain->GetEntry(i);
		trigger_eventTime_list[count_trigger] = eventTime;
		threshold_list[count_trigger] = threshold;
		eventrate_list[count_trigger] = eventRate;
		count_trigger++;
	}
	    

	
}
*/



/*
void fuck(int pol, string pulser_data_name, int date, int start_run, int end_run, int Energy){
 string out_put_root_name = "/media/cyc/For_Linux/CoREAS_Sim/proton/wyn/make_coreas_data_" + std::to_string(Energy*100) + ".root";
 TFile *make_pulser_data = new TFile(out_put_root_name.c_str(), "RECREATE");

 TTree *make_Tree = new TTree("t","TAROGE data");
 
 CreateDataTreeBranches(make_Tree);



  Int_t event_number;
  const char * fname_evt;
  TFile *file;
  TTree *Tree_Muon;


  for(int number = start_run; number<end_run+1; number++){


    //fname_evt =  Form( "/media/cyc/For_Linux/TAROGE4_DATA/data/20211209/run%08d.root",number);
    fname_evt =  Form( "/media/cyc/1p9TB/TAROGE4_DATA/data/%d/run%08d.root", date, number);
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
      //cout<<j<<endl;

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
    
      //cout<<gTrigBits->CountBits()<<endl;
      gTrigBits->SetBitNumber(32, kTRUE);//let it become not all zeros
      //cout<<gTrigBits->CountBits()<<endl;


    for (int i=0;i<1500;i++)
	{
	  pps_x[0][i]= (double)pps_T1H[i]*500./32512;;
	  pps_x[2][i]= (double)pps_T2H[i]*500./32512;;
	  pps_x[4][i]= (double)pps_T3H[i]*500./32512;;
	  pps_x[6][i]= (double)pps_T4H[i]*500./32512;;
	  pps_x[1][i]= (double)pps_T1V[i]*500./32512;;
	  pps_x[3][i]= (double)pps_T2V[i]*500./32512;;
	  pps_x[5][i]= (double)pps_T3V[i]*500./32512;;
	  pps_x[7][i]= (double)pps_T4V[i]*500./32512;;
	  //cout<<pps_x[0][i]<<endl;	  
	}

	for(int ch=0; ch<8; ch++){
		filter_fre(pps_x[ch], 1);
	}
	
 for(int count=0; count<count_sample; count++){
	for(int ch=0; ch<8; ch++){

		for(int k=0; k<1500; k++){
			new_x[ch][k] = pps_x[ch][k];
		}

	}

	for(int ch=0; ch<8; ch++){
		int ant = ch/2;
		double ant_4_factor = 1.;
		// FEE response of all channel already taken it into account
		//if(ant==3){ant_4_factor = 1./sqrt(2.);}
		for(int k=0; k<1500; k++){
			new_x[ch][k] += ant_4_factor*E_waveform[count][ch][k];
		}
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
	
	}// count_sample loop
  }//event loop
  }//run loop

   make_pulser_data->cd();
   make_Tree->Write();
   make_pulser_data->Close();
   delete make_pulser_data;
	
   
}//end of main function
*/


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


/*

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


}

*/

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

