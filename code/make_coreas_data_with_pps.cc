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
#include <string.h>



#include <stdio.h>
#include <dirent.h>
#include <sys/types.h>

using namespace std;
#define PI TMath::Pi()
#define Speed_Of_Light 2.99792458e8

#define HALF_WINDOW_SPACING 75


void GetTChain(string directory, TChain *chaintree); 
void pps_initial(TTree *pps_samples_Tree);
bool CreateDataTreeBranches(TTree* tr);
 


  Short_t triggered_T1H[1500];
  Short_t triggered_T2H[1500];
  Short_t triggered_T3H[1500];
  Short_t triggered_T4H[1500];
  Short_t triggered_T1V[1500];
  Short_t triggered_T2V[1500];
  Short_t triggered_T3V[1500];
  Short_t triggered_T4V[1500];
  
  TBits gTrigBits(64);
  TBits gOverVoltBits(8);
  bool bGood = true;
  double Energy = 0.;
    
    
    
    
  Short_t pps_T1H[1500];
  Short_t pps_T2H[1500];
  Short_t pps_T3H[1500];
  Short_t pps_T4H[1500];
  Short_t pps_T1V[1500];
  Short_t pps_T2V[1500];
  Short_t pps_T3V[1500];
  Short_t pps_T4V[1500];
  double eventTime = 0;
  Int_t gDate = 0;
  Int_t gEventNum = 0;
  Double_t gEventRate = 0.;
  Int_t gThresholdStep = 0;
  
  
  double radius = 0;
  double pattern_phi = 0;
  double ShowerTheta;
  double ShowerPhi;
  Int_t coreas_T1H[1500];
  Int_t coreas_T2H[1500];
  Int_t coreas_T3H[1500];
  Int_t coreas_T4H[1500];
  Int_t coreas_T1V[1500];
  Int_t coreas_T2V[1500];
  Int_t coreas_T3V[1500];
  Int_t coreas_T4V[1500];
  

double loop_make_data(string coreas_root, int pps_entries, TTree *pps_samples_Tree, double scale_factor, TTree *make_tree, double &DSR_radius, int &Triggered_count, int &Total_count);
 
void make_coreas_data_with_pps(int start_run, int final_run, int cr_energy);
 
 
 Double_t new_x[8][1500]={0};//8 antennas
 Double_t pps_x[8][1500]={0};//8 antennas
 Double_t wCoreas_x[8][1500]={0};//8 antennas



  TVector3 *vec_pos = NULL;
  
  struct coreas_sample {
    double radius=-1.;
    double phi=-1.;
    int entry=-1;
    double max_E[2]={0.};// H V
  };


  Short_t new_T1H[1500];
  Short_t new_T2H[1500];
  Short_t new_T3H[1500];
  Short_t new_T4H[1500];
  Short_t new_T1V[1500];
  Short_t new_T2V[1500];
  Short_t new_T3V[1500];
  Short_t new_T4V[1500];
    
int main(int argc, char *argv[]){

	

	int start_run = std::stoi(argv[1]);
	int final_run = std::stoi(argv[2]);
	int cr_energy = std::stoi(argv[3]);
	make_coreas_data_with_pps(start_run, final_run, cr_energy);
	return 1;
}

//TRandom *r3;

//.x make_coreas_data.cc(1, "", 20210408, 18969, 18970)
void make_coreas_data_with_pps(int start_run, int final_run, int cr_energy){

	TFile *pps_samples = new TFile("pps_samples_filtered.root");
 	TTree *pps_samples_Tree = (TTree *) pps_samples->Get("t");

	pps_initial(pps_samples_Tree);
	int pps_entries = pps_samples_Tree->GetEntries();
	
	gOverVoltBits.ResetAllBits(false);


	trigger_sim_initial();

	initial_coreas_coordinate();
	r3 = new TRandom3();
	r3->SetSeed(final_run);
	
	string out_put_ntuple_name = "ntuple_position_factor_with_pps_iron_" + std::to_string(cr_energy) + ".root";
	TFile *f_ntuple = new TFile(out_put_ntuple_name.c_str(),"RECREATE");
	TNtuple *ntuple = new TNtuple("ntuple","t","Energy:ShowerTheta:ShowerPhi:EffectedArea:DSR_radius:Triggered_count:Total_count");


  
  	string out_put_root_name = "make_coreas_data_position_factor_with_pps_iron_" + std::to_string(cr_energy) + ".root";
 	TFile *make_coreas_data = new TFile(out_put_root_name.c_str(), "RECREATE");
 	TTree *make_Tree = new TTree("t","TAROGE data");
 	CreateDataTreeBranches(make_Tree);
 	
	for(int run = start_run; run<=final_run; run++){
		
		const string coreas_root = Form("/mnt/network_hdd/process_coreas_data/iron/Coreas-t4-r%06d.root", run);
		// 0 means exist, -1 means not
		if(access(coreas_root.c_str(), F_OK)!=0)
		{
			cout<<coreas_root<<"  do not exist! skip it"<<endl;
			continue;
		}
		double scale_factor = 1.;
		double Energy = cr_energy;
		
		//this is for 89000 89.7 degrees showers scaling
		//double scale_factor = pow(10, (cr_energy/1000.-18.75));
		//double Energy = cr_energy/10.;
		
		double DSR_radius;
		int Triggered_count, Total_count;
		double EffectedArea = loop_make_data(coreas_root, pps_entries, pps_samples_Tree, scale_factor, make_Tree, DSR_radius, Triggered_count, Total_count);
		cout<<run<<"  "<<Energy<<"   "<<EffectedArea<<endl;
		//DSR_radius:Triggered_count:Total_count
		ntuple->Fill(Energy, ShowerTheta, ShowerPhi, EffectedArea, DSR_radius, Triggered_count, Total_count);
	}
	make_coreas_data->cd();
	make_Tree->Write(0,TObject::kOverwrite);
	make_coreas_data->Close();
	delete make_coreas_data;

	
	
 
	//double EffectedArea = loop_make_data(31009);
	
	
	
	f_ntuple->cd();
	ntuple->Write(0,TObject::kOverwrite);
	f_ntuple->Close();

}

double loop_make_data(string coreas_root, int pps_entries, TTree *pps_samples_Tree, double scale_factor, TTree *make_tree, double &DSR_radius, int &Triggered_count, int &Total_count){	
	//int run = 17015;
	
	TFile *coreas_file = new TFile(coreas_root.c_str());
	
	TTree *EfieldTree = (TTree*) coreas_file->Get("EfieldTree");
	EfieldTree->SetBranchAddress("radius", &radius);
	EfieldTree->SetBranchAddress("phi", &pattern_phi);
	EfieldTree->SetBranchAddress("ShowerTheta", &ShowerTheta);
	EfieldTree->SetBranchAddress("ShowerPhi", &ShowerPhi);
	EfieldTree->SetBranchAddress("T1V", coreas_T1V);
	EfieldTree->SetBranchAddress("T2V", coreas_T2V);
	EfieldTree->SetBranchAddress("T3V", coreas_T3V);
	EfieldTree->SetBranchAddress("T4V", coreas_T4V);
	EfieldTree->SetBranchAddress("T1H", coreas_T1H);
	EfieldTree->SetBranchAddress("T2H", coreas_T2H);
	EfieldTree->SetBranchAddress("T3H", coreas_T3H);
	EfieldTree->SetBranchAddress("T4H", coreas_T4H);
	
	//EfieldTree->SetBranchAddress("vecPos", &vec_pos);
	const int coreas_entries = EfieldTree->GetEntries();
	
	// maximum 120 sample points, 8 angle * 15 radius
	coreas_sample E_sample[200];
	double max_radius=0.;
	
	for(int i=0; i<coreas_entries; i++){
		EfieldTree->GetEntry(i);
		E_sample[i].entry = i;
		E_sample[i].radius = radius;
		E_sample[i].phi = pattern_phi;
		E_sample[i].max_E[0] = (TMath::MaxElement(1500, coreas_T1H) + TMath::MaxElement(1500, coreas_T2H) + TMath::MaxElement(1500, coreas_T3H) + sqrt(2.)*TMath::MaxElement(1500, coreas_T4H)) + fabs((TMath::MinElement(1500, coreas_T1H) + TMath::MinElement(1500, coreas_T2H) + TMath::MinElement(1500, coreas_T3H) + sqrt(2.)*TMath::MinElement(1500, coreas_T4H)));
		
		E_sample[i].max_E[1] = (TMath::MaxElement(1500, coreas_T1V) + TMath::MaxElement(1500, coreas_T2V) + TMath::MaxElement(1500, coreas_T3V) + sqrt(2.)*TMath::MaxElement(1500, coreas_T4V)) + fabs((TMath::MinElement(1500, coreas_T1V) + TMath::MinElement(1500, coreas_T2V) + TMath::MinElement(1500, coreas_T3V) + sqrt(2.)*TMath::MinElement(1500, coreas_T4V)));
		if(radius>max_radius){max_radius=radius;}
	}
	

	

	
	int trigger_count = 0;
	Total_count = 1000;
	for(int sp=0; sp<Total_count; sp++){
		double impact_radius = sqrt(r3->Rndm()) * max_radius;
		double impact_phi = r3->Rndm() * 360;;
		
		double d_square = 1.0e13;
		int s_entry=0;
		int sample_entry[4] = {0};
		for(int i=0; i<coreas_entries; i++){
			double ds = impact_radius*impact_radius + E_sample[i].radius*E_sample[i].radius - 2*impact_radius*E_sample[i].radius*TMath::Cos((impact_phi-E_sample[i].phi)*TMath::DegToRad());
			
			if(ds<d_square){
				d_square = ds;
				sample_entry[0] = i;
				s_entry = E_sample[i].entry;
				}
		}
		
		
		if(impact_radius-E_sample[sample_entry[0]].radius>0){sample_entry[1] = sample_entry[0]+8;}
		else{sample_entry[1] = sample_entry[0]-8;}
		if((impact_phi-E_sample[sample_entry[0]].phi>0)&&(impact_phi-E_sample[sample_entry[0]].phi<46)){
			sample_entry[2] = sample_entry[0]+1;
			sample_entry[3] = sample_entry[1]+1;
		}
		else{
			sample_entry[2] = sample_entry[0]-1;
			sample_entry[3] = sample_entry[1]-1;
		}
		
		double radius_factor[2] = {0.};
		double phi_factor[2] = {0.};
		
		double radius_factor_length = fabs(E_sample[sample_entry[0]].radius-E_sample[sample_entry[1]].radius)+1.0e-12;
		radius_factor[0] = fabs(impact_radius-E_sample[sample_entry[1]].radius)/radius_factor_length;
		radius_factor[1] = fabs(impact_radius-E_sample[sample_entry[0]].radius)/radius_factor_length;
		double phi_factor_length = 45;
		double phi_factor_1 = fabs(impact_phi-E_sample[sample_entry[0]].phi);
		double phi_factor_2 = fabs(impact_phi-E_sample[sample_entry[2]].phi);
		if(phi_factor_1>46){phi_factor_1 = 360-phi_factor_1;}
		if(phi_factor_2>46){phi_factor_2 = 360-phi_factor_2;}
		phi_factor[0] = phi_factor_2/phi_factor_length;
		phi_factor[1] = phi_factor_1/phi_factor_length;
		double distance_factors[4] = {0.};
		distance_factors[0] = radius_factor[0]*phi_factor[0];
		distance_factors[1] = radius_factor[1]*phi_factor[0];
		distance_factors[2] = radius_factor[0]*phi_factor[1];
		distance_factors[3] = radius_factor[1]*phi_factor[1];

		double impact_fields[2] = {0.};
		double position_scale_factors[2] = {0.};
		for(int pol=0; pol<2; pol++){
			for(int i=0; i<4; i++){
				impact_fields[pol] += distance_factors[i]*E_sample[sample_entry[i]].max_E[pol];
				//cout<<"distance factos: "<<distance_factors[i]<<endl;
				//cout<<"E: "<<E_sample[sample_entry[i]].max_E[pol]<<endl;
			}
			position_scale_factors[pol] = impact_fields[pol]/E_sample[sample_entry[0]].max_E[pol];
			//cout<<position_scale_factors[pol]<<endl;
		}
		
		if((sample_entry[0]>coreas_entries-9)||(sample_entry[0]<9)){
			position_scale_factors[0] = 1.;
			position_scale_factors[1] = 1.;
		}
		
		
		
		
			//cout<<s_entry<<endl;
		EfieldTree->GetEntry(s_entry);
		int pps_entry = round(r3->Rndm()*pps_entries);
		pps_samples_Tree->GetEntry(pps_entry);
		bGood = true;
		
		gThresholdStep = 600;
		
		
		if(sp%500==0){
			gTrigBits.ResetAllBits(false);
			for(int i=0; i<1500; i++){
				triggered_T1H[i] = pps_T1H[i];
				triggered_T2H[i] = pps_T2H[i];
				triggered_T3H[i] = pps_T3H[i];
				triggered_T4H[i] = pps_T4H[i];
				triggered_T1V[i] = pps_T1V[i];
				triggered_T2V[i] = pps_T2V[i];
				triggered_T3V[i] = pps_T3V[i];
				triggered_T4V[i] = pps_T4V[i];
			}
				
			make_tree->Fill();
		}
		
		
		double x_vpp[8][4096] = {0.};
		for(int i=0; i<1500; i++){

			x_vpp[0][i] = coreas_T1H[i] * scale_factor * position_scale_factors[0];
			x_vpp[2][i] = coreas_T2H[i] * scale_factor * position_scale_factors[0];
			x_vpp[4][i] = coreas_T3H[i] * scale_factor * position_scale_factors[0];
			x_vpp[6][i] = coreas_T4H[i] * scale_factor * position_scale_factors[0];
			x_vpp[1][i] = coreas_T1V[i] * scale_factor * position_scale_factors[1];
			x_vpp[3][i] = coreas_T2V[i] * scale_factor * position_scale_factors[1];
			x_vpp[5][i] = coreas_T3V[i] * scale_factor * position_scale_factors[1];
			x_vpp[7][i] = coreas_T4V[i] * scale_factor * position_scale_factors[1];
			
		}

		bool need_trigger = false;
		for(int i=0; i<1500; i++){
			if((x_vpp[0][i]>(10*32512./500.))||(x_vpp[1][i]>(10*32512./500.))){
				need_trigger = true;
				break;
			}
			if((x_vpp[0][i]>32512)||(x_vpp[1][i]>32512)){
				bGood = false;
			}
		}
		/*
		for(int i=0; i<1500; i++){
		
			x_vpp[0][i] = (x_vpp[0][i])*500./32512 + r3->Gaus(0.,9.);
			x_vpp[2][i] = (x_vpp[2][i])*500./32512 + r3->Gaus(0.,9.);
			x_vpp[4][i] = (x_vpp[4][i])*500./32512 + r3->Gaus(0.,9.);
			x_vpp[6][i] = (x_vpp[6][i])*500./32512 + r3->Gaus(0.,6.4);
			x_vpp[1][i] = (x_vpp[1][i])*500./32512 + r3->Gaus(0.,9.);
			x_vpp[3][i] = (x_vpp[3][i])*500./32512 + r3->Gaus(0.,9.);
			x_vpp[5][i] = (x_vpp[5][i])*500./32512 + r3->Gaus(0.,9.);
			x_vpp[7][i] = (x_vpp[7][i])*500./32512 + r3->Gaus(0.,6.4);
		}
		*/
		
		for(int i=0; i<1500; i++){
		
			x_vpp[0][i] = (x_vpp[0][i] + pps_T1H[i])*500./32512;
			x_vpp[2][i] = (x_vpp[2][i] + pps_T2H[i])*500./32512;
			x_vpp[4][i] = (x_vpp[4][i] + pps_T3H[i])*500./32512;
			x_vpp[6][i] = (x_vpp[6][i] + pps_T4H[i])*500./32512;
			x_vpp[1][i] = (x_vpp[1][i] + pps_T1V[i])*500./32512;
			x_vpp[3][i] = (x_vpp[3][i] + pps_T2V[i])*500./32512;
			x_vpp[5][i] = (x_vpp[5][i] + pps_T3V[i])*500./32512;
			x_vpp[7][i] = (x_vpp[7][i] + pps_T4V[i])*500./32512;
		}
		
		
		gOverVoltBits.ResetAllBits(false);
		for(int ch=0; ch<8; ch++){
			if(TMath::MaxElement(1500, x_vpp[ch])>500){gOverVoltBits.SetBitNumber(ch, kTRUE);}
			if(TMath::MinElement(1500, x_vpp[ch])<-500){gOverVoltBits.SetBitNumber(ch, kTRUE);}
		}

		//for(int i=0; i<8; i++){
		//	cout<<impact_radius<<"   "<<i<<": "<<*std::max_element(&x_vpp[i][0], &x_vpp[i][0]+1500)<<endl;
		//}
		//double max_T1H = *std::max_element(&x_vpp[0][0], &x_vpp[0][0]+1500);
		//double max_T1V = *std::max_element(&x_vpp[1][0], &x_vpp[1][0]+1500);
		//cout<<max_T1H<<" "<<max_T1V<<endl;
		
		bool is_trigger = false;
		if(need_trigger){
			gTrigBits.ResetAllBits(false);
			if((gThresholdStep<480)||(gThresholdStep>1280)){continue;}
			is_trigger = IsTrigger(x_vpp, gThresholdStep, gTrigBits);
		}

		if(is_trigger){
			trigger_count++;
			//cout<<pps_entry<<endl;
			
			for(int i=0; i<1500; i++){
				triggered_T1H[i] = (Int_t) round(x_vpp[0][i]*32512./500.);
				triggered_T2H[i] = (Int_t) round(x_vpp[2][i]*32512./500.);
				triggered_T3H[i] = (Int_t) round(x_vpp[4][i]*32512./500.);
				triggered_T4H[i] = (Int_t) round(x_vpp[6][i]*32512./500.);
				triggered_T1V[i] = (Int_t) round(x_vpp[1][i]*32512./500.);
				triggered_T2V[i] = (Int_t) round(x_vpp[3][i]*32512./500.);
				triggered_T3V[i] = (Int_t) round(x_vpp[5][i]*32512./500.);
				triggered_T4V[i] = (Int_t) round(x_vpp[7][i]*32512./500.);
			}
				
			make_tree->Fill();
		}
		
		
	}
	
	

	
	
	delete EfieldTree;
	coreas_file->Close();
	delete coreas_file;
	cout<<trigger_count<<"   "<<max_radius<<endl;
	
	DSR_radius = max_radius;
	Triggered_count = trigger_count;
	
	
	return 3.14159*max_radius*max_radius*trigger_count/Total_count;


}

void pps_initial(TTree *pps_samples_Tree){


  
  pps_samples_Tree->SetBranchAddress("date", &gDate);
  pps_samples_Tree->SetBranchAddress("EventNum", &gEventNum);
  pps_samples_Tree->SetBranchAddress("eventTime", &eventTime); //PC timestamp
  pps_samples_Tree->SetBranchAddress("EventRate", &gEventRate);
  pps_samples_Tree->SetBranchAddress("ThresholdStep", &gThresholdStep);
  pps_samples_Tree->SetBranchAddress("T1H", pps_T1H);
  pps_samples_Tree->SetBranchAddress("T2H", pps_T2H);
  pps_samples_Tree->SetBranchAddress("T3H", pps_T3H);
  pps_samples_Tree->SetBranchAddress("T4H", pps_T4H);
  pps_samples_Tree->SetBranchAddress("T1V", pps_T1V);
  pps_samples_Tree->SetBranchAddress("T2V", pps_T2V);
  pps_samples_Tree->SetBranchAddress("T3V", pps_T3V);
  pps_samples_Tree->SetBranchAddress("T4V", pps_T4V);
  

}


bool CreateDataTreeBranches(TTree* tr){

  printf("CreateDataTreeBranches(): %s\t%s\n", tr->GetName(), tr->GetTitle() );
  tr->Branch("Energy", &Energy,"Energy/D");
  tr->Branch("ShowerTheta", &ShowerTheta,"ShowerTheta/D");
  tr->Branch("ShowerPhi", &ShowerPhi,"ShowerPhi/D");
  tr->Branch("date", &gDate, "date/I");
  tr->Branch("event", &gEventNum, "event/I");
  tr->Branch("eventTime", &eventTime, "eventTime/D"); //PC timestamp
  tr->Branch("EventRate", &gEventRate, "EventRate/D");
  tr->Branch("ThresholdStep", &gThresholdStep, "ThresholdStep/I");
  tr->Branch("bGood", &bGood,"bGood/O");   
  tr->Branch("triggerBits", &gTrigBits); //TAROGE-4  TBits, align with ChIDs    Hpol and Vpol
  tr->Branch("overVoltBits", &gOverVoltBits); //TAROGE-4  TBits(kNCh)   overvoltage or not
  
  tr->Branch("radius", &radius, "radius/D");
  tr->Branch("pattern_phi", &pattern_phi, "pattern_phi/D");

  tr->Branch("T1H", triggered_T1H, "T1H[1500]/S");
  tr->Branch("T2H", triggered_T2H, "T2H[1500]/S");
  tr->Branch("T3H", triggered_T3H, "T3H[1500]/S");
  tr->Branch("T4H", triggered_T4H, "T4H[1500]/S");
  tr->Branch("T1V", triggered_T1V, "T1V[1500]/S");
  tr->Branch("T2V", triggered_T2V, "T2V[1500]/S");
  tr->Branch("T3V", triggered_T3V, "T3V[1500]/S");
  tr->Branch("T4V", triggered_T4V, "T4V[1500]/S");
  // /S: 2-byte signed short  /I: 4 byte signed int  /F double
 

 
  return true;
}
