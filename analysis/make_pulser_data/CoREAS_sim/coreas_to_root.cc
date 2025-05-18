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
using namespace std;
#define PI TMath::Pi()
#define Speed_Of_Light 2.99792458e8

#define HALF_WINDOW_SPACING 75
 
void coreas_to_root(const string strCoreasDir);

  
  Int_t T1H[1500];
  Int_t T2H[1500];
  Int_t T3H[1500];
  Int_t T4H[1500];
  Int_t T1V[1500];
  Int_t T2V[1500];
  Int_t T3V[1500];
  Int_t T4V[1500];
  
  Double_t EF_theta[32768];
  Double_t EF_phi[32768];


bool loopCoreasToRoot(int run1, int run2){
  
	get_response();
	
	initial_antenna_direction();
	initial_coreas_coordinate();
	
	
	for(int run=run1; run<=run2; run++){
		//if(run%1000>19){continue;}
		const string strCoreasDir = Form("/media/cyc/SSD_data_disk/Energy_estimation/rootfiles/Coreas-t4-r%06d.root", run);
		// 0 means exist, -1 means not
		if(access(strCoreasDir.c_str(), F_OK)!=0)
		{
			cout<<strCoreasDir<<"  do not exist! skip it"<<endl;
			continue;
		}
		coreas_to_root(strCoreasDir);
		
	}
	
	


	return 1;
}


void coreas_to_root(const string strCoreasDir){


    double ShowerZenithAngle = 0.;
    double ShowerAzimuthAngle = 0.; // shower propagate to; in CoREAS coordinates(x axis = 0°　= magnetic north)
    double _time, _E_north, _E_west, _E_vert;
    vector<double> *vE_north=0;
    vector<double> *vE_west=0;
    vector<double> *vE_vert=0;
    
    double ShowerTheta = 0.;
    double ShowerPhi = 0.;
    double distance_Xmax_m = 0.;

	TFile *CoreasRoot = new TFile(strCoreasDir.c_str(),"UPDATE");
	
	TTree *CRTree = (TTree*) CoreasRoot->Get("CRTree");
	CRTree->SetBranchAddress("zenith_deg", &ShowerZenithAngle);
	CRTree->SetBranchAddress("azimuth_corsika_deg", &ShowerAzimuthAngle);
	CRTree->SetBranchAddress("distance_Xmax_m", &distance_Xmax_m);
	
	CRTree->GetEntry(0);
	
	
	TTree *EfieldTree = (TTree*) CoreasRoot->Get("EfieldTree");
	EfieldTree->SetBranchAddress("vE_north",  &vE_north);
	EfieldTree->SetBranchAddress("vE_west",  &vE_west);
	EfieldTree->SetBranchAddress("vE_vert",  &vE_vert);
	
	auto T1H_branch = EfieldTree->Branch("T1H",T1H,"T1H[1500]/I");
	auto T2H_branch = EfieldTree->Branch("T2H",T2H,"T2H[1500]/I");
	auto T3H_branch = EfieldTree->Branch("T3H",T3H,"T3H[1500]/I");
	auto T4H_branch = EfieldTree->Branch("T4H",T4H,"T4H[1500]/I");
	auto T1V_branch = EfieldTree->Branch("T1V",T1V,"T1V[1500]/I");
	auto T2V_branch = EfieldTree->Branch("T2V",T2V,"T2V[1500]/I");
	auto T3V_branch = EfieldTree->Branch("T3V",T3V,"T3V[1500]/I");
	auto T4V_branch = EfieldTree->Branch("T4V",T4V,"T4V[1500]/I");
	
	auto ETheta_branch = EfieldTree->Branch("ETheta",EF_theta,"EF_theta[32768]/D");
	auto EPhi_branch = EfieldTree->Branch("EPhi",EF_phi,"EF_phi[32768]/D");
	
	auto Theta_branch = EfieldTree->Branch("ShowerTheta", &ShowerTheta,"ShowerTheta/D");
	auto Phi_branch = EfieldTree->Branch("ShowerPhi", &ShowerPhi,"ShowerPhi/D");
	auto Xmax_branch = EfieldTree->Branch("distance_Xmax_m", &distance_Xmax_m,"distance_Xmax_m/D");
	
    
    
    CoreasAngleToStation(ShowerZenithAngle, ShowerAzimuthAngle, ShowerTheta, ShowerPhi);
    StationAngleToAntenna(ShowerTheta, ShowerPhi, antenna_theta, antenna_phi, E_theta_station, E_phi_station);


    double dt_angle[4] = {0};
    for(int k=0; k<3; k++){
    	dt_angle[k] = expectedTimeDiff_angle_with_R(ShowerTheta, ShowerPhi, 50000., coor_ant[6], coor_ant[k*2]);
    } 
	
	double E_theta[4][4096*8] = {0.};
	double E_phi[4][4096*8] = {0.};
	double E_waveform[8][4096*8] = {0.};
	
	for(int entry=0; entry<EfieldTree->GetEntries(); entry++){
		EfieldTree->GetEntry(entry);
	
    	for(int i=0; i<vE_north->size(); ++i){
		
    		_E_north = vE_north->at(i);
    		_E_west = vE_west->at(i);
    		_E_vert = vE_vert->at(i);
    		
    		TVector3 _E_vector = _E_north*coreas_x + _E_west*coreas_y + _E_vert*coreas_z;
    		
    		for(int ant=0; ant<4; ant++){
    			// put waveform in middle, first 400 ns to be 0.
    			E_theta[ant][500+i] = _E_vector * E_theta_station[ant];
    			E_phi[ant][500+i] = _E_vector * E_phi_station[ant];
    		}
    		
  		}
  		
  		for(int i=0; i<4096*8; ++i){
  			EF_theta[i] = E_theta[0][i];
  			EF_phi[i] = E_phi[0][i];
  		}
  		
  		for(int ch=0; ch<8; ch++){
  			int ant = ch/2;
  			int pos_theta = round((antenna_theta[ant] + 90.)/5.);
    		int pos_phi = round(antenna_phi[ant]/5.);
  			coreas_covolution(E_theta[ant], E_phi[ant], E_waveform[ch], pos_theta, pos_phi, ch, dt_angle[ant]);
  		}
  		



		

		for (int i=0;i<1500;i++)
		{
			T1H[i] = (Int_t) round(E_waveform[0][i]*32512./500.);
			T2H[i] = (Int_t) round(E_waveform[2][i]*32512./500.);
			T3H[i] = (Int_t) round(E_waveform[4][i]*32512./500.);
			T4H[i] = (Int_t) round(E_waveform[6][i]*32512./500.);
			T1V[i] = (Int_t) round(E_waveform[1][i]*32512./500.);
			T2V[i] = (Int_t) round(E_waveform[3][i]*32512./500.);
			T3V[i] = (Int_t) round(E_waveform[5][i]*32512./500.);
			T4V[i] = (Int_t) round(E_waveform[7][i]*32512./500.);
	  	  
		}
	T1H_branch->Fill();
	T2H_branch->Fill();
	T3H_branch->Fill();
	T4H_branch->Fill();
	T1V_branch->Fill();
	T2V_branch->Fill();
	T3V_branch->Fill();
	T4V_branch->Fill();
	ETheta_branch->Fill();
	EPhi_branch->Fill();
	Theta_branch->Fill();
	Phi_branch->Fill();
	Xmax_branch->Fill();
		  
 	}// entry loop
 	

EfieldTree->Write(0,TObject::kOverwrite);
CoreasRoot->Close();

   
}//end of main function





