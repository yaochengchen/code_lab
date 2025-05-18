// 20210526, change the antenna coordinates and delay order, from T1H T2H....T3V T4V -->> T1H T1V....T4H T4V
// zero_correct_x is wrong! since antenna 4 has oscillator signal!!
#include "Reconstruction_tchain_new.h"

 
template <typename T>
T** new_Array2D(int row, int col);

template <typename T> 
void delete_Array2D(T **arr, int row, int col);
 
void get_cross_correlation(int pol, int nt_j, double tdoa[2][6], double Xcor[2][6]);
void real_cross_correlation(int pol, int pair, int i, int j);
int filter_time(double filtered_x[1500], double new_x[1500], int Max_i);
int filter_fre(int ch, double &SNR, bool &multi_pulses);
int get_FFT(int ch);
Double_t expectedTimeDiff_angle(double theta, double phi, double *xyz);
Double_t expectedTimeDiff_coord(double *drone_xyz, double *antenna_i, double *antenna_j);
Double_t expectedTimeDiff_angle_with_R(double theta, double phi, double R, double *antenna_i, double *antenna_j);

void filtering(double saSNR[2], int multi_pulses[2]);
void get_simple_Xcor(int pol, double tdoa[2][6], double Xcor[2][6]);
void reconstruction(int pol, double event_theta[2], double event_phi[2], double x_value[2]);
double fix_angle_cross_correlation(int pol, double in_theta, double in_phi);
double multi_pulse_cross_correlation(int pol, double in_theta, double in_phi);
double cross_pol_cross_correlation(int pol, double in_theta, double in_phi);
double spectrum_cut(int pol, int run_pos);

Double_t * fft_back_x = new Double_t[1500];
Double_t ** filtered_x = new_Array2D<double>(8, 1500);
Double_t ** zero_correct_x = new_Array2D<double>(8, 1500);

Double_t pps_ave_power[200][8][199] = {0.};//199 //180-350MHz in 199 bins; one day maximum has 172 runs
Double_t power_frequency[199] = {0.};//180-350MHz in 199 bins;
Double_t trigger_spectrum[8][199] = {0.};

bool Is_saturated();
 
 //std::thread threads[12];

//base center   2021 photogrametry  Tower-4 is back side
double coor_ant[8][3] = {{328914.252, 2697625.953, 709.389},{328914.252, 2697625.953, 709.389},{328919.115, 2697635.05, 707.86},{328919.115, 2697635.05, 707.86},{328923.564, 2697643.311, 706.08},{328923.564, 2697643.311, 706.08},{328936.041, 2697639.685, 704.949},{328936.041, 2697639.685, 704.949}};
double delay_ant[8] = {0,-92.4,53.2,-30.3,84.6,59.5,-1523.1,-743.1};
int peak_position[8] = {0};

 Double_t x[8][1500]={0};//8 antennas
 Double_t new_x[8][1500]={0};//8 antennas
 double remain_x[8][1500]={0};//8 antennas
 double ant_re[8][HALF_WINDOW_WIDTH] = {0};
 double ant_im[8][HALF_WINDOW_WIDTH] = {0};
 //double ant_re[8][2049] = {0};
 //double ant_im[8][2049] = {0};
 double Re_ant[8][751] = {0};
 double Im_ant[8][751] = {0};
  


 double correlation_pair[2][6][32768] = {0};
 
 
double Xcor_tdoa[6][100000]={0};
double Xcor_max[100000][6]={0};
double Sum_x[8] = {0};
double Correlator_re[2][6][HALF_WINDOW_WIDTH] = {0};            
double Correlator_im[2][6][HALF_WINDOW_WIDTH] = {0};


double Expect_Xcor_tdoa[6][100000]={0};

int event_record[100000][2] = {0};// run, entry


int pass_count[10] = {0};




int gRunNum = 0;
int gEventNum = 0;
double eventTime = 0;
int V_bits = 0;
int H_bits = 0;
int iGood = 0;
int is_pps = 0;
double rough_cross_correlation[2] = {0.};
double rough_SNR[2] = {0.};
int multi_pulses[2] = {0};
double reconstruct_cross_correlation[2] = {0.};
double remain_pulse_cross_correlation[2] = {0.};
double other_pol_cross_correlation[2] = {0.};
double pass_trigger_ratio[2] = {0.};
void SetTreeBranches(TTree* tr);
void initial_data();

TRandom3 *r3;


int main(int argc, char *argv[]) {

r3 = new TRandom3();
//cout<<"fuck "<<r3->Gaus(0, 20)<<endl;

string date_name = argv[1];
string run_list = argv[2];

cout<<endl<<date_name<<endl;
cout<<run_list<<endl;

//void Real_CRselect_tchain(string date_name, string run_list){

int which_date = std::stoi(date_name);

TChain *tr_event = new TChain("t");
//TChain *tr_header = new TChain("tHeader");
ifstream inFile(run_list);
	string str;
	int fNFile = 0;
	vector<Long64_t> fvNEntries;
	vector<int> run_number;
	
	int which_run = 0;
	int which_event = 0;

	while(getline(inFile, str)){

		//cout<<"reading: "<< str <<endl;
		if(str.find(".root") != std::string::npos)
		{
			cout<<"here"<<endl;
		    int nPos_1 = str.find_last_of("n");
            int nPos_2 = str.find_last_of(".");
            string run_num="";
            run_num = str.substr(nPos_1+1,nPos_2-nPos_1-1);
            run_number.push_back(std::stoi(run_num));
            
			fNFile +=  tr_event->AddFile(str.c_str());

			fvNEntries.push_back(tr_event->GetEntries());
		}


	}

	inFile.close();
	
	//cout<<fNFile<<endl;
	//for(auto context : fvNEntries){
		//cout<<context<<endl;
		//}
	//for(auto context : run_number){
		//cout<<context<<endl;1.65	*	//}




  const   char * out_filename;

  out_filename =  Form("./candidate_events/%08d.txt", which_date);
  fstream out_put;
  out_put.open(out_filename, ios::out | ios::trunc);// without ios::out   clean the file before write






double t4point_phi = 90 - 117.09;
double t4point_theta = -0.08;
double base_thick = 0.01;

//transfer t4 back side to front side
coor_ant[6][0] += TMath::Cos(t4point_theta*TMath::DegToRad()) * TMath::Cos(t4point_phi*TMath::DegToRad()) * base_thick;
coor_ant[6][1] += TMath::Cos(t4point_theta*TMath::DegToRad()) * TMath::Sin(t4point_phi*TMath::DegToRad()) * base_thick;
coor_ant[6][2] += TMath::Sin(t4point_theta*TMath::DegToRad()) * base_thick;

coor_ant[7][0] += TMath::Cos(t4point_theta*TMath::DegToRad()) * TMath::Cos(t4point_phi*TMath::DegToRad()) * base_thick;
coor_ant[7][1] += TMath::Cos(t4point_theta*TMath::DegToRad()) * TMath::Sin(t4point_phi*TMath::DegToRad()) * base_thick;
coor_ant[7][2] += TMath::Sin(t4point_theta*TMath::DegToRad()) * base_thick;



//set tower-2 to be origin
double set_t2_origin_photo[3] = {328919.115, 2697635.05, 707.86};//same as tower2
for(int i=0; i<8; i++){
for(int j=0; j<3; j++){
coor_ant[i][j] = coor_ant[i][j] - set_t2_origin_photo[j];
cout<<coor_ant[i][j]<<"  ";
}
cout<<endl;
}

  
out_filename =  Form("t_selection_%08d.root", which_date);
TFile *selection_file = new TFile(out_filename,"RECREATE");
TTree *selectionTree = new TTree("selection","selection");

SetTreeBranches(selectionTree);

 
  


TH2F *h2 = new TH2F("h2","h2", 25, -26.2, -23.7, 110, -6, 5);
//TCanvas *c4 = new TCanvas("c4","Canvas Example",1500,1000);

TCanvas *c1 = new TCanvas("c1","Canvas Example",1500,1000);
//c1->Divide(3,2);

TCanvas *c2 = new TCanvas("c2","Canvas Example",1500,1000);
c2->Divide(2,1);

//TCanvas *c3 = new TCanvas("c3","Canvas Example",1500,1000);

//c3->Divide(1,3);

float nt_run;
float nt_event;
float nt_time;
float drn_azimuth;
float drn_elevation;




double event_time[100000] = {0};
double event_theta[100000][2] = {0};
double event_phi[100000][2] = {0};
double x_value[100000][2] = {0};
//double delta_HV_theta[100000] = {0};
//double delta_HV_phi[100000] = {0};



int nt_j = 0;

 
 double Sum_template = 0;
double Re_template[2049] = {0};
double Im_template[2049] = {0};

 

  bool bGood;
  Int_t timeStamp_FPGA; 

  Int_t event_number;

int Event_N = 0;
//double eventTime = 0;
double eventTime_start = 0;

double time_prev = 0;
double time_accum = 0;

double timestamp[100000]={0};

double event[100000]={0};


int pulser_event_count=0;




double cal_timestamp[100000]={0};
double cal_event[100000]={0};
int cal_Event_N = 0;

int flag_check[100000] = {0};

  Double_t t[1500];
  Short_t T1H[1500];
  Short_t T2H[1500];
  Short_t T3H[1500];
  Short_t T4H[1500];
  Short_t T1V[1500];
  Short_t T2V[1500];
  Short_t T3V[1500];
  Short_t T4V[1500];
  TBits* triggerBitsccc = 0;
  



int count_first_cut = 0;




nt_j=0;

      	  //tr_header->SetBranchAddress("runNumber",&gRunNum);
      	  
      	  tr_event->SetBranchAddress("event",&gEventNum);
      	  tr_event->SetBranchAddress("bGood",&bGood);
		  tr_event->SetBranchAddress("timeStamp_FPGA",&timeStamp_FPGA);
		  tr_event->SetBranchAddress("eventTime", &eventTime);
		  tr_event->SetBranchAddress("T1H",T1H);
		  tr_event->SetBranchAddress("T2H",T2H);
		  tr_event->SetBranchAddress("T3H",T3H);
		  tr_event->SetBranchAddress("T4H",T4H);
		  tr_event->SetBranchAddress("T1V",T1V);
		  tr_event->SetBranchAddress("T2V",T2V);
		  tr_event->SetBranchAddress("T3V",T3V);
		  tr_event->SetBranchAddress("T4V",T4V);
		  // H4, H3, H2, H1, V4, V3, V2, V1
		  tr_event->SetBranchAddress("triggerBits",&triggerBitsccc);
		  
 event_number = tr_event->GetEntries();
 
 cout<<event_number<<endl;

 pass_count[0] = event_number;


	double FreqBin = 5./6.;
	for(int k=0; k<199; k++){
		int _k = k+222;
		power_frequency[k] = _k*FreqBin;
	}


	
	int run_pos = 0;
	int num_pps = 0;
	int max = fvNEntries.at(run_pos);
	// initial all value to 1;
	//fill(pps_ave_power[0][0], pps_ave_power[0][0] + 200*8*20, 1);
	// for finding median value as background strength.
	double median[8][20][500] = {0.};
	for(int entry = 0; entry<event_number/100; entry++){
	
		if(entry == max){
			for(int ch=0; ch<8; ch++){
				for(int k=0; k<199; k++){//199
					// take average in dB unit
					//cout<<pps_ave_power[run_pos][ch][k]<<endl;
					//pps_ave_power[run_pos][ch][k] /= ((double) num_pps);
					//cout<<pps_ave_power[run_pos][ch][k]<<endl<<endl;
					sort(median[ch][k], median[ch][k]+num_pps);
					pps_ave_power[run_pos][ch][k] = median[ch][k][num_pps/2];
					//cout<<ch<<"  "<<k<<"   "<<pps_ave_power[run_pos][ch][k]<<"   "<<median[ch][k][num_pps/2]<<endl;
				}
			}
			num_pps = 0;
			run_pos += 1;
			max = fvNEntries.at(run_pos);
			//break;
		}
	
		tr_event->GetEntry(entry);
		if(triggerBitsccc->CountBits() != 0) {continue;}
		
  for (int i=0;i<1500;i++)
	{
	  t[i]=i*0.8;
	  x[0][i]= (float)T1H[i]*500/32512;
	  x[2][i]= (float)T2H[i]*500/32512;
	  x[4][i]= (float)T3H[i]*500/32512;
	  x[6][i]= (float)T4H[i]*500/32512;
	  x[1][i]= (float)T1V[i]*500/32512;
	  x[3][i]= (float)T2V[i]*500/32512;
	  x[5][i]= (float)T3V[i]*500/32512;
	  x[7][i]= (float)T4V[i]*500/32512;	  	  
	}
	
	//cout<<x[0][502]<<"  "<<x[1][502]<<endl;

  if(Is_saturated()) {continue;}
		
	
		for(int ch=0; ch<8; ch++){
	
			get_FFT(ch);
		
			//for(int k=0; k<199; k++){
			//	int _k = k+222;
			//	pps_ave_power[run_pos][ch][k] += Re_ant[ch][_k]*Re_ant[ch][_k] + Im_ant[ch][_k]*Im_ant[ch][_k];
			//}
			for(int k=0; k<199; k++){
				int _k = k + 222;
				double power = Re_ant[ch][_k]*Re_ant[ch][_k] + Im_ant[ch][_k]*Im_ant[ch][_k];

				median[ch][k][num_pps] = power;
				
			}
			
		}
	//selectionTree->Fill();
	num_pps++;
	//cout<<num_pps<<endl;
	//if(num_pps>248){num_pps = 249;}
	}

	// last one need to average out manually;
	for(int ch=0; ch<8; ch++){
		for(int k=0; k<20; k++){//199
		// take average in dB unit
		//pps_ave_power[run_pos][ch][k] /= ((double) num_pps);
		sort(median[ch][k], median[ch][k]+num_pps);
		pps_ave_power[run_pos][ch][k] = median[ch][k][num_pps/2];
		}
	}	
	




	run_pos = 0;
	max = fvNEntries.at(run_pos);
for(int entry = 0; entry<event_number/100; entry++){

	if(entry == max){
		run_pos += 1;
		max = fvNEntries.at(run_pos);
	}

	gRunNum = run_number.at(run_pos);
		
	if(entry%1000==0){cout<<entry<<endl;}// cout<<triggerBitsccc->CountBits()<<endl;}

      //if(triggerBitsccc->CountBits() != 0) {continue;}
      //if(entry%5000){continue;}
      
      tr_event->GetEntry(entry);
      cout<<triggerBitsccc->CountBits()<<endl;
      
      
      //tr_header->GetEntry(entry);
      //cout<<gRunNum<<endl;
      
      int total_trigger_bit = triggerBitsccc->CountBits();
      int V_trigger_bit = triggerBitsccc->CountBits(32);
      int H_trigger_bit = total_trigger_bit - V_trigger_bit;
      
      //cout<<V_trigger_bit<<"    "<<H_trigger_bit<<"    "<<total_trigger_bit<<endl;
      
      V_bits = V_trigger_bit;
      H_bits = H_trigger_bit;
      is_pps = total_trigger_bit == 0;
      //cout<<"is pps?:  "<<is_pps<<endl;
      


if(entry==0) eventTime_start = eventTime;

  for (int i=0;i<1500;i++)
	{
	  t[i]=i*0.8;
	  x[0][i]= (float)T1H[i]*500/32512;
	  x[2][i]= (float)T2H[i]*500/32512;
	  x[4][i]= (float)T3H[i]*500/32512;
	  x[6][i]= (float)T4H[i]*500/32512;
	  x[1][i]= (float)T1V[i]*500/32512;
	  x[3][i]= (float)T2V[i]*500/32512;
	  x[5][i]= (float)T3V[i]*500/32512;
	  x[7][i]= (float)T4V[i]*500/32512;	  	  
	}
	
	//cout<<x[0][502]<<"  "<<x[1][502]<<endl;

  if(Is_saturated()) {continue;}

     event_time[nt_j] = eventTime - eventTime_start;
     
    iGood = bGood;
    initial_data();

     
double tdoa[2][6] = {0};
double Xcor[2][6] = {0};
double saSNR[2] = {0};
int mul_pulses[2] = {0};

filtering(saSNR, mul_pulses);

	rough_SNR[0] = saSNR[0];
	rough_SNR[1] = saSNR[1];
	multi_pulses[0] = mul_pulses[0];
	multi_pulses[1] = mul_pulses[1];

	bool single_strong = false;
	double average_Xcor[2] = {0.};
	for(int pol=0; pol<2; pol++){
		//if((rough_SNR[pol]>3)&&(multi_pulses[pol]<4)){
		if(true){
			get_simple_Xcor(pol, tdoa, Xcor);
			for(int pair = pol; pair <6; pair++){
    			average_Xcor[pol] += Xcor[pol][pair]/6.;
    			}
    		rough_cross_correlation[pol] = average_Xcor[pol];
    		pass_trigger_ratio[pol] = spectrum_cut(pol, run_pos);
    		cout<<pass_trigger_ratio[pol]<<endl;
				
    		//if(average_Xcor[pol]>0.8){
    		if(false){
				get_cross_correlation(pol, nt_j, tdoa, Xcor);
				reconstruction(pol, event_theta[nt_j], event_phi[nt_j], x_value[nt_j]);
					
				reconstruct_cross_correlation[pol] = x_value[nt_j][pol]/6.;
				//cout<<x_value[nt_j][pol]<<endl;
				
				//if(x_value[nt_j][pol]>4.8){
				if(true){
					double multi_pulse_cross_correlator = multi_pulse_cross_correlation(pol, event_theta[nt_j][pol], event_phi[nt_j][pol]);
					remain_pulse_cross_correlation[pol] = multi_pulse_cross_correlator/6.;
					//cout<<remain_pulse_cross_correlation[pol]<<endl<<endl;

					//if(multi_pulse_cross_correlator<3){
					if(true){
						double cross_pol_cross_correlator = cross_pol_cross_correlation(pol, event_theta[nt_j][pol], event_phi[nt_j][pol]);
						other_pol_cross_correlation[pol] = cross_pol_cross_correlator/6.;
						//cout<<pol<<"  cross pol  "<<cross_pol_cross_correlator<<endl;
						if(cross_pol_cross_correlator>3.6){single_strong = true;}
						}
					}
    			}
			}
		}
	


single_strong = false;

if(single_strong) {

//cout<<delta_HV_theta[nt_j]<<"    "<<delta_HV_phi[nt_j]<<endl;


	int run_N = fvNEntries.size();
	int run_pos = 0;
	for( ; run_pos<run_N; run_pos++){
		if(entry < fvNEntries.at(run_pos)) {break;}
	}
	 //cout<<run_pos<<endl;
	 which_run = run_number.at(run_pos);
	 
	 //if(run_pos==0){
	 //	which_event = entry;
	 //	}
	 //else{
	 //	which_event = entry-fvNEntries.at(run_pos-1);
	 //	}

//event_record[nt_j][0] = which_run;
//event_record[nt_j][1] = which_event;

    gRunNum = which_run;
    //gEventNum = which_event;
    cout<<gRunNum<<"   "<<gEventNum<<endl;

//cout<<which_date<<"    "<<event_record[nt_j][0]<<"    "<<event_record[nt_j][1]<<"  "<<event_time[nt_j]<<"    "<<event_theta[0][nt_j]<<"    "<<event_phi[0][nt_j]<<"    "<<event_theta[1][nt_j]<<"    "<<event_phi[1][nt_j]<<"   "<<event_time[nt_j] - event_time[nt_j-1]<<endl;

//nt_j++; 

pass_count[2] += 1;

//cout<<number<<"   "<<entry<<endl;
}

//nt_j++;

selectionTree->Fill();

}// event loop


   selection_file->cd();
   selectionTree->Write("", TObject::kOverwrite);
   selection_file->Close();
   delete selection_file;


cout<<"nt_j:   "<<nt_j<<endl;

for(int i=0; i<nt_j; i++){

	bool record = true;
	for(int j=i-1; j>=0; j--){
	    //caution!!!! cyc, should be 10 instead of 0
		if((event_time[i] - event_time[j]) < 10){
			if((x_value[i][0]>4.2)&&(x_value[j][0]>4.2)&&(fabs(event_theta[0][i]-event_theta[0][j])<2)&&(fabs(event_phi[0][i]-event_phi[0][j])<2)){
				record = false;
				}
			if((x_value[i][1]>4.2)&&(x_value[j][1]>4.2)&&(fabs(event_theta[1][i]-event_theta[1][j])<2)&&(fabs(event_phi[1][i]-event_phi[1][j])<2)){
				record = false;
				}
			}
		else{
			break;
			}
		}
		
	for(int k=i+1; k<nt_j; k++){
	    //caution!!!! cyc, should be 10 instead of 0
		if((event_time[k] - event_time[i]) < 10){
			if((x_value[k][0]>4.2)&&(x_value[i][0]>4.2)&&(fabs(event_theta[0][k]-event_theta[0][i])<2)&&(fabs(event_phi[0][k]-event_phi[0][i])<2)){
				record = false;
				}
			if((x_value[k][1]>4.2)&&(x_value[i][1]>4.2)&&(fabs(event_theta[1][k]-event_theta[1][i])<2)&&(fabs(event_phi[1][k]-event_phi[1][i])<2)){
				record = false;
				}
			}
		else{
			break;
			}
		}

	if(record){
	out_put<<which_date<<"    "<<event_record[i][0]<<"    "<<event_record[i][1]<<"    "<<x_value[i][0]<<"    "<<event_theta[0][i]<<"    "<<event_phi[0][i]<<"    "<<x_value[i][1]<<"    "<<event_theta[1][i]<<"    "<<event_phi[1][i]<<endl;
	pass_count[3] += 1;
		}

	}



cout<<count_first_cut<<"  "<<nt_j<<endl;
cout<<pass_count[0]<<"  "<<pass_count[1]<<"   "<<pass_count[2]<<"   "<<pass_count[3]<<endl;


  out_put.close(); 
    
              }// end of main function
              
              
Double_t * sXcor_time = new Double_t [WINDOW_SPACING];

   int sXcor_N = WINDOW_SPACING;
   Int_t sXcor_size = sXcor_N+1;		
   TVirtualFFT *sXcor_fft_forward = TVirtualFFT::FFT(1, &sXcor_size,"R2C ES K");
   
   TVirtualFFT *sXcor_fft_backward = TVirtualFFT::FFT(1, &sXcor_N, "C2RBACKWARD M K");
   
void filtering(double saSNR[2], int multi_pulses[2]){
	double sSNR[8] = {0};
	bool sMulti_pulses[8] = {false}; 
  
	for(int ant = 0; ant < 8; ant ++){
		Sum_x[ant] = 0;

	
    	int max_i = filter_fre(ant, sSNR[ant], sMulti_pulses[ant]);
    	peak_position[ant] = max_i;
    	filter_time(filtered_x[ant], new_x[ant], max_i);

    //for(int j=0;j<1500;j++){// make waveform to be 100 ns
      //if((j<(max_i-80)) || (j>(max_i+80)))new_x[ant][j] = 0.;
    //}    	

   sXcor_fft_forward->SetPoints(&new_x[ant][max_i-HALF_WINDOW_SPACING]);
   sXcor_fft_forward->Transform();
   sXcor_fft_forward->GetPointsComplex(ant_re[ant],ant_im[ant]);

    	  for(int i = 0; i <WINDOW_SPACING; i++){
          Sum_x[ant] += pow(new_x[ant][max_i-HALF_WINDOW_SPACING+i], 2);
          
          }
   
   }
   
   saSNR[0] = (sSNR[0] + sSNR[2] +sSNR[4] +sSNR[6])/4.;
   saSNR[1] = (sSNR[1] + sSNR[3] +sSNR[5] +sSNR[7])/4.;

   multi_pulses[0] = (sMulti_pulses[0] + sMulti_pulses[2] +sMulti_pulses[4] +sMulti_pulses[6])/4.;
   multi_pulses[1] = (sMulti_pulses[1] + sMulti_pulses[3] +sMulti_pulses[5] +sMulti_pulses[7])/4.;
}


   
void get_simple_Xcor(int pol, double tdoa[2][6], double Xcor[2][6]){

   //double re[8][HALF_WINDOW_WIDTH] = {0};
   //double im[8][HALF_WINDOW_WIDTH] = {0};

	//double Sum_x[8] = {0};
 

   //0 for H-pol 1 for V-pol
   //for(int pol = 0; pol<2; pol++){
      
   	int i = 0;
    int j = 2;
    
   for(int pair = 0; pair <6; pair++){

    switch(pair){
      case 0:
        i = 0+pol; j = 2+pol;
        break;
      case 1:
        i = 0+pol; j = 4+pol;
        break;
      case 2:
        i = 0+pol; j = 6+pol;
        break;
      case 3:
        i = 2+pol; j = 4+pol;
        break;
      case 4:
        i = 2+pol; j = 6+pol;
        break;
      case 5:
        i = 4+pol; j = 6+pol;
        break;

    }
    


        //double Correlator_re[pol][pair][101] = {0};            
        //double Correlator_im[pol][pair][101] = {0};

                  	   for (Int_t k=1; k<HALF_WINDOW_WIDTH; k++){//0 means DC component
            Correlator_re[pol][pair][k] = (ant_re[i][k] * ant_re[j][k] + ant_im[i][k] * ant_im[j][k]); 
                   
            Correlator_im[pol][pair][k] = (ant_re[i][k] * ant_im[j][k] - ant_im[i][k] * ant_re[j][k]);
            //cout<<ant_re[i][k]<<endl;
            
            //if((pol==0)&&(pair==5)&&(k==100)) cout<<Correlator_re[pol][pair][k]<<endl;
                
                                                    }//k  


        //Double_t * Correlator_time = new Double_t [409600];

          sXcor_fft_backward->SetPointsComplex(Correlator_re[pol][pair],Correlator_im[pol][pair]);
          sXcor_fft_backward->Transform();
          sXcor_time = sXcor_fft_backward->GetPointsReal();
          

//if((pol==0)&&(pair==5)) cout<<sXcor_time[0]<<endl;


         double factor = sqrt(Sum_x[i])*sqrt(Sum_x[j])*WINDOW_SPACING;
         double max = 0;
         int mm;
                       for (Int_t k=0; k<WINDOW_SPACING; k++){
                sXcor_time[k] = sXcor_time[k]/factor;
                if (sXcor_time[k] > max){max = sXcor_time[k]; mm = k;}
                                                    }//k 

//cout<<max<<"  "<<mm<<endl;

if(mm>HALF_WINDOW_SPACING-1) {mm -= WINDOW_SPACING;}

//cout<<mm<<endl;

tdoa[pol][pair] = mm*0.8;
Xcor[pol][pair] = max;
//cout<<pol<<"    "<<pair<<"   "<<factor<<"    "<<Xcor[pol][pair]<<"  "<<tdoa[pol][pair]<<endl;

     }// pair loop

    //}// pol loop

}




	
	
Double_t ** Correlator_time = new_Array2D<double>(12, 32768);

int N = 32768;
//TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &N, "C2RBACKWARD M K");
TVirtualFFT *fft_back[12] = {TVirtualFFT::FFT(1, &N, "C2RBACKWARD M K"), TVirtualFFT::FFT(1, &N, "C2RBACKWARD M K"), TVirtualFFT::FFT(1, &N, "C2RBACKWARD M K"), TVirtualFFT::FFT(1, &N, "C2RBACKWARD M K"), 
TVirtualFFT::FFT(1, &N, "C2RBACKWARD M K"), TVirtualFFT::FFT(1, &N, "C2RBACKWARD M K"), TVirtualFFT::FFT(1, &N, "C2RBACKWARD M K"), TVirtualFFT::FFT(1, &N, "C2RBACKWARD M K"), 
TVirtualFFT::FFT(1, &N, "C2RBACKWARD M K"), TVirtualFFT::FFT(1, &N, "C2RBACKWARD M K"), TVirtualFFT::FFT(1, &N, "C2RBACKWARD M K"), TVirtualFFT::FFT(1, &N, "C2RBACKWARD M K")};

// void gen_fft(TVirtualFFT *fft_back[12], int &N){
//     for(int i=0; i<12; i++){
//       fft_back[i] = TVirtualFFT::FFT(1, &N, "C2RBACKWARD M K");
//     }
// }
// gen_fft(fft_back, N);


void foo(int n);
	
//result: correlation_pair[pair][k]; said pair antenna i and j, correlation_pair[pair][k] is the X-correlation value when we move antenna j signal earier k (just shift toward left side).
void get_cross_correlation(int pol, int nt_j, double tdoa[2][6], double Xcor[2][6]){	



   //0 for H-pol 1 for V-pol
   //for(int pol = 0; pol<2; pol++){
      
   	int i = 0;
    int j = 2;
    
   for(int pair = 0; pair <6; pair++){

    switch(pair){
      case 0:
        i = 0+pol; j = 2+pol;
        break;
      case 1:
        i = 0+pol; j = 4+pol;
        break;
      case 2:
        i = 0+pol; j = 6+pol;
        break;
      case 3:
        i = 2+pol; j = 4+pol;
        break;
      case 4:
        i = 2+pol; j = 6+pol;
        break;
      case 5:
        i = 4+pol; j = 6+pol;
        break;

    }


    real_cross_correlation(pol, pair, i, j);
    //threads[pol*6+pair] = std::thread(real_cross_correlation, pol, pair, i, j);
    //std::thread t1(real_cross_correlation, pol, pair, i, j);
    //t1.join();

    //cout<<mm<<endl;

    //tdoa[pol][pair] = mm/32768.*0.8*WINDOW_SPACING;//32768.
    //Xcor[pol][pair] = max;
    //cout<<pol<<"    "<<pair<<"   "<<factor<<"    "<<Xcor[pol][pair]<<"  "<<tdoa[pol][pair]<<endl;

     }// pair loop

    //}// pol loop
// for(int i=0; i<12; i++){
//   threads[i].join();
// }

}
   
void foo(int n){
    //usleep(1000*1000);
    double a = 10;
    cout<<"foo:  "<<n<<endl;
}

void real_cross_correlation(int pol, int pair, int i, int j){

        double Fine_re[16385] = {0};            
        double Fine_im[16385] = {0};

                       for (Int_t k=1; k<HALF_WINDOW_WIDTH; k++){//0 means DC component
            Fine_re[k] = Correlator_re[pol][pair][k]; 
                   
            Fine_im[k] = Correlator_im[pol][pair][k];
            //cout<<ant_re[i][k]<<endl;
            
            //if((pol==0)&&(pair==5)&&(k==100)) cout<<Fine_re[k]<<endl;
                
                                                    }//k 


          int which_fft = pol*6+pair;

          fft_back[which_fft]->SetPointsComplex(Fine_re, Fine_im);
          fft_back[which_fft]->Transform();
          Correlator_time[which_fft] = fft_back[which_fft]->GetPointsReal();


//if((pol==0)&&(pair==5)) cout<<Correlator_time[0]<<endl;


         double factor = sqrt(Sum_x[i])*sqrt(Sum_x[j])*WINDOW_SPACING;
         double max = 0;
         int mm;

                       for (Int_t k=0; k<32768; k++){//32768
                        //cout<<Correlator_time[k]<<endl;
                correlation_pair[pol][pair][k] = Correlator_time[which_fft][k]/factor;
                if (correlation_pair[pol][pair][k] > max){max = correlation_pair[pol][pair][k]; mm = k;}

                                                    }//k 

//cout<<max<<"  "<<mm<<endl;

if(mm>16383) {mm -= 32768;}
//if(mm>149) {mm -= 300;}



}	
	
	


int filter_time(double filtered_x[1500], double new_x[1500], int Max_i){

//int Max_i = 0;
//double Max = 0.;

           //for (int i=0;i<1500;i++)
	//{

	  //if(fabs(x[i])>Max) {Max = fabs(x[i]);  Max_i = i;} 	  

	//}
	//cout<<Max_i<<endl;
	
           for (int i=0;i<1500;i++)
	{

   	       new_x[i] = ((i<Max_i-HALF_WINDOW_SPACING)||(i>Max_i+HALF_WINDOW_SPACING))?0:(0.42-0.5*cos(2*PI*(i-Max_i+HALF_WINDOW_SPACING)/WINDOW_SPACING)+0.08*cos(2*PI*(i-Max_i+HALF_WINDOW_SPACING)/HALF_WINDOW_SPACING))*filtered_x[i];
   	       //new_x[i] = ((i<Max_i-HALF_WINDOW_SPACING)||(i>Max_i+HALF_WINDOW_SPACING))?0:filtered_x[i];
   	       //remain_x[i] = ((i<Max_i-HALF_WINDOW_SPACING)||(i>Max_i+HALF_WINDOW_SPACING))?x[i]:0;
   	       //x[i] = ((i<Max_i-100)||(i>Max_i+100))?0:x[i];
 	                               
	}

           //for (int i=1500;i<4096;i++)
	//{
   	       //new_x[i] = 0;	                               
	//}
	
return Max_i;

}



int filter_N = 1500;
Int_t filter_n_size = filter_N+1;
TVirtualFFT *filter_fft_forward = TVirtualFFT::FFT(1, &filter_n_size,"R2C ES K");
TVirtualFFT *filter_fft_back = TVirtualFFT::FFT(1, &filter_N, "C2RBACKWARD M K");


int get_FFT(int ch){


   filter_fft_forward->SetPoints(x[ch]);
   filter_fft_forward->Transform();
   filter_fft_forward->GetPointsComplex(Re_ant[ch],Im_ant[ch]);



 for (Int_t k=0; k<751; k++){//0 means DC component // from 185-350 MHz
       if((k<222)||(k>420)) {Re_ant[ch][k] = 0.; Im_ant[ch][k] = 0;}
                             }

   double FreqBin = 5/6.;
   int WeakNotchHalfWidth = 3./FreqBin;  //MHz +-
   int StrongNotchHalfWidth = 6./FreqBin;  //MHz +-


  const int NKnownCW = 11;
  double   CWPeak[ NKnownCW ] = {
    165./FreqBin, 187.5/FreqBin, 200./FreqBin, 217.5/FreqBin, 244.167/FreqBin, 249.167/FreqBin, 262./FreqBin, 273.5/FreqBin, 276.65/FreqBin, 299.165/FreqBin, 313.6/FreqBin  //in cal-pulser, not found in forced trig
      //187.5, 272.5 are the strongest
  };
  int CWPeakBin[ NKnownCW ]={0};
  for(int j=0;j<NKnownCW;j++) CWPeakBin[j] = TMath::Nint( CWPeak[j] );
  
  for(int i=0; i<NKnownCW; i++){
//cout<<"(Filter) notch: "<< CWPeak[i]<<endl;
   for(int n= -WeakNotchHalfWidth; n <= WeakNotchHalfWidth; n++){
    //cout<<"Notch: "<< CWPeak[i]+n <<endl;
    Re_ant[ch][ CWPeakBin[i]+n ] *= 0.;
    Im_ant[ch][ CWPeakBin[i]+n ] *= 0.;
   }
  }
}


int filter_fre(int ch, double &SNR, bool &multi_pulses){


   filter_fft_forward->SetPoints(x[ch]);
   filter_fft_forward->Transform();
   filter_fft_forward->GetPointsComplex(Re_ant[ch],Im_ant[ch]);

                             
	Re_ant[ch][0] = 0;
	Im_ant[ch][0] = 0;

 	// filter out oscillator
 	for (Int_t k=0; k<751; k++){//0 means DC component // from 185-350 MHz
       if((k<222)||(k>420)) {Re_ant[ch][k] = 0.; Im_ant[ch][k] = 0;}
                             }	

  	filter_fft_back->SetPointsComplex(Re_ant[ch],Im_ant[ch]);
  	filter_fft_back->Transform();
  	fft_back_x = filter_fft_back->GetPointsReal();


	for(int i=0;i<1500;i++){
		zero_correct_x[ch][i] = fft_back_x[i]/1500.;
  						}





 for (Int_t k=0; k<751; k++){//0 means DC component // from 185-350 MHz
       if((k<222)||(k>420)) {Re_ant[ch][k] = 0.; Im_ant[ch][k] = 0;}
                             }
/*
double FreqBin = 5/6.;
   int WeakNotchHalfWidth = 3./FreqBin;  //MHz +-
   int StrongNotchHalfWidth = 6./FreqBin;  //MHz +-


  const int NKnownCW = 11;
  double   CWPeak[ NKnownCW ] = {
    165./FreqBin, 187.5/FreqBin, 200./FreqBin, 217.5/FreqBin, 244.167/FreqBin, 249.167/FreqBin, 262./FreqBin, 273.5/FreqBin, 276.65/FreqBin, 299.165/FreqBin, 313.6/FreqBin  //in cal-pulser, not found in forced trig
      //187.5, 272.5 are the strongest
  };
  int CWPeakBin[ NKnownCW ]={0};
  for(int j=0;j<NKnownCW;j++) CWPeakBin[j] = TMath::Nint( CWPeak[j] );
  
  for(int i=0; i<NKnownCW; i++){
//cout<<"(Filter) notch: "<< CWPeak[i]<<endl;
   for(int n= -WeakNotchHalfWidth; n <= WeakNotchHalfWidth; n++){
    //cout<<"Notch: "<< CWPeak[i]+n <<endl;
    Re_ant[ch][ CWPeakBin[i]+n ] *= 0.;
    Im_ant[ch][ CWPeakBin[i]+n ] *= 0.;
   }
  }
*/

  filter_fft_back->SetPointsComplex(Re_ant[ch],Im_ant[ch]);
  filter_fft_back->Transform();
  fft_back_x = filter_fft_back->GetPointsReal();



int Max_i = 0;
double Max = 0.;

for(int i=0;i<1500;i++){
	filtered_x[ch][i] = fft_back_x[i]/1500.;
	//cout<<filtered_x[ch][i]<<endl;
  						}

for(int i=200;i<1300;i++){
	if(fabs(filtered_x[ch][i])>Max) {
		Max = fabs(filtered_x[ch][i]);  Max_i = i; //cout<<"i: "<<i<<"  max:  "<<Max<<endl;
		}
  	}

double power_sum = 0.;
int count_over_30 = 0;
int half_pulse_width = 40;
/* 						
for(int i=200;i<1300;i++){
	if((i<Max_i-half_pulse_width)||(i>Max_i+half_pulse_width)){
	    double square_x = pow(filtered_x[ch][i],2);
		power_sum += square_x;
		//if(square_x>400) {count_over_30++;}
		}
						}
power_sum = sqrt(power_sum/(1100.-2.*half_pulse_width-1.));
SNR = Max/power_sum;
*/

SNR = Max;

//
double pulse_threshold = 0.6*Max;
if(pulse_threshold<30.){pulse_threshold = 30.;}

for(int i=200;i<1300;i++){
	if((i<Max_i-half_pulse_width)||(i>Max_i+half_pulse_width)){
		if(fabs(filtered_x[ch][i])>pulse_threshold) {count_over_30++;}
		}
						}
						
multi_pulses = false;
if(count_over_30>=5){multi_pulses=true;}
  						
return Max_i;

}

//theta:-pi/2.~pi/2;  phi:0~2*pi, x-axis direction is phi = 0;
//xyz: antenna2 - antenna1
//return: relative time between signal arrive antenna2 and antenna1, if it is minus, means signal arrive antenna2 first.
Double_t expectedTimeDiff_angle(double theta, double phi, double *xyz)
         {
   double n[3] = {0};

   phi = phi*TMath::DegToRad();
   theta = theta*TMath::DegToRad();
   
   n[0] = (-1)*TMath::Cos(theta)*TMath::Cos(phi);
   n[1] = (-1)*TMath::Cos(theta)*TMath::Sin(phi);
   n[2] = (-1)*TMath::Sin(theta);

   double dt = (n[0]*xyz[0] + n[1]*xyz[1] + n[2]*xyz[2])/(Speed_Of_Light);


   // s to ns
   return dt/1e-9;
         }


Double_t expectedTimeDiff_angle_with_R(double theta, double phi, double R, double *antenna_i, double *antenna_j)
         {
//cout<<"output:   "<<R<<endl;
   theta = theta * TMath::DegToRad();
   phi = phi * TMath::DegToRad();
   double drone_xyz[3] = {0};
   //R=10000000;
   drone_xyz[0] = TMath::Cos(theta)*TMath::Cos(phi)*R;
   drone_xyz[1] = TMath::Cos(theta)*TMath::Sin(phi)*R;
   drone_xyz[2] = TMath::Sin(theta)*R;
   return expectedTimeDiff_coord(drone_xyz, antenna_i, antenna_j);
         }
         
                  
         
Double_t expectedTimeDiff_coord(double *drone_xyz, double *antenna_i, double *antenna_j)
         {
         double sum_i = 0;
         double sum_j = 0;
    for(int k=0; k<3; k++){
        sum_i += pow((drone_xyz[k] - antenna_i[k]), 2);
        sum_j += pow((drone_xyz[k] - antenna_j[k]), 2);
        }

   double dt = (sqrt(sum_j) - sqrt(sum_i))/(Speed_Of_Light);
   //cout<<dt/1e-9<<"wtf?   "<<endl;


   // s to ns
   return dt/1.0e-9;
         }

bool Is_saturated()
     {
      bool is_saturated = false;
      for(int i=0; i<8; i++){
        for(int j=0; j<1500; j++){
          if(x[i][j] > 490) {is_saturated = true;}
        }
      }
      return is_saturated;
     }
         
         
         

template <typename T>
T** new_Array2D(int row, int col)  
{  
    int size = sizeof(T);  
    int point_size = sizeof(T*);  
    //先申请内存，其中sizeof(T*) * row表示存放row个行指针   
    T **arr = (T **) malloc(point_size * row + size * row * col);  
    if (arr != NULL)  
    {     
        T *head = (T*)((long)arr + point_size * row);  
        for (int i = 0; i < row; ++i)  
        {  
            arr[i] =  (T*)((long)head + i * col * size);  
            for (int j = 0; j < col; ++j)  
                new (&arr[i][j]) T{}; // {} means initial to 0
        }  
    }  
    return (T**)arr;  
}  
//释放二维数组   
template <typename T>  
void delete_Array2D(T **arr, int row, int col)  
{  
    for (int i = 0; i < row; ++i)  
        for (int j = 0; j < col; ++j)  
            arr[i][j].~T();  
    if (arr != NULL)  
        free((void**)arr);  
}


void reconstruction(int pol, double event_theta[2], double event_phi[2], double x_value[2]){

double R = 10000;

double inter_max[2] = {0};
   //0 for H-pol 1 for V-pol
   //for(int pol = 0; pol<2; pol++){
      
   	int i = 0;
    int j = 2;
    
double inter_map[56][100] = {0};//theta: -6~50; phi: -76~24

     double max[3] = {0};
     int max_theta[3] = {0};
     int max_phi[3] = {0};
  
   for(int pair = 0; pair <6; pair++){

  //threads[pol*6+pair].join();

    switch(pair){
      case 0:
        i = 0+pol; j = 2+pol;
        break;
      case 1:
        i = 0+pol; j = 4+pol;
        break;
      case 2:
        i = 0+pol; j = 6+pol;
        break;
      case 3:
        i = 2+pol; j = 4+pol;
        break;
      case 4:
        i = 2+pol; j = 6+pol;
        break;
      case 5:
        i = 4+pol; j = 6+pol;
        break;

    }



    for(int i_theta=0; i_theta<56; i_theta++){
      for(int i_phi=0; i_phi<100; i_phi++){
          double theta = i_theta - 1.;//cyc
          double phi = -76 + i_phi;
          double delta_t = expectedTimeDiff_angle_with_R(theta, phi, R, coor_ant[i], coor_ant[j]) + (-delay_ant[j] + delay_ant[i])/1000.;
          
          delta_t = delta_t - (peak_position[j] - peak_position[i])*0.8;
          delta_t = delta_t/(0.8*WINDOW_SPACING/32768.);
          if(delta_t<0) delta_t += 32768;
          if(delta_t<0||delta_t>32768) continue; //cout<<"delta_t error!  "<<delta_t<<endl;

          //cout<<delta_t<<endl;
          //cout<<correlation_pair[pol][pair][(int) delta_t]<<endl;
          inter_map[i_theta][i_phi] += correlation_pair[pol][pair][(int) round(delta_t)];
          //inter_map[theta+90][phi+90] += theta+phi;
          //h2->SetBinContent(theta+90, phi+90, correlation_pair[pol][pair][(int) delta_t]);
        }
      }
      
    if(pair == 5) {


     for (int ci=0; ci<56; ci++){
      for (int cj=0; cj<100; cj++){
          if(inter_map[ci][cj] > max[0]){
          		if((abs(ci-max_theta[0])<=2)&&(abs(cj-max_phi[0])<=1)){
          			max[0] = inter_map[ci][cj]; max_theta[0] = ci; max_phi[0] = cj;
          		}
          		else if((abs(ci-max_theta[1])<=2)&&(abs(cj-max_phi[1])<=1)){
          			max[1] = max[0]; max_theta[1] = max_theta[0]; max_phi[1] = max_phi[0];
          			max[0] = inter_map[ci][cj]; max_theta[0] = ci; max_phi[0] = cj;
          		}
          		else{
          			max[2] = max[1]; max_theta[2] = max_theta[1]; max_phi[2] = max_phi[1];
          			max[1] = max[0]; max_theta[1] = max_theta[0]; max_phi[1] = max_phi[0];
          			max[0] = inter_map[ci][cj]; max_theta[0] = ci; max_phi[0] = cj;
          		}
          		//cout<<a<<"   "<<max_theta[0]<<"   "<<max_theta[1]<<"    "<<max_phi[0]<<"   "<<max_phi[1]<<endl;
          }
     }
   }

     }
     
     	

	}//pair loop





//////////////////////////////////////////////////////////////////////
//fine search

double fine_inter_map[3][40][20] = {0};
int fine_max_theta = 0;
int fine_max_phi = 0;
double fine_max = 0;
int max_w = 0;

   for(int pair = 0; pair <6; pair++){

    switch(pair){
      case 0:
        i = 0+pol; j = 2+pol;
        break;
      case 1:
        i = 0+pol; j = 4+pol;
        break;
      case 2:
        i = 0+pol; j = 6+pol;
        break;
      case 3:
        i = 2+pol; j = 4+pol;
        break;
      case 4:
        i = 2+pol; j = 6+pol;
        break;
      case 5:
        i = 4+pol; j = 6+pol;
        break;

    }



	for(int w = 0; w < 3; w++){
    	for(int i_theta=0; i_theta<40; i_theta++){
      	for(int i_phi=0; i_phi<20; i_phi++){
          double theta = -2 + i_theta/10. + max_theta[w]-1;//cyc
          double phi = -1 + i_phi/10. + max_phi[w]-76;
          double delta_t = expectedTimeDiff_angle_with_R(theta, phi, R, coor_ant[i], coor_ant[j]) + (-delay_ant[j] + delay_ant[i])/1000.;
          
          delta_t = delta_t - (peak_position[j] - peak_position[i])*0.8;
          delta_t = delta_t/(0.8*WINDOW_SPACING/32768.);
          if(delta_t<0) delta_t += 32768;
          if(delta_t<0||delta_t>32768) continue; 
          fine_inter_map[w][i_theta][i_phi] += correlation_pair[pol][pair][(int) round(delta_t)];
        }
      }
      
    if(pair == 5) {


     for (int ci=0; ci<40; ci++){
      for (int cj=0; cj<20; cj++){
          if(fine_inter_map[w][ci][cj] > fine_max){
          			fine_max = fine_inter_map[w][ci][cj]; fine_max_theta = ci; fine_max_phi = cj; max_w = w;
          }
       }
     }
     }
     }//w loop  
     	

	}//pair loop
	
     inter_max[pol] = fine_max;
     
     event_theta[pol] = -2 + fine_max_theta/10. + max_theta[max_w]-1;//cyc
     event_phi[pol] = -1 + fine_max_phi/10. + max_phi[max_w]-76;
     x_value[pol] = inter_max[pol];

	
	//}//pol loop

}






//result: correlation_pair[pair][k]; said pair antenna i and j, correlation_pair[pair][k] is the X-correlation value when we move antenna j signal earier k (just shift toward left side).
double fix_angle_cross_correlation(int pol, double in_theta, double in_phi){	

   double cross_correlator = 0;
   double R = 10000;
   
   double sum_factor[8] = {0.};
   for(int ant=0; ant<4; ant++){
   		int ch = 2*ant + pol;
   		for(int k=0; k<1500; k++){
   			sum_factor[ch] += remain_x[ch][k] * remain_x[ch][k];
   		}
   		sum_factor[ch] = sqrt(sum_factor[ch]);
   }


   //0 for H-pol 1 for V-pol
   //for(int pol = 0; pol<2; pol++){
      
   	int i = 0;
    int j = 2;
    
   for(int pair = 0; pair <6; pair++){

    switch(pair){
      case 0:
        i = 0+pol; j = 2+pol;
        break;
      case 1:
        i = 0+pol; j = 4+pol;
        break;
      case 2:
        i = 0+pol; j = 6+pol;
        break;
      case 3:
        i = 2+pol; j = 4+pol;
        break;
      case 4:
        i = 2+pol; j = 6+pol;
        break;
      case 5:
        i = 4+pol; j = 6+pol;
        break;

    } 	
    	
    	double delta_t = expectedTimeDiff_angle_with_R(in_theta, in_phi, R, coor_ant[i], coor_ant[j]) + (-delay_ant[j] + delay_ant[i])/1000.;
    	int shift_points = round(delta_t/0.8);//+r3->Integer(100);
    	//cout<<"shift_points:   "<<shift_points<<endl;
    	
    	double sum = 0.;
    	for(int k=0; k<1500; k++){
    		int l = k + shift_points;
    		if(l<0||l>=1500){continue;}
   			sum += remain_x[i][k] * remain_x[j][l];
   		}
   		double factor = sum_factor[i] * sum_factor[j];
   		
   		cross_correlator += sum/factor;
   		
   }// pair loop
   //}// pol loop
   
   return cross_correlator;
   
}
   		
    	
double multi_pulse_cross_correlation(int pol, double in_theta, double in_phi){
   	
	for(int ant=0; ant<4; ant++){
		int ch = 2*ant + pol;
		//peak_position[ch]=502;
		//cout<<peak_position[ch]<<endl;
		for(int k =0; k<1500; k++){
			remain_x[ch][k] = ((k<peak_position[ch]-HALF_WINDOW_SPACING)||(k>peak_position[ch]+HALF_WINDOW_SPACING))?zero_correct_x[ch][k]:0;//filtered_x[ch][k]:0;
			}
		}
	double cross_correlation = fix_angle_cross_correlation(pol, in_theta, in_phi);
	
	cout<<"cross_correlation:  "<<cross_correlation<<endl;
	
	return cross_correlation;
}

  
// if pol = V-pol, it gives cross correlation of H-pol at V-pol reconstructed direction;	
double cross_pol_cross_correlation(int pol, double in_theta, double in_phi){
   	
	int cross_pol = (pol+1)%2;
	
	for(int ant=0; ant<4; ant++){
		int old_ch = 2*ant + pol;
		int ch = 2*ant + cross_pol;
		for(int k =0; k<1500; k++){
			// if pol = V-pol, modify H-pol's data with V-pol's peak position.
			remain_x[ch][k] = ((k<peak_position[old_ch]-HALF_WINDOW_SPACING)||(k>peak_position[old_ch]+HALF_WINDOW_SPACING))?0:filtered_x[ch][k];
			}
		}
	double cross_correlation = fix_angle_cross_correlation(cross_pol, in_theta, in_phi);
	
	return cross_correlation;
}   	



double spectrum_cut(int pol, int run_pos){

	int total_bin = 0;
	int total_pass = 0;
	int _k = 0;
	int k_pps = 0;
	for(int ant=0; ant<4; ant++){
		int ch = ant*2 + pol;
		for(int k=0; k<20; k++){
			
			_k = k+22;
			k_pps = k;
			//cout<<"k: "<<k<<"  k_pps: "<<k_pps<<endl;
			//if(k>20) break;
			
			//double power = (Re_ant[ch][_k]*Re_ant[ch][_k] + Im_ant[ch][_k]*Im_ant[ch][_k]);
			double power = (ant_re[ch][_k]*ant_re[ch][_k] + ant_im[ch][_k]*ant_im[ch][_k]);
			double pps_power = pps_ave_power[run_pos][ch][k_pps];
			
			//cout<<pps_power<<"   "<<power<<endl;
			
			//cout<<pps_power<<"   "<<power<<endl;
			//trigger_spectrum[ch][k] = -1000;
			if((pps_power == 0)||(power == 0)){
				trigger_spectrum[ch][k] = -1000;
				//cout<<pps_power<<"   "<<power<<endl;
				//continue;
			}
			
			else{
				
			    // 2.*(T*1.0e06)*(im_full[l]*im_full[l]+re_full[l]*re_full[l])/(n_size*n_size*50.*1000.);
			    // cancel out other terms, only factor of 10 left;
			    //power = 10*log10(10*power);
				//double different = power - pps_power;
				power = 10*power;
				double different = 10*log10(power/pps_power);
				//cout<<pps_power<<"   "<<power<<endl;
				trigger_spectrum[ch][k] = different;
				total_bin++;
				//
				//3dB higher
				if(different>3){total_pass++;}
				
				//if(different>10){cout<<k<<" different:  "<<different<<endl;}
			}
			//cout<<"total_bin: "<<total_bin<<endl;
			
		}
		//cout<<"trigger_spectrum:  "<<trigger_spectrum[ch][10]<<endl;
	}
	cout<<"total bins, total pss:  "<<total_bin<<"  "<<total_pass<<endl;
	return total_pass/((double) total_bin);	

}


void SetTreeBranches(TTree* tr){
	
	tr->Branch("run",&gRunNum, "run/I");
	tr->Branch("event",&gEventNum, "event/I");
	tr->Branch("is_pps", &is_pps, "is_pps/I");
	tr->Branch("eventTime",&eventTime, "eventTime/D"); //PC timestamp
	tr->Branch("iGood",&iGood, "iGood/I");
	tr->Branch("V_bits",&V_bits, "V_bits/I");
	tr->Branch("H_bits",&H_bits, "H_bits/I");
	tr->Branch("rough_cross_correlation", rough_cross_correlation, "rough_cross_correlation[2]/D");
	tr->Branch("SNR", rough_SNR, "SNR[2]/D");
	tr->Branch("trigger_spectrum", trigger_spectrum, "trigger_spectrum[8][199]/D");
	tr->Branch("pass_trigger_ratio", pass_trigger_ratio, "pass_trigger_ratio[2]/D");
	tr->Branch("power_frequency", power_frequency, "power_frequency[199]/D");
	tr->Branch("multi_pulses", multi_pulses, "multi_pulses[2]/I");
	tr->Branch("reconstruct_cross_correlation", reconstruct_cross_correlation, "reconstruct_cross_correlation[2]/D");
	tr->Branch("remain_pulse_cross_correlation", remain_pulse_cross_correlation, "remain_pulse_cross_correlation[2]/D");
	tr->Branch("other_pol_cross_correlation", other_pol_cross_correlation, "other_pol_cross_correlation[2]/D");
	

}

void initial_data(){
	for(int pol=0; pol<2; pol++){
		rough_cross_correlation[pol] = -1000;
		reconstruct_cross_correlation[pol] = -1000;
		remain_pulse_cross_correlation[pol] = -1000;
		other_pol_cross_correlation[pol] = -1000;
	}
}
