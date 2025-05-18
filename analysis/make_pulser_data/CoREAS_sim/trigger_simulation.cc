#include "./read_pulse.cc"

void filter_fre(double x[1500], bool save_100MHz);

void trigger_simulation(){



	TRandom *r3 = new TRandom3();


	TFile* f0 = new TFile("./Trigger_sim/efficiency_vs_strength.root");
  
  	TGraph * gr_eff_vs_strength[1000][8][8] = {NULL};
	
	for(int step=300; step<1300; step++){
		for(int board=1; board<9; board++){
			for(int band=1; band<9; band++){
				char *txtname = Form("gr_Step%d_Board%d_Band%d", step, board, band); 
				gr_eff_vs_strength[step-300][board-1][band-1] = (TGraph*)f0->Get(txtname); 
			}
		}
	}
	


    
  
  bool bGood;
  Int_t timeStamp_FPGA;
  double eventTime = 0;
  Double_t t[1500];
  Short_t T1H[1500];
  Short_t T2H[1500];
  Short_t T3H[1500];
  Short_t T4H[1500];
  Short_t T1V[1500];
  Short_t T2V[1500];
  Short_t T3V[1500];
  
  Short_t T4V[1500];
  
  TBits * gTrigBits = NULL;
  TBits * gOverVoltBits = NULL;
  

		  
		  
		  
  		  TChain *pps_event = new TChain("t");
  		  pps_event->AddFile("/home/cyc/software/TAROGE-4_analysis/make_pulser_data/CoREAS_Sim/make_pulser_data_57000.root");//
  		  //pps_event->AddFile("/media/cyc/For_Linux/TAROGE4_DATA/data/20201218/run00015891.root");
          pps_event->SetBranchAddress("bGood",&bGood);
		  pps_event->SetBranchAddress("timeStamp_FPGA",&timeStamp_FPGA);
		  pps_event->SetBranchAddress("eventTime", &eventTime);
  		  pps_event->SetBranchAddress("triggerBits", &gTrigBits);
  		  pps_event->SetBranchAddress("overVoltBits", &gOverVoltBits);
		  pps_event->SetBranchAddress("T1H",T1H);
		  pps_event->SetBranchAddress("T2H",T2H);
		  pps_event->SetBranchAddress("T3H",T3H);
		  pps_event->SetBranchAddress("T4H",T4H);
		  pps_event->SetBranchAddress("T1V",T1V);
		  pps_event->SetBranchAddress("T2V",T2V);
		  pps_event->SetBranchAddress("T3V",T3V);
		  pps_event->SetBranchAddress("T4V",T4V);


	
	
	
	read_pulse();
	
	
	
	double x[8][4096] = {0.};
	double template_x[8][4096] = {0.};

		  
 int event_number = pps_event->GetEntries();
 //cout<<"event_number:  "<<event_number<<endl;

 int pps_count = 0;
 for (int j=0;j<event_number;j++){

      pps_event->GetEntry(j);


		for(int i=0;i<1500;i++)
			{
	  		x[0][i]= (float)T1H[i]*500/32512.;
	  		x[2][i]= (float)T2H[i]*500/32512.;
	  		x[4][i]= (float)T3H[i]*500/32512.;
	  		x[6][i]= (float)T4H[i]*500/32512.;
	  		x[1][i]= (float)T1V[i]*500/32512.;
	  		x[3][i]= (float)T2V[i]*500/32512.;
	  		x[5][i]= (float)T3V[i]*500/32512.;
	  		x[7][i]= (float)T4V[i]*500/32512.;	  	  
			}      
      			
		for(int i=0;i<8;i++)
			{
			filter_fre(x[i], 0);				
			}
      		

		
		
		bool is_trigger = false;
		for(int pol=0; pol<2; pol++){
			int ant_count = 0;
			for(int board=0;board<4;board++){
				double sub_band_power[8] = {0.};
				get_sub_band_strength(x[board*2+pol], sub_band_power);
				
				//cout<<endl;
				int band_count = 0;
				for(int band=1; band<9; band++){
					//cout<<sub_band_power[band-1]<<"  ";
					double vpp = sqrt(sub_band_power[band-1]/template_sub_band_power[band-1]) * template_v_pp;
					//cout<<board<<"  "<<band<<"  "<<vpp<<endl;
					double efficiency = gr_eff_vs_strength[350][board+pol*4][band-1]->Eval(vpp);
					double random_number = r3->Rndm();
					bool band_trigger = false;
					if(random_number<efficiency){
						band_trigger = true;
					}
					// cyc caution! probability of trigger by CW 
					random_number = r3->Rndm();
					if(random_number<0.3){
						//band_trigger = true;
					}
					if(band_trigger == true){
						band_count ++;
					}
				}
				if(band_count>=6){
					ant_count++;
				}
				//cout<<band_count<<"  ";
			}
			if(ant_count>=4){
				//cout<<"trigger!!"<<endl;
				is_trigger = true;
			}
		}
		
		if(is_trigger){pps_count++;}
		

	} // j loop
		cout<<pps_count<<"  "<<event_number<<endl;
	




}






int N = 1500;
   TVirtualFFT *fft_forward = TVirtualFFT::FFT(1, &N,"R2C");
      
   TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &N, "C2RBACKWARD M K");
   
void filter_fre(double x[1500], bool save_100MHz){


   
   	
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
