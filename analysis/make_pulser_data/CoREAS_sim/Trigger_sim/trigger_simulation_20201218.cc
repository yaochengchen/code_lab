#include "../read_pulse.cc"

void filter_fre(double x[1500], bool save_100MHz);

void trigger_simulation(){



	TRandom *r3 = new TRandom3();


	TFile* f0 = new TFile("efficiency_vs_strength.root");
  
  	TGraph * gr_eff_vs_strength[1000][8][8] = {NULL};
	
	for(int step=300; step<1300; step++){
		for(int board=1; board<9; board++){
			for(int band=1; band<9; band++){
				char *txtname = Form("gr_Step%d_Board%d_Band%d", step, board, band); 
				gr_eff_vs_strength[step-300][board-1][band-1] = (TGraph*)f0->Get(txtname); 
			}
		}
	}
	
	


TFile *ntuple_file = new TFile("./ntuple_20201218.root", "update");

TNtuple * ntuple = (TNtuple*) ntuple_file->Get("ntuple");

float nt_run;
float nt_event;
float nt_time;
float drn_azimuth;
float drn_elevation;
ntuple->SetBranchAddress("run",&nt_run);
ntuple->SetBranchAddress("event",&nt_event);
ntuple->SetBranchAddress("time",&nt_time);
ntuple->SetBranchAddress("drn_azimuth",&drn_azimuth);
ntuple->SetBranchAddress("drn_elevation",&drn_elevation);

float rel_x;
float rel_y;
float rel_z;
ntuple->SetBranchAddress("rel_x",&rel_x);
ntuple->SetBranchAddress("rel_y",&rel_y);
ntuple->SetBranchAddress("rel_z",&rel_z);

float V_peak;
float flag;
ntuple->SetBranchAddress("V_p",&V_peak);
ntuple->SetBranchAddress("flag",&flag);

float distance;
float attenuation;
//ntuple->SetBranchAddress("distance",&distance);
ntuple->SetBranchAddress("attenuation",&attenuation);

float trig_eff;
auto trig_eff_branch = ntuple->Branch("trig_eff",&trig_eff,"trig_eff/F");
float V_peak_4 = 0.;
auto V_p_branch = ntuple->Branch("V_p_4",&V_peak_4,"V_p_4/F");
	
const Int_t nt_nentries = (Int_t)ntuple->GetEntries() - 1;


    
  
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
  
  		  TChain *tr_event = new TChain("t");
  		  vector<Long64_t> fvNEntries;
  		  tr_event->AddFile("/media/cyc/For_Linux/TAROGE4_DATA/data/20201218/run00015888.root");
		  fvNEntries.push_back(tr_event->GetEntries());
  		  tr_event->AddFile("/media/cyc/For_Linux/TAROGE4_DATA/data/20201218/run00015889.root"); // /20201218/run00015889
  		  fvNEntries.push_back(tr_event->GetEntries());
  		  
  
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
		  
		  
		  
  		  TChain *pps_event = new TChain("t");
  		  pps_event->AddFile("/media/cyc/For_Linux/TAROGE4_DATA/data/20201218/run00015890.root");//
  		  pps_event->AddFile("/media/cyc/For_Linux/TAROGE4_DATA/data/20201218/run00015891.root");
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
	
	double eff_list[10000] = {0.};
	double V_p_list[10000] = {0.};
	int pulser_count = 0;
	for(int m=0; m<nt_nentries; m++){
      	ntuple->GetEntry(m);
      	
      	distance = sqrt(rel_x*rel_x+rel_y*rel_y+rel_z*rel_z);
      	if(nt_time<962||nt_time>1586){continue;}//198 405  
      	//if((attenuation<72)||(attenuation>74)){continue;}
      	//if(attenuation<76||attenuation>78){continue;}
      	
      	pulser_count++;
      	
      	
      	int t_event_number = (nt_run-15888)*fvNEntries[0] + nt_event; //15888
      	tr_event->GetEntry(t_event_number);

		for(int i=0;i<1500;i++)
			{
	  		t[i]=i*0.8;
	  		template_x[0][i]= (float)T1H[i]*500/32512.;
	  		template_x[2][i]= (float)T2H[i]*500/32512.;
	  		template_x[4][i]= (float)T3H[i]*500/32512.;
	  		template_x[6][i]= (float)T4H[i]*500/32512.;
	  		template_x[1][i]= (float)T1V[i]*500/32512.;
	  		template_x[3][i]= (float)T2V[i]*500/32512.;
	  		template_x[5][i]= (float)T3V[i]*500/32512.;
	  		template_x[7][i]= (float)T4V[i]*500/32512.;	  	  
			}
			
		double max[8] = {0.};
		double max_i[8] = {0.};
		
		for(int i=0;i<8;i++)
			{
			//filter_fre(template_x[i], 0);
			for(int k=0;k<1500;k++)
				{if(fabs(template_x[i][k])>max[i]){max_i[i] = k; max[i] = fabs(template_x[i][k]);}
				}
			}
		
		cout<<max[0]<<"   "<<max[2]<<"   "<<max[4]<<"   "<<max[6]<<endl;
		V_peak_4 = max[6];
		V_p_list[m] = V_peak_4;
		V_p_branch->Fill();
		//if(V_p_list[m]>100) {cout<<V_p_list[m]<<endl;}
		//if(max[0]<70){continue;}

int trigger_count = 0;
int pps_count = 0;
		  
 int event_number = pps_event->GetEntries();
 //cout<<"event_number:  "<<event_number<<endl;

 for(int j=0;j<event_number;j+=50){
		
      pps_event->GetEntry(j);
      //cout<<gTrigBits->CountBits()<<endl;
      if(gTrigBits->CountBits() != 0) {continue;}
      //cout<<j<<endl;
      //pps_count++; continue;
      //cout<<j<<endl;
		
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
			for(int k=0;k<1500;k++){
				if((k<(max_i[i]-31))||(k>(max_i[i]+31))){
				//	x[i][k] = 0.;
					}
				else{
					x[i][k] += template_x[i][k];
					}
				}
				
			}
      		
		
		
		int ant_count = 0;
		for(int board=0;board<4;board++){
				double sub_band_power[8] = {0.};
				get_sub_band_strength(x[board*2], sub_band_power);
				
				//cout<<endl;
				int band_count = 0;
				for(int band=1; band<9; band++){
					//cout<<sub_band_power[band-1]<<"  ";
					double vpp = sqrt(sub_band_power[band-1]/template_sub_band_power[band-1]) * template_v_pp;
					//cout<<board<<"  "<<band<<"  "<<vpp<<endl;
					double efficiency = gr_eff_vs_strength[460][board][band-1]->Eval(vpp);
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
				trigger_count++;
			}
			else{
				//cout<<"no trigger!!"<<endl;
			}
			//cout<<endl<<ant_count<<endl;
			pps_count++;
	} // j loop
	
	trig_eff = trigger_count/18.;
	cout<<attenuation<<"   "<<drn_elevation<<"   "<<"trigger_eff:  "<<trig_eff<<endl;
	eff_list[m] = trig_eff;
	
	trig_eff_branch->Fill();
	
	cout<<pps_count<<"  "<<trigger_count<<"   "<<pulser_count<<endl;
	} // m loop
	
	
	cout<<"fuck"<<endl;
	

	
	ntuple->Write("", TObject::kOverwrite); 
	ntuple_file->Close();
	
	//cout<<pps_count<<"  "<<trigger_count<<"   "<<pulser_count<<endl;


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
