#include "../read_pulse.cc"

void filter_fre(double x[1500], bool save_100MHz);

void trigger_simulation(){




TFile *ntuple_file = new TFile("./ntuple.root");

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

const Int_t nt_nentries = (Int_t)ntuple->GetEntries() - 1;

TH1F *wyn = new TH1F("a", "a", 40, 66, 95);
	
	int flag_count = 0;
	double average_V_peak = 0;
	double V_list[10000] = {0.};
	double E_list[10000] = {0.};
	int count_list = 0;
	for(int m=0; m<nt_nentries; m++){
		//cout<<m<<"  "<<nt_nentries<<endl;
      	ntuple->GetEntry(m);
      	distance = sqrt(rel_x*rel_x+rel_y*rel_y+rel_z*rel_z);
      	
      	if(nt_time<198||nt_time>405||distance<2200e-15){continue;}
      	//if(((int) nt_run)!=16946){continue;}
      	//cout<<flag<<endl;
      	if(flag == 1){
      		cout<<flag_count<<"   "<<average_V_peak/flag_count<<endl;
      		V_list[count_list] = average_V_peak/flag_count;
      		E_list[count_list] = flag_count/5.;
      		count_list ++;
      		
      		average_V_peak = attenuation;//V_peak;
      		flag_count = 1;
      	}
      	else{
      		average_V_peak += attenuation;//V_peak;
      		flag_count++;
      	}
      	
      	wyn->Fill(attenuation);

	}
	//cout<<V_list[1000]<<endl;
	
	TGraph * cyc = new TGraph(count_list-2, &V_list[1], &E_list[1]);
	cyc->Draw("AP");
	cyc->SetMarkerStyle(22);
	
	wyn->Draw();
    
  
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
  
    	const char * res_fname =  Form("/media/cyc/For_Linux/TAROGE4_DATA/data/20190925/run%08d.root", 2112);
    	  TFile * res_file = new TFile(res_fname);
    	  TTree * Tree_res = (TTree*) res_file->Get("t");
          Tree_res->SetBranchAddress("bGood",&bGood);
		  Tree_res->SetBranchAddress("timeStamp_FPGA",&timeStamp_FPGA);
		  Tree_res->SetBranchAddress("eventTime", &eventTime);
		  Tree_res->SetBranchAddress("T1H",T1H);
		  Tree_res->SetBranchAddress("T2H",T2H);
		  Tree_res->SetBranchAddress("T3H",T3H);
		  Tree_res->SetBranchAddress("T4H",T4H);
		  Tree_res->SetBranchAddress("T1V",T1V);
		  Tree_res->SetBranchAddress("T2V",T2V);
		  Tree_res->SetBranchAddress("T3V",T3V);
		  Tree_res->SetBranchAddress("T4V",T4V);
		  



	TRandom *r3 = new TRandom3();


	TFile* f0 = new TFile("efficiency_vs_strength.root");
  
  	TGraph * gr_eff_vs_strength[500][8][8] = {NULL};
	
	for(int step=400; step<900; step++){
		for(int board=1; board<9; board++){
			for(int band=1; band<9; band++){
				char *txtname = Form("gr_Step%d_Board%d_Band%d", step, board, band); 
				gr_eff_vs_strength[step-400][board-1][band-1] = (TGraph*)f0->Get(txtname); 
			}
		}
	}
	
	
	gr_eff_vs_strength[100][0][0]->Draw();	
	for(int board=1; board<9; board++){
		for(int band=1; band<9; band++){
			gr_eff_vs_strength[100][board-1][band-1]->Draw("same");
			gr_eff_vs_strength[100][board-1][band-1]->SetLineColor(8*(board-1)+(band-1));
			gr_eff_vs_strength[400][board-1][band-1]->Draw("same");
			gr_eff_vs_strength[400][board-1][band-1]->SetLineColor(8*(board-1)+(band-1));
		}
	}
	
	
	read_pulse();
	
	
	
	double x[8][4096] = {0.};
	double template_x[8][4096] = {0.};
	
	double eff_list[10000] = {0.};
	double V_p_list[10000] = {0.};
	int pulser_count = 0;
	int pps_count = 0;
	for(int m=0; m<nt_nentries; m++){
      	ntuple->GetEntry(m);
      	if(((int) nt_run)!=2112){continue;}
      	distance = sqrt(rel_x*rel_x+rel_y*rel_y+rel_z*rel_z);
      	if(nt_time<198||nt_time>405||distance<2200e-15){continue;}//198 405
      	//if((attenuation<72)||(attenuation>74)){continue;}
      	//if(attenuation<76||attenuation>78){continue;}
      	
      	pulser_count++;
      	
      	
      	Tree_res->GetEntry(nt_event);

		for(int i=0;i<1500;i++)
			{
	  		t[i]=i*0.8;
	  		template_x[0][i]= (float)T1H[i]*500/32512;
	  		template_x[2][i]= (float)T2H[i]*500/32512;
	  		template_x[4][i]= (float)T3H[i]*500/32512;
	  		template_x[6][i]= (float)T4H[i]*500/32512;
	  		template_x[1][i]= (float)T1V[i]*500/32512;
	  		template_x[3][i]= (float)T2V[i]*500/32512;
	  		template_x[5][i]= (float)T3V[i]*500/32512;
	  		template_x[7][i]= (float)T4V[i]*500/32512;	  	  
			}
			
		double max[8] = {0.};
		double max_i[8] = {0.};
		for(int i=0;i<8;i++)
			{
			filter_fre(template_x[i], 0);
			for(int k=0;k<1500;k++)
				{if(fabs(template_x[i][k])>max[i]){max_i[i] = k; max[i] = fabs(template_x[i][k]);}
				}
			}
		//cout<<max[0]<<"   "<<max[2]<<"   "<<max[4]<<"   "<<max[6]<<endl;
		double V_peak_4 = max[6];
		V_p_list[m] = V_peak_4;
		if(V_p_list[m]>100) {cout<<V_p_list[m]<<endl;}
		//if(max[0]<70){continue;}

int trigger_count = 0;

//pps data 20210408, 18969, 18970
	int start_run = 2114;
	int end_run = 2115;
	int date = 20190925;
  for(int number = start_run; number<end_run+1; number++){


    //fname_evt =  Form( "/media/cyc/For_Linux/TAROGE4_DATA/data/20211209/run%08d.root",number);
    const char * fname_evt =  Form("/media/cyc/For_Linux/TAROGE4_DATA/data/%d/run%08d.root", date, number);
    //Form( "/media/cyc/1p9TB/TAROGE4_DATA/data/%d/run%08d.root", date, number);
    TFile *file = new TFile(fname_evt);
    TTree *Tree_Muon = (TTree*) file->Get("t");
          Tree_Muon->SetBranchAddress("bGood",&bGood);
		  Tree_Muon->SetBranchAddress("timeStamp_FPGA",&timeStamp_FPGA);
		  Tree_Muon->SetBranchAddress("eventTime", &eventTime);
  		  Tree_Muon->SetBranchAddress("triggerBits", &gTrigBits);
  		  Tree_Muon->SetBranchAddress("overVoltBits", &gOverVoltBits);
		  Tree_Muon->SetBranchAddress("T1H",T1H);
		  Tree_Muon->SetBranchAddress("T2H",T2H);
		  Tree_Muon->SetBranchAddress("T3H",T3H);
		  Tree_Muon->SetBranchAddress("T4H",T4H);
		  Tree_Muon->SetBranchAddress("T1V",T1V);
		  Tree_Muon->SetBranchAddress("T2V",T2V);
		  Tree_Muon->SetBranchAddress("T3V",T3V);
		  Tree_Muon->SetBranchAddress("T4V",T4V);
		  
 int event_number = Tree_Muon->GetEntries();

 for (int j=0;j<event_number;j+=1){

      Tree_Muon->GetEntry(j);
      if(gTrigBits->CountBits() != 0) {continue;}
      //cout<<j<<endl;

		for(int i=0;i<1500;i++)
			{
	  		x[0][i]= (float)T1H[i]*500/32512;
	  		x[2][i]= (float)T2H[i]*500/32512;
	  		x[4][i]= (float)T3H[i]*500/32512;
	  		x[6][i]= (float)T4H[i]*500/32512;
	  		x[1][i]= (float)T1V[i]*500/32512;
	  		x[3][i]= (float)T2V[i]*500/32512;
	  		x[5][i]= (float)T3V[i]*500/32512;
	  		x[7][i]= (float)T4V[i]*500/32512;	  	  
			}      
      			
		for(int i=0;i<8;i++)
			{
			filter_fre(x[i], 0);
			for(int k=0;k<1500;k++){
				if((k<(max_i[i]-31))||(k>(max_i[i]+31))){
					x[i][k] = 0.;
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
					double efficiency = gr_eff_vs_strength[134][board][band-1]->Eval(vpp);
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
				if(band_count>=5){
					ant_count++;
				}
				//cout<<band_count<<"  ";
			}
			if(ant_count>=3){
				//cout<<"trigger!!"<<endl;
				trigger_count++;
			}
			else{
				//cout<<"no trigger!!"<<endl;
			}
			//cout<<endl<<ant_count<<endl;
			pps_count++;
	} // j loop
	file->Close();
	delete file;
	} // number loop
	
	trig_eff = trigger_count/448.;
	cout<<attenuation<<"   "<<"trigger_eff:  "<<trig_eff<<endl;
	eff_list[m] = trig_eff;
	
	//cout<<pps_count<<"  "<<trigger_count<<"   "<<pulser_count<<endl;
	} // m loop
	
	
	ntuple_file->Close();
	//delete ntuple_file;
	cout<<"fuck"<<endl;
	
	ntuple_file = new TFile("./ntuple.root", "update");
	ntuple = (TNtuple*) ntuple_file->Get("ntuple");
	auto trig_eff_branch = ntuple->Branch("trig_eff",&trig_eff,"trig_eff/F");
	float V_peak_4 = 0.;
	auto V_p_branch = ntuple->Branch("V_p_4",&V_peak_4,"V_p_4/F");
	for(int m=0; m<nt_nentries; m++){
      	ntuple->GetEntry(m);
      	trig_eff = eff_list[m];
      	V_peak_4 = V_p_list[m];
		trig_eff_branch->Fill();
		V_p_branch->Fill();
	}
	
	ntuple->Write("", TObject::kOverwrite); 
	ntuple_file->Close();
	
	//cout<<pps_count<<"  "<<trigger_count<<"   "<<pulser_count<<endl;
      		
	/*
	double x[4096] = {0.};
	double snr_list[160] = {0.};
	double efficiency_list[160] = {0.};
	for(int strength = 100; strength<260; strength++){
		
		for(int i=0; i<1500; i++){
			x[i] = strength*signal_strength[i]/template_v_pp;
		}
		//double sub_band_power[8] = {0.};
		//get_sub_band_strength(x, sub_band_power);
		//cout<<sqrt(sub_band_power[7]/template_sub_band_power[7])<<endl;
		//cout<<sqrt(sub_band_power[5]/template_sub_band_power[5])<<endl<<endl;
		
		int sys_trigger = 0;
		for(int i=0;i<100; i++){
			int ant_count = 0;
			for(int board=1; board<5; board++){
				//double noise = (r3->Rndm()-0.5)*20.;
				
				double antenna_4_factor = 1;
				if(board==4){
					antenna_4_factor = 1./sqrt(2.);
				}
				double board_x[4096] = {0};
				for(int i=512; i<575; i++){
					
					board_x[i] = antenna_4_factor*(x[i] + (r3->Rndm()-0.5)*20.);
				}
				double sub_band_power[8] = {0.};
				get_sub_band_strength(board_x, sub_band_power);
				
				
				int band_count = 0;
				for(int band=1; band<9; band++){
					double vpp = sqrt(sub_band_power[band-1]/template_sub_band_power[band-1]) * template_v_pp;
					//cout<<band<<"  "<<vpp<<endl;
					double efficiency = gr_eff_vs_strength[390][board-1][band-1]->Eval(vpp);
					double random_number = r3->Rndm();
					bool band_trigger = false;
					if(random_number<efficiency){
						band_trigger = true;
					}
					// cyc caution! probability of trigger by CW 
					random_number = r3->Rndm();
					if(random_number<0.3){
						band_trigger = true;
					}
					if(band_trigger == true){
						band_count ++;
					}
				}
				if(band_count>=6){
					ant_count++;
				}
			}
			if(ant_count>=4){
				sys_trigger++;
			}
		}
		double sys_efficiency = sys_trigger/100.;
		cout<<strength/2.<<"   "<<100*sys_efficiency<<endl;
		snr_list[strength-40] = strength/2./10.;
		efficiency_list[strength-40] = 100*sys_efficiency;
	}
	
	TGraph *cyc = new TGraph(160, snr_list, efficiency_list);
	cyc->Draw();
	*/

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
