
int pol = 1;

void tree_selection_simulation(){

	Double_t trigger_spectrum[8][20] = {0.};
	double pass_trigger_ratio[2] = {0.};
	double eventTime;
	int is_pps = 0;
	int V_bits = 0;
	int H_bits = 0;
	double SNR[8] = {0.};
	double pulse_power[8] = {0.};
    double reconstruct_cross_correlation[2] = {0.};
    double rough_cross_correlation[2] = {0.};
    double remain_pulse_cross_correlation[2] = {0.};
    double other_pol_cross_correlation[2] = {0.};
    
    double reconstructed_theta[2] = {0.};
    double reconstructed_phi[2] = {0.};
    double TDOA[6] = {0.};
	
	//double pps_ratio[1000000] = {0.};
	//double pulser_ratio[1000000] = {0.};
	//double RF_ratio[1000000] = {0.};
	//double pps_snr[1000000] = {0.};
	//double pulser_snr[1000000] = {0.};
	//double RF_snr[1000000] = {0.};
	int pps_count = 0;
	int RF_count = 0;
	int pulser_count = 0;
	int gRunNum=0;
	int gEventNum = 0;
	
	double start_time = 0;
	
	double ShowerTheta;
	double ShowerPhi;
	
	double theta_list[29], phi_list[19];
	
	// pass_trigger_ratio > 0.7, rough_cross_correlation>0.7; 70 mV signal, 26./551. pass, 94 mV signal, 506./551. pass, 
	// why using this? because the pps data added in simulation waveform is random in time, so pass_trigger_ratio is not very precise since it needs the evaluation of the environment noise level a run.
	//auto fa2 = new TF1("fa2","(20./551.)*x-1374./551.",60,80);
	auto fa2 = new TF1("fa2","0.043557169*x-2.5662432",60,80);
	


	const  char * root_filename =  "/media/cyc/Data_disk/Linux/TAROGE4_PROCESSED_DATA/240822_t_selection_20250230.root";
	TFile *RF_file = new TFile(root_filename);
    TTree *Tree_RF = (TTree*) RF_file->Get("selection");
    cout<<"here2"<<endl;
    Tree_RF->SetBranchAddress("eventTime", &eventTime); //PC timestamp
    Tree_RF->SetBranchAddress("is_pps", &is_pps);
    Tree_RF->SetBranchAddress("V_bits", &V_bits);
    Tree_RF->SetBranchAddress("H_bits", &H_bits);
    Tree_RF->SetBranchAddress("run", &gRunNum);
    Tree_RF->SetBranchAddress("event", &gEventNum);
    
	Tree_RF->SetBranchAddress("theta", reconstructed_theta);
	Tree_RF->SetBranchAddress("phi", reconstructed_phi);
	Tree_RF->SetBranchAddress("TDOA", TDOA);
	
    Tree_RF->SetBranchAddress("pulse_power", pulse_power);
    Tree_RF->SetBranchAddress("SNR", SNR);
    Tree_RF->SetBranchAddress("rough_cross_correlation", rough_cross_correlation);
	Tree_RF->SetBranchAddress("reconstruct_cross_correlation", reconstruct_cross_correlation);
	Tree_RF->SetBranchAddress("remain_pulse_cross_correlation", remain_pulse_cross_correlation);
	Tree_RF->SetBranchAddress("other_pol_cross_correlation", other_pol_cross_correlation);
    Tree_RF->SetBranchAddress("trigger_spectrum", trigger_spectrum);
    Tree_RF->SetBranchAddress("pass_trigger_ratio", pass_trigger_ratio);
    int Entries = Tree_RF->GetEntries();
    
    Tree_RF->SetBranchAddress("ShowerTheta", &ShowerTheta);
    Tree_RF->SetBranchAddress("ShowerPhi", &ShowerPhi);
    
    
    
    double total_count[8][29][19] = {0};
    double pass_count[8][29][19] = {0};
    double analysis_efficiency[8][29][19] = {0};
    double analysis_efficiency_pt[8][19][29] = {0};
    for(int entry=0; entry<Entries; entry++){
    	//cout<<entry<<endl;
    	Tree_RF->GetEntry(entry);
    	
    	if(is_pps!=0){continue;}
    	
    	if(gRunNum>2000){gRunNum /= 10;}
    	int e = (gRunNum-1750)/25;
    	if(e==8){e=7;}//for 1950
    	//cout<<gRunNum<<endl;
    	double st = round(ShowerTheta);
    	//cout<<st<<endl;

    	int i_theta;
    	if(st==50) i_theta = 0;
    	if(st==48) i_theta = 1;
    	if(st==46) i_theta = 2;
    	if(st==44) i_theta = 3;
    	if(st==42) i_theta = 4;
    	if(st==40) i_theta = 5;
    	if(st==38) i_theta = 6;
    	if(st==36) i_theta = 7;
    	if(st==34) i_theta = 8;
    	if(st==32) i_theta = 9;
    	if(st==30) i_theta = 10;
    	if(st==28) i_theta = 11;
    	if(st==26) i_theta = 12;
    	if(st==25) i_theta = 13;
    	if(st==24) i_theta = 14;
    	if(st==22) i_theta = 15;
    	if(st==20) i_theta = 16;
    	if(st==18) i_theta = 17;
    	if(st==16) i_theta = 18;
    	if(st==15) i_theta = 19;
    	if(st==14) i_theta = 20;
    	if(st==12) i_theta = 21;
    	if(st==10) i_theta = 22;
    	if(st==8) i_theta = 23;
    	if(st==6) i_theta = 24;
    	if(st==5) i_theta = 25;
    	if(st==4) i_theta = 26;
    	if(st==2) i_theta = 27;
    	if(st==0) i_theta = 28;
    	int i_phi = round((ShowerPhi + 67.2845)/5.);
    	//if(i_theta>28||i_theta<0)cout<<st<<"   "<<i_theta<<endl;
    	//if(i_phi>18||i_phi<0)cout<<ShowerPhi<<"   "<<i_phi<<endl;
    	//if(e>8||e<0)cout<<e<<endl;
    	theta_list[i_theta] = 90 - ShowerTheta;
    	phi_list[i_phi] = ShowerPhi;
    	
    	
    	total_count[e][i_theta][i_phi] += 1;
    	//cout<<i_theta<<"   "<<i_phi<<"  "<<e<<"   "<<total_count[e][i_theta][i_phi]<<endl;
   	

    	if(H_bits>V_bits){pol = 0;}
    	else{pol = 1;}
    	
		double highest_power = 0;
		double lowest_power = 1.0e12;
		pulse_power[6+pol] = 2*pulse_power[6+pol];
		SNR[6+pol] = sqrt(2)*SNR[6+pol];
		double average_peak = 0.;
		for(int ant=0; ant<4; ant++){
			if(pulse_power[ant+pol*4]>highest_power){highest_power=pulse_power[ant+pol*4];}
			if(pulse_power[ant+pol*4]<lowest_power){lowest_power=pulse_power[ant+pol*4];}
			average_peak += SNR[ant+pol*4]/4.;
			//cout<<pulse_power[2*ant+pol]<<endl;
		}
		for(int ch=0; ch<8; ch++){
			
			//cout<<pulse_power[ch]<<endl;
		}
		//cout<<endl;
		//cout<<lowest_power<<"   "<<highest_power<<endl;
		//cout<<rough_cross_correlation[pol]<<endl;
    	
    	//cout<<reconstruct_cross_correlation[pol]<<endl;
    	//&&lowest_power/highest_power>-10.0
    	//if(pass_trigger_ratio[pol]>0.8&&reconstruct_cross_correlation[pol]>0.8&&fabs(ShowerTheta-reconstructed_theta[pol])<3){
    	//if(reconstruct_cross_correlation[pol]>0.75&&fabs(ShowerTheta-reconstructed_theta[pol])<3){
    	//if(SNR[0]>80){
    	//if(true){
    	//if(reconstruct_cross_correlation[pol]>0.75&&fabs(TDOA[1])<60&&TDOA[5]<-17){

		//	pass_count[e][i_theta][i_phi] ++;
 
    					
		//		}
	
    	//if(e==2&&i_theta==2) cout<<i_theta<<"  "<<i_phi<<"   "<<pass_count[e][i_theta][i_phi]<<"   "<<total_count[e][i_theta][i_phi]<<endl;	
		
		
		if(rough_cross_correlation[pol]>0.7){
			double average_v_peak = (SNR[0] + SNR[1] + SNR[2] + sqrt(2.)*SNR[3])/4.;
			double trigger_probability = fa2->Eval(average_v_peak);
			if(trigger_probability<0){trigger_probability=0;}
			if(trigger_probability>1){trigger_probability=1;}
			pass_count[e][i_theta][i_phi] += trigger_probability;
		}

    }
    
    for(int i_theta = 2; i_theta<3; i_theta++){
		for(int i_phi = 0; i_phi<19; i_phi++){
    		//cout<<i_theta<<"  "<<i_phi<<"   "<<pass_count[2][i_theta][i_phi]<<"   "<<total_count[2][i_theta][i_phi]<<endl;
    	}
    }
    
    
    TFile *f_ntuple = new TFile("ntuple_analysis_efficiency_iron.root","RECREATE");
	TNtuple *ntuple = new TNtuple("ntuple","t","e:i_theta:i_phi:efficiency");
    
    double analysis_efficiency_t[8][29] = {0.};
    double t_pass[8][29] = {0.};
    double t_count[8][29] = {0.};
    double analysis_efficiency_p[8][19] = {0.};
    double p_pass[8][19] = {0.};
    double p_count[8][19] = {0.};
    for(int e=0; e<8; e++){
    	for(int i_theta = 0; i_theta<29; i_theta++){
			for(int i_phi = 0; i_phi<19; i_phi++){
				t_pass[e][i_theta] += pass_count[e][i_theta][i_phi];
				p_pass[e][i_phi] += pass_count[e][i_theta][i_phi];
				t_count[e][i_theta] += total_count[e][i_theta][i_phi];
				p_count[e][i_phi] += total_count[e][i_theta][i_phi];
				
				if(total_count[e][i_theta][i_phi]>0) {analysis_efficiency[e][i_theta][i_phi] = pass_count[e][i_theta][i_phi]/total_count[e][i_theta][i_phi];}
				else{analysis_efficiency[e][i_theta][i_phi] = 0;}
				//cout<<analysis_efficiency[e][i_theta][i_phi]<<endl;
				analysis_efficiency_pt[e][i_phi][i_theta] = analysis_efficiency[e][i_theta][i_phi];
				
				//if(e==1&&i_theta==3) cout<<i_theta<<"  "<<i_phi<<"   "<<pass_count[e][i_theta][i_phi]<<"   "<<total_count[e][i_theta][i_phi]<<"  "<<analysis_efficiency_pt[e][i_phi][i_theta]<<endl;
				
				ntuple->Fill(e, i_theta, i_phi, analysis_efficiency[e][i_theta][i_phi]);
			}
		}
		
	}
	
	
	for(int e=7; e>=0; e--){
		for(int i_theta = 0; i_theta<29; i_theta++){
			if(t_count[e][i_theta]>0){
				analysis_efficiency_t[e][i_theta] = t_pass[e][i_theta]/t_count[e][i_theta];
			}
				
			//if(e==1&&i_theta==3) 
			//cout<<i_theta<<"   "<<t_pass[e][i_theta]<<"   "<<t_count[e][i_theta]<<"  "<<analysis_efficiency_t[e][i_theta]<<endl;
		}
		TGraph *cyc = new TGraph(29, theta_list, analysis_efficiency_t[e]);
		if(e==7) cyc->Draw("APL");
		else cyc->Draw("PL same");
		cyc->GetYaxis()->SetRangeUser(0,1);
		cyc->SetLineColor(e+1);
		if(e==4) cyc->SetLineColor(e+3);
		cyc->SetLineWidth(2);
	}
	//TGraph *cyc = new TGraph(19, phi_list, analysis_efficiency[1][4]);
	//cyc->Draw("APL");
	
	
	f_ntuple->cd();
	f_ntuple->Write(0, TObject::kOverwrite);
	f_ntuple->Close();


    
}
