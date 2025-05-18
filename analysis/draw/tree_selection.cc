
int pol = 1;

void tree_selection(){

	Double_t trigger_spectrum[8][20] = {0.};
	double pass_trigger_ratio[2] = {0.};
	double eventTime;
	int is_pps = 0;
	int V_bits = 0;
	int H_bits = 0;
	double SNR[8] = {0.};
	double pulse_power[8] = {0.};
    double reconstruct_cross_correlation[2] = {0.};
    double remain_pulse_cross_correlation[2] = {0.};
    double other_pol_cross_correlation[2] = {0.};
    
    double reconstructed_theta[2] = {0.};
    double reconstructed_phi[2] = {0.};
    double TDOA[6] = {0.};
	
	double pps_ratio[100000] = {0.};
	double pulser_ratio[100000] = {0.};
	double RF_ratio[100000] = {0.};
	double pps_snr[100000] = {0.};
	double pulser_snr[100000] = {0.};
	double RF_snr[100000] = {0.};
	int pps_count = 0;
	int RF_count = 0;
	int pulser_count = 0;
	int gRunNum=0;
	int gEventNum = 0;
	
	double start_time = 0;
	
	cout<<"here1"<<endl;

    const  char * out_filename;
    int which_date = 20211012;
    out_filename =  Form("./candidate_events/%08d.txt", which_date);
    fstream out_put;
    out_put.open(out_filename, ios::out | ios::trunc);// without ios::out   clean the file before write
  
	const  char * root_filename =  Form("../0709t_selection_%08d.root", which_date);
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
	Tree_RF->SetBranchAddress("reconstruct_cross_correlation", reconstruct_cross_correlation);
	Tree_RF->SetBranchAddress("remain_pulse_cross_correlation", remain_pulse_cross_correlation);
	Tree_RF->SetBranchAddress("other_pol_cross_correlation", other_pol_cross_correlation);
    Tree_RF->SetBranchAddress("trigger_spectrum", trigger_spectrum);
    Tree_RF->SetBranchAddress("pass_trigger_ratio", pass_trigger_ratio);
    int Entries = Tree_RF->GetEntries();
    
    TH1F * RF_trigger = new TH1F("RF_trigger", "RF_trigger", 100, -50, 50);
    TH1F * pps_trigger = new TH1F("pps_trigger", "pps_trigger", 100, -50, 50);
    for(int entry=0; entry<Entries; entry++){
    	//cout<<entry<<endl;
    	Tree_RF->GetEntry(entry);
    	
    	if(entry==0) {start_time = eventTime; cout<<start_time<<endl;}

       const time_t tm_start = (const time_t) eventTime;
       tm *tm_gmt = gmtime(&tm_start);
       int hour = tm_gmt->tm_hour+8;
       if(hour>=24){hour-=24;}
           	
    	if(H_bits>V_bits){pol = 0;}
    	else{pol = 1;}
    	
    	// boresight -26 degree
    	reconstructed_phi[pol] = reconstructed_phi[pol] + 26;
    	if(reconstructed_phi[pol]>180) {reconstructed_phi[pol] -= 360;}
    	    	
    	//selection->Draw("run:SNR[1]", "is_pps==0&&reconstruct_cross_correlation[1]>0.7&&remain_pulse_cross_correlation[1]<0.5&&other_pol_cross_correlation[1]>0.6&&pass_trigger_ratio[1]>0.8")

		double time = eventTime - start_time;
    	
		double highest_power = 0;
		double lowest_power = 1.0e12;
		pulse_power[6+pol] = 2*pulse_power[6+pol];
		SNR[6+pol] = sqrt(2)*SNR[6+pol];
		double average_peak = 0.;
		for(int ant=0; ant<4; ant++){
			if(pulse_power[2*ant+pol]>highest_power){highest_power=pulse_power[2*ant+pol];}
			if(pulse_power[2*ant+pol]<lowest_power){lowest_power=pulse_power[2*ant+pol];}
			average_peak += SNR[2*ant+pol]/4.;
		}
		
		//cout<<lowest_power<<"   "<<highest_power<<endl;
    	
    	//&&hour>=0&&hour<6
    	//&&remain_pulse_cross_correlation[pol]<0.5&&TDOA[1]>-50&&pulse_power[6+pol]/pulse_power[2+pol]<1
    	if(pass_trigger_ratio[pol]>0.6&&reconstruct_cross_correlation[pol]>0.7&&remain_pulse_cross_correlation[pol]<0.5&&fabs(TDOA[1])<50&&TDOA[5]<-17&&(hour>=-22||hour<=60)&&time>18000&&time<21600){
    				//cout<<pol<<endl;&&reconstructed_phi[pol]>-6.5&&reconstructed_phi[pol]<-6
    				//&&(hour>=22||hour<=6)
    			
    				out_put<<which_date<<"    "<<gRunNum<<"    "<<gEventNum<<"    "<<reconstruct_cross_correlation[0]<<"    	"<<reconstructed_theta[0]<<"    "<<reconstructed_phi[0]<<"    "<<reconstruct_cross_correlation[1]<<"    "<<reconstructed_theta[1]<<"    	"<<reconstructed_phi[1]<<endl;


    				// 1% of data
    				if(true){
    					if(is_pps==1){
    						//V-pol
    						pps_ratio[pps_count] = pass_trigger_ratio[pol];
    						pps_snr[pps_count] = SNR[4+pol];
    						//cout<<pps_snr[pps_count]<<endl;
    						pps_count++;
    					}
    					if(is_pps==0){
    						//V-pol
    						RF_ratio[RF_count] = pass_trigger_ratio[pol];
    						RF_snr[RF_count] = SNR[4+pol];
    						RF_count++;
    					}
    				}
    					
				}
			//}
		//}
    				
    	
    	for(int ant=0; ant<4; ant++){
    		int i = 2*ant + pol;
    		for(int k=1; k<20; k++){
    			if(is_pps==1){
    				pps_trigger->Fill(trigger_spectrum[i][k]);
    			}
    			if(is_pps==0){
    				//if(SNR[1]>100)
    				RF_trigger->Fill(trigger_spectrum[i][k]);
    			}
    			
    		}
    	}

    }
    Double_t factor = 1.;
	RF_trigger->Scale(factor/RF_trigger->GetEntries());
	RF_trigger->SetLineColor(2);
    RF_trigger->Draw();
	pps_trigger->Scale(factor/pps_trigger->GetEntries());
	pps_trigger->SetLineColor(3);
    pps_trigger->Draw("same");
    

    
    //cout<<pps_count<<"  "<<RF_count<<endl;
    //cout<<pps_snr[10]<<"   "<<pps_ratio[10]<<endl;
    TGraph * pps_pass = new TGraph(pps_count, pps_snr, pps_ratio);
    TGraph * RF_pass = new TGraph(RF_count, RF_snr, RF_ratio);

    double draw_point[2][2]={{0, 500}, {0, 1}};
    TGraph * draw_graph = new TGraph(2, draw_point[0], draw_point[1]);
    draw_graph->Draw("AP");
    draw_graph->SetMarkerColor(0);

    RF_pass->SetMarkerColor(2);
    RF_pass->SetMarkerStyle(4);  
    RF_pass->SetMarkerSize(0.5);    
    RF_pass->Draw("P same");
  
    pps_pass->Draw("P same");
    pps_pass->SetMarkerColor(3);
    pps_pass->SetMarkerStyle(4);  
    pps_pass->SetMarkerSize(0.5);

    
}
