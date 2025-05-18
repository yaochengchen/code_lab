


void draw_tree(){

	Double_t trigger_spectrum[8][20] = {0.};
	double pass_trigger_ratio[2] = {0.};
	int is_pps = 0;
	double SNR[2] = {0.};
	int V_bits = 0;
	int H_bits = 0;
		
	double pps_ratio[100000] = {0.};
	double pulser_ratio[10000] = {0.};
	//double RF_ratio[1000000] = {0.};
	double *RF_ratio = new double[10000000];
	double pps_snr[100000] = {0.};
	double pulser_snr[10000] = {0.};
	//double RF_snr[1000000] = {0.};
	double *RF_snr = new double[10000000];
	int pps_count = 0;
	int RF_count = 0;
	int pulser_count = 0;

TH1F * barson_RF = new TH1F("RF_trigger", "RF_trigger", 77, 0, 76);
TH1F * barson_pps = new TH1F("RF_trigger", "RF_trigger", 77, 0, 76);
TH1F * barson_pulser = new TH1F("RF_trigger", "RF_trigger", 77, 0, 76);
	
	TFile *RF_file = new TFile("../0702t_selection_20210411.root");
    TTree *Tree_RF = (TTree*) RF_file->Get("selection");
    Tree_RF->SetBranchAddress("is_pps", &is_pps);
    Tree_RF->SetBranchAddress("V_bits", &V_bits);
    Tree_RF->SetBranchAddress("H_bits", &H_bits);
    Tree_RF->SetBranchAddress("SNR", SNR);
    Tree_RF->SetBranchAddress("trigger_spectrum", trigger_spectrum);
    Tree_RF->SetBranchAddress("pass_trigger_ratio", pass_trigger_ratio);
    int Entries = Tree_RF->GetEntries();
    
    TH1F * RF_trigger = new TH1F("RF_trigger", "RF_trigger", 100, -50, 50);
    TH1F * pps_trigger = new TH1F("pps_trigger", "pps_trigger", 100, -50, 50);
    for(int entry=0; entry<1*Entries/10.; entry++){
    	Tree_RF->GetEntry(entry);
    	
    	int pol = 1;
    	if(H_bits>V_bits){pol=0;}
    	
    	int cyc = 0;
    	double threshold = 10*log10(pow(SNR[pol]/60., 2)) + 3;
    	if(threshold<3) {threshold=3;}
    	//int total_count = 0;
    	//int pass_count = 0;
    	for(int ant=0; ant<4; ant++){
    		int i = 2*ant + pol;
    		for(int k=1; k<20; k++){
    			if(is_pps==1){
    				pps_trigger->Fill(trigger_spectrum[i][k]);
    				//if(trigger_spectrum[i][k]>threshold){cyc++;}
    			}
    			if(is_pps==0){
    				//if(SNR[1]>100)
    				RF_trigger->Fill(trigger_spectrum[i][k]);
    				//if(trigger_spectrum[i][k]>threshold){cyc++;}
    			}
    			
    			//if(trigger_spectrum[i][k] != 0){total_count++;}
    			//trigger_spectrum[i][k]
    			
    		}
    	}
    	
    	if(is_pps==1){
    		//cout<<H_bits<<"   "<<V_bits<<"   "<<pol<<endl;
    		//V-pol
    		//pol = 0;
    		pps_ratio[pps_count] = pass_trigger_ratio[pol]*80;
    		//cout<<pps_ratio[pps_count]<<endl;
    		barson_pps->Fill(pps_ratio[pps_count]);
    		pps_snr[pps_count] = SNR[pol];
    		//cout<<pps_snr[pps_count]<<endl;
    		pps_count++;
    	}
    	if(is_pps==0){
    		//V-pol
    		RF_ratio[RF_count] = pass_trigger_ratio[pol]*80;
    		barson_RF->Fill(RF_ratio[RF_count]);
    		RF_snr[RF_count] = SNR[pol];
    		RF_count++;
    	}
    }
    Double_t factor = 1.;
	RF_trigger->Scale(factor/RF_trigger->GetEntries());
	RF_trigger->SetLineColor(2);
    RF_trigger->Draw();
	pps_trigger->Scale(factor/pps_trigger->GetEntries());
	pps_trigger->SetLineColor(3);
    pps_trigger->Draw("same");
    
    
    
	TFile *pulser_file = new TFile("../0625t_selection_20240102.root");
    TTree *Tree_pulser = (TTree*) pulser_file->Get("selection");
    Tree_pulser->SetBranchAddress("is_pps", &is_pps);
    Tree_pulser->SetBranchAddress("SNR", SNR);
    Tree_pulser->SetBranchAddress("V_bits", &V_bits);
    Tree_pulser->SetBranchAddress("H_bits", &H_bits);
    Tree_pulser->SetBranchAddress("trigger_spectrum", trigger_spectrum);
    Tree_pulser->SetBranchAddress("pass_trigger_ratio", pass_trigger_ratio);
    Entries = Tree_pulser->GetEntries();
    
    TH1F * pulser_trigger = new TH1F("pulser_trigger", "pulser_trigger", 100, -50, 50);
    TH1F * pulser_pps_trigger = new TH1F("pulser_pps_trigger", "pulser_pps_trigger", 100, -50, 50);
    for(int entry=0; entry<Entries; entry++){
    	Tree_pulser->GetEntry(entry);
    	
    	int pol = 1;
    	if(H_bits>V_bits){pol=0;}
    	
    	for(int ant=0; ant<4; ant++){
    		//this is a H-pol pulser
    		int i = 2*ant + 0;
    		for(int k=1; k<20; k++){
    			if(is_pps==1){
    				pulser_pps_trigger->Fill(trigger_spectrum[i][k]);
    			}
    			if(is_pps==0){
    				pulser_trigger->Fill(trigger_spectrum[i][k]);
    			}
    			
    		}
    	}
    	if(is_pps==0){
    		//H-pol
    		pulser_ratio[pulser_count] = pass_trigger_ratio[pol]*80;
    		barson_pulser->Fill(pulser_ratio[pulser_count]);
    		pulser_snr[pulser_count] = SNR[pol];
    		//cout<<SNR[0]<<endl;
    		pulser_count++;
    	}
    }

    factor = 1.;
	pulser_trigger->Scale(factor/pulser_trigger->GetEntries());
	pulser_trigger->SetLineColor(4);
    pulser_trigger->Draw("same");
	pulser_pps_trigger->Scale(factor/pulser_pps_trigger->GetEntries());
	pulser_pps_trigger->SetLineColor(5);
    //pulser_pps_trigger->Draw("same");
    
    //cout<<pps_count<<"  "<<RF_count<<endl;
    //cout<<pps_snr[10]<<"   "<<pps_ratio[10]<<endl;
    TGraph * pps_pass = new TGraph(pps_count, pps_snr, pps_ratio);
    TGraph * RF_pass = new TGraph(RF_count, RF_snr, RF_ratio);
    TGraph * pulser_pass = new TGraph(pulser_count, pulser_snr, pulser_ratio);

    double draw_point[2][2]={{0, 500}, {0, 76}};
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
    
    pulser_pass->Draw("P same");
    pulser_pass->SetMarkerColor(4);
    pulser_pass->SetMarkerStyle(4);  
    pulser_pass->SetMarkerSize(0.5);

/*
TCanvas *c1 = new TCanvas("c1","Canvas Example",1500,1000);
c1->SetLogy();
c1->cd();
barson_RF->Draw("");
barson_RF->SetLineColor(2);
barson_pulser->Draw("same");
barson_pulser->SetLineColor(3);
barson_pps->Draw("same");
barson_pps->SetLineColor(4);
*/


    
}
