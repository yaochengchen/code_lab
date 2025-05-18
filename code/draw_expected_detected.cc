

void draw_expected_detected(){



double BinEdge_theta_list[7] = {40, 50.33, 59.29, 67.48, 75.20, 82.66, 90};
double BinCenter_theta_list[6] = {0.};

	for(int i=0; i<6; i++){
		BinCenter_theta_list[i] = (BinEdge_theta_list[i] + BinEdge_theta_list[i+1])/2.;
	}
	
	

double candidates_angle[55][2]={{6.7, 20.8}, {16.6, -52.9}, {12.4, -1.3}, {6.1, -7.6}, {7.4, -31.6}, {13.3, -5.8}, {24.6, -19}, {2.61, -22.56}, {25.9, -31.3}, {14, -11.5}, {29, -13.9}, {8.7, -11.8}, {25.7, -31}, {7.4, -31.6}, {4.6, -16.69}, {25.3, -4.1}, {28.3, -7.6}, {9.7 -30.7}, {28.7, -12.1}, {9.3, -59.5}, {27.7, -44.5}, {34.4, -42.7}, {22, -62.5}, {10, -15.1}, {21, -16.3}, {11.7, -12.4}, {27.7, -43}, {0.5, -56.65}, {31.6, -23.8}, {15.4, 0.2}, {9.6, 1.1}, {2.1, 18.9}, {18.3, -54.7}, {28.3, -50.2}, {27.3, -1.9}, {2.9, -24.86}, {51, 7.4}, {19.3, -45.4}, {30.6, -12.1}, {12, -1.3}, {27, -7.9}, {17.7, -29.2}, {11, -37.6}, {9.3, -10.9}, {9.7, -24.7}, {15.7, -49.9}, {41.6, -32.2}, {30.7, -11.2}, {3.1, -31.1}, {31.7, 13.1}, {41.6, 21.8}, {31.4, 5.6}, {18.3, -43.9}, {21.6, 12.8}, {18.7, 11.9}};


	double candidates_theta[55]={0};//{0,0,4,3,3};
	double candidates_phi[55]={0.};
	//double candidates_theta_error[5]={0,0,2,1.732,1.732};
	double all_one[55];
	for(int i=0; i<55; i++){
		candidates_phi[i] = candidates_angle[i][1];//-26;
		candidates_theta[i] = 90-candidates_angle[i][0];
		all_one[i] = 1;
	}
	
	/*
	double phi_number[19] = {0.};
	double phi_number_error[19] = {0.};
	for(int k=0; k<19; k++){
		for(int i=0; i<55; i++){
			if((candidates_phi[i]>(phi_list[k]-2.5))&&(candidates_phi[i]<=(phi_list[k]+2.5))){
				phi_number[k] ++;
			}
		}
		phi_number_error[k] = sqrt(phi_number[k]);
	}
	*/
	
	
	double theta_number[6] = {0.};
	double theta_number_error[6] = {0.};
	
	for(int i=0; i<55; i++){
		for(int j=0; j<6; j++){
			if((candidates_theta[i]>=BinEdge_theta_list[j])&&(candidates_theta[i]<BinEdge_theta_list[j+1])){
				theta_number[j] ++;
			}
		}
		
	}
		
	for(int k=0; k<6; k++){
		theta_number_error[k] = sqrt(theta_number[k]);
		//cout<<theta_number[k]<<endl;
	}
	
	// value, lower error, upper error
	double _sim_E_number_proton[8][3] = {{0.00357411, 0.0028027, 0.0032987}, {3.22523, 2.32338, 0.967263}, {17.4937, 2.356, 2.13259}, {15.9009, 1.60489, 1.34227}, {9.63429, 0.707777, 0.691813}, {4.85426, 0.322734, 0.320265}, {2.36363, 0.142801, 0.142304}, {1.54225, 0.0878298, 0.0876577}};
	
	double _sim_E_number_iron[8][3] = {{0.00390955, 0.00324253,  0.00343924}, {1.53169, 0.875435,  0.493969}, {12.1149, 1.53669,  1.51047}, {12.0996, 1.35266,  1.10667}, {8.13694, 0.631109,  0.622381}, {4.22384, 0.292098,  0.287454}, {2.22029, 0.137298,  0.136791}, {1.56164, 0.0891545,  0.0889985}};
	
	double _sim_theta_number_proton[6][3] = {{2.6882, 0.526643, 0.452674}, {5.44584, 0.895225, 0.718818}, {9.32936, 1.37159, 1.05793}, {14.1728, 2.18723, 1.41834}, {14.779, 1.56956, 1.32226}, {8.6026, 0.99796, 0.717437}};
	double _sim_theta_number_iron[6][3] = {{2.29286, 0.311305, 0.309072}, {4.66892, 0.6256, 0.544351}, {7.59952, 0.911171, 0.770605}, {11.5478, 1.35376, 1.2049}, {9.92264, 1.10977, 0.932074}, {5.86102, 0.606076, 0.489171}};
	
	double sim_E_number_proton[3][8] = {0.};
	double sim_E_number_iron[3][8] = {0.};
	double sim_theta_number_proton[3][6] = {0.};
	double sim_theta_number_iron[3][6] = {0.};
	
	for(int i=0; i<8; i++){
		for(int j=0; j<3; j++){
			sim_E_number_proton[j][i] = _sim_E_number_proton[i][j];
			sim_E_number_iron[j][i] = _sim_E_number_iron[i][j];
			if(i<6){
				sim_theta_number_proton[j][i] = _sim_theta_number_proton[i][j];
				sim_theta_number_iron[j][i] = _sim_theta_number_iron[i][j];
			}
		}
	}
	
	
	TLegend *leg_R = new TLegend(0.1,0.7,0.25,0.9);
	
	TGraphAsymmErrors *ysm = new TGraphAsymmErrors(6,  BinCenter_theta_list, sim_theta_number_proton[0], 0, 0, sim_theta_number_proton[1], sim_theta_number_proton[2]);
	ysm->SetMarkerStyle(8);
	ysm->SetMarkerColor(2);
	ysm->SetLineColor(2);
	ysm->SetFillColor(2);
	ysm->SetFillStyle(3005);
	//ysm->Draw("AL 3");
	ysm->GetXaxis()->SetTitle("Zenith (#circ)");
	//ysm->GetXaxis()->SetTitle("log_{10}(E) eV");
	//ysm->GetYaxis()->SetTitle("Event Rate [N/year]");
	ysm->GetYaxis()->SetTitle("Number of events");
	//leg_R->AddEntry(ysm, "Primary proton", "l");
	
TGraphAsymmErrors *hxj = new TGraphAsymmErrors(6,  BinCenter_theta_list, sim_theta_number_iron[0], 0, 0, sim_theta_number_iron[1], sim_theta_number_iron[2]);
	hxj->SetMarkerStyle(8);
	hxj->SetMarkerColor(3);
	hxj->SetLineColor(3);
	hxj->SetFillColor(3);
	hxj->SetFillStyle(3005);
	//hxj->Draw("L 3");
	//leg_R->AddEntry(hxj, "Primary iron", "l");
	
	
	TGraphErrors *lmx = new TGraphErrors(6, BinCenter_theta_list, theta_number, 0, theta_number_error);
	//lmx->Draw("P same");
	lmx->SetMarkerStyle(8);
	lmx->SetMarkerSize(1);
	//leg_R->AddEntry(lmx, "Observed", "l");
	
	
	
	
	
	double data_Bin_E_list[8] = {17.375, 17.625, 17.875, 18.125, 18.375, 18.625, 18.875, 19.25};
	double Error_data_Bin_E_list[8] = {0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.25};
	
	TGraphAsymmErrors *syt = new TGraphAsymmErrors(8,  data_Bin_E_list, sim_E_number_proton[0], Error_data_Bin_E_list, Error_data_Bin_E_list, sim_E_number_proton[1], sim_E_number_proton[2]);
	syt->SetMarkerStyle(8);
	syt->SetMarkerColor(2);
	syt->SetLineColor(2);
	syt->SetFillColor(2);
	syt->SetFillStyle(3005);
	//syt->Draw("AL 3");
	//syt->GetXaxis()->SetTitle("Zenith (#circ)");
	syt->GetXaxis()->SetTitle("log_{10}(E) eV");
	//syt->GetYaxis()->SetTitle("Event Rate [N/year]");
	syt->GetYaxis()->SetTitle("Number of events");
	//leg_R->AddEntry(syt, "Simulation", "l");	
	

	TGraphAsymmErrors *lyf = new TGraphAsymmErrors(8,  data_Bin_E_list, sim_E_number_iron[0], Error_data_Bin_E_list, Error_data_Bin_E_list, sim_E_number_iron[1], sim_E_number_iron[2]);
	lyf->SetMarkerStyle(8);
	lyf->SetMarkerColor(2);
	lyf->SetLineColor(2);
	lyf->SetFillColor(2);
	lyf->SetFillStyle(3005);
	lyf->Draw("AL 3");
	//lyf->GetXaxis()->SetTitle("Zenith (#circ)");
	lyf->GetXaxis()->SetTitle("log_{10}(E) eV");
	//lyf->GetYaxis()->SetTitle("Event Rate [N/year]");
	lyf->GetYaxis()->SetTitle("Number of events");
	leg_R->AddEntry(lyf, "Simulation", "l");
		
	
	double iron_E_number[8] = {1, 10, 16, 11, 6, 6, 3, 2};
	double iron_E_number_error[8] = {0.};
	for(int k=0; k<8; k++){
		iron_E_number_error[k] = sqrt(iron_E_number[k]);
		//cout<<E_number_error[k]<<endl;
	}

	
	TGraphErrors *iron_cyq = new TGraphErrors(8, data_Bin_E_list, iron_E_number, Error_data_Bin_E_list, iron_E_number_error);
	iron_cyq->Draw("P same");
	iron_cyq->SetMarkerStyle(8);
	iron_cyq->SetMarkerSize(1);
	iron_cyq->SetLineColor(1);
	leg_R->AddEntry(iron_cyq, "Observed", "l");



	double proton_E_number[8] = {3, 8, 14, 16, 5, 3, 2, 4};
	double proton_E_number_error[8] = {0.};
	for(int k=0; k<8; k++){
		proton_E_number_error[k] = sqrt(proton_E_number[k]);
		//cout<<E_number_error[k]<<endl;
	}
	
	TGraphErrors *proton_cyq = new TGraphErrors(8, data_Bin_E_list, proton_E_number, Error_data_Bin_E_list, proton_E_number_error);
	//proton_cyq->Draw("P same");
	proton_cyq->SetMarkerStyle(8);
	proton_cyq->SetMarkerSize(1); 
	proton_cyq->SetLineColor(1);
	//leg_R->AddEntry(proton_cyq, "Observed", "l");
	
	
	
	
	leg_R->Draw();


}
