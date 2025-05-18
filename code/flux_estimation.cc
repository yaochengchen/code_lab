


	// value, lower error, upper error
	double _sim_E_number_proton[8][3] = {{0.00357411, 0.0028027, 0.0032987}, {3.22523, 2.32338, 0.967263}, {17.4937, 2.356, 2.13259}, {15.9009, 1.60489, 1.34227}, {9.63429, 0.707777, 0.691813}, {4.85426, 0.322734, 0.320265}, {2.36363, 0.142801, 0.142304}, {1.54225, 0.0878298, 0.0876577}};
	
	double _sim_E_number_iron[8][3] = {{0.00390955, 0.00324253,  0.00343924}, {1.53169, 0.875435,  0.493969}, {12.1149, 1.53669,  1.51047}, {12.0996, 1.35266,  1.10667}, {8.13694, 0.631109,  0.622381}, {4.22384, 0.292098,  0.287454}, {2.22029, 0.137298,  0.136791}, {1.56164, 0.0891545,  0.0889985}};
	
double proton_E_number[7] = {11, 14, 16, 5, 3, 2, 4};

double iron_E_number[7] = {11, 16, 11, 6, 6, 3, 2};
	
	

void flux_estimation(){




	double sim_E_number_proton[3][7] = {0.};
	double sim_E_number_iron[3][7] = {0.};
	
	for(int i=1; i<7; i++){
		for(int j=0; j<3; j++){
			sim_E_number_proton[j][i] = _sim_E_number_proton[i+1][j];
			sim_E_number_iron[j][i] = _sim_E_number_iron[i+1][j];
		}
	}
	
	for(int j=0; j<3; j++){
			sim_E_number_proton[j][0] = _sim_E_number_proton[0][j] + _sim_E_number_proton[1][j];
			sim_E_number_iron[j][0] = _sim_E_number_iron[0][j] + _sim_E_number_iron[1][j];
	}

double data_Bin_E_list[7] = {17.5, 17.875, 18.125, 18.375, 18.625, 18.875, 19.25};

double Error_data_Bin_E_list[7] = {0.25, 0.125, 0.125, 0.125, 0.125, 0.125, 0.25};





	TFile *f = new TFile("/media/cyc/For_Linux/CoREAS_Sim/proton/Auger5045-combined-spectrum-data-2019.root", "read"); //Dir_Sim

	TGraphAsymmErrors *gr_EJ = (TGraphAsymmErrors*) f->Get("gr_EJ");
	TGraphAsymmErrors *gr_logJ = (TGraphAsymmErrors*) f->Get("gr_logJ");
	
	
	TCanvas *c = new TCanvas("c_cr", "cr spec", 0, 0, 1500, 2*800);
	//c->Divide(1,2);
	//c->cd(1);	gPad->SetGrid(1,1);
	//gPad->SetLogy(1);	
	//gr_EJ->Draw("ap");

	c->cd();	gPad->SetGrid(1,1);
	gr_logJ->Draw("apl");
	
	int n = gr_logJ->GetN();
	
	double *x = gr_logJ->GetX();
	double *y = gr_logJ->GetY();
	double new_x[1000] = {0.};
	double new_y[1000] = {0.};
				
	for(int i=0; i<n; i++){
		new_x[i] = x[i];
		new_y[i] = log10(pow(10, y[i]) *365.*24.*60.*60.*1000.*1000.);
	}
	
	
	
	TLegend *leg_R = new TLegend(0.7,0.7,0.9,0.9);
	
	
	TGraph * new_auger = new TGraph(n, new_x, new_y);
	new_auger->Draw("AP");
	new_auger->SetMarkerStyle(33);
	new_auger->SetMarkerSize(1.5);
	new_auger->SetMarkerColor(4);
	//new_auger->SetLineWidth(2);
	
	new_auger->GetXaxis()->SetTitle("log_{10}(E) eV");
	new_auger->GetYaxis()->SetTitle("log_{10}(flux[eV^{-1}km^{-2}yr^{-1}sr^{-1}])");
	leg_R->AddEntry(new_auger, "Auger(2019)", "p");
	
	
	
	
	
	double flux_proton[7] = {0.};
	double flux_error_proton[2][7] = {0.};

	
	for(int i=0; i<7; i++){
		flux_proton[i] = log10(pow(10, gr_logJ->Eval(data_Bin_E_list[i]))*proton_E_number[i]/sim_E_number_proton[0][i]*365.*24.*60.*60.*1000.*1000.);
	}
	
	for(int i=0; i<7; i++){
		flux_error_proton[0][i] = flux_proton[i] - log10(pow(10, gr_logJ->Eval(data_Bin_E_list[i]))*(proton_E_number[i]-sqrt(proton_E_number[i]))/sim_E_number_proton[0][i]*365.*24.*60.*60.*1000.*1000.);
		
		flux_error_proton[1][i] = log10(pow(10, gr_logJ->Eval(data_Bin_E_list[i]))*(proton_E_number[i]+sqrt(proton_E_number[i]))/sim_E_number_proton[0][i]*365.*24.*60.*60.*1000.*1000.) - flux_proton[i];
	}
	

	gStyle->SetEndErrorSize(10);
	TGraphAsymmErrors *cyc_proton = new TGraphAsymmErrors(7, data_Bin_E_list, flux_proton, Error_data_Bin_E_list, Error_data_Bin_E_list, flux_error_proton[0], flux_error_proton[1]);
	//cyc_proton->Draw("P same");
	//cyc_proton->SetMarkerSize(0.5);
	cyc_proton->SetMarkerColor(2);
	cyc_proton->SetMarkerStyle(21);
	cyc_proton->SetLineWidth(3);
	//leg_R->AddEntry(cyc_proton, "TAROGE-4 (this work)", "p");
	
	
	
	double flux_iron[7] = {0.};
	double flux_error_iron[2][7] = {0.};

	
	for(int i=0; i<7; i++){
		flux_iron[i] = log10(pow(10, gr_logJ->Eval(data_Bin_E_list[i]))*iron_E_number[i]/sim_E_number_iron[0][i]*365.*24.*60.*60.*1000.*1000.);
	}
	
	for(int i=0; i<7; i++){
		flux_error_iron[0][i] = flux_iron[i] - log10(pow(10, gr_logJ->Eval(data_Bin_E_list[i]))*(iron_E_number[i]-sqrt(iron_E_number[i]))/sim_E_number_iron[0][i]*365.*24.*60.*60.*1000.*1000.);
		
		flux_error_iron[1][i] = log10(pow(10, gr_logJ->Eval(data_Bin_E_list[i]))*(iron_E_number[i]+sqrt(iron_E_number[i]))/sim_E_number_iron[0][i]*365.*24.*60.*60.*1000.*1000.) - flux_iron[i];
	}
	

	gStyle->SetEndErrorSize(10);
	TGraphAsymmErrors *cyc_iron = new TGraphAsymmErrors(7, data_Bin_E_list, flux_iron, Error_data_Bin_E_list, Error_data_Bin_E_list, flux_error_iron[0], flux_error_iron[1]);
	cyc_iron->Draw("P same");
	//cyc_iron->SetMarkerSize(0.5);
	cyc_iron->SetMarkerColor(2);
	cyc_iron->SetMarkerStyle(21);
	cyc_iron->SetLineWidth(3);	
	leg_R->AddEntry(cyc_iron, "TAROGE-4 (this work)", "p");
	
	
	
	
	
	
	
	

	
double arianna[5][2][1] = 
{{17.81313868613139, -15.957605985037407},
{17.88686131386861, -15.957605985037407},
{17.81386861313868, -15.670822942643390},
{17.74087591240876, -15.953449709060681},
{17.81313868613139, -16.398171238570242}};

double tm[5][2][1] = 
{{17.976642335766424, -15.924355777223608},
{18.148905109489053, -15.924355777223607},
{17.975912408759125, -15.733167082294264},
{17.803649635036496, -15.924355777223608},
{17.975912408759125, -16.539484621778886}};

double anita[5][2][1] = 
{{18.46204379562044, -17.852867830423940},
{18.51897810218978, -17.852867830423940},
{18.46204379562044, -17.615960099750623},
{18.39781021897810, -17.852867830423940},
{18.46204379562044, -18.044056525353284}};


double error_arianna[2][2][1] = {0.};// x->y, lower->high
error_arianna[0][0][0] = arianna[2][0][0] - arianna[3][0][0];
error_arianna[0][1][0] = arianna[1][0][0] - arianna[2][0][0];
error_arianna[1][0][0] = arianna[0][1][0] - arianna[4][1][0];
error_arianna[1][1][0] = arianna[2][1][0] - arianna[0][1][0];
cout<<error_arianna[1][1][0]<<endl;


TGraphAsymmErrors *f_arianna = new TGraphAsymmErrors(1, arianna[0][0], arianna[0][1], error_arianna[0][0], error_arianna[0][1], error_arianna[1][0], error_arianna[1][1]);
f_arianna->Draw("PL same");
f_arianna->SetLineColor(3);
f_arianna->SetLineWidth(3);
f_arianna->SetMarkerColor(3);
f_arianna->SetMarkerStyle(8);
leg_R->AddEntry(f_arianna, "ARIANNA (2017)", "p");

double error_tm[2][2][1] = {0.};// x->y, lower->high
error_tm[0][0][0] = tm[2][0][0] - tm[3][0][0];
error_tm[0][1][0] = tm[1][0][0] - tm[2][0][0];
error_tm[1][0][0] = tm[0][1][0] - tm[4][1][0];
error_tm[1][1][0] = tm[2][1][0] - tm[0][1][0];
cout<<error_tm[1][1][0]<<endl;


TGraphAsymmErrors *f_tm = new TGraphAsymmErrors(1, tm[0][0], tm[0][1], error_tm[0][0], error_tm[0][1], error_tm[1][0], error_tm[1][1]);
f_tm->Draw("PL same");
f_tm->SetLineColor(41);
f_tm->SetLineWidth(3);
f_tm->SetMarkerColor(41);
f_tm->SetMarkerStyle(21);
leg_R->AddEntry(f_tm, "TAROGE-M (2023)", "p");

double error_anita[2][2][1] = {0.};// x->y, lower->high
error_anita[0][0][0] = anita[2][0][0] - anita[3][0][0];
error_anita[0][1][0] = anita[1][0][0] - anita[2][0][0];
error_anita[1][0][0] = anita[0][1][0] - anita[4][1][0];
error_anita[1][1][0] = anita[2][1][0] - anita[0][1][0];
cout<<error_anita[1][1][0]<<endl;


TGraphAsymmErrors *f_anita = new TGraphAsymmErrors(1, anita[0][0], anita[0][1], error_anita[0][0], error_anita[0][1], error_anita[1][0], error_anita[1][1]);
f_anita->Draw("PL same");
f_anita->SetLineColor(6);
f_anita->SetLineWidth(3);
f_anita->SetMarkerColor(6);
f_anita->SetMarkerStyle(23);
leg_R->AddEntry(f_anita, "ANITA (2015)", "p");
	
leg_R->Draw();

}
