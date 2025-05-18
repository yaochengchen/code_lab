/*
17.625   2.74048   8
17.875   18.7915   16
18.125   17.3889   14
18.375   10.8685   6
18.625   6.24959   4
*/

void flux_estimation(){
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
	
	double flux[2][6] = {0.};
	double flux_error[2][6] = {0.};
	flux[0][0] = 17.625;
	flux[0][1] = 17.875;
	flux[0][2] = 18.125;
	flux[0][3] = 18.375;
	flux[0][4] = 18.625;
	flux[0][5] = 19.125;
	
	flux[1][0] = log10(pow(10, gr_logJ->Eval(17.625))*8/2.74048);
	flux[1][1] = log10(pow(10, gr_logJ->Eval(17.875))*16/18.7915);
	flux[1][2] = log10(pow(10, gr_logJ->Eval(18.125))*14/17.3889);
	flux[1][3] = log10(pow(10, gr_logJ->Eval(18.375))*6/10.8685);
	flux[1][4] = log10(pow(10, gr_logJ->Eval(18.625))*4/6.24959);
	
	flux_error[0][0] = 0.125;
	flux_error[0][1] = 0.125;
	flux_error[0][2] = 0.125;
	flux_error[0][3] = 0.125;
	flux_error[0][4] = 0.125;
	flux_error[0][5] = 0.125;
	
	flux_error[1][0] = log10(pow(10, gr_logJ->Eval(17.625))*(8+sqrt(8))/2.74048) - flux[1][0];
	flux_error[1][1] = log10(pow(10, gr_logJ->Eval(17.875))*(16+sqrt(16))/18.7915) - flux[1][1];
	flux_error[1][2] = log10(pow(10, gr_logJ->Eval(18.125))*(14+sqrt(14))/17.3889) - flux[1][2];
	flux_error[1][3] = log10(pow(10, gr_logJ->Eval(18.375))*(6+sqrt(6))/10.8685) - flux[1][3];
	flux_error[1][4] = log10(pow(10, gr_logJ->Eval(18.625))*(4+sqrt(4))/6.24959) - flux[1][4];
	
	for(int i=0; i<5; i++){
		cout<<flux[1][i]<<"   "<<flux_error[1][i]<<endl;
	}
	gStyle->SetEndErrorSize(10);
	TGraphErrors *cyc = new TGraphErrors(5, flux[0], flux[1], flux_error[0], flux_error[1]);
	cyc->Draw("same P");
	//cyc->SetMarkerSize(0.5);
	cyc->SetMarkerColor(2);
	cyc->SetMarkerStyle(21);
	cyc->SetLineWidth(5);
	
}
