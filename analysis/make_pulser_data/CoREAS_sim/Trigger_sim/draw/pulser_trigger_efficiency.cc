void pulser_trigger_efficiency(){
	double attenuation[9] = {69, 71, 73., 75., 77., 79., 81., 83., 85};
	double pulser_V_p[9] = {0.};
	double pulser_count[9] = {583, 579, 580., 557., 35., 25., 10., 1, 0.}; //ntuple->Draw("(3*V_p+V_p*sqrt(2))/4.","time>198&&time<405&&(time_stamp<0.1||time_stamp>0.9||(time_stamp>0.4&&time_stamp<0.7))&&attenuation>72&&attenuation<74&&trig_eff>0.85")
	double pulser_eff[9] = {0.};
	
	double sim_eff[10] = {1., 1., 0.9688, 0.85, 0.6171, 0.3484, 0.1862, 0.02455, 0., 0.};
	double sim_V_p[10] = {0.};
	
	double err_pulser[9] = {0.};
	double err_sim[7] = {0.};
	double err_pulser_V[7] = {0.};
	double err_sim_V[7] = {0.};
	
	
	for(int i=0; i<9; i++){
		pulser_V_p[i] = sqrt(1.) * 76.21/sqrt(pow(10, 2*(i-2)/10.)); // 2dB gap, 73 dB is 55.61 average V_p; sqrt(2.) is because of splitter in T4 system
		//cout<<(76.21/sqrt(pow(10, 2*i/10.)) + 76.21/sqrt(pow(10, 2*(i+1)/10.)))/2.<<endl;
		err_pulser_V[i] = 9.514/sqrt(pow(10, 2*(i-2)/10.))/1.; // 73 dB rms is 5.171
		
		double pulser_sent = (5*18*(1.-0.022*14));// 5 pulses per second, 18 second, roughly 19Hz trigger rate, dead time for each event is 0.024s
		if(i==2){pulser_sent = (4*19*(1.-0.022*14));}// only 4 pulses for first attenuation, but it has 19 seconds.
		if(i<4) pulser_sent = 584.;
		
		err_pulser[i] = sqrt(pulser_count[i])/pulser_sent;
		
		pulser_eff[i] = pulser_count[i]/pulser_sent;
		
		//cout<<pulser_eff[i]<<"   "<<err_pulser[i]<<"   "<<pulser_sent<<endl;
		
		if(pulser_eff[i]>1){pulser_eff[i] = 1;}
		
		err_sim_V[i] = err_pulser_V[i];
		err_sim[i] = sqrt((pulser_count[i]*448)*sim_eff[i])/(pulser_count[i]*448);

							}
	
	
	  for(int i=0; i<7; i++){
	  	  sim_V_p[i+2] = pulser_V_p[i+2];
	  	  }
	  	  sim_V_p[1] = sqrt(1.) * 76.21*sqrt(pow(10, 2*1/10.));
	  	  sim_V_p[0] = sqrt(1.) * 76.21*sqrt(pow(10, 2*2/10.));
							



	//TGraph * gr1 = new TGraphErrors(7, pulser_V_p, sim_eff, err_sim_V, err_sim);
	TGraph * gr1 = new TGraph(10, sim_V_p, sim_eff);

	gr1->Draw("AL");
	gr1->SetLineColor(3);
	gr1->SetLineWidth(2);
	gr1->SetMarkerStyle(8);
	gr1->SetMarkerColor(1);
	gr1->GetYaxis()->SetRangeUser(0,1.01);
	gr1->GetXaxis()->SetTitle("V_{peak} (mV)");
	gr1->GetYaxis()->SetTitle("Efficiency");
	
	

	TGraph * gr = new TGraphErrors(9, pulser_V_p, pulser_eff, err_pulser_V, err_pulser);

	gr->Draw("P same");
	gr->SetLineColor(2);
	gr->SetLineWidth(2);
	gr->SetMarkerStyle(21);
	gr->SetMarkerColor(2);
	gr->GetYaxis()->SetRangeUser(0,1.01);
	
	TLegend *leg_R = new TLegend(0.1,0.75,0.35,0.9);
	leg_R->AddEntry(gr,"Pulser","p");
	leg_R->AddEntry(gr1,"simulation","l");

	leg_R->Draw();
					
}
