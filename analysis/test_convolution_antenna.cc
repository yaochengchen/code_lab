#include "./make_pulser_data/CoREAS_Sim/utils.cc"

void convolute(double t_E_theta[4096*8], double x[1500], int pos_theta, int pos_phi, int ch){
		
  	FFT_forward_coreas_theta->SetPoints(t_E_theta);
	FFT_forward_coreas_theta->Transform();
	FFT_forward_coreas_theta->GetPointsComplex(theta_re_coreas, theta_im_coreas);
	
	for(int k=0; k<2048*8+1; k++){
		if((k<500)||(k>2000)) {
			theta_re_coreas[k] = 0;
			theta_im_coreas[k] = 0;
		}
		else{
			TComplex f_E_theta(theta_re_coreas[k], theta_im_coreas[k]);
			TComplex A_theta = TComplex(ANT_H_theta_re[pos_theta][pos_phi][k], ANT_H_theta_im[pos_theta][pos_phi][k]);
			//TComplex A_theta = TComplex(0.5, 0);
		
			TComplex A_waveform = f_E_theta*A_theta;
			
			TComplex FEE(FEE_re[ch][k], FEE_im[ch][k]);
			TComplex result = A_waveform*FEE;
			FFT_re_coreas[k] = result.Re();
			FFT_im_coreas[k] = result.Im();
			
			//FFT_re_coreas[k] = A_waveform.Re();
			//FFT_im_coreas[k] = A_waveform.Im();
		}
		
	}

	FFT_backward_coreas->SetPointsComplex(FFT_re_coreas, FFT_im_coreas);
	FFT_backward_coreas->Transform();
	FFT_time_coreas = FFT_backward_coreas->GetPointsReal();
	for(int i=0; i<1500; i++){
		x[i] = FFT_time_coreas[i*8]/(4096.*8.);// *1000: V to mV
		
	}

}


void test_convolution_antenna(){

	get_response();

	double fine_t[4096*8];
	double t[1500];
	double x[3][4096*8] = {0};
	for(int i=0; i<4096*8; i++){
		fine_t[i] = 0.1*i;
	}
	
	for(int i=0; i<1500; i++){
		t[i] = 0.8*i;
	}
	
	for(int i=5000; i<=5030; i++){
		x[0][i] = (i-5000)*4.;
	}
	
	for(int i=5060; i>5030; i--){
		x[0][i] = (5060-i)*4.;
	}
	
	double a[2] = {0, 1500};
	double b[2] = {-200, 200};
	TGraph *abc = new TGraph(2, a, b);
	abc->Draw("AP");
	
	TGraph *cyc = new TGraph(4096*8, fine_t, x[0]);
	cyc->Draw("PL same");
	
	
	convolute(x[0], x[1], 18, 0, 3);
	
	TGraph *wyn = new TGraph(1500, t, x[1]);
	wyn->Draw("PL same");
	
	convolute(x[0], x[2], 18, 0, 5);
	
	TGraph *hxj = new TGraph(1500, t, x[2]);
	hxj->Draw("PL same");
	hxj->SetLineColor(2);


}






