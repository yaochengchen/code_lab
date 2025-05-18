#include "./read_pulse.cc"

void filter_fre(double x[1500], bool save_100MHz);

TGraph * gr_eff_vs_strength[1000][8][8] = {NULL};
TRandom *r3;

void trigger_sim_initial(){
	read_pulse();
	r3 = new TRandom3();

	TFile* f0 = new TFile("./Trigger_sim/efficiency_vs_strength.root");

	for(int step=300; step<1300; step++){
		for(int board=1; board<9; board++){
			for(int band=1; band<9; band++){
				char *txtname = Form("gr_Step%d_Board%d_Band%d", step, board, band); 
				gr_eff_vs_strength[step-300][board-1][band-1] = (TGraph*)f0->Get(txtname); 
			}
		}
	}
}

bool IsTrigger(double x[8][4096], int threshold_step, TBits& gTrigBits){

		bool is_trigger = false;
		for(int pol=0; pol<2; pol++){
			int ant_count = 0;
			for(int board=0;board<4;board++){
				double sub_band_power[8] = {0.};
				get_sub_band_strength(x[board*2+pol], sub_band_power);
				//cout<<endl;
				int band_count = 0;
				for(int band=1; band<9; band++){
					//cout<<sub_band_power[band-1]<<"  ";
					double vpp = sqrt(sub_band_power[band-1]/template_sub_band_power[band-1]) * template_v_pp;
					//cout<<board<<"  "<<band<<"  "<<vpp<<endl;
					double efficiency = gr_eff_vs_strength[threshold_step-300][board+pol*4][band-1]->Eval(vpp);
					double random_number = r3->Rndm();
					bool band_trigger = false;
					if(random_number<efficiency){
						band_trigger = true;
					}
					//H-pol 0; V-pol 32
					if(band_trigger == true){
						band_count ++;
						gTrigBits.SetBitNumber(pol*32+board*8+band-1, kTRUE);
					}
				}
				if(band_count>=6){
					ant_count++;
				}
				//cout<<band_count<<"  ";
			}
			if(ant_count>=4){
				//cout<<"trigger!!"<<endl;
				is_trigger = true;
			}
		}
		
		return is_trigger;
		
}







