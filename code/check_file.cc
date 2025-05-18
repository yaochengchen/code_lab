/*
void check_file(){
	for(int run = 10000; run<90000; run++){
		const char* dir_cr = Form("../iron/%d/", run); 
		const string strCoreasDir = Form("SIM%06d_coreas/", 1);
		const string wf_files = dir_cr + strCoreasDir;
		// 0 means exist, -1 means not
		if(access(dir_cr, F_OK)==0){
			if(access(wf_files.c_str(), F_OK)!=0)
			{
				cout<<wf_files<<"  do not exist! skip it"<<endl;
				//return 0;
			}
		}
	}
}
*/


void check_file(){
	for(int run = 1; run<63; run++){
		int true_run = run*100000+1;
		const char* dir_cr = Form("../iron/%d/", true_run); 
		const string strCoreasDir = Form("SIM%06d_coreas/", 1);
		const string wf_files = dir_cr + strCoreasDir;
		// 0 means exist, -1 means not
		//if(access(dir_cr, F_OK)==0){
			if(access(wf_files.c_str(), F_OK)!=0)
			{
				cout<<wf_files<<"  do not exist! skip it"<<endl;
				//return 0;
			}
			else{
				cout<<run<<endl;
			}
		//}
	}
}
