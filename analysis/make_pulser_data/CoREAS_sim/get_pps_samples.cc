#include <stdio.h>
#include <dirent.h>
#include <sys/types.h>
#include "utils.cc"
#include "folder_helper.cc"

using namespace std;

void GetTChain(string directory, TChain *chaintree); 

void get_pps_samples(){


  Short_t T1H[1500];
  Short_t T2H[1500];
  Short_t T3H[1500];
  Short_t T4H[1500];
  Short_t T1V[1500];
  Short_t T2V[1500];
  Short_t T3V[1500];
  Short_t T4V[1500];
  TBits* triggerBits = 0;
  double eventTime = 0;
  Int_t gDate = 0;
  Int_t gEventNum = 0;
  Double_t gEventRate = 0.;
  Int_t gThresholdStep = 0;


	TFile *pps_samples = new TFile("pps_samples_cyc.root","RECREATE");
	
 	TTree *pps_samples_Tree = new TTree("t","pps samples");

  pps_samples_Tree->Branch("date", &gDate, "date/I");
  pps_samples_Tree->Branch("EventNum", &gEventNum, "EventNum/I");
  
  pps_samples_Tree->Branch("eventTime", &eventTime, "eventTime/D"); //PC timestamp
  
  pps_samples_Tree->Branch("EventRate", &gEventRate, "EventRate/D");
  pps_samples_Tree->Branch("ThresholdStep", &gThresholdStep, "ThresholdStep/I");

  pps_samples_Tree->Branch("T1H", &T1H, "T1H[1500]/S");
  pps_samples_Tree->Branch("T2H", &T2H, "T2H[1500]/S");
  pps_samples_Tree->Branch("T3H", &T3H, "T3H[1500]/S");
  pps_samples_Tree->Branch("T4H", &T4H, "T4H[1500]/S");
  pps_samples_Tree->Branch("T1V", &T1V, "T1V[1500]/S");
  pps_samples_Tree->Branch("T2V", &T2V, "T2V[1500]/S");
  pps_samples_Tree->Branch("T3V", &T3V, "T3V[1500]/S");
  pps_samples_Tree->Branch("T4V", &T4V, "T4V[1500]/S");
	


	TRandom *r3 = new TRandom3();
	
	vector<int> RootFileDates;
	GetRootFileDates("/media/cyc/1p9TB/TAROGE4_DATA/data/", RootFileDates);
	GetRootFileDates("/media/cyc/For_Linux/TAROGE4_DATA/data/", RootFileDates);
  
    
	for (int i=0; i<RootFileDates.size(); ++i){
		
		int date = RootFileDates[i];
		if(date<20210900){continue;}
		if(date>20220600){continue;}
		
		if((date>=20220528)){continue;}
		
		if((date==20211110)||(date==20211111)||(date==20211222)||(date==20211223)||(date==20220211)||(date==20220212)){continue;}
		
		string root_path = "/media/cyc/1p9TB/TAROGE4_DATA/data/";
		if(date>20220631){root_path = "/media/cyc/For_Linux/TAROGE4_DATA/data/";}
		
		string trigger_path = "/media/cyc/For_Linux/TAROGE4_DATA/trigger/";
		
		
		gDate = date;
		cout<<gDate<<endl;

		TChain *pps_chain = new TChain("t");
		string pps_directory = root_path + to_string(gDate) + "/";
		GetTChain(pps_directory, pps_chain);
		pps_chain->SetBranchAddress("eventTime", &eventTime);
		pps_chain->SetBranchAddress("triggerBits",&triggerBits);
    	const int pps_entries = pps_chain->GetEntries();
    	//cout<<pps_entries<<endl;
    	
    	if(pps_entries<2000){continue;}

		int count_pps=0;
		double *pps_eventTime_list = new double[pps_entries];
		int *entry_list = new int[pps_entries];
		
		for(int i=0; i<pps_entries; i++){
			pps_chain->GetEntry(i);
			if(triggerBits->CountBits()==0){
				entry_list[count_pps] = i;
				pps_eventTime_list[count_pps] = eventTime;
				count_pps++;
			}
		}
		
	
	
		  pps_chain->SetBranchAddress("T1H",T1H);
		  pps_chain->SetBranchAddress("T2H",T2H);
		  pps_chain->SetBranchAddress("T3H",T3H);
		  pps_chain->SetBranchAddress("T4H",T4H);
		  pps_chain->SetBranchAddress("T1V",T1V);
		  pps_chain->SetBranchAddress("T2V",T2V);
		  pps_chain->SetBranchAddress("T3V",T3V);
		  pps_chain->SetBranchAddress("T4V",T4V);
		  
		  
		  


		double eventRate = 0;
		double threshold = 0;
	
		TChain *trigger_chain = new TChain("t");
		string trigger_directory = trigger_path + to_string(gDate) + "/";
		GetTChain(trigger_directory, trigger_chain);
		const int trigger_entries = trigger_chain->GetEntries();
		//cout<<trigger_entries<<endl;
		trigger_chain->SetBranchAddress("eventTime", &eventTime);
		trigger_chain->SetBranchAddress("threshold",&threshold);
		trigger_chain->SetBranchAddress("eventRate",&eventRate);
    	

		int count_trigger=0;
		double trigger_eventTime_list[trigger_entries];
		double threshold_list[trigger_entries];
		double eventrate_list[trigger_entries];
		for(int i=0; i<trigger_entries; i++){
			trigger_chain->GetEntry(i);
			trigger_eventTime_list[count_trigger] = eventTime;
			threshold_list[count_trigger] = threshold;
			eventrate_list[count_trigger] = eventRate;
			count_trigger++;
		}
	
		TGraph *threshold_vs_time = new TGraph(count_trigger, trigger_eventTime_list, threshold_list);
		TGraph *eventrate_vs_time = new TGraph(count_trigger, trigger_eventTime_list, eventrate_list);
	
	
		for(int i=0; i<1000; i++){
			if(i%100==0) {cout<<i<<endl;}
			int random_pps = round(r3->Rndm()*count_pps);
			int random_entry = entry_list[random_pps];
			double random_time = pps_eventTime_list[random_pps];
			pps_chain->GetEntry(random_entry);
			if(eventTime != random_time){
				cout<<"Error!! eventTime should be the same!"<<endl;
				cout<<eventTime<<"   "<<random_time<<endl;
				continue;
			}

			gThresholdStep = round(threshold_vs_time->Eval(random_time));
			gEventRate = eventrate_vs_time->Eval(random_time);
			
			if(gThresholdStep<400){continue;}
			if(gEventRate<4){continue;}
		
			gEventNum = random_entry;
			cout<<gThresholdStep<<"   "<<gEventRate<<endl;
		
			pps_samples_Tree->Fill();
		}
	
		delete threshold_vs_time;
		delete eventrate_vs_time;
		//pps_chain->Close();
    	delete pps_chain;
    	//trigger_chain->Close();
    	delete trigger_chain;
    	delete [] pps_eventTime_list;
    	delete [] entry_list;
    }
	
	
	    

	pps_samples->cd();
    pps_samples_Tree->Write();
    pps_samples->Close();
    
	
}




void GetTChain(string directory, TChain *chaintree){

    DIR *dirp = NULL;
    struct dirent *dir_entry = NULL;
    char namebuf[100] = {0};
    if((dirp = opendir(directory.c_str())) == NULL){
        cout<<"Opendir fail!  "<<directory<<endl;
        return -1;
    }
    vector<string> files;
    while((dir_entry = readdir(dirp)) != NULL){
        string root_file_name = directory + dir_entry->d_name;
    	//int nPos = root_file_name.find_last_of(".root");
    	//cout<<root_file_name<<endl;
    	//cout<<nPos<<endl;
    	if(endsWith(root_file_name, ".root")){
    		files.push_back(root_file_name);
    	}
    	
    }
   
    closedir(dirp);
    
    sort(files.begin(), files.end());
    for(vector<string>::iterator it = files.begin(); it != files.end(); it++)
    	chaintree->AddFile((*it).c_str());
}
