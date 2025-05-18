void test_tchain(){

	TChain *tr_event = new TChain("t");
	TChain *tr_config = new TChain("tHeader");
	
	
	ifstream inFile("run_list.txt");
	string str;
	int fNFile = 0;
	vector<Long64_t> fvNEntries;
	vector<int> run_number;

	while(getline(inFile, str)){

		cout<<"reading: "<< str <<endl;
		if( str.find(".root") != std::string::npos)
		{
		    int nPos_1 = str.find_last_of("n");
            int nPos_2 = str.find_last_of(".");
            string run_num="";
            run_num = str.substr(nPos_1+1,nPos_2-nPos_1-1);
            run_number.push_back(std::stoi(run_num));
            
			fNFile +=  tr_event->AddFile( str.c_str() );
			tr_config->AddFile( str.c_str() );

			fvNEntries.push_back( tr_event->GetEntries() );
		}


	}

	inFile.close();
	
	for(auto context : fvNEntries){
		cout<<context<<endl;
		}
	for(auto context : run_number){
		cout<<context<<endl;
		}
	cout<<endl;
	cout<<fNFile<<endl;
	cout<<tr_config->GetEntries()<<endl;
	cout<<tr_event->GetEntries()<<endl;
	
	int run_N = fvNEntries.size();
	int i = 0;
	int tchain_event = 10100;
	int run_entry = 0;
	int which_run = 0;
	int which_event = 0;
	for( ; i<run_N; i++){
		if(tchain_event < fvNEntries.at(i)) {break;}
	}
	cout<<i<<endl;
	 which_run = run_number.at(i);
	 
	 if(i==0){
	 	which_event = tchain_event;
	 	}
	 else{
	 	which_event = tchain_event-fvNEntries.at(i-1);
	 	}
	 	
	 cout<<which_run<<"  "<<which_event<<endl;
}



