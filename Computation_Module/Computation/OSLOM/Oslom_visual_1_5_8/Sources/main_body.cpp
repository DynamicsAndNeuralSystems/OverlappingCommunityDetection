







int main(int argc, char * argv[]) {
		
		
	if(argc<2) {
		program_statement(argv[0]);
		return -1;
	}

	
	
	if(paras._set_(argc, argv)==false)
		return -1;

	
	paras.print();
	
	if(paras.seed_random==-1)
		srand_file();
	else
		srand5(paras.seed_random);
		
	string netfile=paras.file1;
	clomnet luca(netfile, paras.threshold);
	
	
	//luca.draw("dd");

	
	
	char direc_char[1000];
	cast_string_to_char(paras.file1, direc_char);
	char char_to_use[1000];
	sprintf(char_to_use, "nohup mkdir %s_oslo_files > mk_hup", direc_char);
	int sy=system(char_to_use);
	sprintf(char_to_use, "nohup rm -r %s_oslo_files/* > mk_hup", direc_char);
	sy=system(char_to_use);

	
	cout<<"output files will be written in directory: "<<direc_char<<"_oslo_files"<<endl;
	
	
		
	
	
	
	cout<<"clomnet:: "<<luca.size()<<" nodes and "<<luca.stubs()<<" stubs;\t average degree = "<<luca.stubs()/luca.size()<<endl;
	
	if(luca.size()==0) {
		cerr<<"networkfile not found"<<endl;
		return -1;
		
	}
	LOG_TABLE._set_(cast_int(luca.stubs()));
	
	
	


	
	sprintf(char_to_use, "%s_oslo_files/every_group", direc_char);
	if(paras.every_group)
		fileout.open(char_to_use);
	
	
	
	if(paras.fast==false && paras.veryfast==false) {


		int soft_partitions_written=0;
			
					
		sprintf(char_to_use, "%s_oslo_files/partitions_level_0", direc_char);
		string partitions_level_0_string(char_to_use);
		ofstream plz_out(char_to_use);


		
		
		
		// ------------------------ collecting modules -----------------------------------------
		
		external_program_to_call(netfile, luca, plz_out, soft_partitions_written);
		
		time_t t=time(NULL);
		paras.initial_time=t;
		
		for(int i=0; i<paras.Or; i++) {
			
			
			t=time(NULL);
			int this_time=t;
			++this_time;
			
			
			if(paras.time_seconds<0 || this_time - paras.initial_time < paras.time_seconds)
				cout<<"RUN #"<<i<<endl;
			
			luca.gather(paras.gather_runs);
			
			if(luca.Mcoll.size()>0) {
				luca.post_processing(luca.Mcoll, true);
				luca.print_modules(true, plz_out);				// not homeless nodes
				soft_partitions_written++;
			}
		
		}
		
		
		
		if(paras.value) {
			
			luca.Mcoll.clear();
			luca.hint(paras.file2);
			
			if(luca.Mcoll.size()>0) {
				luca.post_processing(luca.Mcoll, true);
				luca.print_modules(true, plz_out);				// not homeless nodes
				soft_partitions_written++;
			}
		
		}
		
		
		if(paras.value_load) {
		
			luca.Mcoll.clear();
			luca.load(paras.file_load);

			if(luca.Mcoll.size()>0) {
				luca.post_processing(luca.Mcoll, true);
				luca.print_modules(true, plz_out);				// not homeless nodes
				soft_partitions_written++;
			}
			
			

		}


		
		
		
	
	
	

		// ------------------------ collecting modules -----------------------------------------

		
		
		cout<<"pruning all the modules collected. Partitions found: "<<soft_partitions_written<<endl;
		luca.Mcoll.clear();
		luca.load(partitions_level_0_string);
		if(soft_partitions_written>1)
			luca.post_processing(luca.Mcoll, true);
		
		
		
		cout<<"checking homeless nodes"<<endl;
		luca.try_to_assign_homeless(luca.Mcoll, false);
		
		
		
		
		
	
	} else {
		
		luca.fast_gather_and_process(false);
		
		int sy=system("mv louvain_fast.part louvain_fast_0.part");
		
		if(paras.veryfast==false) {
			luca.post_processing(luca.Mcoll, true);
			cout<<"checking homeless nodes"<<endl;
			luca.try_to_assign_homeless(luca.Mcoll, false);
		
		}
		
		
	
	}
	
	
	
	if(paras.homeless_anyway) {
				
		deque<int> homel;
		luca.Mcoll.homeless(homel);
	
		int before_procedure = homel.size();
		
		while(homel.size()>0) {
		
			cout<<"assigning homeless nodes: option -all. Homeless at this point: "<<before_procedure<<endl;
			
			luca.try_to_assign_homeless(luca.Mcoll, true);
			luca.Mcoll.homeless(homel);
			if(homel.size()>= before_procedure)
				break;
			before_procedure=homel.size();
		
		}
	
	}
	else {
	
		deque<set<int> > memberships_;				// this only refers to modules
		deque<deque<int> > modules_;		
		map<int, double> module_bs_;					// it mpas the module id into the b-score
		
		memberships_=luca.Mcoll.memberships;
		modules_=luca.Mcoll.modules;
		module_bs_=luca.Mcoll.module_bs;
		
		deque<int> homel;
		luca.Mcoll.homeless(homel);
	
		int before_procedure = homel.size();
		
		while(homel.size()>0) {
		
			//cout<<"assigning homeless nodes: option -all. Homeless at this point: "<<before_procedure<<endl;
			luca.try_to_assign_homeless(luca.Mcoll, true);
			luca.Mcoll.homeless(homel);
			if(homel.size()>= before_procedure)
				break;
			before_procedure=homel.size();
		
		}
		
		
		sprintf(char_to_use, "%s_oslo_files/tp_without_singletons", direc_char);	
		luca.print_modules(false, string(char_to_use));				// homeless nodes printed
			
		luca.Mcoll.memberships=memberships_;
		luca.Mcoll.modules=modules_;
		luca.Mcoll.module_bs=module_bs_;

	
	}
		
		
	
	
	
	sprintf(char_to_use, "%s_oslo_files/tp", direc_char);
	
	cout<<"writing final solution in file "<<char_to_use<<endl;
	luca.print_modules(false, string(char_to_use));				// homeless nodes printed
	sprintf(char_to_use, "cp %s_oslo_files/tp .", direc_char);
	int csy=system(char_to_use);
	
	sprintf(char_to_use, "%s_oslo_files/statistics_level_0.dat", direc_char);
	ofstream stout(char_to_use);
	luca.print_statistics(stout);
	

	
	
	if(paras.hierarchy)
		hierarchy(luca, direc_char);	
	
	
		
	
	
	return 0;
	
	

}



