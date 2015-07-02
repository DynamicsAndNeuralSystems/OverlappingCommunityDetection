


// to insert a new parameter there are four steps:
// 1- define it in the class
// 2- initialize it
// 3- set the flag
// 4- set the parameter




void general_program_statement(char * b) {
	
	
	cout<<"USAGE: "<<b<<" -f network.dat -uw(-w)"<<endl<<endl;
	cout<<"-uw must be used if you want to use the unweighted null model; -w otherwise."<<endl;
	cout<<"network.dat is the list of edges. Please look at ReadMe.pdf for more details."<<endl<<endl;
	
	cout<<"***************************************************************************************************************************************************"<<endl;
	cout<<"RECOMMENDED OPTIONS"<<endl;
	
	
	
	//cout<<"\n[-noh]: it doesn't compute higher hierarchical levels. This option can be used if you have a big network and you want to start to see the first level partition."<<endl;
	
	cout<<"\n[-time hours]:\t\t\tthe maximum time for the seaching step (for each hierarchical level, postprocesses excluded) RECOMMENDED!"<<endl;

	cout<<"[-infomap runs]:\t\tcalls infomap and uses its output as a starting point. runs is the number of times you want to call infomap."<<endl;
	cout<<"[-copra runs]:\t\t\tsame as above using copra."<<endl;
	cout<<"[-louvain runs]:\t\tsame as above using louvain method."<<endl;
	cout<<"[-fast]:\t\t\tcalls only louvain method for a faster execution."<<endl;
	cout<<"[-veryfast]:\t\t\tcalls only louvain method and skips the post-processing."<<endl;
	
	
	cout<<"\n\nEXAMPLE:"<<endl;
	cout<<b<<" -f example.dat -uw -time 0.005 -infomap 3 -copra 2 -louvain 1 "<<endl;
	
	
	
	//cout<<"[-Opt optimization_flag]:\tdecides how much time to spend for the module search and how many programs to call, according to the accuracy requested (higher the optimization_flag, more accurate the results).";
	//cout<<"\n\t\t\t\tFor instance, use -Opt 0 for a fast (less accurate) execution, -Opt 24 for an accurate (slower) execution."<<endl;
	
	
	
	
	cout<<"\n\n\n";
	cout<<"***************************************************************************************************************************************************"<<endl;
	cout<<"ADVANCED OPTIONS"<<endl;	
	
	cout<<"\n[-t T]:\t\t\t\tsets the threshold equal to T, default value is 0.1"<<endl;
	cout<<"[-all]:\t\t\t\tno singletons. If you don't use this flag, the program generally finds a number of nodes which are not assigned to any module.\n\t\t\t\tIf you don't like these homeless guys, set this option: ";
	cout<<"the program will assign each node with at least one not homeless neighbor. This only applies to the lowest hierarchical level."<<endl;
	cout<<"[-r R]:\t\t\t\tsets the number of runs for the first hierarchical level, bigger this value, more accurate the output (of course, it takes more). Default value is 10."<<endl;
	cout<<"[-hr R]:\t\t\tsets the number of runs  for higher hierarchical levels. Default value is 100 (the method should be faster since the aggregated network is usually much smaller)."<<endl;
	cout<<"[-seed m]:\t\t\tsets m equal to the seed for the random number generator. (instead of reading from time_seed.dat)"<<endl; 
	cout<<"[-hint filename]:\t\ttakes a partition from filename. The file is expected to have the nodes belonging to the same cluster on the same line."<<endl;
	cout<<"[-load filename]:\t\ttakes modules from a tp file you already got in a previous run."<<endl;
	cout<<"[-e \"./my_program -options- NETx \" \"community_output\"]:\t\t\tthis option can be used to analyze the network with another program and automatically use the modules given by that."<<endl;

	
	cout<<"Please look at ReadMe.pdf for a detailed explanation."<<endl;
	cout<<"***************************************************************************************************************************************************"<<endl;



	cout<<"\n\n\n";
	cout<<"***************************************************************************************************************************************************"<<endl;
	cout<<"OUTPUT FILES"<<endl<<endl;
	
	
	cout<<"The program will create a directory called \"[network.dat]_oslo_files\". If the directory is not empty it will cleared, so be careful if you want to save some previous output files.\n"<<endl;
	cout<<"All the files will be written in this directory. "<<endl;
	cout<<"The first level partition will be written in a file called \"tp\", the next ";
	cout<<"hierchical network will be recorded as \"net1\", "<<endl;
	cout<<"the second level partition will be called \"tp1\" and so on."<<endl;
	cout<<"For convenience, the first level partition will be also written in a file called \"tp\" located in the same folder where the program is."<<endl;
	cout<<"***************************************************************************************************************************************************"<<endl;
	
	
	cout<<endl<<endl;
	cout<<"PLEASE LOOK AT ReadMe.pdf for more details. Thanks!"<<endl<<endl<<endl;
	
	
	
}




void error_statement(char * b) {

	cerr<<"\n\n************************************************************"<<endl;
	cerr<<"ERROR while reading parameters from command line... Please read program instructions or type: \n"<<b<<endl;
	cerr<<"************************************************************"<<endl;


}






class Parameters {
	
	
public:
	
	
	
	
	Parameters();
	~Parameters(){};
	void print();
	bool _set_(int argc, char * argv[]);
	
	
	string file1;
	string file2;
	string file_load;
	
	
	//*******************************************************
	
	bool weighted;
	bool hierarchy;
	bool homeless_anyway;

	double threshold;
	double part_rate;
	double part_rate_search;


	int Or;
	int hier_gather_runs;

	double coverage_inclusion_module_collection;
	double coverage_percentage_minimality;
	double coverage_percentage_fusion;		
	double check_inter_p;


	int setic_stopper_global;
	int setic_stopper_local;
	int CUT_Off;
	int bscore_loops;	
	int bscore_loops_check;
	int iterative_stopper;
	double gather_max_per100;
	double hierarchy_convergence;

	int gather_runs;
	int sub_gather_runs;



	int maxborder_nodes;
	double maxbg_ordinary;	
	double maxbg_condensed;
	double max_self_loops;





	bool list;
	bool value;
	bool value_load;
	bool print_cbs;
	bool every_group;
	
	deque<string> to_run;
	deque<string> to_run_part;
	
	
	int seed_random;
	int infomap_runs;
	int copra_runs;
	int louvain_runs;

	
	
	int time_seconds;
	int initial_time;
	
	int optimization_flag;
	bool fast;
	bool veryfast;
	bool hint_flag;
	
	
	
private:
	
	map<string, int> command_flags;
	bool set_flag_and_number(string program_name, int & argct, int & number_to_set, int argc, char * argv[]);
	bool set_flag_and_number(string program_name, int & argct, double & number_to_set, int argc, char * argv[]);
	
	
};



bool Parameters::set_flag_and_number(string program_name, int & argct, int & number_to_set, int argc, char * argv[]) {



	
	argct++;
	if(argct==argc) {
		
		cout<<"you didn't set the number of "<<program_name<<endl;
		
		error_statement(argv[0]);
		return false;
	}
	
	string tt=argv[argct];
	double ttt;
	if(cast_string_to_double(tt, ttt)==false) {
	
		cout<<"you didn't set the number of "<<program_name<<endl;

	
		error_statement(argv[0]);
		return false;
	}
	
	number_to_set=cast_int(ttt);
	
	if(number_to_set<0) {
	
		cout<<" the number of "<<program_name<<" must be positive"<<endl;

		error_statement(argv[0]);
		return false;
	}
	
	
	
	return true;
	
	
	
}




bool Parameters::set_flag_and_number(string program_name, int & argct, double & number_to_set, int argc, char * argv[]) {



	
	argct++;
	if(argct==argc) {
		
		cout<<"you didn't set the number of "<<program_name<<endl;
		
		error_statement(argv[0]);
		return false;
	}
	
	string tt=argv[argct];
	double ttt;
	if(cast_string_to_double(tt, ttt)==false) {
	
		cout<<"you didn't set the number of "<<program_name<<endl;

	
		error_statement(argv[0]);
		return false;
	}
	
	number_to_set=ttt;
	
	if(number_to_set<0) {
	
		cout<<" the number of "<<program_name<<" must be positive"<<endl;

		error_statement(argv[0]);
		return false;
	}
	
	
	
	return true;
	
	
	
}




Parameters::Parameters() {
	
	//**************************************************************************
	
	seed_random=-1;
	
	
	weighted=false;
	hierarchy=true;
	homeless_anyway=false;

	threshold= 0.1;					// this is the P-value for the significance of the module



	
	part_rate=0.5;							// this is the threshold value of the clean-up procedure. if you decrease this, you get more noise
	part_rate_search=0.05;					// this is the threshold value of the clean-up procedure (in the search mode)
	
	
		
	
	Or=10;							// this is the number of global runs in the gather function		(first level)
	hier_gather_runs=100;			// this is the number of global runs in the gather function		(higher level)
	
	
	
	coverage_inclusion_module_collection=0.49999;			// this is used to see if two modules are higly similar in processing the clusters (big_module)
	coverage_percentage_minimality=0.7;						// if you increase this value the program tends to find bigger clusters
	coverage_percentage_fusion=0.8;							// this is used to see when fusing clusters how much is left
	check_inter_p=0.333;									// this parameter is a check parameter for the fusion of quite similar clusters
	
	
	
	// all these are parameters for the convergence of the search process - or the clean up procedure - all this could be improved 
	setic_stopper_global=10;
	setic_stopper_local=10;
	CUT_Off=200;
	bscore_loops=3;
	bscore_loops_check=300;
	iterative_stopper=10;							// this is to prevent the iterative procedure to last too long. this can happen in case of strong backbones (just an idea, not sure)
	gather_max_per100= 0.91;						// if the user wants to use the option time, a single run stops when this percentage of the network has been reached
	hierarchy_convergence=0.05;						// this parameter is used to stop the hierarchical process when not enough modules are found
	
	gather_runs=1;					// this is the number of runs in the gather function (parameter for the convergence)
	sub_gather_runs=2;				// this is the number of runs in the gather sub function (parameter for the convergence)		

	
	// these is stuff to speed up the code in looking for "good" neighbors
	maxborder_nodes=50;
	maxbg_ordinary=0.5;
	maxbg_condensed=0.4;
	max_self_loops=1.;
	
	
	
	//	these are some flags to read input files
	list=false;
	value=false;
	value_load=false;
	print_cbs=false;
	every_group=true;
		
	//	these are some flags to read input files
	infomap_runs=0;
	copra_runs=0;
	louvain_runs=0;
	optimization_flag=-1;
	
	time_seconds=-1;
	
	veryfast=false;
	fast=false;
	hint_flag=false;
	
	//**************************************************************************
	
	command_flags.insert(make_pair("-f", 1));
	command_flags.insert(make_pair("-hint", 2));
	command_flags.insert(make_pair("-load", 11));
	command_flags.insert(make_pair("-t", 3));
	command_flags.insert(make_pair("-r", 13));
	command_flags.insert(make_pair("-hr", 17));
	command_flags.insert(make_pair("-w", 14));
	command_flags.insert(make_pair("-uw", 15));
	command_flags.insert(make_pair("-e", 16));
	command_flags.insert(make_pair("-all", 9));
	command_flags.insert(make_pair("-seed", 18));
	command_flags.insert(make_pair("-infomap", 19));
	command_flags.insert(make_pair("-copra", 20));
	command_flags.insert(make_pair("-louvain", 21));
	command_flags.insert(make_pair("-time", 22));
	command_flags.insert(make_pair("-fast", 23));
	command_flags.insert(make_pair("-veryfast", 27));
	
	
	//command_flags.insert(make_pair("-Opt", 24));
	//command_flags.insert(make_pair("-moses", 26));
	//command_flags.insert(make_pair("-noh", 10));
	//command_flags.insert(make_pair("-v", 4));
	//command_flags.insert(make_pair("-r", 5));
	//command_flags.insert(make_pair("-list", 7));
	//command_flags.insert(make_pair("-maxborder", 8));
	//command_flags.insert(make_pair("-home", 9));
	//command_flags.insert(make_pair("-no-every", 12));

	
	
}



void Parameters::print() {
	
	
	cout<<"**************************************"<<endl;
	cout<<"Threshold:\t\t\t"<<threshold<<endl;
	cout<<"Network file:\t\t\t"<<file1<<endl;
	
	
	if(weighted)
		cout<<"Weighted: yes"<<endl;
	else
		cout<<"Weighted: no"<<endl;
		
	
	if(optimization_flag>=0)
		cout<<"Optimization flag: "<<optimization_flag<<endl;
	

	if(time_seconds>=0)
		cout<<"Maximum time for the module search: "<<time_seconds<<" secs"<<endl;
	
		
	for(int i=0; i<to_run.size(); i++)
		cout<<"String to run: ["<<to_run[i]<<"]\t\t\t\t\t\tModule file: ["<<to_run_part[i]<<"]"<<endl;
	

	if(fast)
		cout<<"-fast option selected"<<endl;
	if(veryfast)
		cout<<"-veryfast option selected"<<endl;

		
	if(value)
		cout<<"Hint from file:\t\t\t"<<file2<<endl;
	if(value_load)
		cout<<"tp-file:\t\t\t"<<file_load<<endl;

	cout<<"First Level Run:\t\t\t"<<Or<<endl;
	cout<<"Higher Level Run:\t\t\t"<<hier_gather_runs<<endl;

	cout<<"**************************************"<<endl<<endl;
	
	
	
	
	
}






bool Parameters::_set_(int argc, char * argv[]) {
	
	int argct = 0;
	string temp;
	
	if (argc <= 1) { // if no arguments, return error_statement about program usage.
		
		error_statement(argv[0]);
		return false;
	}
	
	
	bool f_set=false;
	bool set_weighted=false;

	
	
	while (++argct < argc) {			// input file name
	
		
		
		cout<<"setting "<<argv[argct]<<endl;
		
		temp = argv[argct];
		map<string, int>::iterator itf=command_flags.find(temp);
		
		if(itf==command_flags.end()) {
			error_statement(argv[0]);
			return false;
		}
		
		int vp=itf->second;
		if(vp==1) {
			
			argct++;
			if(argct==argc) {
				error_statement(argv[0]);
				return false;
			}
			
			file1=argv[argct];
			f_set=true;
			
		}
		else if(vp==2) {
			
			argct++;
			if(argct==argc) {
				error_statement(argv[0]);
				return false;
			}
			
			file2=argv[argct];
			value=true;
			
		}
		else if(vp==11) {
			
			argct++;
			if(argct==argc) {
				error_statement(argv[0]);
				return false;
			}
			
			file_load=argv[argct];
			value_load=true;
			
		}
		else if(vp==3) {
			
			argct++;
			if(argct==argc) {
				error_statement(argv[0]);
				return false;
			}
			
			string tt=argv[argct];
			if(cast_string_to_double(tt, threshold)==false) {
				error_statement(argv[0]);
				return false;
			}
			
		}
		else if(vp==4) {
			
			print_cbs=true;
			
		}
		else if(vp==12) {
			
			every_group=true;
			
		}
		else if(vp==14) {
			
			weighted=true;
			set_weighted=true;
			
		}else if(vp==15) {
			
			weighted=false;
			set_weighted=true;
			
		}
		else if(vp==9) {
			
			homeless_anyway=true;
		}
		else if(vp==23) {
			
			fast=true;
		}
		else if(vp==27) {
			
			veryfast=true;
		}
		else if(vp==7) {
			
			list=true;
			
		}
		else if(vp==10) {
			
			hierarchy=false;
			
		}
		else if(vp==5) {
			
			argct++;
			if(argct==argc) {
				error_statement(argv[0]);
				return false;
			}
			
			string tt=argv[argct];
			double ttt;
			if(cast_string_to_double(tt, ttt)==false) {
				error_statement(argv[0]);
				return false;
			}
			
			gather_runs=cast_int(ttt);
			
		}
		else if(vp==6) {
			
			argct++;
			if(argct==argc) {
				error_statement(argv[0]);
				return false;
			}
			
			string tt=argv[argct];
			if(cast_string_to_double(tt, coverage_percentage_minimality)==false) {
				error_statement(argv[0]);
				return false;
			}
			
			
		}
		else if(vp==8) {
			
			argct++;
			if(argct==argc) {
				error_statement(argv[0]);
				return false;
			}
			
			string tt=argv[argct];
			double ttt;
			if(cast_string_to_double(tt, ttt)==false) {
				error_statement(argv[0]);
				return false;
			}
			
			maxborder_nodes=cast_int(ttt);
			
		}
		
		else if(vp==13) {
			
			argct++;
			if(argct==argc) {
				error_statement(argv[0]);
				return false;
			}
			
			string tt=argv[argct];
			double ttt;
			if(cast_string_to_double(tt, ttt)==false) {
				error_statement(argv[0]);
				return false;
			}
			
			Or=cast_int(ttt);
			
		}
		else if(vp==18) {
			
			argct++;
			if(argct==argc) {
				
				cout<<"you didn't set the seed for the random number generator"<<endl;
				
				error_statement(argv[0]);
				return false;
			}
			
			string tt=argv[argct];
			double ttt;
			if(cast_string_to_double(tt, ttt)==false) {
			
				cout<<"you didn't set the seed for the random number generator"<<endl;

			
				error_statement(argv[0]);
				return false;
			}
			
			seed_random=cast_int(ttt);
			
			if(seed_random<=0) {
			
				cout<<"the seed for the random number generator must be positive"<<endl;

				error_statement(argv[0]);
				return false;
			}
			
			
			
			
			
		}
		else if(vp==19) {
			
			if(set_flag_and_number("runs for infomap", argct, infomap_runs, argc,  argv)==false)
				return false;
		}
		else if(vp==20) {
			
			if(set_flag_and_number("runs for copra", argct, copra_runs, argc,  argv)==false)
				return false;
		}
		else if(vp==21) {
			
			if(set_flag_and_number("runs for louvain method", argct, louvain_runs, argc,  argv)==false)
				return false;
		}
		else if(vp==22) {
			
			
			double hours;
			if(set_flag_and_number("hours (option -time)", argct, hours, argc,  argv)==false)
				return false;
			
			time_seconds=cast_int(hours*3600);
			
		}
		else if(vp==24) {
			
			
			if(set_flag_and_number("optimization flag (option -Opt)", argct, optimization_flag, argc,  argv)==false)
				return false;
		}
		else if(vp==17) {
			
			argct++;
			if(argct==argc) {
				error_statement(argv[0]);
				return false;
			}
			
			string tt=argv[argct];
			double ttt;
			if(cast_string_to_double(tt, ttt)==false) {
				error_statement(argv[0]);
				return false;
			}
			
			hier_gather_runs=cast_int(ttt);
			
		}
		else if(vp==16) {
			
			argct++;
			
			if(argct==argc) {
				
				cout<<" string with \"my_program\" is missing"<<endl;
				
				error_statement(argv[0]);
				return false;
			}
			
			string tr1=argv[argct];
			argct++;
			
			if(argct==argc) {
			
				cout<<" string with the output file of \"my_program\" is missing"<<endl;

				error_statement(argv[0]);
				return false;
			}
			string tr2=argv[argct];
			
			to_run.push_back(tr1);
			to_run_part.push_back(tr2);
			
			
		}
		else {
			
			error_statement(argv[0]);
			return false;
		}
		
		
		cout<<" ... done"<<endl;
		
	}
	
	if(f_set==false) {
		
		cerr<<"\n\n************************************************************"<<endl;
		cout<<"ERROR: you didn't set the file with the network.  Please read program instructions or type: \n"<<argv[0]<<endl;
		cerr<<"************************************************************"<<endl;
		
		return false;
		
	}
	
	if(set_weighted==false) {	
		
		
		cerr<<"\n\n************************************************************"<<endl;
		cout<<"ERROR: you didn't set the option -w (weighted network) or -uw (unweighted network).  Please read program instructions or type: \n"<<argv[0]<<endl;
		cerr<<"************************************************************"<<endl;
		
		return false;
		
	}
	
	
	if(weighted==false) {
	
		maxbg_ordinary=0.2;
		maxbg_condensed=0.03;
	
	}
	
	
	
	if(optimization_flag>=0) {
	
		if(optimization_flag==0) {
		
			if(time_seconds<0)
				time_seconds=0;
			
			infomap_runs=1;
			copra_runs=1;
			louvain_runs=1;

		}

		else if(optimization_flag==1) {
		
			if(time_seconds<0)
				time_seconds=3600;
			
			infomap_runs=1;
			copra_runs=1;
			louvain_runs=1;

		}
		
		else if(optimization_flag==2) {
		
			if(time_seconds<0)
				time_seconds=3600 * (optimization_flag);
			
			infomap_runs=3;
			copra_runs=3;
			louvain_runs=1;

		}
	
	}

	
	
			
	for(int i=0; i<infomap_runs; i++) {
		
		char number_r[1000];
		
		//cout<<"************** "<<string(argv[0])<<endl;
		string pros(argv[0]);
		
		if(pros=="./oslom_undir")
			sprintf(number_r, "./infomap_undir_script NETx %d %d", irand(10000000), 1);
		else
			sprintf(number_r, "./infomap_dir_script NETx %d %d", irand(10000000), 1);
		
		string sr(number_r);
		
		//cout<<"here "<<endl;
		to_run.push_back(sr);
		to_run_part.push_back("infomap.part");
	
	}
	
	
	for(int i=0; i<copra_runs; i++) {
		
		char number_r[1000];
		sprintf(number_r, "java -cp copra.jar COPRA NETx -v 5 -w");			
		string sr(number_r);
		
		
		to_run.push_back(sr);
		to_run_part.push_back("clusters-NETx");
	
	}

		
	for(int i=0; i<louvain_runs; i++) {
		
		char number_r[1000];
		sprintf(number_r, "./louvain_script -f NETx");			
		string sr(number_r);
		
		
		to_run.push_back(sr);
		to_run_part.push_back("louvain.part");
	
	}
	
	
	
		
	
	
	
	time_t t=time(NULL);
	initial_time=t;
	
	
	if(time_seconds>=0 && to_run.size()==0) {
		
		cout<<"\nWARNING::: you set option -time and you didn't use any option for external programs."<<endl;
		cout<<"THIS is not very recommendable because if the program runs out of time some pieces of the network could be not explored"<<endl;
		cout<<"I would use options: -infomap, -copra, -louvain or at least one of them"<<endl<<endl;
	
	
	}
	
	
	if(fast)
		bscore_loops_check=125;
	
	if(veryfast)
		bscore_loops_check=15;

	
		
	
	return true;
}




