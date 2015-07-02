


/************** this function is to replace the keyword NETx */

bool manipulate_string(string s, string netfile, string & outs) {



	outs.clear();
	string find_s("NETx");
	bool netx_found=false;
	
	
	int jj=0;
	string tmp;
		
	for (int j=0; j<s.size(); j++) {
		
		
		
		tmp.push_back(s[j]);

		
		
		if(s[j]==find_s[jj]) {
			
			jj++;
			
			if(jj==find_s.size()) {
			
				netx_found=true;
				tmp.clear();
				outs+= netfile;
				
			}
		
		
		}
		else {
			jj=0;
			outs+=tmp;
			tmp.clear();
		}
					
	
	
	
	
	
	}
	
	
	
	
	
	
	return netx_found;





}




/* this function is to call a different program */
void external_program_to_call(string network_file, clomnet & matteo, ostream & plz_out, int & soft_partitions_written) {



	for(int ei=0; ei<paras.to_run.size(); ei++) {
			
		
			
		
		string output_string;
		
		
		
		if(manipulate_string(paras.to_run[ei], network_file, output_string)==false) {
			
			cout<<"In string: "<<paras.to_run[ei]<<", keyword NETx was missing. The string cannot be run"<<endl;
		
		} else {
			
			
			
			
			string to_run_p;
			manipulate_string(paras.to_run_part[ei], network_file, to_run_p);

		
			char exec_this[2000];
			cast_string_to_char(output_string, exec_this);
			
			cout<<"running "<<exec_this<<endl;
			int sy=system(exec_this);
			
		
			matteo.Mcoll.clear();
			matteo.hint(to_run_p);
			
			
			if(matteo.Mcoll.size()>0) {
				matteo.post_processing(matteo.Mcoll, true);
				matteo.print_modules(true, plz_out);				// not homeless nodes	
				soft_partitions_written++;
			}
			
		}
	
	}



}



void main_collection(int level, clomnet & matteo, char * direc_char) {




	char char_to_use[1000];
	sprintf(char_to_use, "%s_oslo_files/net%d", direc_char, level);
	string sl2(char_to_use);

	
	int soft_partitions_written=0;

	if(paras.fast==false && paras.veryfast==false) {
			
		sprintf(char_to_use, "%s_oslo_files/partitions_level_%d", direc_char, level);
		string partitions_level_string(char_to_use);
		ofstream plz_out(char_to_use);

		
		sprintf(char_to_use, "%s_oslo_files/every_group_%d", direc_char, level);
		if(paras.every_group) {
			
			fileout.close();
			fileout.open(char_to_use);
	
		}
		
		
		cout<<"copying "<<sl2<<" to the main folder"<<endl;
		
			
		char char_to_copy[1000];
		sprintf(char_to_copy, "cp %s_oslo_files/net%d oslo_network_h", direc_char, level);
		int csy=system(char_to_copy);
		
		// ------------------------ collecting modules -----------------------------------------
		
		
		external_program_to_call("oslo_network_h", matteo, plz_out, soft_partitions_written);
		
		time_t t=time(NULL);
		paras.initial_time=t;
		
		for(int i=0; i<paras.hier_gather_runs; i++) {
			
			t=time(NULL);
			int this_time=t;
			++this_time;

			if(paras.time_seconds<0 || this_time - paras.initial_time < paras.time_seconds)
				cout<<"RUN #"<<i<<endl;
			
			matteo.gather(paras.gather_runs);
			
			if(matteo.Mcoll.size()>0) {
				matteo.post_processing(matteo.Mcoll, true);
				matteo.print_modules(true, plz_out);				// not homeless nodes	
				soft_partitions_written++;
			}
		
		}

		
		
		// ------------------------ collecting modules -----------------------------------------

		
		
		cout<<"pruning all the modules collected. Partitions found: "<<soft_partitions_written<<endl;
		matteo.Mcoll.clear();
		matteo.load(partitions_level_string);
		if(soft_partitions_written>1)
			matteo.post_processing(matteo.Mcoll, true);
		
		cout<<"checking homeless nodes"<<endl;
		matteo.try_to_assign_homeless(matteo.Mcoll, false);


	}
	else {
	
	
		matteo.fast_gather_and_process(false);
		
		if(paras.veryfast==false) {
		
			matteo.post_processing(matteo.Mcoll, true);
			cout<<"checking homeless nodes"<<endl;
			matteo.try_to_assign_homeless(matteo.Mcoll, false);
			
		
		}
			
		char mvl[1000];
		sprintf(mvl, "mv louvain_fast.part louvain_fast_%d.part", level);
		int sy=system(mvl);

		
	
	}
	



}















int hierarchy(clomnet & _luca_, char * direc_char) {
	

	map<int, set<int> > idmodule_idnodes;
	
	map<int, map<int, pair<int, double> > > neigh_weight_s;		// this maps the module id into the neighbor module ids and weights
	_luca_.module_id_first(idmodule_idnodes);
	
	
	
	_luca_.set_upper_network(neigh_weight_s, _luca_.Mcoll);
	
	int level=0;
	int previuous_dim=_luca_.size();
	
	while(true) {
	
	
			// upper level
			
			
			clomnet matteo(paras.threshold, neigh_weight_s);


			cout<<"clomnet level:: "<< ++level <<" ********************************** "<<endl;
			cout<<matteo.size()<<" nodes and "<<matteo.stubs()<<" edges;\t average degree = "<<matteo.stubs()/matteo.size()<<endl;
			
			char char_to_use[1000];
			sprintf(char_to_use, "%s_oslo_files/net%d", direc_char, level);
			string sl2(char_to_use);
			matteo.draw(sl2);
			
			
			if(neigh_weight_s.size()>=(1.-paras.hierarchy_convergence)*previuous_dim || previuous_dim<=2) {
				
				// this is the exiting strategy but it's not the best possible. It could be better to compare modules found at this level with those ones at the lower level (to be checked)
				cout<<"hierarchies done ********* "<<endl;
				break;
			}

			
			main_collection(level, matteo, direc_char);
			
			
			if(matteo.effective_groups()==0) {
				cout<<"hierarchies done ********* "<<endl;
				break;
			}
			
			
			
			
			sprintf(char_to_use, "%s_oslo_files/short_tp%d", direc_char, level);
			string sl(char_to_use);
			

			matteo.print_modules(false, sl);		// not homeless nodes
			
			previuous_dim=matteo.size();
			
			
			
			
			matteo.set_upper_network(neigh_weight_s, matteo.Mcoll);
			
			
			// ************************* printing stuff **************************************
			
			map<int, set<int> > idmodule_idnodes_here;
			matteo.module_id_later(idmodule_idnodes, idmodule_idnodes_here);
			
			sprintf(char_to_use, "%s_oslo_files/tp%d", direc_char, level);
			ofstream outt(char_to_use);
			for(map<int, set<int> >::iterator itm = idmodule_idnodes_here.begin(); itm!=idmodule_idnodes_here.end(); itm++) {

				outt<<"#module "<<itm->first<<" size: "<<idmodule_idnodes_here[itm->first].size()<<" bs: "<<matteo.Mcoll.module_bs[itm->first]<<endl;
				set<int> & a= itm->second;
				for(set<int>:: iterator its = a.begin(); its!=a.end(); its++ )
					outt<<*its<<" ";
				outt<<endl;
			
			}
			
			
			idmodule_idnodes= idmodule_idnodes_here;
			_luca_.Mcoll.clear();
			string stat_to_load(char_to_use);
			_luca_.load(stat_to_load);
			sprintf(char_to_use, "%s_oslo_files/statistics_level_%d.dat", direc_char, level);
			ofstream stout(char_to_use);
			_luca_.print_statistics(stout);
			
			// ************************* printing stuff **************************************
		
		
		
		
		}
	
	
	
	return 0;



}


