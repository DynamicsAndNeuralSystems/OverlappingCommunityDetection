




//#include "clomnet_local.h"
#define maximum_level 2




void set_changendi_cum(deque<double> & changendi_cum, double tstrength, int dim, int cutoff) {
	
	
	
	int flat_until= cast_int(2*tstrength/dim * 3);
	flat_until=min(dim/2, flat_until);
	
	int max_p= max(cutoff, flat_until);		// this is something which might be optimized
	max_p=min(dim/2, max_p);
	
	
	
	powerlaw(max_p, flat_until+1, 3, changendi_cum);
	deque<double> distr;
	distribution_from_cumulative(changendi_cum, distr);
	double ac=1;
	
	if(distr.size()>0)
		ac=distr[0];
	
	
	for(int i=0; i<flat_until; i++)
		distr.push_front(ac);
	
	
	
	normalize_one(distr);
	
	
	cumulative_from_distribution(changendi_cum, distr);
	
	


}


class clomnet : public clomnet_local {

	
	public:
		
		clomnet(deque<deque<int> > & b, deque<deque<pair<int, double> > > & c, deque<int> & d, double cl);
		clomnet(string a, double cl);
		clomnet(double cl, map<int, map<int, pair<int, double> > > & A);
		~clomnet(){};
		
		int gather(int );
		int hint(string);
		void print_modules(bool not_homeless, string tp);
		void print_modules(bool not_homeless, ostream &);
		int fill_homel();
		int effective_groups();
		
		int load(string filename);
		
		void module_id_first(map<int, set<int> > & A);
		void module_id_later(map<int, set<int> > & A, map<int, set<int> > & B);
		
				
		int LEVEL;

		int matteo_pairs();
		//int try_to_assign_homeless();
	


		int checklevel(module_collection &);
		int post_processing();
		int post_processing(module_collection & mall, bool post_p_print);
	
		int check_ex_pieces(module_collection & mall);
		
		
		// ********** things to be used in gather_ **********
		module_collection Mcoll;
		// ********** things to be used in gather_ **********

		bool fast_minimality_bg(deque<int> & group, double bs, ne_modules &);
		int fast_gather_and_process(bool);

		void print_statistics(ostream & outt);

	
	private:		
		
		double clean_t;

		
		int check_intersection();
		int check_intersection(deque<int> & to_check);
		int fusion_intersection(set<pair<int, int> > & pairs_to_check);
		inline bool fusion_intersection(int ai1, int ai2, deque<int> & new_insertions);
		
		//int check_pairs(multimap<int, int> & A);
		int upper_fusion(set<pair<int, int> > & p1);
		double fusion(deque<int> & a1, deque<int> & a2, deque<int> & group);
		bool fusion_module_its_subs(deque<int> & a, deque<deque<int> > & its_submodules);

		/*****************************/
		
		
		
		
		
		/*****************************/
		inline int single_run(set<int> & uncovered, ne_modules & already_found, int max_stopper, bool suubb, bool fshell);

		int single_run_gather(ne_modules & already_found, int max_stopper, bool);
		int gather(int number_of_runs, bool f);
		bool minimality_bg(deque<int> & group, double bs, set<int> & uncovered, bool g_print, ne_modules & already_found);
		bool minimality_bg(deque<int> & group, double bs, set<int> & uncovered, bool g_print);
		deque<double> setic_changendi_cum;
		/*****************************/		
		int check_too_big(deque<deque<int> > & ten, int rowi, int max_size);
		
};





clomnet::clomnet(double cl, map<int, map<int, pair<int, double> > > & A) : clomnet_local(A) { 

	clean_t=cl;
	Mcoll._set_(dim, double(oneM)/dim);
	LEVEL=0;
	
}




clomnet::clomnet(deque<deque<int> > & b, deque<deque<pair<int, double> > > & c, deque<int> & d, double cl) : clomnet_local(b, c, d) { 

	clean_t=cl;
	Mcoll._set_(dim, double(oneM)/dim);
	LEVEL=0;
	
}



clomnet::clomnet(string a, double cl) : clomnet_local(a) { 
	
	clean_t=cl; 
	Mcoll._set_(dim, double(oneM)/dim);
	LEVEL=0;
	
}



int clomnet::gather(int number_of_runs, bool suubb) {
	
	Mcoll.clear();


	
	if(LEVEL>maximum_level)
		return 0;
	
	//if(paras.hint_flag)
		//cout<<"hint flag"<<endl;
	
	
	if(paras.hint_flag==false && paras.time_seconds>=0) {
		
		time_t t=time(NULL);
		int this_time=t;
		
		++this_time;
	
		if(this_time - paras.initial_time > paras.time_seconds)
			return 0;
	
	}
	
	
	
	//cout<<"gathering"<<endl;
	
		
	if(dim==0) {
		
		cout<<"network empty"<<endl;
		cout<<"STOP"<<endl;
		//cherr();
		return -1;
	}
	
	
	
	setic_changendi_cum.clear();
	set_changendi_cum(setic_changendi_cum, double(oneM)/2., dim, paras.CUT_Off);
	
	
	
	
	ne_modules already_found(dim);
	

	
	if(suubb==false)
		cout<<"starting to collect modules: "<<endl;
	
	int stopper=0;
	while(true) {
		
		
		if(stopper>=number_of_runs)
			break;
		
		
		/*if(stopper==0)
			Mcoll.average_degree=(5.*oneM)/dim;
		else
			Mcoll.average_degree=(2.*oneM)/dim;*/
		
		if(suubb==false)
			cout<<"\n\n\n*********************** starting a new search on a network of "<<dim<<" nodes "<<already_found.size()<<" clusters to be tested found"<<endl;
		
		//if(stopper%UM_num==0 && suubb==false)
		if(stopper==0 && suubb==false)
			single_run_gather(already_found, dim/3, suubb);
		else if(stopper==0 && suubb==true)
			single_run_gather(already_found, dim/10, suubb);
		else
			single_run_gather(already_found, paras.setic_stopper_global, suubb);
		
		
		stopper++;
		
		
	}
	
	
	
	return 0;

}


int clomnet::post_processing() {

	checklevel(Mcoll);
	try_to_assign_homeless(Mcoll, false);
	
	return 0;
	
}


int clomnet::post_processing(module_collection & mall, bool post_p_print) {
	
	
	if(post_p_print)
		cout<<"checking overlaps and inclusions between groups"<<endl;
	checklevel(mall);
	
	
	
	/*if(post_p_print)
		cout<<"checking homeless nodes"<<endl;
	
	try_to_assign_homeless(mall, false);*/
	
	
	
	return 0;
	
}
















//*
int clomnet::gather(int a) {

	return gather(a, false);
	
}//*/





inline int clomnet::single_run(set<int> & uncovered, ne_modules & already_found, int max_stopper, bool suubb, bool fshell) {
	
	
	int initial_dim=uncovered.size();
	
	int stopper=0;
	int groups_found=0;


	while(!uncovered.empty()) {

		
		
		int A=random_from_set(uncovered);
		//cout<<"A: "<<A<<endl;
		
		deque<int> group;
		double bs= setic_simann(A, group, clean_t, setic_changendi_cum, 1, paras.setic_stopper_local, suubb);
		
		
		/*cout<<"GROUP:::"<<endl;
		print_id(group, cout);*/
		
		//cherr();
		
		//if(group.size()>0  && already_found.insert(group)==true ) {
		if(group.size()>0) {
		
			
					
			deque<int> grc;
			bs=b_clean_up_search_iterative(group, grc, paras.bscore_loops_check);
			
			
			if(grc.size()>0 && minimality_bg(grc, bs, uncovered, !suubb, already_found))
				stopper=0;
			else
				stopper++;
		}
		else {
			
			stopper++;
			for(int j=0; j<group.size(); j++)
				uncovered.erase(group[j]);
		
		}
		uncovered.erase(A);
		
		
		if(stopper>=max_stopper)
			break;
		
		if(fshell==false) {
			
			if(suubb==false) {
				cout<<"main-percentage covered "<<1 - double(uncovered.size())/initial_dim<<"\t\tgood clusters so far:\t"<<Mcoll.size()<<"\ttotal trial cluster: "<<already_found.size()<<
				"\t\t\t--- steps since the last good "<<stopper<<endl;
				
				
				if(paras.time_seconds>=0) {
					
					if(1 - double(uncovered.size())/initial_dim > paras.gather_max_per100)
						return 0;
					
					
					time_t t=time(NULL);
					int this_time=t;
					
					++this_time;
					
					if(this_time - paras.initial_time > paras.time_seconds)
						return 0;
				
				}

				
					
			}
			
			/*else if(dim>min_print_s && paras.print_cbs && ran4()<0.01)
				cout<<"sub-percentage covered "<<1 - double(uncovered.size())/initial_dim<<" good clusters so far:\t"<<Mcoll.size()<<"\ttotal trial cluster: "
				<<already_found.size()<<"\t\t\t---stopper "<<stopper<<" dim: "<<dim<<endl;*/
			//cout<<"percentage covered "<<1 - double(uncovered.size())/dim<<" good clusters so far:\t"<<Mcoll.size()<<"\ttotal trial cluster: "<<already_found.size()<<"\t\t\t---stopper "<<stopper<<endl;
		}
			
		
	
	}
	
	
	return 0;

}





int clomnet::single_run_gather(ne_modules & already_found, int max_stopper, bool suubb) {

	
	
	
	// ************* this is done for one single run
	
	set<int> uncovered;
	for(int i=0; i<dim; i++)		// I think this should be more important: this is where we choose the initial seeds
		uncovered.insert(i);
	
	
	
	single_run(uncovered, already_found, max_stopper, suubb, false);
	
	
	return 0;
	
}





/*
int clomnet::first_shell(ne_modules & already_found, bool suubb) {

	
	
	// ************* this is done for one single run
	
	
	deque<int> homel;
	Mcoll.homeless(homel);

	
	set<int> uncovered;
	for(int ii=0; ii<homel.size(); ii++) if(vertices[homel[ii]]->links->size()>5 && vertices[homel[ii]]->links->size()<dim-1) {
	
	
		set<int> G;
		for(int i=0; i<vertices[homel[ii]]->links->size(); i++)
			G.insert(vertices[homel[ii]]->links->l[i]);
			
		G.insert(homel[ii]);
			
			
		if(kin_m(G)<2.5*G.size())
			uncovered.insert(homel[ii]);
	
	}
	
	if(suubb==false && uncovered.size()>0)
		cout<<"\n\n\n***********************"<<endl<<"checking homeless nodes: "<<uncovered.size()<<" nodes to check"<<endl<<endl;

		
	
	
	
	single_run(uncovered, already_found, homel.size(), suubb, true);
	
	
	/*Mcoll.homeless(homel);
	cout<<"still_homeless: "<<homel.size()<<endl;* /
	
	return 0;
	
	



}
*/




/*
int clomnet::trying(ne_modules & already_found, deque<int> & group) {


	set<int> uncovered;
			
	if(already_found.insert(group)==true ) {
	

		deque<int> grc;
		double bs=b_clean_up_check(group, grc, paras.bscore_loops_check);
		
		//if(paras.print_cbs)
		//cout<<"FIRST SHELL - GROUP OF "<<grc.size()<<" nodes, out of "<<dim<<endl;
				
		if(grc.size()>0) {
			minimality_bg(grc, bs, uncovered);
			return 0;
		}
		
				
	}


	return -1;

}
*/






int clomnet::fast_gather_and_process(bool suubb) {
	
		
	if(LEVEL>maximum_level)
		return 0;

	
	deque<deque<int> >  ten;
	draw_consecutive("netconsec.dat", "labelconsec.dat");
	
	
	cout<<"calling louvain method: "<<dim<<endl;
	int sy= system("./convert -i netconsec.dat -o network.bin -w");
	sy=system("./community network.bin -l -1 -w > network.tree");
	sy = system("./hierarchy network.tree -l 1 > leve1.dat");
		
	get_partition_from_file_list("leve1.dat", ten);
	cout<<ten.size()<<" modules found"<<endl;
	
			
	for(int i=0; i<ten.size(); i++)
		check_too_big(ten, i, 10000);
	
	
	
	
	if(suubb==false) {
		
		cout<<ten.size()<<" modules found on a network of "<<dim<<" nodes"<<endl;
		ofstream outt("louvain_fast.part");
		print_id(ten, outt);
	
	}
	
	set<int> uncovered;
	ne_modules already_found(dim);

	for(int i=0; i<ten.size(); i++) {
		
		
		deque<int> grc;
		
		if(suubb==false && (i%100==0))
			cout<<"processed modules: "<<double(i)/ten.size()<<"\t good modules found "<<Mcoll.size()<<"\t outf of "<<i<<endl;		
		
		
		if(suubb==false)
			cout<<i<<". trying group of "<<ten[i].size()<<" nodes on a net of size: "<<dim<<endl;
			
	
		double bs=b_clean_up_search(ten[i], grc, paras.bscore_loops_check);
		
		
		if(suubb==false)
			cout<<i<<". NEW GROUP OF "<<grc.size()<<" nodes "<<endl;

		if(grc.size()>0) {
			
			if(grc.size()<100)
				minimality_bg(grc, bs, uncovered, false, already_found);
			else
				fast_minimality_bg(grc, bs, already_found);
		}
	
	}

	
	Mcoll.compute_inclusions(false);
	return 0;
	

}



int clomnet::check_too_big(deque<deque<int> > & ten, int rowi, int max_size) {

	
	if(ten[rowi].size()>max_size) {
	
		
		deque<deque<int> > link_per_node;
		deque<deque<pair<int, double> > > weights_per_node;
		set_subgraph(ten[rowi], link_per_node, weights_per_node);
		clomnet giovanni(link_per_node, weights_per_node, ten[rowi], clean_t);
		giovanni.draw_consecutive("netconsec.dat", "labelconsec.dat");


		int sy= system("./convert -i netconsec.dat -o network.bin -w");
		sy=system("./community network.bin -l -1 -w > network.tree");
		sy = system("./hierarchy network.tree -l 1 > leve1.dat");
		
		deque<deque<int> > tens;
		
		get_partition_from_file_list("leve1.dat", tens);
		
		if(tens.size()>1) {
			
			cout<<"module of "<<ten[rowi].size()<<" nodes split in "<<tens.size()<<" modules"<<endl;
			
			ten[rowi].clear();
			for(int i=0; i<tens.size(); i++) {
				
				giovanni.deque_id(tens[i]);
				ten.push_back(tens[i]);
				
			}
		}
		
		
	}
	
	
	
	return 0;

}

int clomnet::hint(string filename) {

	
	paras.hint_flag=true;
	
	cout<<"getting partition from file: "<<filename<<endl;
	deque<deque<int> > A;
	
	if(paras.list==false)
		get_partition_from_file(filename, A);
	else	
		get_partition_from_file_list(filename, A);
	
	
	
	//printm(A);
	if(translate(A)==-1)
		return -1;
	
	cout<<A.size()<<" groups found"<<endl;
	//ofstream ttr("translated");
	//printm(A, ttr);
	//cherr();
	
	ne_modules already_found(dim);

	set<int> uncovered;
	for(int ii=0; ii<A.size(); ii++) /*if(A[ii].size()==474)*/ {

		deque<int> group;
		
		cout<<"processing group number "<<ii<<" size: "<<A[ii].size()<<endl;//" kin: "<<double(kin_m(A[ii]))/ktot_m(A[ii])<<endl;

		double bcu=b_clean_up_search(A[ii], group, paras.bscore_loops_check);
		
		if(group.size()>0) {
		
			cout<<"**** GROUP OF **** "<<group.size()<<" nodes"<<endl;
			//cherr();
			minimality_bg(group, bcu, uncovered, false, already_found);
			
		}
		else
			cout<<"bad group"<<endl;
		
		
		//*/
	}
	
	paras.hint_flag=false;
	
	return 0;
	
}



int clomnet::load(string filename) {

	
	// this function is to read a file in the tp-format
	
	
	cout<<"getting partition from tp-file: "<<filename<<endl;
	deque<double> bss;
	deque<deque<int> > A;
	
	get_partition_from_file_tp_format(filename, A, bss);
		
	if(translate(A)==-1)
		return -1;
	
	cout<<A.size()<<" groups found"<<endl;
	
	/*char b[1000];
	cast_string_to_char(filename, b);
	ifstream gin(b);
	string s;
	
	deque<double> bss;
	while(getline(gin, s)) {
	
		deque<string>  v;
		separate_strings(s, v);
		if(s.size()>0 && (v[0]=="#module" || v[0]=="#group") )
			bss.push_back(cast_string_to_double(v[5]));
		
		
	}*/
	
	cout<<bss.size()<<" bss found"<<endl;
	
	//time_t t1=time(NULL);
	
	
	for(int ii=0; ii<A.size(); ii++) {

		//cout<<"inserting group number "<<ii<<" size: "<<A[ii].size()<<endl;
		Mcoll.insert(A[ii], bss[ii]);
		
		
	}
	
	
	//checklevel();

	//time_t t2=time(NULL);
	//cout<<"included "<<t2-t1<<" secs"<<endl;

	return 0;
	
}



/*
int clomnet::load(string filename) {

	
	// this function is to read a file in the tp-format
	
	
	cout<<"getting partition from tp-file: "<<filename<<endl;
	deque<deque<int> > A;
	
	get_partition_from_file_list(filename, A);
		
	if(translate(A)==-1)
		return -1;
	
	cout<<A.size()<<" groups found"<<endl;
	
	char b[1000];
	cast_string_to_char(filename, b);
	ifstream gin(b);
	string s;
	
	deque<double> bss;
	while(getline(gin, s)) {
	
		deque<string>  v;
		separate_strings(s, v);
		if(s.size()>0 && (v[0]=="#module" || v[0]=="#group") )
			bss.push_back(cast_string_to_double(v[5]));
		
		
	}
	
	cout<<bss.size()<<" bss found"<<endl;
	
	//time_t t1=time(NULL);
	
	
	for(int ii=0; ii<min(int(A.size()), 500000); ii++) {

		cout<<"inserting group number "<<ii<<" size: "<<A[ii].size()<<endl;
		//Mcoll.insert(A[ii], bss[ii]);
		
		
		
		deque<deque<int> > link_per_node;
		deque<deque<pair<int, double> > > weights_per_node;
		set_subgraph(A[ii], link_per_node, weights_per_node);
		clomnet giovanni(link_per_node, weights_per_node, A[ii], clean_t);
		
		print_id(A[ii], cout);
		print_ri(A[ii]);
		giovanni.draw("loaded");
		cherr();
		
		
		
		
		
	}

	//time_t t2=time(NULL);
	//cout<<"included "<<t2-t1<<" secs"<<endl;

	return 0;
	
}

*/

bool clomnet::minimality_bg(deque<int> & group, double bs, set<int> & uncovered, bool g_print) {

	ne_modules already_found(0);
		
	return minimality_bg(group, bs, uncovered, g_print, already_found);


}






bool clomnet::minimality_bg(deque<int> & group, double bs, set<int> & uncovered, bool g_print, ne_modules & already_found) {

	if(already_found.memberships.size()!=0 && already_found.insert(group)==false) {
		
		
		return false;
	
	
	}
	
		
	//   ******************** pieces **********************
	//   ******************** pieces **********************
	//   ******************** pieces **********************
	//   ******************** pieces **********************
	
	//if(paras.print_cbs)
		//cout<<"minimality..."<<endl;
	
	
	
	deque<int> smaller;
	deque<int> leftp;
	if (Mcoll.pieces(group, leftp, smaller)==true) {
		
			
				
		for(int j=0; j<group.size(); j++)
			uncovered.erase(group[j]);
		
		deque<int> cg2_0;
		double bs__= b_clean_up_search(leftp, cg2_0, paras.bscore_loops);
		deque<int> cg2;
		bs__= b_clean_up_check(cg2_0, cg2, paras.bscore_loops_check);
		
		if(cg2.size()>0) {
			
			group=cg2;
			bs=bs__;
		}
		else		
			return false;
	
	
	}
	
	//   ******************** pieces **********************
	//   ******************** pieces **********************
	//   ******************** pieces **********************
	//   ******************** pieces **********************
	//   ******************** pieces **********************
	
	
	//if(paras.print_cbs)
		//cout<<"pieces done..."<<endl;
	
	
	deque<deque<int> > subM;
	deque<double> bss;
	
	
	bool will_be_inserted=false;
	
	{	//******************  giovanni stuff   ******************
		
		deque<deque<int> > link_per_node;
		deque<deque<pair<int, double> > > weights_per_node;
		set_subgraph(group, link_per_node, weights_per_node);
		clomnet giovanni(link_per_node, weights_per_node, group, clean_t);
		
		/*if(id_of(group[0])==2) {
			print_id(group, cout);
			//cherr();
		}*/
		
		
		//if(paras.print_cbs)
			//cout<<"subgraph:: "<<giovanni.size()<<" nodes and "<<giovanni.stubs()<<" stubs;\t average degree = "<<giovanni.stubs()/giovanni.size()<<endl;
		/*
		prints(group);
		print_id(group, cout);
		giovanni.draw("giov");
		cherr();
		//*/
		//giovanni.gather_sub(subM, bss);

		//if(paras.print_cbs)
			//cout<<"**************   start sub "<<endl;
		
		giovanni.LEVEL= LEVEL+1;
		
		giovanni.gather(paras.sub_gather_runs, true);
		giovanni.post_processing();
		giovanni.Mcoll.set_partition(subM, bss);
		
		//print_id(subM, cout);
		/*
		if(id_of(group[0])==2) {
			print_id(group, cout);
			//cherr();
		}*/
		
		//if(paras.print_cbs)
			//cout<<"end sub **************"<<endl;
		
		
		//*******************   left   ****************
		//*******************   left   ****************
		//*******************   left   ****************
		
		
		if(subM.size()>0) {		// this is done because it can happen that a good module is too big to be not random
	
			deque<int> left;
			giovanni.Mcoll.homeless(left);	
			if(left.size()<giovanni.dim) {
				
				
				
				
				deque<int> cg2;
				double bs__=giovanni.b_clean_up_check(left, cg2, paras.bscore_loops_check);
				if(cg2.size()>0) {
					
					subM.push_back(cg2);
					bss.push_back(bs__);
					
				}
			}
		}

		
		//*******************   left   ****************
		//*******************   left   ****************
		//*******************   left   ****************
		
		
		
		
		
		
		
		for(int i=0; i<subM.size(); i++) {
			
			giovanni.deque_id(subM[i]);
			
			deque<int> grbe;
			bss[i]=b_clean_up_check(subM[i], grbe,  paras.bscore_loops_check);
			subM[i]=grbe;
			
			
			// so now you know these modules are cosy and cleaned
		
		}
		
		
		
		int note=0;
		for(int i=0; i<subM.size(); i++) if(subM[i].size()==0)
			note++;
		
		if(subM.size()<=note +1) {
			
			//subM.clear();
			will_be_inserted=true;
		
		}
		
		

				

		
		
	}   //******************  giovanni stuff   ******************
	
	
	/*
	//cout<<"-------------------------"<<endl;
	cout<<subM.size()<<" subgraphs found: "<<endl;
	for(int i=0; i<subM.size(); i++) {
		cout<<subM[i].size()<<endl;
		//print_true(subM[i]);
	}//*/

	
	set<int> a;
	for(int i=0; i<subM.size(); i++) for(int j=0; j<subM[i].size(); j++)
		a.insert(subM[i][j]);
	
		
	//cout<<"sub coverage: "<<double(a.size())/group.size()<<endl;
	
	
	

	
	bool something = false;
	if(a.size()>paras.coverage_percentage_minimality*group.size()) {
		
		
		
		//if(paras.print_cbs)
			//cout<<"GROUP ERASED"<<endl;
		
		
		for(int i=0; i<subM.size(); i++) if(subM[i].size()>0) {
			
			//cout<<"trying some sugges"<<endl;
			
			if(minimality_bg(subM[i], bss[i], uncovered, g_print, already_found))
				something=true;
		}
	
		return something;

		
	}
	else
		will_be_inserted=true;
		
	
	if(will_be_inserted) {
	

		//if(paras.print_cbs)
			//cout<<"GROUP ACCEPTED"<<endl;

		
		for(int j=0; j<group.size(); j++)
			uncovered.erase(group[j]);
		
		
		
		int new_name;
		if(group.size()> 1 && Mcoll.insert(group, bs, new_name)) {
			
			
			
			if(g_print && paras.every_group) {
				
				
				
				fileout<<"#module "<<new_name<<" size: "<<Mcoll.modules[new_name].size()<<" bs: "<<Mcoll.module_bs[new_name]<<"; network of "<<dim<<" nodes"<<endl;
				print_id(group, fileout);
			
			}
			
			
			something=true;
			//Mcoll.insert(new_name, subM);
			
		}
		
		
		return something;
		
		
	
	}
	
	
	

}










void clomnet::print_modules(bool not_homeless, string tp) {
	
	
	Mcoll.arrange_module_id();		// this function is to have all the modules in consecutive order starting from 0
	fill_homel();
	cout<<"******** module_collection ******** "<<effective_groups()<<" modules. writing... "<<endl;

			
	deque<int> netlabs;
	for(int i=0; i<dim; i++)
		netlabs.push_back(id_of(i));
	
	
	
	char b[1000];
	cast_string_to_char(tp, b);
	ofstream out1(b);
	
	Mcoll.print(out1, netlabs, not_homeless);
	cout<<"DONE   ****************************"<<endl;
	



}




void clomnet::print_modules(bool not_homeless, ostream & outt) {
	
	
	Mcoll.arrange_module_id();		// this function is to have all the modules in consecutive order starting from 0
	fill_homel();
	cout<<"******** module_collection ******** "<<effective_groups()<<" modules. writing... "<<endl;

			
	deque<int> netlabs;
	for(int i=0; i<dim; i++)
		netlabs.push_back(id_of(i));
	
	
	
	
	Mcoll.print(outt, netlabs, not_homeless);
	cout<<"DONE   ****************************"<<endl;
	



}




bool clomnet::fusion_module_its_subs(deque<int> & A, deque<deque<int> > & its_submodules) {

	
	// this is done for a group that you know it is a good group
	
	
	// return true if A won against its submodules
	// ******************************************
	
	

	
	//*
	
	//cout<<"sub modules to check"<<its_submodules.size()<<endl;
	
	set<int> all_the_groups;
	for(int i=0; i<its_submodules.size(); i++) {
	
		//print_id(its_submodules[i], cout);
		
		for(int j=0; j<its_submodules[i].size(); j++)
			all_the_groups.insert(its_submodules[i][j]);
	
	}
	for(int i=0; i<A.size(); i++)
		all_the_groups.insert(A[i]);
	
	deque<int> group;
	set_to_deque(all_the_groups, group);
	
	//*/
	
	
	
		
	
	{	//******************  giovanni stuff   ******************
	
		
		
		//cout<<"ok- here we are"<<endl;
		//print_id(all_the_groups, cout);
		
		deque<deque<int> > link_per_node;
		deque<deque<pair<int, double> > > weights_per_node;
		set_subgraph(group, link_per_node, weights_per_node);
		clomnet giovanni(link_per_node, weights_per_node, group, clean_t);
		
				
		giovanni.translate(its_submodules);
		
		//------------------------------------ cleaning up submodules --------------------------
		
		//cout<<"-> modules: "<<its_submodules.size()<<endl;
		//printm(its_submodules);
		
		for(int i=0; i<its_submodules.size(); i++)
			giovanni.Mcoll.insert(its_submodules[i], 1e-3);
		
		
		
		giovanni.Mcoll.compute_inclusions(false);
		giovanni.Mcoll.set_partition(its_submodules);
		giovanni.Mcoll.clear();
		
		
		
		//cout<<"-> modules: "<<its_submodules.size()<<endl;
		//printm(its_submodules);

		
		//------------------------------------ cleaning up submodules --------------------------
		
		for(int i=0; i<its_submodules.size(); i++) {
			
			deque<int> grbe;
			giovanni.b_clean_up_check(its_submodules[i], grbe,  paras.bscore_loops_check);
			its_submodules[i]=grbe;
		
		}
		
		
		
		int note=0;
		for(int i=0; i<its_submodules.size(); i++) if(its_submodules[i].size()==0)
			note++;
		
		if(its_submodules.size()<=note +1)
			return true;		
		
		
		
			
		set<int> a;
		for(int i=0; i<its_submodules.size(); i++) for(int j=0; j<its_submodules[i].size(); j++)
			a.insert(its_submodules[i][j]);

		
		
		if(a.size()>paras.coverage_percentage_minimality*group.size())
			return false;
		else
			return true;

		
		
	}   //******************  giovanni stuff   ******************
	
	
	

}




int clomnet::check_ex_pieces(module_collection & mall) {

	
	// this function si to check if ex_pieces can win again the small groups it contains
	//cout<<"checking ex-pieces"<<endl;
	deque<int> this_will_be_inserted;
	
	
	for(int i=0; i<mall.ex_pieces.size(); i++) {
		
		deque<int> not_covered;
		deque<int> smaller;
		deque<deque<int> > its_submodules;

		mall.pieces(mall.ex_pieces[i], not_covered, smaller);
		for(int j=0; j<smaller.size(); j++)
			its_submodules.push_back(mall.modules[smaller[j]]);

		
		//cout<<"submodules"<<endl;
		//printm(its_submodules);
		
	
		if(fusion_module_its_subs(mall.ex_pieces[i], its_submodules)) {
			
			//cout<<"fusion_ won"<<endl;
			this_will_be_inserted.push_back(i);
			
		} //else
			//cout<<"fusion_ lost"<<endl;
	
	}
	
			
	
	//prints(this_will_be_inserted);
	
	//cout<<"size1: "<<mall.size()<<endl;
	//cout<<"inserting again ex_almost_equal and this_will_be_inserted "<<endl;

	
	for(int i=0; i<this_will_be_inserted.size(); i++)
		mall.insert(mall.ex_pieces[this_will_be_inserted[i]], 1e-3);
	
	
	for(int i=0; i<mall.ex_almost_equal.size(); i++)
		mall.insert(mall.ex_almost_equal[i], 1e-3);
	
	//cout<<"size2: "<<mall.size()<<endl;
	
	mall.erase_included(false);
	
		
	
	
	
	
	//cherr();

	return 0;
	
}




int clomnet::checklevel(module_collection & mall) {

	
	
	
	if(effective_groups()==0)
		return 0;
	
	//cout<<"mall compute_inclusions... "<<dim<<endl;
	
	
	
	
		
	mall.ex_pieces.clear();
	mall.ex_almost_equal.clear();
	
	//cout<<"mall compute_inclusions... "<<dim<<endl;
	mall.compute_inclusions(true);
	check_ex_pieces(mall);
	
		
	
	//cout<<"mall check_intersection... "<<dim<<endl;
	check_intersection();
	mall.compute_inclusions(false);
		
	
	
	
	
	
	
	
	if(matteo_pairs()>0) {
		
		mall.compute_inclusions(false);
		
		if(matteo_pairs()>0)
			mall.compute_inclusions(false);	
	
	}
	
	return 0;

}

//*


int clomnet::check_intersection() {



	deque<int> to_check;
	for(map<int, double>::iterator itM = Mcoll.module_bs.begin(); itM!=Mcoll.module_bs.end(); itM++)
		to_check.push_back(itM->first);
	
	
	

	return check_intersection(to_check);
	
}




int clomnet::check_intersection(deque<int> & to_check) {

	
	
	
	// this function is done for a very peculiar problem of the method.
	//cout<<"check intersection"<<endl;
	
	
	set<pair<int, int> > pairs_to_check;
	
	
	for(deque<int>::iterator itM = to_check.begin(); itM!=to_check.end(); itM++) if( Mcoll.module_bs.find(*itM) != Mcoll.module_bs.end() )  {
		
		deque<int> & c= Mcoll.modules[*itM];
		
		/*
		cout<<"ccccccccc"<<endl;
		prints(c);
		
		//*/
		
		
		
		map<int, int> com_ol;		// it maps the index of the modules into the overlap (overlap=numeber of overlapping nodes)
	
		for(int i=0; i<c.size(); i++)
			for(set<int>:: iterator itj=Mcoll.memberships[c[i]].begin(); itj!=Mcoll.memberships[c[i]].end(); itj++)
				int_histogram(*itj, com_ol);
	
	
	
	
		for(map<int, int> :: iterator cit= com_ol.begin(); cit!=com_ol.end(); cit++) if( cit->first != *itM ) {
				
				
			/*cout<<"overlaps "<<endl;
			prints(c);
			prints(Mcoll.modules[cit->first]);
			cout<<double(cit->second) / min(Mcoll.modules[cit->first].size(), c.size())<<" <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;*/
			
			if(  double(cit->second) / min(Mcoll.modules[cit->first].size(), c.size()) > paras.check_inter_p  ) {		// they have a few nodes in common
			
				// this means that these guys have quite some overlap in common ...
				pairs_to_check.insert( make_pair( min(*itM, cit->first) , max(*itM, cit->first) ) );
				//pairs_to_check.insert( make_pair(itm->first, cit->first) );
			
			}
		
		}
		
	}
	
	
	
	/*if(fusion_intersection(c, Mcoll.modules[cit->first], itM->first, cit->first))
					return check_intersection();*/
					
					
	return fusion_intersection(pairs_to_check);
	
	
}





int clomnet::fusion_intersection(set<pair<int, int> > & pairs_to_check) {

	/*
	if(paras.print_cbs) {	
		cout<<"Intersection check"<<endl;
	}
	//*/
	
	
	
	
	
	deque<int> new_insertions;
	
	
	for(set<pair<int, int> >::iterator ith= pairs_to_check.begin(); ith!=pairs_to_check.end(); ith++ ) if(ith->first < ith->second ) 
	if(Mcoll.module_bs.find(ith->first)!=Mcoll.module_bs.end())	if(Mcoll.module_bs.find(ith->second)!=Mcoll.module_bs.end()) {
	
	
	
		//		first, you need to check if both the modules in the pair are still in mcoll
		
		
		deque<int> & a1= Mcoll.modules[ith->first];
		deque<int> & a2= Mcoll.modules[ith->second];
		int min_s=min(a1.size(), a2.size());
		
		
		
		deque<int> group_intsec;
		set_intersection(a1.begin(), a1.end(), a2.begin(), a2.end(), back_inserter(group_intsec));



		//		if they are, you need to check if they are not almost equal.


		if(double(group_intsec.size()) / min_s >= paras.coverage_inclusion_module_collection) {
		
			int em=ith->first;
			if(a1.size() < a2.size())
				em=ith->second;
			else if(a1.size()==a2.size() && Mcoll.module_bs[ith->first]>Mcoll.module_bs[ith->second])
				em=ith->second;
		
			Mcoll.erase(em);
		
		
		}		
		else
			fusion_intersection(ith->first, ith->second, new_insertions);
		
		
	
	}
	

	if(new_insertions.size()>0)
		return check_intersection(new_insertions);
	
	return 0;




}




inline bool clomnet::fusion_intersection(int ai1, int ai2, deque<int> & new_insertions) {

	
	
	
	deque<int> & a1= Mcoll.modules[ai1];
	deque<int> & a2= Mcoll.modules[ai2];
	
	deque<int> group;
	set_union(a1.begin(), a1.end(), a2.begin(), a2.end(), back_inserter(group));
	



	if(group.size()!=dim) {	//******************  giovanni stuff   ******************
		
		
		
		
		deque<deque<int> > link_per_node;
		deque<deque<pair<int, double> > > weights_per_node;
		set_subgraph(group, link_per_node, weights_per_node);
		clomnet giovanni(link_per_node, weights_per_node, group, clean_t);
		
		deque<deque<int> > A;
		A.push_back(a1);
		A.push_back(a2);
		
		giovanni.translate(A);
		
		
		/*
		cout<<"giovanni dim "<<giovanni.dim<<endl;
		giovanni.print_id(A[0], cout);
		giovanni.print_id(A[1], cout);
		//*/
		
		deque<int> grc1;
		double bs=giovanni.b_clean_up_check(A[0], grc1, paras.bscore_loops_check);
		deque<int> grc2;
		bs=giovanni.b_clean_up_check(A[1], grc2, paras.bscore_loops_check);
		
				
		
		if(grc1.size()==0 || grc2.size()==0) {		// if one or the other is random you can take the fusion_ (if it's good)
			
			
			//cout<<"died"<<endl;
			
			deque<int> group_2;
			bs=b_clean_up_check(group, group_2, paras.bscore_loops_check);
			
			
			if(group_2.size()>paras.coverage_percentage_fusion * group.size()) {
				
				
				//cout<<group_2.size()<<" "<<group.size()<<endl;
				
				Mcoll.erase(ai1);
				Mcoll.erase(ai2);
				
				
				set<int> uncovered;
				minimality_bg(group_2, bs, uncovered, false);

				
				return true;
			}
			else
				return false;
		
		
		}
		
		//return false;
		
		//cout<<"alive"<<endl;
		giovanni.deque_id(grc1);
		giovanni.deque_id(grc2);
		
		
		//prints(grc1);
		//prints(grc2);
		
		
		deque<int> cg1;
		double bs__1=b_clean_up_search(grc1, cg1, paras.bscore_loops_check);		
			
		deque<int> cg2;
		double bs__2=b_clean_up_search(grc2, cg2, paras.bscore_loops_check);

		
			
		
		
		deque<int> inters;
		set_intersection(cg1.begin(), cg1.end(), cg2.begin(), cg2.end(), back_inserter(inters));
		
		deque<int> unions;
		set_union(cg1.begin(), cg1.end(), cg2.begin(), cg2.end(), back_inserter(unions));


		//cout<<"new inter "<<double(inters.size())/min(cg1.size(), cg2.size())<<endl;
		
		
		if(double(inters.size())/min(cg1.size(), cg2.size())<paras.check_inter_p && cg1.size()>0 && cg2.size()>0 && (unions.size()>paras.coverage_percentage_fusion * group.size())) {
			
			
				
			Mcoll.erase(ai1);
			Mcoll.erase(ai2);
			
			int newi;
			
			
			Mcoll.insert(cg1, bs__1, newi);
			new_insertions.push_back(newi);
			
			Mcoll.insert(cg2, bs__2, newi);
			new_insertions.push_back(newi);
			
			
			
			/*print_id(cg1, cout);
			print_id(cg2, cout);*/
			
			return true;
			
		}




	}

	return false;
	
	
	
}






int clomnet::matteo_pairs() {
	
	
	map<int, map<int, pair<int, double> > > neigh_weight_s;		// this maps the module id into the neighbor module ids and weights
	set_upper_network(neigh_weight_s, Mcoll);
	multimap<int, int> A;
	
	
	if(neigh_weight_s.size()>0) {	// upper level stuff
		
		clomnet matteo(paras.threshold, neigh_weight_s);
		matteo.check_pairs(A, clean_t);
		
	}
	
	
	
	set<pair<int, int> > p1;
	
	for(multimap<int, int>:: iterator itmm=A.begin(); itmm!=A.end(); itmm++) {
		p1.insert( make_pair( min(itmm->first, itmm->second) , max(itmm->first, itmm->second) ) );
	}
		
	
	int something_done=upper_fusion(p1);
	
	
	return something_done;


}




int map_game(map<int, int> & new_names, int to_search) {

	
	int to_search_first= to_search;
	map<int, int> :: iterator itf= new_names.find(to_search);
	
	bool entered=false;
	
	while(itf!=new_names.end()) {
		
		//cout<<"replacing: "<<to_search<<" with "<<itf->second<<endl;
		
		to_search=itf->second;
		itf= new_names.find(to_search);
		entered=true;
	
	}
	
	if(entered)
		new_names[to_search_first]=to_search;
	
	
	return to_search;
	

}







int clomnet::upper_fusion(set<pair<int, int> > & p1) {
	
	
	int something_done=-1;
	
		
	
	
	
	map<int, int> new_names;
	
	

		
	for(set<pair<int, int> >::iterator ith= p1.begin(); ith!=p1.end(); ith++ ) if(ith->first < ith->second ) {
	
		
		int a1_name=ith->first;
		
		
		
		
		if(Mcoll.module_bs.find(a1_name)==Mcoll.module_bs.end()) {	

			a1_name=map_game(new_names, a1_name);
		
		}
		
		int a2_name=ith->second;
		
		
		if(Mcoll.module_bs.find(a2_name)==Mcoll.module_bs.end()) {
			
			a2_name=map_game(new_names, a2_name);

			
		
		}
		
	
		if(a1_name!=-1 && a2_name!=-1) {
		
			
			deque<int> & a1= Mcoll.modules[a1_name];
			deque<int> & a2= Mcoll.modules[a2_name];
			{

			
				deque<int> a12;
				
				
				double fbs=fusion(a1, a2, a12);
				
				
				
				if(fbs<1) {
					
					something_done=1;
					
									
					Mcoll.erase(a1_name);
					Mcoll.erase(a2_name);
					
					int newi;
					
					Mcoll.insert(a12, fbs, newi);
					
					//cout<<"fused: "<<a1_name<<" "<<a2_name<<" ->  "<<newi<<endl;
					
					new_names[a1_name]=newi;
					new_names[a2_name]=newi;
					
					
			
					
				}
				
				
								
				
			}
		}
	}
	
	

	return something_done;



}


double clomnet::fusion(deque<int> & a1, deque<int> & a2, deque<int> & group) {

	
		
	//cout<<"Fusion_ "<<a1.size()<<" "<<a2.size()<<endl;


	/*
	cout<<"-************"<<endl;	
	print_id(a1, cout);
	print_id(a2, cout);
	cout<<"-************"<<endl;
	cherr();
	//*/
	
	

	
	
	group.clear();
	set_union(a1.begin(), a1.end(), a2.begin(), a2.end(), back_inserter(group));
	
	
	
	{	//******************  giovanni stuff   ******************
		
		deque<deque<int> > link_per_node;
		deque<deque<pair<int, double> > > weights_per_node;
		set_subgraph(group, link_per_node, weights_per_node);
		clomnet giovanni(link_per_node, weights_per_node, group, clean_t);
		
		deque<deque<int> > A;
		A.push_back(a1);
		A.push_back(a2);
		
		giovanni.translate(A);
		
		
		/*
		cout<<"giovanni dim "<<giovanni.dim<<endl;
		print_id(A[0], cout);
		print_id(A[1], cout);
		//*/
		
		deque<int> grc1;
		double bs=giovanni.b_clean_up_check(A[0], grc1, paras.bscore_loops_check);
		deque<int> grc2;
		bs=giovanni.b_clean_up_check(A[1], grc2, paras.bscore_loops_check);
		
		if(grc1.size()==0 || grc2.size()==0) {		// if one or the other is random you can take the fusion_ (if it's good)
			
			//cout<<"both random"<<endl;
			grc1.clear();
			bs=b_clean_up_check(group, grc1, paras.bscore_loops_check);
			if(grc1.size()>paras.coverage_percentage_fusion * group.size()) {
				group=grc1;
				return bs;
			}
			else {
				return 1.1;
			}
		}
				
	}


	return 1.1;
	
}











/*
int clomnet::one_core(deque<int> & thisrun, deque<int> & si, deque<int> & nextrun, set<int> & ultimate) {
	
	
	
	for(int i=0; i<thisrun.size(); i++) if(ultimate.find(thisrun[i])==ultimate.end()) {
			
		ultimate.insert(thisrun[i]);
		
		for(int j=0; j<vertices[thisrun[i]]->links->size(); j++) {
			
			
			int neigh=vertices[thisrun[i]]->links->l[j];
			si[neigh]--;
			
			if(si[vertices[thisrun[i]]->links->l[j]]==1)
				nextrun.push_back(neigh);
				
		}
	}
	
	
	
	
	return 0;

}*/






/*
int clomnet::set_one_core(deque<deque<int> > & M, deque<int> & homel) {
	
	
	M.clear();
	/*map<int, int> A;
	get_id_label(A);* /
	
	
		
		
		
	deque<int> si(dim);
	for(int i=0; i<dim; i++)
		si[i]=-1;
		
	
	deque<int> thisrun;
	for(int i=0; i<homel.size(); i++) {
	
		si[homel[i]]=cast_int(vertices[homel[i]]->links->size());
		
		if(si[homel[i]]==1)
			thisrun.push_back(homel[i]);
	}
	
	
	//si[A[9]]=-1;
	
	
	set<int> ultimate;
	int sisi=0;
	
	while(thisrun.size()>0 ) {
		
		
		
		/*cout<<"thisrun: "<<++sisi<<endl;
		print_id(thisrun, cout);* /
		
		deque<int> nextrun;
		one_core(thisrun, si, nextrun, ultimate);
		
		thisrun=nextrun;
		
	}

	
	deque<int> prop;
	set_to_deque(ultimate, prop);
	
	deque<deque<int> > link_per_node;
	deque<deque<pair<int, double> > > weights_per_node;
	set_subgraph(prop, link_per_node, weights_per_node);
	static_network giovanni;
	giovanni.set_graph(link_per_node, weights_per_node, prop);
	
	giovanni.set_connected_components(M);
		
	for(int j=0; j<M.size(); j++)
		giovanni.deque_id(M[j]);
	

	return 0;

}*/


/*
int clomnet::check_props() {

	return check_props(false);

}*/



/*
int clomnet::check_props(bool anyway) {


	deque<int> homel;
	Mcoll.homeless(homel);

	
	
	deque<deque<int> >  M;
	set_one_core(M, homel);
	
	map<int, deque<int> > to_check;			// module - homeless nodes added to that
	

	
	set<int> homes;
	deque_to_set(homel, homes);
	
	bool something=false;

	
	
	for(int i=0; i<M.size(); i++) {
			
			
		/*cout<<"props ***************** "<<i<<endl;
		print_id(M[i], cout);* /
		
		set<int> connected_nodes;		
		
		//**************************************************************************
		for(int j=0; j<M[i].size(); j++) if(vertices[M[i][j]]->links->size()>1) {
			
			for(int k=0; k<vertices[M[i][j]]->links->size(); k++) {

				if(homes.find(vertices[M[i][j]]->links->l[k])==homes.end())
					connected_nodes.insert(vertices[M[i][j]]->links->l[k]);
				
			}
		}
		
		//**************************************************************************
		
		
		
		for(set<int>:: iterator its=connected_nodes.begin(); its!=connected_nodes.end(); its++) {
			
			
			
			//cout<<"connected_nodes "<<*its<<endl;
			
			
			for(set<int>:: iterator itk=Mcoll.memberships[*its].begin(); itk!=Mcoll.memberships[*its].end(); itk++) {
				
				int belongs_to= *itk;
				
				if(to_check.find(belongs_to)==to_check.end()) {
					deque<int> void_d;
					to_check[belongs_to]=void_d;
				}
				
				deque<int> & thid = to_check[belongs_to];
				for(int h=0; h<M[i].size(); h++)
					thid.push_back(M[i][h]);

				
			
			}
		
		}
		
		//**************************************************************************
		
		
		
		
		
		if(connected_nodes.size()==0) {
		
		
			if(anyway) {

				something=true;
				Mcoll.insert(M[i], ran4()+paras.threshold);
			}
			else {
				
				
				
				//cout<<"no friends!"<<endl;
				
				
				deque<int> grbe;
				double bs=b_clean_up_check(M[i], grbe,  paras.bscore_loops_check);
				
				//cout<<"union_deque "<<itm->first<<" size: "<<grbe.size()<<endl;
				
				
				if(grbe.size()>1) {
					
					
					//cout<<"inserted with no friends "<<grbe.size()<<endl;
					something=true;
					Mcoll.insert(grbe, bs);
					
				}
				
			}
		
		
		}
		//**************************************************************************
	
	
	}
	
	
	// **** try the groups with the homeless //******************
	
	for(map<int, deque<int> > :: iterator itm= to_check.begin();  itm!=to_check.end(); itm++) {
	
		
		
		deque<int> union_deque=Mcoll.modules[itm->first];
		for(int i=0; i<itm->second.size(); i++)
			union_deque.push_back(itm->second[i]);
		
		/*
		cout<<"union_deque "<<itm->first<<" size: "<<union_deque.size()<<endl;
		print_id(union_deque, cout);//* /
		
		
		if(anyway) {

			something=true;
			Mcoll.insert(union_deque, ran4()+paras.threshold);
		}
		else {
		
			deque<int> grbe;
			double bs=b_clean_up_check(union_deque, grbe,  paras.bscore_loops_check);
			
			//cout<<"union_deque "<<itm->first<<" size: "<<grbe.size()<<endl;
			
			
			if(grbe.size()>1) {
				
				something=true;
				Mcoll.insert(grbe, bs);
				
			}
			
		}
	
	}
	
	
	if(something)
		Mcoll.compute_inclusions();
	


	return 0;


}
*/

int clomnet::fill_homel() {

	Mcoll.fill_gaps();

}



void clomnet::module_id_first(map<int, set<int> > & A) {
	
	
	// A maps the module index into the node labels
	
	
	deque<int> netlabs;
	for(int i=0; i<dim; i++)
		netlabs.push_back(id_of(i));
	
	
	Mcoll.id_nodes(A, netlabs);


}


void put_into_set(set<int> & a, set<int> & b) {

	
	for(set<int>:: iterator its = a.begin(); its!=a.end(); its++ )
		b.insert(*its);
}


void clomnet::module_id_later(map<int, set<int> > & A, map<int, set<int> > & B) {
	
	
	module_id_first(B);
	
	for(map<int, set<int> >::iterator itm = B.begin(); itm!=B.end(); itm++) {
		
		
		set<int> union_set;
		set<int> & a= itm->second;
		for(set<int>:: iterator its = a.begin(); its!=a.end(); its++ )
			put_into_set(A[*its], union_set);
		
		B[itm->first]=union_set;
	
	}


}






int clomnet::effective_groups() {

	
		
	int nmod=0;
	for(map<int, double >::iterator itm = Mcoll.module_bs.begin(); itm!=Mcoll.module_bs.end(); itm++) if(Mcoll.modules[itm->first].size() > 1)		
		nmod++;

	return nmod;



}




bool clomnet::fast_minimality_bg(deque<int> & group, double bs, ne_modules & already_found) {

	
	
	deque<deque<int> > subM;
	deque<double> bss;
	
	
	bool will_be_inserted=false;
	
	{	//******************  giovanni stuff   ******************
		
		deque<deque<int> > link_per_node;
		deque<deque<pair<int, double> > > weights_per_node;
		set_subgraph(group, link_per_node, weights_per_node);
		clomnet giovanni(link_per_node, weights_per_node, group, clean_t);
		
				
		giovanni.LEVEL= LEVEL+1;
		
		//cout<<"giovanni calls louvain "<<giovanni.dim<<endl;
		giovanni.fast_gather_and_process(true);
		giovanni.Mcoll.set_partition(subM, bss);
		
			
		
		if(subM.size()>0) {		// this is done because it can happen that a good module is too big to be not random
	
			deque<int> left;
			giovanni.Mcoll.homeless(left);	
			if(left.size()<giovanni.dim) {
				
				deque<int> cg2;
				double bs__=giovanni.b_clean_up_check(left, cg2, paras.bscore_loops_check);
				if(cg2.size()>0) {
					
					subM.push_back(cg2);
					bss.push_back(bs__);
					
				}
			}
		}

		
		//*******************   left   ****************
		//*******************   left   ****************
		//*******************   left   ****************
		
		
		
		
		
		
		
		for(int i=0; i<subM.size(); i++) {
			
			giovanni.deque_id(subM[i]);
			
			deque<int> grbe;
			bss[i]=b_clean_up_check(subM[i], grbe,  paras.bscore_loops_check);
			subM[i]=grbe;
			
			
			// so now you know these modules are cosy and cleaned
		
		}
		
		
		
		int note=0;
		for(int i=0; i<subM.size(); i++) if(subM[i].size()==0)
			note++;
		
		if(subM.size()<=note +1) {
			
			//subM.clear();
			will_be_inserted=true;
		
		}
		
		

				

		
		
	}   //******************  giovanni stuff   ******************
	
	

	
	set<int> a;
	for(int i=0; i<subM.size(); i++) for(int j=0; j<subM[i].size(); j++)
		a.insert(subM[i][j]);
	
		
	//cout<<"sub coverage: "<<double(a.size())/group.size()<<endl;
	
	
	

	
	bool something = false;
	if(a.size()>paras.coverage_percentage_minimality*group.size()) {
		
		
		
			
		for(int i=0; i<subM.size(); i++) if(subM[i].size()>0) {
			
			//cout<<"trying some sugges"<<endl;
			
			if(subM[i].size()<100) {
				set<int> uncovered;
				if(minimality_bg(subM[i], bs, uncovered, false, already_found))
					something=true;
			} else {
				if(fast_minimality_bg(subM[i], bs, already_found))
					something=true;
			}
			
		}
	
		return something;

		
	}
	else
		will_be_inserted=true;
		
	
	if(will_be_inserted) {
	
	
		int new_name;
		if(group.size()> 1 && Mcoll.insert(group, bs, new_name))
			something=true;
		
		
		return something;		
	
	}
	
	
}




void clomnet::print_statistics(ostream & outt) {
	
	
	int nmod=0;
	for(map<int, double >::iterator itm = Mcoll.module_bs.begin(); itm!=Mcoll.module_bs.end(); itm++) if(Mcoll.modules[itm->first].size() > 1)		
		nmod++;
	
	
	int cov=0;
	for(int i=0; i<Mcoll.memberships.size(); i++)
		cov+=Mcoll.memberships[i].size();
	
	
	//deque<int> degree_of_homeless;
	deque<int> homel;
	Mcoll.homeless(homel);
	/*for(int i=0; i<homel.size(); i++)
		degree_of_homeless.push_back(vertices[homel[i]]->stub_number);*/
		
	
	outt<<"number of modules: "<<nmod<<endl;
	outt<<"number of covered nodes: "<<dim - homel.size()<<" fraction of homeless nodes: "<<double(homel.size())/dim<<endl;
	outt<<"average number of memberships of covered nodes: "<<double(cov)/(dim - homel.size())<<endl;
	//outt<<"average degree of homeless nodes: "<<average_func(degree_of_homeless)<<" dev: "<<sqrt(variance_func(degree_of_homeless))<<endl;
	outt<<"average community size: "<<double(cov)/nmod<<endl;

}




