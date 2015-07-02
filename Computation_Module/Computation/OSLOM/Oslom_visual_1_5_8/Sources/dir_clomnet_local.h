


// neigh update could be optimized better




inline double log_zero(double a) {

	if(a<=0)
		return -1e20;
	else
		return log(a);

}




class clomnet_local : public static_network {
	
	
	public:
		
		

		clomnet_local(deque<deque<int> > & b, deque<deque<pair<int, double> > > & c, deque<int> & d): static_network() { set_graph(b,c,d); set_maxbord();  };
		clomnet_local(string a): static_network() { set_graph(a); set_maxbord();  };
		clomnet_local(map<int, map<int, pair<int, double> > > & A) : static_network() { set_graph(A);  set_maxbord();  };
		~clomnet_local(){};
		
		
		
		double setic_simann(int , deque<int> & , double, const deque<double> & , int , int, bool);
		//double CBS(deque<int> & );
		double refinement_single(deque<int> & );
		double refinement_single_only_c(deque<int> & gr_cleaned, set<int> & ONLY_AMONG_THIS);
		double b_clean_up_search_iterative(deque<int> _c_, deque<int> & gr_cleaned, int l);
		double b_clean_up_check(const deque<int> & _c_, deque<int> & gr_cleaned, int l);
		double b_clean_up_search(const deque<int> & _c_, deque<int> & gr_cleaned, int l);
		//int check_all();

		/*
		void print_ri(deque<int> & aa) {
		
			for(int i=0; i<aa.size(); i++)
				cout<<RI[id_of(aa[i])]<<" ";
			cout<<endl;
		
		
		};//*/
		
		/*
		int this_color(int a) {
			
			int nn=0;
			for(int i=0; i<dim; i++)
				if(RI[id_of(i)]==a)
					nn++;
			
			
			return nn;
			
		
		
		};*/
		
		int set_upper_network(map<int, map<int, pair<int, double> > > & neigh_weight_f, module_collection & module_coll);
		int try_to_assign_homeless(module_collection & module_coll, bool);
		int check_pairs(multimap<int, int> & A, double clean_t);
	
	private:
		
		//double b_clean_up_once_yes_no(const deque<int> & _c_, deque<int> & gr_cleaned, int loops);
		void include_prate(map<int, int> part_rate_firsts, deque<int> & gr_cleaned, int times, double );
		
		double CBS(int &);
		double ord_border_score(deque<pair<double, int> > & xs, int & until);
		
		int setic_clean_up(int , int, const double &, bool);
		
	
		int set_cgroup_and_neighs(const deque<int> &);
		
		bool insert_the_best(double & cmax, double threshold, int & benode);
		bool insert_the_best(double);
		bool erase_the_worst(double threshold, int & wnode, set<int> & erased);
		bool erase_the_worst(double);
		bool erase_the_worst(double threshold, int & wnode);

		void insertion(int changendi, deque<int> & inserted, bool);
		int set_maxbord();
		
		
		
		void set_bg();
		double CBS();
		

		
		deque<int> core_p;
		double lcb_temp;
		double lowest_cs;
		int gstopper;
		double temp;
		double  max_r_bord;
		int maxb_nodes;
		double b_clean_up_once(const deque<int> & _c_, deque<int> & gr_cleaned, int loops, bool only_c);
		
		void insert_cgroup(int );
		void erase_cgroup(int );
		void erase_cgroup(int wnode, set<int> & erased);

		int try_to_assign_homeless_help(module_collection & module_coll, map<int, deque<int> > & to_check);

		
		
		// ************* things to update *************************
																//*
		weighted_tabdeg cgroup;									//*
		weighted_tabdeg neighs;									//*
																//*
		int kin_cgroup;											//*
		int ktot_cgroup_in;										//*
		int ktot_cgroup_out;									//*
																//*
																//*
		// ************* things to update *************************

		
		

};





int clomnet_local::set_maxbord() {
	
	
	// here there should be the function with arrange the weights (basically it has to divide each weight by the average expected...)
	
	
	
	
	
	
	
	//cout<<"new network of "<<dim<<endl;
	
	
	//	this function is used because if there are condensed nodes the program gets stuck
	
	
	
	double max_degin_out=0;
	
	if(dim<2) {
		
		max_r_bord=paras.maxbg_ordinary;
		return -1;
	}
		
	
	for(int i=0; i<dim; i++) {
		
		if(vertices[i]->instub_number * vertices[i]->outstub_number > max_degin_out)
			max_degin_out=vertices[i]->instub_number * vertices[i]->outstub_number;
		
		
	}
	
	
	//cout<<"*************** "<<max_deg<<" "<<2*tstrength<<" "<<max_deg*max_deg / (2*tstrength)<<endl;
	
	if(max_degin_out / double(oneM) > paras.max_self_loops) {		/* network has at least one condensed nodes ...*/
		
		max_r_bord=min(paras.maxbg_condensed, paras.maxbg_condensed * 100. /double(dim));				/* this is the maximum r you can have in the border  */
		maxb_nodes=paras.maxborder_nodes/5;
		
		//if(paras.print_cbs)
			//cout<<"condensed nodes: max_deg: "<<max_deg<<" expected self loops "<<max_deg*max_deg / (twoM)<<" max_r_bord "<<max_r_bord<<endl;

		
	} else if(max_degin_out / double(oneM) > paras.max_self_loops/2.) {
		
		max_r_bord=paras.maxbg_ordinary/2;
		maxb_nodes=paras.maxborder_nodes/2;
		//if(paras.print_cbs)
			//cout<<dim<<"-> almost condensed nodes: max_deg: "<<max_deg<<" expected self loops "<<max_deg*max_deg / double(twoM) <<" max_r_bord "<<max_r_bord<<endl;
	
	
	} else {
	
		max_r_bord=paras.maxbg_ordinary;
		maxb_nodes=paras.maxborder_nodes;
	}
	
	return 0;


}




/*
double clomnet_local::CBS(deque<int> & A) {


	set_cgroup_and_neighs(A);
	return CBS();

}
*/


/*

double clomnet_local::CBS() {


	int wnode, kw, degree_ref;
	return CBS(wnode, kw, degree_ref);



}*/








double clomnet_local::CBS(int & wnode) {

	
	// here kin_cgroup is the internal degree of C not C_minus_w! because you don't know the worst internal node yet!
	// wnode is the worst node, kw the worst internal degree, degree_ref its degree
	
	

		
	//------------------------------ compute the worst internal node
	
	int Nstar= dim - cgroup.size();
	int nn= neighs.size();
	int kout_g_in = ktot_cgroup_in - kin_cgroup;
	int kout_g_out = ktot_cgroup_out - kin_cgroup;
	int tmin= oneM - ktot_cgroup_in;
	int tmout= oneM - ktot_cgroup_out;
	
	double wf;
	cgroup.worst_node(wnode, wf, kout_g_in, kout_g_out, Nstar, nn, tmin, tmout);
	
	//------------------------------ compute the worst internal node
	
	
	
	if(wnode==-1)
		return 1;

	
	return 0;

}



bool clomnet_local::insert_the_best(double threshold) {

	
	int benode;
	double cmax;
	return insert_the_best(cmax, threshold, benode);


}


bool clomnet_local::insert_the_best(double & cmax, double threshold, int & benode) {
	
	// this function is to look for the best node in neigh and to insert it if its cscore it is lower than threshold
	
	
	//cout<<"\n\n\n***********************************inserting the best..."<<endl;
	/*cout<<"neighs..."<<endl;
	neighs.print_nodes(cout);
	cout<<"cgroup..."<<endl;
	cgroup.print_nodes(cout);
	*/

	
	
	int Nstar= dim - cgroup.size();
	int nn= neighs.size();
	int kout_g_in = ktot_cgroup_in - kin_cgroup;
	int kout_g_out = ktot_cgroup_out - kin_cgroup;
	int tmin= oneM - ktot_cgroup_in;
	int tmout= oneM - ktot_cgroup_out;
	
	//cout<<"nstar "<<Nstar<<" "<<nn<<" "<<kout_g<<" "<<tm<<endl;
	
	
	
	double lowest_r=1;
	neighs.best_node(benode, lowest_r, kout_g_in, kout_g_out, Nstar, nn, tmin, tmout);
	
	if(benode==-1)
		return false;
	
	cmax= 1 - pow(1 - lowest_r, nn);


	if(cmax>threshold)
		return false;
	
	insert_cgroup(benode);


	//cherr();

	
	return true;

}

bool clomnet_local::erase_the_worst(double threshold) {

	
	
	int wnode;
	bool b= erase_the_worst(threshold, wnode);

	//cout<<"worst node: "<<id_of(wnode)<<" "<<kw<<" "<<degree_ref<<" "<<cmax<<endl;
	
	return b;
}




bool clomnet_local::erase_the_worst(double threshold, int & wnode) {
	
	// this function is to look for the worst node in cgroup and to erase it if its cscore it is bigger than threshold
	
	
	CBS(wnode);

	
	if(cgroup.size()==0) {
		
		//cerr<<"impossible to erase::::empty!!"<<endl;
		return false;
	
	
	}

		
		
	erase_cgroup(wnode);

		
	return true;

}


bool clomnet_local::erase_the_worst(double threshold, int & wnode, set<int> & erased) {
	
	// this function is to look for the worst node in cgroup and to erase it if its cscore it is bigger than threshold
	
	
	CBS(wnode);
	if(cgroup.size()==0)
		return false;

		
			
	erase_cgroup(wnode, erased);

		
	return true;

}



void clomnet_local::insertion(int changendi, deque<int> & inserted, bool suubb) {


	for(int i=0; i<changendi; i++) {
		int benode;
		double cwn;
		insert_the_best(cwn, 1.1, benode);
		if(benode!=-1)
			inserted.push_back(benode);
	}
	
	
	
	while(true) {
		
		int benode;
		double cmax;
		if(insert_the_best(cmax, paras.threshold, benode)) {
			
			inserted.push_back(benode);
			
			
		} else
			break;
	
	}
	
	/*
	if(suubb==false && inserted.size()>dim*0.85 && dim>100) {
		
		
		//cout<<"network of dim: "<<dim<<endl;
		cout<<"WARNING: if the algorithm is taking long there might be problems!"<<endl; 
	
	}
	//*/
	

}







void clomnet_local::set_bg() {
	
	
	if(lcb_temp<lowest_cs) {
			
			lowest_cs=lcb_temp;
			gstopper=0;
			
			cgroup.set_deque(core_p);

	}
	
}










int clomnet_local::set_cgroup_and_neighs(const deque<int> & G) {

	
	kin_cgroup=0;
	ktot_cgroup_in=0;
	ktot_cgroup_out=0;
	cgroup.clear();
	neighs.clear();
	for(int i=0; i<G.size(); i++)
		insert_cgroup(G[i]);
	
	
	return 0;
	

}







double clomnet_local::refinement_single_only_c(deque<int> & gr_cleaned, set<int> & ONLY_AMONG_THIS) {
	
	
	
	// gr_cleaned is the group with the border included
	// border_group is just the set of nodes in the border if there are more than one
	// this function is going to work including only the node in the original group
	
	

	gr_cleaned.clear();
	//border_group.clear();
	
	
	if(cgroup.size()==dim) {
		
		//cerr<<"all the network"<<endl;
		return 1;
	
	}
	
	
	/*neighs.external_xs(ktot_cgroup - kin_cgroup, twoM - ktot_cgroup, dim - cgroup.size(), neighs.size(), logrs_lab, maxb_nodes/double(dim-cgroup.size()), max_r_bord, ONLY_AMONG_THIS);
	
	int until;
	double bs= ord_border_score(logrs_lab, until);
	*/
	
	
	double b_sc=1.1;
	
	
	while(cgroup.size()>1) {
		
		int wnode, until;
		bool b= erase_the_worst(-1, wnode);
		deque<pair<double, int> > logrs_lab;
		neighs.external_xs(ktot_cgroup_in - kin_cgroup, ktot_cgroup_out - kin_cgroup,
							oneM - ktot_cgroup_in, oneM - ktot_cgroup_out, 
							dim - cgroup.size(), neighs.size(), logrs_lab, maxb_nodes/double(dim-cgroup.size()), max_r_bord, ONLY_AMONG_THIS);

		/*
		cout<<"worst node:  only "<<id_of(wnode)<<endl;
		//cherr();
		cout<<"************  external nodes...."<<endl;
		deque<int> HH;
		cgroup.set_deque(HH);
		for(int i=0; i<logrs_lab.size(); i++) {
			
			cout<<"node: "<<id_of(logrs_lab[i].second)<<" fit  "<<logrs_lab[i].first<<" links: "<<vertices[logrs_lab[i].second]->kplus_m(HH)<<" mtlw "<<vertices[logrs_lab[i].second]->kplus_w(HH)<<endl;
		
		
		}
		cout<<"^_^_^_^_^_^_^_^_  external nodes...."<<endl;
		//*/
		
		double bs= ord_border_score(logrs_lab, until);
		int bsize= logrs_lab.size()-until;
		
		if(bsize>0) {
						
			
			cgroup.set_deque(gr_cleaned);
			for(int j=until; j<logrs_lab.size(); j++) {
				gr_cleaned.push_back(logrs_lab[j].second);
				/*if(bsize>1)
					border_group.push_back(logrs_lab[j].second);*/
			}
			b_sc=bs;
			
			break;
		
		
			
			
		}
	}
	
	
		
	
		
	return b_sc;



}




double clomnet_local::refinement_single(deque<int> & gr_cleaned) {
	
	
	// gr_cleaned is the group with the border included
	// border_group is just the set of nodes in the border if there are more than one
	
	
	
	
	//cout<<"group: "<<cgroup.size()<<endl;
	

	gr_cleaned.clear();
	//border_group.clear();
	
	
	if(cgroup.size()==dim) {
		
		//cerr<<"all the network"<<endl;
		return 1;
	
	}
	
	
	
	//cgroup.print_nodes(cout);
	
	//cout<<"refinement_single:: "<<ktot_cgroup_in<<" "<<ktot_cgroup_out<<" "<<kin_cgroup<<endl;
	//cherr();
	
	
	
	
	
	deque<pair<double, int> > logrs_lab;
	neighs.external_xs(ktot_cgroup_in - kin_cgroup, ktot_cgroup_out - kin_cgroup,
						oneM - ktot_cgroup_in, oneM - ktot_cgroup_out, 
						dim - cgroup.size(), neighs.size(), logrs_lab, maxb_nodes/double(dim-cgroup.size()), max_r_bord);
	
	
	
	
	/*
	deque<int> HH;
	cgroup.set_deque(HH);
	print_id(HH, cout);
	cout<<"************  external nodes...."<<endl;
	for(int i=0; i<logrs_lab.size(); i++) {
		
		cout<<"node: "<<id_of(logrs_lab[i].second)<<" fit  "<<logrs_lab[i].first<<" links: "
		<<(vertices[logrs_lab[i].second]->kplus_m(HH)).first<<" "<<vertices[logrs_lab[i].second]->kplus_m(HH).second
		<<" mtlw "<<vertices[logrs_lab[i].second]->kplus_w(HH).first<<" "<<vertices[logrs_lab[i].second]->kplus_w(HH).second<<endl;
		
		
	}
	cout<<"^_^_^_^_^_^_^_^_  external nodes...."<<endl;
	//*/
	int until;
	double bs= ord_border_score(logrs_lab, until);
	
	
	
	
	
	for(int j=until; j<logrs_lab.size(); j++)
		insert_cgroup(logrs_lab[j].second);

	//check_all();
	// now, it starts the dynamical update - to erase bad nodes and include good ones.
	
	//cout<<"starting dynamical update ... "<<endl;
	
	
	double b_sc=1.1;
	
	
	while(cgroup.size()>1) {
		
		int wnode;
		bool b= erase_the_worst(-1, wnode);
		
		
		
		/*cout<<"**************************************************"<<endl;
		cout<<"**************************************************"<<endl;
		cout<<"**************************************************"<<endl;

		cout<<"WORST node: "<<wnode<<endl;
		
		cgroup.print_nodes(cout);
		neighs.print_nodes(cout);

		cout<<"ktot_cgroup_in-out "<<ktot_cgroup_in<<" "<<ktot_cgroup_out<<" kout_in-out"<<
		ktot_cgroup_in - kin_cgroup<<" "<<ktot_cgroup_out - kin_cgroup<<endl;

		cout<<"**************************************************"<<endl;
		cout<<"**************************************************"<<endl;
		cout<<"**************************************************"<<endl;*/
		
		neighs.external_xs(ktot_cgroup_in - kin_cgroup, ktot_cgroup_out - kin_cgroup,						
							oneM - ktot_cgroup_in, oneM - ktot_cgroup_out, 						
							dim - cgroup.size(), neighs.size(), logrs_lab, maxb_nodes/double(dim-cgroup.size()), max_r_bord);
		
		
		
		
		/*
		cout<<"************  external nodes...."<<endl;
		cgroup.set_deque(HH);

		
		for(int i=0; i<logrs_lab.size(); i++) {
		
			cout<<"node: "<<id_of(logrs_lab[i].second)<<" fit  "<<logrs_lab[i].first<<" links: "
			<<vertices[logrs_lab[i].second]->kplus_m(HH).first<<" "<<vertices[logrs_lab[i].second]->kplus_m(HH).first
			<<" mtlw "<<vertices[logrs_lab[i].second]->kplus_w(HH).first<<" "<<vertices[logrs_lab[i].second]->kplus_w(HH).second<<endl;
		
		
		}
		cout<<"^_^_^_^_^_^_^_^_  external nodes...."<<endl;//*/
		
		
		//check_all();

		
		
		bs= ord_border_score(logrs_lab, until);
		int bsize= logrs_lab.size()-until;
		
		if(bsize>0) {
			
			
			
			cgroup.set_deque(gr_cleaned);
			for(int j=until; j<logrs_lab.size(); j++) {
				gr_cleaned.push_back(logrs_lab[j].second);
				/*if(bsize>1)
					border_group.push_back(logrs_lab[j].second);*/
			}
			b_sc=bs;
			
			//
			//cout<<"group, bsize, all "<<cgroup.size()<<" "<<bsize<<" "<<gr_cleaned.size()<<endl;

			
			break;
		
		
			
			
		}
	}
	
	
	
		
	
	
	return b_sc;



}




double clomnet_local::ord_border_score(deque<pair<double, int> > & xs, int & until) {
	
	
	// this function says until where one could include the border
	// N are the nodes out fo the cluster,
	// th is the threshold
	
	until=xs.size();
	
	int Nstar;
	if(paras.weighted)
		Nstar=neighs.size();
	else
		Nstar=dim-cgroup.size();

	double critical_xi= -log(1-paras.threshold)/fitted_exponent(dim - cgroup.size());
	//cout<<"critical xi "<<critical_xi<<" "<<cgroup.size()<<" - xs.size "<<xs.size()<<" "<<neighs.size()<<" fitted exp: "<<fitted_exponent(dim - cgroup.size())<<" log ... "<<log(1-paras.threshold)<<endl;
	
	
	int pos=Nstar;
	double min_xi=1;
	for(int j=xs.size()-1; j>=0; j--) {
		
		
		double min_=right_tail_order_exponential(Nstar, pos, xs[j].first);
		//cout<<"fitness, score: "<<xs[j].first<<" "<<min_<<" position: "<<pos<<" label: "<<id_of(xs[j].second)<<endl;

		--pos;
				
		if(min_<critical_xi) {
			until=j;
			min_xi=min_;
		
			break;

		}
		
		
	}
	
	
	return pron_min_exp(dim - cgroup.size(), min_xi);

}





double clomnet_local::b_clean_up_search_iterative(deque<int> _c_, deque<int> & gr_cleaned, int l) {



	//*
	{
	
		// check *******************************
		deque<int> inserted;
		set_cgroup_and_neighs(_c_);
		insertion(0, inserted, true);
	
	
		if(double(cgroup.size())>dim*0.85)
			return b_clean_up_search(_c_, gr_cleaned, l);
				
		
		// check *******************************
	
	}//*/



	double bs= b_clean_up_search(_c_, gr_cleaned, l);
	//return bs;
	//cout<<"----------------------------"<<endl;
	
	
	int stopp=0;
	
	while(gr_cleaned.size()>_c_.size()) {
		
		_c_=gr_cleaned;
		bs = b_clean_up_search(_c_, gr_cleaned, l);
		++stopp;
		
		
		//cout<<_c_.size()<<" "<<gr_cleaned.size()<<endl;
		/*print_id(_c_, cout);
		print_ri(_c_);*/
		if(stopp==paras.iterative_stopper)
			break;
	
	}
	
	//cout<<"----------------------------"<<endl;
	
	return bs;
	
}



double clomnet_local::b_clean_up_search(const deque<int> & _c_, deque<int> & gr_cleaned, int l) {

	
	
	deque<int> first;
	//cout<<"b_clean_up_once... loops "<<l<<endl;
	double bs=b_clean_up_once(_c_, gr_cleaned, l, false);

	if(bs<1) {
	
		first=gr_cleaned;
		bs=b_clean_up_once(first, gr_cleaned, l, true);
	}
	
	
	
	return bs;
	
}




/*double clomnet_local::b_clean_up_check(const deque<int> & _c_, deque<int> & gr_cleaned, int l) {

	
	double bs=b_clean_up_once(_c_, gr_cleaned, l, false);
	deque<int> first=gr_cleaned;
			

	if(bs<1)
		bs=b_clean_up_once(first, gr_cleaned, l, true);
		
	
	
	
	return bs;
	
}
*/



double clomnet_local::b_clean_up_check(const deque<int> & _c_, deque<int> & gr_cleaned, int l) {

	// double check!

	b_clean_up_once(_c_, gr_cleaned, l, true);
	
	
	
	
	deque<int> first=gr_cleaned;
	
	/* uncomment this for a triple check!
	b_clean_up_once(first, gr_cleaned, l, true);
	first=gr_cleaned;
	//*/
	
	return b_clean_up_once(first, gr_cleaned, l, true);
	

}




double clomnet_local::b_clean_up_once(const deque<int> & _c_, deque<int> & gr_cleaned, int loops, bool only_c) {
		
	
	
	double participation_rate;		// this is the number of time a node must appear in first to be part of gr_cleaned
	
	if(only_c==true)
		participation_rate=paras.part_rate;
	else
		participation_rate=paras.part_rate_search;	
	
	
	gr_cleaned.clear();
	
	//cout<<"SIZE: "<<_c_.size()<<endl;
	
	if(_c_.size()<=1)
		return 1;
	
	
	
	deque<double> bss;
	map<int, int> part_rate_firsts;
	
	
	weighted_tabdeg previous_tab1;
	weighted_tabdeg previous_tab2;
	
	set_cgroup_and_neighs(_c_);
	previous_tab1._set_(cgroup);
	previous_tab2._set_(neighs);
	
	int kin_cgroup_prev=kin_cgroup;
	int ktot_cgroup_prev_in=ktot_cgroup_in;
	int ktot_cgroup_prev_out=ktot_cgroup_out;
	
	set<int> ONLY_AMONG_THIS;
	deque_to_set(_c_, ONLY_AMONG_THIS);
	
	
	int ii=0;
		

	
	while(true) {
	
		ii++;
		//cout<<"****************************************************************loops: "<<ii<<" l: "<<loops<<endl;
	
		
		
		
		double bs=1;

		deque<int> first;
		
		if(only_c==false) {
			
			cgroup._set_(previous_tab1);
			neighs._set_(previous_tab2);
			kin_cgroup=kin_cgroup_prev;
			ktot_cgroup_in=ktot_cgroup_prev_in;
			ktot_cgroup_out=ktot_cgroup_prev_out;
			
			

			
			bs=refinement_single(first);
		
		}
		else {
		
			cgroup._set_(previous_tab1);
			neighs._set_(previous_tab2);
			kin_cgroup=kin_cgroup_prev;
			ktot_cgroup_in=ktot_cgroup_prev_in;
			ktot_cgroup_out=ktot_cgroup_prev_out;
			

			bs=refinement_single_only_c(first, ONLY_AMONG_THIS);
				
		}
		
		

		if(bs<1) {
			
			bss.push_back(bs);
			for(int j=0; j<first.size(); j++)
				int_histogram(first[j], part_rate_firsts);
		
		}
		
				
		if(ii%50==0) {
			
			bool stopnow=true;
			
			

			if( bss.size() - 0.5 * ii < 3 * sqrt(bss.size() * (1 - double(bss.size())/ii)))
				stopnow=false;
			else if(sqrt(variance_func(bss)/bss.size())>0.01)
				stopnow=false;
			else if(only_c==true) for(map<int, int>::iterator itm= part_rate_firsts.begin(); itm!=part_rate_firsts.end(); itm++ ) {
				if(fabs(itm->second - participation_rate*bss.size() ) < sqrt(double(itm->second))) {
					
					stopnow=false;
					break;
				}
				
			}
			
			if(stopnow) {
				//cout<<"BREAKER::: "<<ii<<endl;
				break;
			}
			
			
		
		}
		
		

		
		if(ii==loops)
			break;
		
		
	}
	
	
	
	
	if(double(bss.size())>0.5*ii) {				// it has to find a group below threshold at least half loops
				
		include_prate(part_rate_firsts, gr_cleaned, bss.size(), participation_rate);		
		return average_func(bss);
	}	
	
	
	
	
	return 1;
	
	
	
}


void clomnet_local::include_prate(map<int, int> part_rate_firsts, deque<int> & gr_cleaned, int times, double participation_rate_par) {
	
		
	for(map<int, int>::iterator itm=part_rate_firsts.begin(); itm!=part_rate_firsts.end(); itm++) 
		if(itm->second>participation_rate_par*times) {
			gr_cleaned.push_back(itm->first);
		}
		
}








double clomnet_local::setic_simann(int A, deque<int> & ic, double clean_t, const deque<double> & changendi_cum, int min_p, int maxstopper, bool suubb) {
	
	// A is the initial node
	// ic is where the ic will be set if any
	// the others are parameters
	
	
		
	core_p.clear();
	
	kin_cgroup=0;
	ktot_cgroup_in=0;
	ktot_cgroup_out=0;
	cgroup.clear();
	neighs.clear();
	
	lowest_cs=0;
	lcb_temp=0;

	
	//************************** PARAMETERS for the simulated annealing which I can keep constant (I believe)
	
	
	temp=0.1;
	//if(first_shell)
		//temp=1e-7;
	
	
	double tstep=0.999;
	
	gstopper=0;
	
	
	
	double trigger=log(paras.threshold);
	
	
	
	//************************** PARAMETERS for the simulated annealing...
	
	
	
	//-------------------------  Initialization
	
	insert_cgroup(A);
	int noterasable=A;
	
	
	
	
	bool not_triggered= true;

	while(true) {
		
		
		
		if(not_triggered && lowest_cs<trigger) {
			
			not_triggered=false;
			maxstopper=paras.setic_stopper_local/2;
		
		}
				
		temp*=tstep;


		int changendi=lower_bound(changendi_cum.begin(), changendi_cum.end(), ran4()) - changendi_cum.begin() + min_p;
		changendi=min(changendi, neighs.size());		
		
		//cout<<"changendi: "<<changendi<<endl;
		
		
		if(setic_clean_up(noterasable, changendi, clean_t, suubb)<0)
			return 1;
		
		
		
		set_bg();
		
				
		
		if(gstopper==maxstopper)	
			break;
		
		gstopper++;
		
		//check_all();
		

	}

	
	double bcu;
	//if(first_shell==false) {
	
	bcu=b_clean_up_search_iterative(core_p, ic, 100);
	
		
	return bcu;

}





int clomnet_local::setic_clean_up(int noterasable, int changendi, const  double & clean_t, bool suubb) {
	
		
	deque<int> inserted;
	deque<int> inserted2;
	deque<int> erased;
	insertion(changendi, inserted, suubb);
	
	
		
	

	
	
	
	
	if(cgroup.size()==dim) {
			
		return -1;
	
	}


		
	deque<pair<double, int> > logrs_lab;
	//cout<<"neighs: "<<neighs.size()<<" "<<ktot_cgroup<<" "<<cast_int(2*tstrength) - ktot_cgroup<<endl;
	neighs.external_xs(ktot_cgroup_in - kin_cgroup, ktot_cgroup_out - kin_cgroup,						
						oneM - ktot_cgroup_in, oneM - ktot_cgroup_out, 
						dim - cgroup.size(), neighs.size(), logrs_lab, maxb_nodes/double(dim-cgroup.size()), max_r_bord);
	
	
	
	int until;
	double bs= ord_border_score(logrs_lab, until);
	
	
	//cout<<"this ord_border_score done... "<<until<<" "<<logrs_lab.size()<<endl;
	
	for(int j=until; j<logrs_lab.size(); j++) {
		insert_cgroup(logrs_lab[j].second);
		inserted.push_back(logrs_lab[j].second);
	}

	double bs_=1;
	while(cgroup.size()>1) {
		
		int wnode;
		bool b= erase_the_worst(-1, wnode);
		erased.push_back(wnode);
		
		
		neighs.external_xs(ktot_cgroup_in - kin_cgroup, ktot_cgroup_out - kin_cgroup,						
							oneM - ktot_cgroup_in, oneM - ktot_cgroup_out, 						
							dim - cgroup.size(), neighs.size(), logrs_lab, maxb_nodes/double(dim-cgroup.size()), max_r_bord);
		
			
		bs= ord_border_score(logrs_lab, until);
		int bsize= logrs_lab.size()-until;
		
		
		
		if(bsize>0) {
			
			int average_t=10;
			deque<double> bss_good;
			deque<double> bss_all;
					
			for(int k=0; k<average_t; k++) {
			
			
				neighs.external_xs(ktot_cgroup_in - kin_cgroup, ktot_cgroup_out - kin_cgroup,						
									oneM - ktot_cgroup_in, oneM - ktot_cgroup_out,
									dim - cgroup.size(), neighs.size(), logrs_lab, maxb_nodes/double(dim-cgroup.size()), max_r_bord);
				
				double hh=ord_border_score(logrs_lab, until);
				if(hh<clean_t) {
					bss_good.push_back(hh);
				}
				
				
				bss_all.push_back(hh);
			}
			
		
			
			
			if(bss_good.size()> average_t/2) {
				bs=average_func(bss_good);
			} else {
				bs=average_func(bss_all);
			}
			
			


			for(int j=until; j<logrs_lab.size(); j++) {
				insert_cgroup(logrs_lab[j].second);
				inserted2.push_back(logrs_lab[j].second);
			}
			
		
			bs_= bs;
			break;
		}
	
	}
	
	bool accept_the_move=true;
	
	if(cgroup.is_internal(noterasable)==false)
		accept_the_move=false;
	else {
	


		double newcs=log_zero(bs_);
		double diff= lcb_temp - newcs;
		
		double bol= exp(diff/temp);		// I want the difference to be positive
		
		if(ran4()<bol)
			accept_the_move=true;
		else 
			accept_the_move=false;
		
	
	}
	
	if(accept_the_move) {
		//cout<<"move accepted- newscore: "<<bs_<<" oldscore: "<<exp(lcb_temp)<<endl;
		lcb_temp=log_zero(bs_);
		//if(paras.print_cbs)
			//cout<<"move accepted: "<<bs_<<endl;
	}
	else {
		
		//cout<<"move refused- newscore: "<<bs_<<" oldscore: "<<exp(lcb_temp)<<endl;
		
		
		for(int i=0; i<inserted2.size(); i++)
			erase_cgroup(inserted2[i]);
		
		
		for(int i=0; i<erased.size(); i++)
			insert_cgroup(erased[i]);

	
		
		
		for(int i=0; i<inserted.size(); i++) if(cgroup.is_internal(inserted[i]))
			erase_cgroup(inserted[i]);
		
		//cout<<" inserted "<<endl;
		//print_id(inserted, cout);
		
	}

	
		
	return 0;

}


//*

/*
int clomnet_local::check_all() {
	

	deque<int> A;
	cgroup.set_deque(A);
	
	
	//cout<<"checking this "<<kin_cgroup<<endl;
	//prints(A);
	
	int kinc=0;
	int kinc2=0;
	int ktotcin=0;
	int ktotcout=0;
	
	for(map<int, facts> :: iterator itm =cgroup.lab_facts.begin(); itm!=cgroup.lab_facts.end(); itm++) {
	
		cherr(itm->second.indegree - vertices[itm->first]->instub_number, 1e-6);
		cherr(itm->second.outdegree - vertices[itm->first]->outstub_number, 1e-6);
		
		
		//cout<<"checking group node: "<<itm->first<<" "<<itm->second.indegree<<" "<<vertices[itm->first]->instub_number<<" "<<itm->second.internal_indegree<<endl;
		
		
		cherr(itm->second.internal_indegree - (vertices[itm->first]->kplus_m(A)).first, 1e-6);
		cherr(itm->second.internal_outdegree - (vertices[itm->first]->kplus_m(A)).second, 1e-6);


		cherr(itm->second.minus_log_total_wrin - (vertices[itm->first]->kplus_w(A)).first, 1e-6);
		cherr(itm->second.minus_log_total_wrout - (vertices[itm->first]->kplus_w(A)).second, 1e-6);
		
		kinc+=itm->second.internal_indegree;
		kinc2+=itm->second.internal_outdegree;
		ktotcin+=itm->second.indegree;
		ktotcout+=itm->second.outdegree;
	
	}
	
	//cout<<"group passed"<<endl;
	
	
	for(map<int, facts> :: iterator itm =neighs.lab_facts.begin(); itm!=neighs.lab_facts.end(); itm++) {



		cherr(itm->second.indegree - vertices[itm->first]->instub_number, 1e-6);
		cherr(itm->second.outdegree - vertices[itm->first]->outstub_number, 1e-6);
		

		//cout<<"checking node: "<<itm->first<<" "<<itm->second.indegree<<" "<<vertices[itm->first]->instub_number<<" "<<itm->second.internal_indegree<<endl;

		
		
		cherr(itm->second.internal_indegree - (vertices[itm->first]->kplus_m(A)).first, 1e-6);
		cherr(itm->second.internal_outdegree - (vertices[itm->first]->kplus_m(A)).second, 1e-6);

		/*cout<<"checking node2: "<<itm->first<<" "<<itm->second.indegree<<" "<<vertices[itm->first]->instub_number<<" "<<itm->second.internal_indegree<<endl;

		cout<<"checking nod3: "<<itm->first<<" "<<itm->second.minus_log_total_wrin<<" "<<(vertices[itm->first]->kplus_m(A)).first<<" "<<endl;
		cout<<"checking nod3: "<<itm->first<<" "<<itm->second.minus_log_total_wrout<<" "<<(vertices[itm->first]->kplus_m(A)).second<<" ***  "<<endl;
		
		/*cout<<"wwwwwwwww "<<vertices[2]->outlinks->w[0].second<<" "<<vertices[2]->outlinks->l[0]<<endl;
		cout<<"wwwwwwwww "<<vertices[2]->outlinks->w[1].second<<" "<<vertices[2]->outlinks->l[1]<<endl;
		cout<<"wwwwwwwww "<<vertices[2]->outlinks->w[2].second<<" "<<vertices[2]->outlinks->l[2]<<endl;* /
		
		cherr(itm->second.minus_log_total_wrin - (vertices[itm->first]->kplus_w(A)).first, 1e-6);
		cherr(itm->second.minus_log_total_wrout - (vertices[itm->first]->kplus_w(A)).second, 1e-6);
		
		
	}
	
	
	//cout<<"neighs passed"<<endl;

	
	set<int> AA;
	deque_to_set(A, AA);
	
	cherr(kin_cgroup- kin_m(A), 1e-6);
	cherr(kin_cgroup- kin_m(AA), 1e-6);

	cherr(ktot_cgroup_in - ktot_m(A).first, 1e-6);
	cherr(ktot_cgroup_in - ktot_m(AA).first, 1e-6);


	cherr(ktot_cgroup_out - ktot_m(A).second, 1e-6);
	cherr(ktot_cgroup_out - ktot_m(AA).second, 1e-6);


	cherr(kin_cgroup-kinc, 1e-6);
	cherr(kin_cgroup-kinc2, 1e-6);
	cherr(ktot_cgroup_in-ktotcin, 1e-6);
	cherr(ktot_cgroup_out-ktotcout, 1e-6);
	
	//cout<<"check passed"<<endl;
	return 0;






}
//*/








int clomnet_local::set_upper_network(map<int, map<int, pair<int, double> > > & neigh_weight_f, module_collection & module_coll) {
	
	
	
	// loop on all the edges of the network...

	neigh_weight_f.clear();

	if(module_coll.size()==0)
		return 0;
	
	
	
	map<int, map<int, pair<double, double> > >neigh_weight_s;


	
	for(map<int, double>::iterator its = module_coll.module_bs.begin(); its!=module_coll.module_bs.end(); its++) {

		
		//cout<<"NAMES:: "<<m_name<<endl;
		map<int, pair<double, double> > neigh_weight;
		map<int, pair<int, double> > ooo;
		neigh_weight_s.insert(make_pair(its->first, neigh_weight));
		neigh_weight_f.insert(make_pair(its->first, ooo));
		
	}
	
	
	
	for(int i=0; i<dim; i++) {
		
		set<int> & mem1= module_coll.memberships[i];
		
		for(int j=0; j<vertices[i]->outlinks->size(); j++) {
			
			int & neigh = vertices[i]->outlinks->l[j];
			set<int> & mem2= module_coll.memberships[neigh];
			
			double denominator = mem1.size() * mem2.size();
			// I add a link between all different modules
			
			
			//cout<<"denomi "<<denominator<<endl;
			
			
			//**************************************************************************************************
			if(paras.weighted) {
				for(set<int>:: iterator itk=mem1.begin(); itk!=mem1.end(); itk++) for(set<int>:: iterator itkk=mem2.begin(); itkk!=mem2.end(); itkk++) if(*itk!=*itkk) {
				
					int_histogram(*itkk, neigh_weight_s[*itk], double(vertices[i]->outlinks->w[j].first)/denominator, vertices[i]->out_original_weights[j]/denominator);
				}
			} else { 
				
				
				
				for(set<int>:: iterator itk=mem1.begin(); itk!=mem1.end(); itk++) for(set<int>:: iterator itkk=mem2.begin(); itkk!=mem2.end(); itkk++) if(*itk!=*itkk) {
				
					int_histogram(*itkk, neigh_weight_s[*itk], double(vertices[i]->outlinks->w[j].first)/denominator, 1);
			
				}

			}
			//**************************************************************************************************

		
		}
		
		
	
	}
	
	
	
	for(map<int, map<int, pair<double, double> > >:: iterator itm= neigh_weight_s.begin(); itm!= neigh_weight_s.end(); itm++) {
		
		
		
		for(map<int, pair<double, double> >::iterator itm_ = itm->second.begin(); itm_!=itm->second.end(); itm_++) {
			
			
			neigh_weight_f[itm->first].insert(make_pair(itm_->first, make_pair( cast_int(itm_->second.first), itm_->second.second)));
			
			
		
		}
	
	}
	
	
	
		
	
	
	
	return 0;
	
}







int clomnet_local::try_to_assign_homeless(module_collection & module_coll, bool anyway) {
	
	
	


	map<int, deque<int> > to_check;
	
	
	
	
	int before_procedure = try_to_assign_homeless_help(module_coll, to_check);
	
	bool something=false;
	
	for(map<int, deque<int> > :: iterator itm= to_check.begin();  itm!=to_check.end(); itm++) {
	
		deque<int> union_deque=module_coll.modules[itm->first];
		
		/*
		if(union_deque.size()==0) {
			
			cout<<"what the hell??? "<<endl;
		
		}*/
		
		
		for(int i=0; i<itm->second.size(); i++)
			union_deque.push_back(itm->second[i]);
		
		/*
		cout<<"union_deque "<<itm->first<<" size: "<<union_deque.size()<<endl;
		print_id(union_deque, cout);
		print_ri(union_deque);
		//*/
		
		
		if(anyway) {

			something=true;
			module_coll.insert(union_deque, ran4()+paras.threshold);
		}
		else {
		
			deque<int> grbe;
			double bs=b_clean_up_check(union_deque, grbe,  paras.bscore_loops_check);
			
			//cout<<"union_deque after "<<itm->first<<" size: "<<grbe.size()<<endl;
			
			
			if(grbe.size()>1) {
				
				something=true;
				module_coll.insert(grbe, bs);
				
			}
			
		}
	
	}
	
	
	if(something)
		module_coll.compute_inclusions(false);
	




	return before_procedure;

}









int clomnet_local::check_pairs(multimap<int, int> & A, double clean_t) {
	
	
	
	// this function needs to be checked ....
	
	A.clear();

	for(int i=0; i<dim; i++) {
		
		for(int j=0; j<vertices[i]->outlinks->size(); j++) {
			
			int & neigh= vertices[i]->outlinks->l[j];
			//cout<<id_of(i)<<" "<<id_of(neigh)<<" kin: "<<vertices[i]->links->w[j].first<<" mtlw "<<vertices[i]->links->w[j].second<<endl;
			
			
			int x= vertices[i]->outlinks->w[j].first;
			int k1= vertices[i]->outstub_number;
			int k2= vertices[vertices[i]->outlinks->l[j]]->instub_number;
			int tm = oneM - vertices[vertices[i]->outlinks->l[j]]->outstub_number;
			//double win =  vertices[i]->outlinks->w[j].second;
			
			/*if(id_of(i)==1)
				cherr();*/
			//double r= LOG_TABLE.cum_hyper_right(x, k2, tm, k1);
			//cout<<"ok"<<endl;
			
			
			double r= compute_global_fitness_randomized_short(x, k2, tm, k1, vertices[i]->outlinks->w[j].second);
			
			double c= 1. - pow(1.-r, dim-1);
			//cout<<"->check_pairs: "<<i<<" "<<neigh<<"; x: "<<x<<" "<<k1<<" "<<k2<<" "<<vertices[i]->outlinks->w[j].second<<" r: "<<r<<" "<<c<<endl;
			
			if(c<clean_t) {
				//cout<<"check_pairs: "<<id_of(i)<<" "<<id_of(vertices[i]->links->l[j])<<" "<<c<<endl;
				//cherr();
				A.insert(make_pair(id_of(i), id_of(vertices[i]->outlinks->l[j])));
			}
		
		
		}
	
	
	}

	


	return 0;


}








void clomnet_local::erase_cgroup(int wnode) {
	
	
	
	/*cout<<"erasing erase_cgroup "<<wnode<<endl;
	
	
	
	cout<<"before erasing "<<endl;
	
	cgroup.print_nodes(cout);
	neighs.print_nodes(cout);
	
	cout<<" ********************* ... erasing ************** "<<endl;*/
	
	
	
	
	
	
	map<int, facts>::iterator itm= cgroup.lab_facts.find(wnode);
	if(itm!=cgroup.lab_facts.end()) {

		int kpin = itm->second.internal_indegree;
		int ktin = itm->second.indegree;
		int kpout = itm->second.internal_outdegree;
		int ktout = itm->second.outdegree;
		double mtlwin= itm->second.minus_log_total_wrin;
		double mtlwout= itm->second.minus_log_total_wrout;


		kin_cgroup-=  kpin + kpout;
		ktot_cgroup_out-=  ktout;
		ktot_cgroup_in-=  ktin;
		
		
		int kout_g_in= ktot_cgroup_in - kin_cgroup;
		int kout_g_out= ktot_cgroup_out - kin_cgroup;
		
		int tmin= oneM - ktot_cgroup_in;
		int tmout= oneM - ktot_cgroup_out;

		
		
		double fi= compute_global_fitness_randomized(kpin, kout_g_in, kpout, kout_g_out, tmin, tmout, ktin, ktout, mtlwin, mtlwout, neighs.size()+1, dim - cgroup.size() +1);
		neighs.edinsert(wnode, kpin, kpout, ktin, ktout, mtlwin, mtlwout, fi);
		
		
		cgroup.erase(wnode);


		deque<int> tobe;
		pair <int, pair<int, double> > OPA;
	
	
		
		
		for(int i=0; i<vertices[wnode]->inlinks->size(); i++) {
			
			
			
			OPA = vertices[wnode]->outlinks->posweightof(vertices[wnode]->inlinks->l[i]);

			
			if(cgroup.update_group(vertices[wnode]->inlinks->l[i], -vertices[wnode]->inlinks->w[i].first, -OPA.second.first, -vertices[wnode]->inlinks->w[i].second, -OPA.second.second, 
						dim - cgroup.size(), neighs.size(), kout_g_in, kout_g_out, tmin, tmout, vertices[vertices[wnode]->inlinks->l[i]]->instub_number, vertices[vertices[wnode]->inlinks->l[i]]->outstub_number,tobe)==false)
				neighs.update_neighs(vertices[wnode]->inlinks->l[i], -vertices[wnode]->inlinks->w[i].first, -OPA.second.first, -vertices[wnode]->inlinks->w[i].second, -OPA.second.second, 
						dim - cgroup.size(), kout_g_in, kout_g_out, tmin, tmout, vertices[vertices[wnode]->inlinks->l[i]]->instub_number, vertices[vertices[wnode]->inlinks->l[i]]->outstub_number);
		
		}
		
		
		
		for(int i=0; i<vertices[wnode]->outlinks->size(); i++) {
		
		
					
			OPA = vertices[wnode]->inlinks->posweightof(vertices[wnode]->outlinks->l[i]);
		
			if(OPA.first==-1) if(cgroup.update_group(vertices[wnode]->outlinks->l[i], 0, -vertices[wnode]->outlinks->w[i].first, 0, -vertices[wnode]->outlinks->w[i].second,
						dim - cgroup.size(), neighs.size(), kout_g_in, kout_g_out, tmin, tmout, vertices[vertices[wnode]->outlinks->l[i]]->instub_number, vertices[vertices[wnode]->outlinks->l[i]]->outstub_number,tobe)==false)
				neighs.update_neighs(vertices[wnode]->outlinks->l[i], 0, -vertices[wnode]->outlinks->w[i].first, 0, -vertices[wnode]->outlinks->w[i].second,
						dim - cgroup.size(), kout_g_in, kout_g_out, tmin, tmout, vertices[vertices[wnode]->outlinks->l[i]]->instub_number, vertices[vertices[wnode]->outlinks->l[i]]->outstub_number);
		
		}

		
		
		/*cout<<"END erasing "<<endl;
		
		cgroup.print_nodes(cout);
		neighs.print_nodes(cout);
		
		cout<<" ********************* END erasing ************** "<<endl;*/
		

		
		
		/*cout<<"************* erase_cgroup check"<<endl;
		
		cout<<"TTTTTTTTTTTTTooooooooooooBBBBBBBBBBBBBBeeeeeeeeeeeee"<<endl;
		
		prints(tobe);*/
		
		
		for(int i=0; i<tobe.size(); i++)
			erase_cgroup(tobe[i]);

		
		////check_all();
		
		
	}

}


void clomnet_local::erase_cgroup(int wnode, set<int> & erased) {
	
	
	
	/*cout<<"erasing erase_cgroup "<<wnode<<endl;
	
	
	
	cout<<"before erasing "<<endl;
	
	cgroup.print_nodes(cout);
	neighs.print_nodes(cout);
	
	cout<<" ********************* ... erasing ************** "<<endl;*/
	
	
	
	
	
	
	map<int, facts>::iterator itm= cgroup.lab_facts.find(wnode);
	if(itm!=cgroup.lab_facts.end()) {

		int kpin = itm->second.internal_indegree;
		int ktin = itm->second.indegree;
		int kpout = itm->second.internal_outdegree;
		int ktout = itm->second.outdegree;
		double mtlwin= itm->second.minus_log_total_wrin;
		double mtlwout= itm->second.minus_log_total_wrout;


		kin_cgroup-=  kpin + kpout;
		ktot_cgroup_out-=  ktout;
		ktot_cgroup_in-=  ktin;
		
		
		int kout_g_in= ktot_cgroup_in - kin_cgroup;
		int kout_g_out= ktot_cgroup_out - kin_cgroup;
		
		int tmin= oneM - ktot_cgroup_in;
		int tmout= oneM - ktot_cgroup_out;

		
		
		double fi= compute_global_fitness_randomized(kpin, kout_g_in, kpout, kout_g_out, tmin, tmout, ktin, ktout, mtlwin, mtlwout, neighs.size()+1, dim - cgroup.size() +1);
		neighs.edinsert(wnode, kpin, kpout, ktin, ktout, mtlwin, mtlwout, fi);
		
		
		cgroup.erase(wnode);
		erased.insert(wnode);
		
		
		
		deque<int> tobe;
		pair <int, pair<int, double> > OPA;
	
	
		
		
		for(int i=0; i<vertices[wnode]->inlinks->size(); i++) {
			
			
			
			OPA = vertices[wnode]->outlinks->posweightof(vertices[wnode]->inlinks->l[i]);

			
			if(cgroup.update_group(vertices[wnode]->inlinks->l[i], -vertices[wnode]->inlinks->w[i].first, -OPA.second.first, -vertices[wnode]->inlinks->w[i].second, -OPA.second.second, 
						dim - cgroup.size(), neighs.size(), kout_g_in, kout_g_out, tmin, tmout, vertices[vertices[wnode]->inlinks->l[i]]->instub_number, vertices[vertices[wnode]->inlinks->l[i]]->outstub_number,tobe)==false)
				neighs.update_neighs(vertices[wnode]->inlinks->l[i], -vertices[wnode]->inlinks->w[i].first, -OPA.second.first, -vertices[wnode]->inlinks->w[i].second, -OPA.second.second, 
						dim - cgroup.size(), kout_g_in, kout_g_out, tmin, tmout, vertices[vertices[wnode]->inlinks->l[i]]->instub_number, vertices[vertices[wnode]->inlinks->l[i]]->outstub_number);
		
		}
		
		
		
		for(int i=0; i<vertices[wnode]->outlinks->size(); i++) {
		
		
					
			OPA = vertices[wnode]->inlinks->posweightof(vertices[wnode]->outlinks->l[i]);
		
			if(OPA.first==-1) if(cgroup.update_group(vertices[wnode]->outlinks->l[i], 0, -vertices[wnode]->outlinks->w[i].first, 0, -vertices[wnode]->outlinks->w[i].second,
						dim - cgroup.size(), neighs.size(), kout_g_in, kout_g_out, tmin, tmout, vertices[vertices[wnode]->outlinks->l[i]]->instub_number, vertices[vertices[wnode]->outlinks->l[i]]->outstub_number,tobe)==false)
				neighs.update_neighs(vertices[wnode]->outlinks->l[i], 0, -vertices[wnode]->outlinks->w[i].first, 0, -vertices[wnode]->outlinks->w[i].second,
						dim - cgroup.size(), kout_g_in, kout_g_out, tmin, tmout, vertices[vertices[wnode]->outlinks->l[i]]->instub_number, vertices[vertices[wnode]->outlinks->l[i]]->outstub_number);
		
		}

		
		
		
		
	
		
		
		/*cout<<"END erasing "<<endl;
		
		cgroup.print_nodes(cout);
		neighs.print_nodes(cout);
		
		cout<<" ********************* END erasing ************** "<<endl;*/
		

		
		
		/*cout<<"************* erase_cgroup check"<<endl;
		
		cout<<"TTTTTTTTTTTTTooooooooooooBBBBBBBBBBBBBBeeeeeeeeeeeee"<<endl;
		
		prints(tobe);*/
		
		
		for(int i=0; i<tobe.size(); i++)
			erase_cgroup(tobe[i], erased);

		
		////check_all();
		
		
	}

}


void clomnet_local::insert_cgroup(int wnode) {

		
		
	// this function is to insert benode into cgroup  updating all the system, neighs - kin_cgroup - ktot_cgroup
	
	/*
	*/
	/*
	cout<<"INSERTING insert_cgroup "<<wnode<<endl;
	
	
	
	cout<<"before inserting "<<endl;
	
	cgroup.print_nodes(cout);
	neighs.print_nodes(cout);
	*/
	
	
	int kpin, ktin, kpout, ktout;
	double mtlwin, mtlwout;
	
	{
	
		map<int, facts>::iterator itm= neighs.lab_facts.find(wnode);
		if(itm!=neighs.lab_facts.end()) {
			
			
			kpin = itm->second.internal_indegree;
			ktin = itm->second.indegree;
			
			kpout = itm->second.internal_outdegree;
			ktout = itm->second.outdegree;

			mtlwin= itm->second.minus_log_total_wrin;
			mtlwout= itm->second.minus_log_total_wrout;
		
		
		
		} else {
			
			kpin=0;
			kpout=0;
			ktin=vertices[wnode]->instub_number;
			ktout=vertices[wnode]->outstub_number;
			mtlwin=0;
			mtlwout=0;
			
		}
		
	}	
	

	
	int kout_g_in= ktot_cgroup_in - kin_cgroup;
	int kout_g_out= ktot_cgroup_out - kin_cgroup;	
	int tmin= oneM - ktot_cgroup_in;
	int tmout= oneM - ktot_cgroup_out;
	
	//cout<<"node: "<<wnode<<endl;
	double fi= compute_global_fitness_randomized(kpin, kout_g_in, kpout, kout_g_out, tmin, tmout, ktin, ktout, mtlwin, mtlwout, neighs.size(), dim - cgroup.size());
	
	kin_cgroup+=  kpin + kpout;
	ktot_cgroup_in+= ktin;
	ktot_cgroup_out+= ktout;
	kout_g_in= ktot_cgroup_in - kin_cgroup;
	kout_g_out= ktot_cgroup_out - kin_cgroup;	
	tmin= oneM - ktot_cgroup_in;
	tmout= oneM - ktot_cgroup_out;
	
	
	
	//cout<<"INSERTING::: stuff to update "<<kin_cgroup<<" "<<ktot_cgroup_in<<" "<<tmin<<" "<<kout_g_in<<endl;
	
	cgroup.edinsert(wnode, kpin, kpout, ktin, ktout, mtlwin, mtlwout, fi);
	neighs.erase(wnode);
	
	
	
	deque<int> tobe;
	pair <int, pair<int, double> > OPA;
	
	for(int i=0; i<vertices[wnode]->inlinks->size(); i++) {
	
		OPA = vertices[wnode]->outlinks->posweightof(vertices[wnode]->inlinks->l[i]);
		
		if(cgroup.update_group(vertices[wnode]->inlinks->l[i], vertices[wnode]->inlinks->w[i].first, OPA.second.first, vertices[wnode]->inlinks->w[i].second,  OPA.second.second, 
					dim - cgroup.size(), neighs.size(), kout_g_in, kout_g_out, tmin, tmout, vertices[vertices[wnode]->inlinks->l[i]]->instub_number, vertices[vertices[wnode]->inlinks->l[i]]->outstub_number,tobe)==false)
			neighs.update_neighs(vertices[wnode]->inlinks->l[i], vertices[wnode]->inlinks->w[i].first, OPA.second.first, vertices[wnode]->inlinks->w[i].second, OPA.second.first, 
					dim - cgroup.size(), kout_g_in, kout_g_out, tmin, tmout, vertices[vertices[wnode]->inlinks->l[i]]->instub_number, vertices[vertices[wnode]->inlinks->l[i]]->outstub_number);
	
	}
	
	
	
	for(int i=0; i<vertices[wnode]->outlinks->size(); i++) {
		
		
		//cout<<"???? "<<vertices[wnode]->outlinks->l[i]<<" "<<vertices[wnode]->outlinks->w[i].first<<endl;
		
		OPA = vertices[wnode]->inlinks->posweightof(vertices[wnode]->outlinks->l[i]);
		
		if(OPA.first==-1) if(cgroup.update_group(vertices[wnode]->outlinks->l[i], 0, vertices[wnode]->outlinks->w[i].first, 0, vertices[wnode]->outlinks->w[i].second,
					dim - cgroup.size(), neighs.size(), kout_g_in, kout_g_out, tmin, tmout, vertices[vertices[wnode]->outlinks->l[i]]->instub_number, vertices[vertices[wnode]->outlinks->l[i]]->outstub_number,tobe)==false)
			neighs.update_neighs(vertices[wnode]->outlinks->l[i], 0, vertices[wnode]->outlinks->w[i].first, 0, vertices[wnode]->outlinks->w[i].second,
					dim - cgroup.size(), kout_g_in, kout_g_out, tmin, tmout, vertices[vertices[wnode]->outlinks->l[i]]->instub_number, vertices[vertices[wnode]->outlinks->l[i]]->outstub_number);
	
	}

		
	
	/*
	cout<<"END inserting "<<endl;
	
	cgroup.print_nodes(cout);
	neighs.print_nodes(cout);
	
	cout<<" ********************* END inserting ************** "<<endl;//*/
	//cherr();
	
	
	//cout<<"------pork----"<<endl;
	//check_all();
	//cherr();
	
}





int clomnet_local::try_to_assign_homeless_help(module_collection & module_coll, map<int, deque<int> > & to_check) {
	
	
	
	to_check.clear();
	module_coll.put_gaps();
	
	
	
	
		
	
	//if(paras.print_cbs)
		//cout<<"checking homeless nodes "<<endl;
	
	
	deque<int> homel;
	module_coll.homeless(homel);
	
	int before_procedure= homel.size();
	
	
	if(homel.size()==0)
		return before_procedure;
	
	/*cout<<"homel"<<endl;
	print_id(homel, cout);*/
	
	
	set<int> called;						// modules connected to homeless nodes
	map<int, set<int> > homel_module;		// maps the homeless node with the modules it's connected to
	
	
	
	for(int i=0; i<homel.size(); i++) {
	
		
		
		set<int> thish;
		for(int j=0; j<vertices[homel[i]]->inlinks->size(); j++) {
		
			int & neigh=vertices[homel[i]]->inlinks->l[j];
			
			for(set<int>:: iterator itk=module_coll.memberships[neigh].begin(); itk!=module_coll.memberships[neigh].end(); itk++) {
			
				called.insert(*itk);
				thish.insert(*itk);
			}
		}
		
		
		for(int j=0; j<vertices[homel[i]]->outlinks->size(); j++) {
		
			int & neigh=vertices[homel[i]]->outlinks->l[j];
			
			for(set<int>:: iterator itk=module_coll.memberships[neigh].begin(); itk!=module_coll.memberships[neigh].end(); itk++) {
			
				called.insert(*itk);
				thish.insert(*itk);
			}
		}

		
		
		
		if(thish.size()>0)
			homel_module[homel[i]]=thish;
	
	
	
	}
	
	
	
	
	map<int, int> module_kin;
	map<int, int> module_ktotin;
	map<int, int> module_ktotout;
	
	for(set<int>:: iterator its=called.begin(); its!=called.end(); its++ ) {
		
		module_kin[*its]=cast_int(kin_m(module_coll.modules[*its]));
		module_ktotin[*its]=cast_int(ktot_m(module_coll.modules[*its]).first);
		module_ktotout[*its]=cast_int(ktot_m(module_coll.modules[*its]).second);
	}
	
	
	
	for(map<int, set<int> > :: iterator itm= homel_module.begin();  itm!=homel_module.end(); itm++) {
	
		
		double cmin=1.1;
		int belongs_to=-1;
		
		//cout<<"homeless node: "<<id_of(itm->first)<<endl;
		
		for(set<int>:: iterator its= itm->second.begin(); its!=itm->second.end(); its++) {
			
			int kin_node_in=cast_int(vertices[itm->first]->kplus_m(module_coll.modules[*its]).first);
			int kin_node_out=cast_int(vertices[itm->first]->kplus_m(module_coll.modules[*its]).second);
			
						
			int kout_g_in= module_ktotin[*its] - module_kin[*its];
			int tmin= cast_int(oneM - module_ktotin[*its]);
			double kinw_in= vertices[itm->first]->kplus_w(module_coll.modules[*its]).first;
			
			int kout_g_out= module_ktotout[*its] - module_kin[*its];
			int tmout= cast_int(oneM - module_ktotout[*its]);
			double kinw_out= vertices[itm->first]->kplus_w(module_coll.modules[*its]).second;

			
			
			
			double rh1= compute_global_fitness_randomized_short(kin_node_in, kout_g_out, tmin, vertices[itm->first]->instub_number, kinw_in);
			double rh2= compute_global_fitness_randomized_short(kin_node_out, kout_g_in, tmout, vertices[itm->first]->outstub_number, kinw_out);
						
			double rh=min(rh1, rh2);
			
			//double cs=  1 - pow(1 - rh, dim - module_coll.modules[*its].size());
			//cout<<"rh: "<<rh<<" ..."<<endl;
			/*cout<<"module: "<<*its<<" kin: "<<module_kin[*its]<<"  ktot: "<<module_ktot[*its]<<" kin h "<<kin_node<<endl;
			print_ri(module_coll.modules[*its]);*/

			if(rh<cmin) {
				
				cmin=rh;
				belongs_to=*its;
			
			}
			
		}
		
		
		if(belongs_to!=-1) {
			
			if(to_check.find(belongs_to)==to_check.end()) {
				deque<int> void_d;
				to_check[belongs_to]=void_d;
			}
			
			to_check[belongs_to].push_back(itm->first);
		
		}
		
		
		
		//if(paras.print_cbs)
			//cout<<"homeless node: "<<id_of(itm->first)<<" belongs_to "<<belongs_to<<" cmin... "<<cmin<<endl;
			//cherr();
		
		
	}
	
	
	//if(paras.print_cbs)
		//cout<<"homeless node: "<<homel.size()<<" try_to_assign: "<<homel_module.size()<<" modules to check: "<<to_check.size()<<endl;

	
	



	return before_procedure;

}

















