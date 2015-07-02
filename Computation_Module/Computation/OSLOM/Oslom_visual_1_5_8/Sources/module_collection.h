#include "./standard_package/standard_include.cpp"















class module_collection {		

	// so all the labels refers to the index in modules (or submodules)

	
	public:
	
		module_collection(int d, double kave);
		module_collection(){};
		~module_collection(){};
		
		
		//bool truu(clomnet & a) {cout<<a.size()<<endl; };
		
		bool clear();
		//bool take_from(module_collection & Mcoll);
		
		bool insert(deque<int> & c, double bs, int & new_name);
		bool insert(deque<int> & c, double bs);
		//bool insert(int, deque<deque<int> > & submodules);
		int coverage();
		void fill_gaps();
		void put_gaps();

		int size() { return module_bs.size(); };
		
		void homeless(deque<int> & h);
		bool compute_inclusions(bool );
		bool pieces(const deque<int> & a, deque<int> & not_covered, deque<int> & smaller);
			
		void print(ostream & outt, deque<int> & netlabels, bool);
		void id_nodes(map<int, set<int> > & A, deque<int> & netlabels);
		
		
		
		void set_partition(deque<deque<int> > & A);
		void set_partition(deque<deque<int> > & A, deque<double> & b);
		
		void _set_(int dim, double kave);
		bool erase(int);
		bool erase_included(bool );
		bool erase_pieces(bool );

		
		void arrange_module_id();
		
		//***********************************************************
		//*****************			DATA		 ********************
		//***********************************************************
		
		
		deque<set<int> > memberships;				// this only refers to modules
		deque<deque<int> > modules;		
		map<int, double> module_bs;					// it mpas the module id into the b-score
		
		deque<deque<int> > ex_pieces;						// in this data-struct I have memory of the pieces I erased for further check (list of id's)
		deque<deque<int> > ex_almost_equal;				// in this data-struct I have memory of the groups (which was not pieces) that I erased for further check (list of modules)
		
		//map<int, deque<int> > module_submodules;	// this maps the module id in its submodule ids
		//***********************************************************
		//***********************************************************
		//***********************************************************

		
		
		//double average_degree;

	




	private:
		
		//deque<deque<int> > submodules;
		
		bool check_already(const deque<int> & c);
		bool big_module(const deque<int> & c, deque<int> & smaller);
		bool big_module(const deque<int> & c, set<int> & smaller, int module_id);
		bool pieces(int module_id, map<int, set<int> > & module_smaller_id);
		bool almost_equal(int module_id, deque<int> & smaller);
		bool all_pieces(map<int, set<int> > & module_smaller_id, set<int> & check_this, bool);
		
		bool erase_first_shell(map<int, deque<int> > & erase_net, bool);
		//double tot_mems;		// total number of memberships of all the node
		
		
		
		
	
};




module_collection::module_collection(int dim, double kave) {

	_set_(dim, kave);
}


void module_collection::_set_(int dim, double kave) {

	
	set<int> first;
	for(int i=0; i<dim; i++)
		memberships.push_back(first);
	
	
	//average_degree=kave;
	//tot_mems=0;
	


}



/*
bool module_collection::take_from(module_collection & Mcoll) {
	
	
	
	//cout<<memberships.size()<<" ... "<<Mcoll.memberships.size()<<endl;
	
	
	for(map<int, deque<int> >::iterator itm = Mcoll.module_submodules.begin(); itm!=Mcoll.module_submodules.end(); itm++) if(Mcoll.modules[itm->first].size() > 1) {		
		deque<int> & module_n= Mcoll.modules[itm->first];
		insert(module_n, Mcoll.module_bs[itm->first]);
	}


	//cout<<"doe"<<endl;
	

	return true;


}*/

bool module_collection::insert(deque<int> & c, double bs) {


	int new_name;
	return insert(c, bs, new_name);
	
}




bool module_collection::insert(deque<int> & c, double bs, int & new_name) {

	
	if(bs==0)
		bs=ran4() * 1e-10;
	
	
	
	sort(c.begin(), c.end());
	new_name=-1;
	
	if(check_already(c)==true) {
	
		
		//cout<<"??"<<endl;
		
		new_name=modules.size();
		for(int i=0; i<c.size(); i++)
			memberships[c[i]].insert(new_name);	
		
		
		modules.push_back(c);
		//deque<int> void_;
		//module_submodules.insert(make_pair(new_name, void_));
		
		module_bs[new_name]=bs;
		
		//tot_mems+=c.size();
		
		
		
		/*
		if(tot_mems>average_degree*memberships.size()) {		// there are so many overlaps that it would be worth to check the inclusions (to speed up the code)
			
			
			//cout<<"total memberships: "<<tot_mems<<" computing inclusions"<<endl;
			compute_inclusions();
			tot_mems=0;
		
		
		}//*/
				
		//cout<<"accepted ... erase "<<endl;
		//prints(c);
		
		
		return true;
		
	}
	
	//cout<<"not accepted: "<<endl;
	//prints(c);
	
	return false;
	
}



/*
bool module_collection::insert(int module_id, deque<deque<int> > & submodules_matrix) {
	

	if(module_bs.find(module_id)==module_bs.end())		// it only erases modules, not subs
		return false;

	deque<int> new_label_for_subs;
	
	
	for(int i=0; i<submodules_matrix.size(); i++) if(submodules_matrix[i].size()>0)  {
		
		sort(submodules_matrix[i].begin(), submodules_matrix[i].end());
		
		int new_name=submodules.size();
		submodules.push_back(submodules_matrix[i]);
		
		new_label_for_subs.push_back(new_name);
	}
	
	
	deque<int> & oldl= module_submodules[module_id];
	for(int i=0; i<new_label_for_subs.size(); i++)
		oldl.push_back(new_label_for_subs[i]);


	return true;
	
}*/







/*
int module_collection::insert(deque<int> & c, module_collection & mc) {



	deque<int> new_label_for_subs;
	for(map<int, deque<int> >::iterator itm=mc.module_submodules.begin(); itm!=mc.module_submodules.end(); itm++) {
		
		int new_name=submodules.size();
		submodules.push_back(mc.modules[itm->first]);
		
		new_label_for_subs.push_back(new_name);


	}
	
	
	int new_name=modules.size();
	for(int i=0; i<c.size(); i++)
		memberships[c[i]].push_back(new_name);	
	
	
	modules.push_back(c);	
	
	module_submodules.insert(make_pair(new_name, new_label_for_subs));




}
*/





//*
void module_collection::fill_gaps() {

	for(int i=0; i<memberships.size(); i++)
		if(memberships[i].size()==0) {
			
			deque<int> new_d;
			new_d.push_back(i);
			insert(new_d, 1.);
			
		
		}

}



void module_collection::put_gaps() {
		
	deque<int> to_erase;
	
	
	for(int i=0; i<modules.size(); i++) {
		
		if(modules[i].size()==1)
			to_erase.push_back(i);
	}
	
	
	for(int i=0; i<to_erase.size(); i++)
		erase(to_erase[i]);
	
}



//*/



void module_collection::homeless(deque<int> & h) {
	
	h.clear();
	
	for(int i=0; i<memberships.size(); i++)
		if(memberships[i].size()<1)
			h.push_back(i);
	
	for(int i=0; i<modules.size(); i++) {
		
		if(modules[i].size()==1)
			h.push_back(modules[i][0]);
	
	}
	
	
	sort(h.begin(), h.end());
		

}


int module_collection::coverage() {
	
	
	// this function returns the number of nodes which are covered by at least one module
	
	int cov=0;
	for(int i=0; i<memberships.size(); i++)
		if(memberships[i].size()>0)
			cov++;
	
	
	
	return cov;


}






bool module_collection::erase(int a) {
	
	// it erases module a 
	
	/*cout<<"i'm erasing this stuff... "<<a<<endl;
	if(a==1)
		cherr();*/
	
	if(module_bs.find(a)==module_bs.end())		// it only erases not empty modules
		return false;
	
	deque<int> & nodes_a = modules[a];
	
	for(int i=0; i<nodes_a.size(); i++)
		memberships[nodes_a[i]].erase(a); 
	
			
	/*deque<int> & submodule_labels= module_submodules[a];
	for(int j=0; j<submodule_labels.size(); j++)
		submodules[submodule_labels[j]].clear();*/

	
	//module_submodules.erase(a);
	modules[a].clear();
	module_bs.erase(a);
	
	
	
	return true;

}






bool module_collection::clear() {
	

	int D=memberships.size();
	
	modules.clear();
	memberships.clear();
	module_bs.clear();
	ex_pieces.clear();
	ex_almost_equal.clear();

	
	
	
	
	_set_(D, 0);
		
	return true;

}









void module_collection::id_nodes(map<int, set<int> > & A, deque<int> & netlabels) {

	// A maps the module index into the node labels
	
	
	A.clear();
	for(map<int, double>::iterator itm = module_bs.begin(); itm!=module_bs.end(); itm++) {
		
		deque<int> & module_nodes= modules[itm->first];
		
		set<int> labseq;
		for(int i=0; i<module_nodes.size(); i++) {
			labseq.insert(netlabels[module_nodes[i]]);
		}
		
		
		
		A[itm->first]=labseq;
	
	
	
	}
}




void module_collection::print(ostream & outt, deque<int> & netlabels, bool not_homeless) {
	
	
	// maybe I should do a check to see if some subgroups are too similar
	
	// netlabels[i] is the id of node i
	
	
	//printm(modules);
	//printm(memberships);
	
	
		
	deque<int> homel;
	homeless(homel);
	
	
	int nmod=0;
	for(map<int, double >::iterator itm = module_bs.begin(); itm!=module_bs.end(); itm++) if(not_homeless==false || modules[itm->first].size() > 1) {
		
		
		nmod++;
		
		
		deque<int> & module_nodes= modules[itm->first];
		outt<<"#module "<<itm->first<<" size: "<<modules[itm->first].size()<<" bs: "<<module_bs[itm->first]<<endl;
		
		deque<int> labseq;
		for(int i=0; i<module_nodes.size(); i++) {
			labseq.push_back(netlabels[module_nodes[i]]);
		}
		
		sort(labseq.begin(), labseq.end());
		
		for(int i=0; i<labseq.size(); i++) {
			outt<<labseq[i]<<" ";
		}
		
		//outt<<endl<<"\t# ";
		/*
		for(int i=0; i<labseq.size(); i++) {
			outt<<RI[labseq[i]]<<" ";
		}
		//*/

		outt<<endl;
		
		/*
		if(itm->second.size()>0) {
			
			deque<int> & submodule_labels= itm->second;
			
			for(int j=0; j<submodule_labels.size(); j++) {
				
				deque<int> & ds= submodules[submodule_labels[j]];
				
				deque<int> subseq;
				for(int k=0; k<ds.size(); k++)
					subseq.push_back(netlabels[ds[k]]);
				
				sort(subseq.begin(), subseq.end());
				
				
				if(includes(labseq.begin(), labseq.end(), subseq.begin(), subseq.end())) {
					outt<<"#\t[";
					for(int k=0; k<subseq.size(); k++)
						outt<<subseq[k]<<" ";
					outt<<"\t]"<<endl;
				}
			}
		
		
		}*/
		

	}
	
	
	/*
	outt<<"#homeless: size: "<<homel.size()<<" nodes: ";
	deque<int> hms;
	for(int k=0; k<homel.size(); k++)
		hms.push_back(netlabels[homel[k]]);
				
	sort(hms.begin(), hms.end());
				
				
	for(int k=0; k<hms.size(); k++)
		outt<<hms[k]<<" ";
	outt<<endl;*/
	
	

}



void module_collection::set_partition(deque<deque<int> > & A) {


	A.clear();

	for(map<int, double >::iterator itm = module_bs.begin(); itm!=module_bs.end(); itm++) if(modules[itm->first].size()>1)
		A.push_back(modules[itm->first]);
		





}


void module_collection::set_partition(deque<deque<int> > & A, deque<double> & b) {


	A.clear();
	b.clear();
	
	for(map<int, double>::iterator itm = module_bs.begin(); itm!=module_bs.end(); itm++) if(modules[itm->first].size()>1){
		A.push_back(modules[itm->first]);
		b.push_back(module_bs[itm->first]);
	}





}











bool module_collection::big_module(const deque<int> & c, deque<int> & smaller) {

	
	
	
	smaller.clear();
	
	map<int, int> com_ol;		// it maps the index of the modules into the overlap (overlap=numeber of overlapping nodes)
	
	for(int i=0; i<c.size(); i++) {
		
		for(set<int>:: iterator itj=memberships[c[i]].begin(); itj!=memberships[c[i]].end(); itj++)
			int_histogram(*itj, com_ol);
	
	}
	
	
	
	for(map<int, int>::iterator itm=com_ol.begin(); itm!=com_ol.end(); itm++) if(modules[itm->first].size()<= c.size())
		if(double(itm->second) / modules[itm->first].size() >= paras.coverage_inclusion_module_collection)
			smaller.push_back(itm->first);
	
	if(smaller.size()==0)
		return false;
	
	
	
	return true;


}



bool module_collection::big_module(const deque<int> & c, set<int> & smaller, int module_id) {
	
	
	
	
	// this function puts in smaller all the module ids contained in c, which has module id = module_id 
	
	
	smaller.clear();
	
	map<int, int> com_ol;		// it maps the index of the modules into the overlap (overlap=numeber of overlapping nodes)
	
	for(int i=0; i<c.size(); i++) {
		
		for(set<int>:: iterator itj=memberships[c[i]].begin(); itj!=memberships[c[i]].end(); itj++)
			int_histogram(*itj, com_ol);
	
	}
	
	
	for(map<int, int>::iterator itm=com_ol.begin(); itm!=com_ol.end(); itm++) if(modules[itm->first].size()<= c.size()  && itm->first!=module_id)
		if(double(itm->second) / modules[itm->first].size() >= paras.coverage_inclusion_module_collection)
			smaller.insert(itm->first);
	
	if(smaller.size()==0)
		return false;
	
	
	
	return true;


}




bool module_collection::check_already(const deque<int> & c) {

	// returns false if the module is already present
	
	
	
	map<int, int> com_ol;		// it maps the index of the modules into the overlap (overlap=numeber of overlapping nodes)
	
	for(int i=0; i<c.size(); i++) {
		
		for(set<int>:: iterator itj=memberships[c[i]].begin(); itj!=memberships[c[i]].end(); itj++)
			int_histogram(*itj, com_ol);
		
	
	}
	
	
	
	for(map<int, int>::iterator itm=com_ol.begin(); itm!=com_ol.end(); itm++) {
		
		if(itm->second==c.size() && itm->second==modules[itm->first].size())
			return false;
			
	}
	
		
	return true;



}









bool module_collection::pieces(const deque<int> & a, deque<int> & not_covered, deque<int> & smaller) {

	
	// return true if a is the union of more than two different pieces
	
	not_covered.clear();

	smaller.clear();
	
	if(big_module(a, smaller)==false)
		return false;
		
	set<int> aa;
	
	//cout<<"smaller.size(): "<<smaller.size()<<endl;
	for(int i=0; i<smaller.size(); i++) {
		
		deque<int> & bb = modules[smaller[i]];
		
		//cout<<"bb size: "<<bb.size()<<endl;
		for(int j=0; j<bb.size(); j++)
			aa.insert(bb[j]);	
	
	}
	
	
	//if(paras.print_cbs)
		//cout<<"pieces: "<<aa.size()<<" vs "<<paras.coverage_percentage_minimality * a.size()<<endl;
	
	
	if(aa.size() < paras.coverage_percentage_minimality * a.size())
		return false;
	
	
	
	sort(smaller.begin(), smaller.end());
	for(int i=0; i<smaller.size(); i++) {
		
	
		
		
		deque<int> smaller2;
		big_module(modules[smaller[i]], smaller2);
		sort(smaller2.begin(), smaller2.end());
		deque<int> cc;
		set_intersection(smaller.begin(), smaller.end(), smaller2.begin(), smaller2.end(), back_inserter(cc));
		
		//cout<<"intersection: "<<cc.size()<<" vs "<<smaller.size()<<endl;
		
		if(cc.size()==smaller.size())
			return false;				// return false if some submodule contains all the other submodules (in this case it is not a union of different groups)
	
	
	}
	
	
	for(int i=0; i<a.size(); i++) if(aa.find(a[i])==aa.end())
		not_covered.push_back(a[i]);
	
	
	return true;





}







inline bool module_collection::all_pieces(map<int, set<int> > & module_smaller_id, set<int> & check_this, bool ex_flag) {



	//cout<<"all pieces before"<<endl;
	//printm(modules);
	



	set<int> erase_modules;
	
	//cout<<"check_this.size() "<<check_this.size()<<endl;
	
	for(set<int>::iterator its= check_this.begin();  its!=check_this.end(); its++) {
		
		if(pieces(*its, module_smaller_id))		// if *its is the union of different modules
			erase_modules.insert(*its);
		
		//cout<<*its<<endl;
		
	}
	
	
	//cout<<"erase modules pieces--> "<<erase_modules.size()<<endl;
		
	bool again=false;
	set<int> to_check_again; 

	for(set<int>::iterator its=erase_modules.begin(); its!=erase_modules.end(); its++) {
	
		if(module_bs.find(*its)!=module_bs.end()) {
		
			deque<int> & c= modules[*its];

			for(int i=0; i<c.size(); i++) for(set<int>:: iterator itj=memberships[c[i]].begin(); itj!=memberships[c[i]].end(); itj++)
				to_check_again.insert(*itj);
				
			
			/*set<int> dhelp(module_smaller_id[*its]);
			
			//************************
			cout<<"these are the pieces"<<endl;
			for(set<int>:: iterator itts= dhelp.begin(); itts!= dhelp.end(); itts++)
				prints(modules[*itts]);
			
			//****************************
			//*/
			
			
			if(ex_flag)
				ex_pieces.push_back(modules[*its]);
			
			erase(*its);
			
			module_smaller_id.erase(*its);
			again=true;
		}
	}
	
	
	check_this.clear();
	
	for(set<int>::iterator its=to_check_again.begin(); its!=to_check_again.end(); its++) if(module_bs.find(*its)!=module_bs.end()) {
		
		set<int> left;
		big_module(modules[*its], left, *its);
		module_smaller_id[*its]=left;
		check_this.insert(*its);
	
	}
	
	
	//cout<<"modules as it is"<<endl;
	//printm(modules);
	
	return true;
	

}


bool module_collection::erase_pieces(bool ex_flag) {
	
	map<int, set<int> > module_smaller_id;		// this maps the module_id into the modules which it "contains"
	
	//time_t t1=time(NULL);

	for(map<int, double>::iterator itm= module_bs.begin();  itm!=module_bs.end(); itm++) {
		
		set<int> included;
		big_module(modules[itm->first], included, itm->first);
		module_smaller_id.insert(make_pair(itm->first, included));
		
		
	}

	//time_t t2=time(NULL);
	//cout<<"all smallers "<<t2-t1<<" secs"<<endl;
	
	set<int> initial;
	for(map<int, set<int> >::iterator itm= module_smaller_id.begin();  itm!=module_smaller_id.end(); itm++) {
		if(itm->second.size()>0)
			initial.insert(itm->first);
	}
	
	
	
	while(true) {
	
		all_pieces(module_smaller_id, initial, ex_flag);		// this function is to erase all the modules which are unions of different modules
		
		if(initial.size()==0)
			break;
	
	}
	
	
}




bool module_collection::erase_included(bool ex_flag) {




	map<int, deque<int> > erase_net;
	for(map<int, double >::iterator itm= module_bs.begin();  itm!=module_bs.end(); itm++) {
	

		deque<int> smaller;
		almost_equal(itm->first, smaller);
		erase_net[itm->first]=smaller;
	
	
	}
	
	
	while(true) {
		
		if(erase_first_shell(erase_net, ex_flag)==false)
			break;

	}


}



bool module_collection::compute_inclusions(bool ex_flag) {


	
	erase_pieces(ex_flag);
	if(ex_flag==false)
		erase_included(ex_flag);
	
	
	return true;

}



bool module_collection::erase_first_shell(map<int, deque<int> > & erase_net, bool ex_flag) {
	
	

	
	
	
	
	bool again=false;
	set<int> roots;
	
	for(map<int, double >::iterator itm= module_bs.begin();  itm!=module_bs.end(); itm++)
		roots.insert(itm->first);

	
	
	
	for(map<int, deque<int> >::iterator itm= erase_net.begin();  itm!=erase_net.end(); itm++) {
	

		deque<int> & smaller=itm->second;
		
		for(int i=0; i<smaller.size(); i++)
			roots.erase(smaller[i]);
	
	}
	
	
	//cout<<"roots:"<<endl;
	//prints(roots);
	


	for(set<int>::iterator its=roots.begin(); its!=roots.end(); its++) {
		
		deque<int> & smaller=erase_net[*its];
		//cout<<"root: "<<*its<<endl;
		//cout<<"smaller: ";
		//prints(smaller);
		
		
		for(int i=0; i<smaller.size(); i++) {
		
			if(module_bs.find(smaller[i])!=module_bs.end()) {
				
				
				if(ex_flag)
					ex_almost_equal.push_back(modules[smaller[i]]);
				
				erase(smaller[i]);
				erase_net.erase(smaller[i]);
				again=true;
			
			}
		}
		
	}

	
	return again;

}



bool module_collection::almost_equal(int module_id, deque<int> & smaller) {




	// c is the module you want to know about
	// smaller is set to contain the module ids contained by module_id
	
	
	
	deque<int> & c= modules[module_id];
	
	map<int, int> com_ol;		// it maps the index of the modules into the overlap (overlap=numeber of overlapping nodes)
	
	for(int i=0; i<c.size(); i++) {
		
		for(set<int>:: iterator itj=memberships[c[i]].begin(); itj!=memberships[c[i]].end(); itj++)
			int_histogram(*itj, com_ol);
	
	}
	
	
			
	for(map<int, int>::iterator itm=com_ol.begin(); itm!=com_ol.end(); itm++) if(itm->first!=module_id) {
		
		
		int other_size= modules[itm->first].size();
		int min_s=min(int(c.size()), other_size);
			
			
		//if(paras.print_cbs)
			//cout<<"overlap percentage between: "<<itm->first<<" and "<<module_id<<" "<<double(itm->second) / min_s<<endl;
			
			
		if(double(itm->second) / min_s >= paras.coverage_inclusion_module_collection) {
				

				
			int em=itm->first;
			if(c.size() > other_size)
				smaller.push_back(em);
			
			
			else if(c.size()==other_size && module_bs[module_id]<module_bs[itm->first])
				smaller.push_back(em);

			
			//cout<<"module_id itm->second "<<module_id<<" "<<itm->first<<" "<<c.size()<<" "<<modules[itm->first].size()<<" "<<module_bs[module_id];
			//cout<<" "<<module_bs[itm->first]<<" em: "<<em<<endl;
				
			
		}
		
			
	}
	
	
	return true;

}


bool module_collection::pieces(int module_id, map<int, set<int> > & module_smaller_id) {

	
	// return true if a is the union of more than two different pieces
	
	const deque<int> & a= modules[module_id];

	const set<int> & smaller=module_smaller_id[module_id];
	if(smaller.size()==0)
		return false;
		
	set<int> aa;
	
	/*cout<<"module size: "<<a.size()<<endl;
	prints(a);
	cout<<"-----------------------"<<endl;*/
	
	
	//cout<<"smaller.size(): "<<smaller.size()<<endl;
	for(set<int>:: iterator iti=smaller.begin(); iti!=smaller.end(); iti++) {
		
		deque<int> & bb = modules[*iti];
		
		/*cout<<"bb size: "<<bb.size()<<endl;
		prints(bb);*/
		for(int j=0; j<bb.size(); j++)
			aa.insert(bb[j]);	
	
	}
	
	
	
	//cout<<"pieces: "<<aa.size()<<" vs "<<paras.coverage_percentage_minimality * a.size()<<endl;
	
	
	if(aa.size() < paras.coverage_percentage_minimality * a.size())
		return false;
	
	
	
	for(set<int>:: iterator iti=smaller.begin(); iti!=smaller.end(); iti++) {
		
		const set<int> & smaller2 =module_smaller_id[*iti];
		
		deque<int> cc;
		set_intersection(smaller.begin(), smaller.end(), smaller2.begin(), smaller2.end(), back_inserter(cc));
		
		//cout<<"intersection: "<<cc.size()<<" vs "<<smaller.size()<<endl;
		
		if(cc.size()+1==smaller.size())
			return false;
	
	
	}
	
	

	

	
	
	return true;





}



void module_collection::arrange_module_id() {
	
	
	
	/*cout<<"before"<<endl;
	printm(memberships);
	printm(modules);
	prints(module_bs);*/
		
	map<int, int> from_old_index_to_new;


	{
		
		deque<deque<int> > modules2;
		map<int, double> module_bs2;
	
	
		for(map<int, double> :: iterator itm= module_bs.begin(); itm!=module_bs.end(); itm++) {
			from_old_index_to_new.insert(make_pair(itm->first, from_old_index_to_new.size()));
			modules2.push_back(modules[itm->first]);
			module_bs2[from_old_index_to_new.size()-1]=itm->second;
		}
		
		
		modules=modules2;
		module_bs=module_bs2;
	}
	
	/*cout<<"from_old_index_to_new"<<endl;
	prints(from_old_index_to_new);*/
	
	
	for(int i=0; i<memberships.size(); i++) {
		
		set<int> first;
		for(set<int>:: iterator its= memberships[i].begin(); its!=memberships[i].end(); its++)
			first.insert(from_old_index_to_new[*its]);
		
		memberships[i]=first;
	}
	
	/*cout<<"after"<<endl;
	printm(memberships);
	printm(modules);
	prints(module_bs);*/

}





