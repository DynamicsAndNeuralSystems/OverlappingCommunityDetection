
// external_xs can be implemented better

#include "./standard_package/standard_include.cpp"
//#include "border_score.h"



//#include <gsl/gsl_cdf.h>
//#include <gsl/gsl_randist.h>
//#include <gsl/gsl_sf_gamma.h>
//#define bin_vs_hyper 0.2







//#include <gsl/gsl_sf_erf.h>
# define sqrt_two 1.41421356237
# define num_up_to 5









double right_error_function(double x) {
	
	
	return 0.5*erfc(x/sqrt_two);

}


double log_together(double minus_log_total, int number) {
	
	
	if(number==0)
		return 0.5;
	
	// this function will be better! (depending on number...)
	
	
	if(number<11) {
	
		double fa=1;
		double zum=1;
		
		for(int i=1; i<number; i++) {
		
			fa *= minus_log_total / i;
			zum+=fa;
		
		}
	
		
		return zum*exp(-minus_log_total);
	
	
	
	}
	
	
	
	double mu=number;
	return right_error_function((minus_log_total-mu)/sqrt(mu));

}




inline double right_tail_order_exponential(int N, int pos, double x) {


	// N is the total number of variables, pos is from 1 to N. 1 is the smallest.
	// the routine computes the probality c_pos=  p(X_pos >= -x)
	
	
	double y= 1 +x;
	
	if(y>1)
		return 0;
	
	//return gsl_cdf_binomial_P(pos-1, y, N);
	//cout<<pos<<" pos-N "<<N<<endl;
	//cout<<"pos: "<<pos<<" N: "<<N<<" log p: "<<log(y)<<" "<< LOG_TABLE.cum_binomial_left(pos, N, y)<<endl;
	
	return LOG_TABLE.cum_binomial_left(pos, N, y);

}





double fitted_exponent(int N) {



	double l=log(double(N));

	
	if(N>100)
		return 4.2*l -8.5;

	if(N>30)
		return 3.5*l -5.5;

	if(N>7)
		return 2.5*l -2;

	if(N>1)
		return 1.3*l +0.1;

	
	return 1;



}



inline double pron_min_exp(int N, double xi) {


	// this should return the probability that the minimum of the quantiles p(c_min<=xi)
	
	//cout<<"-> "<<fitted_exponent(N)<<endl;
	
	return 1 - exp(-fitted_exponent(N)*xi);


}


/*
double compute_r_hyper_not_calibrated(int kin_node, int kout_g, int tm, int degree_node) {

	double r=1;
	if(kin_node>0)
		r=gsl_cdf_hypergeometric_Q(kin_node-1, kout_g, tm - kout_g, degree_node);
	
	//cout<<"hyp: "<<kin_node-1<<" "<<kout_g<<" "<<tm-kout_g<<" "<<degree_node<<endl;
	if(r<0)
		r=0;
	
	return r;

}
*/

/*
inline double compute_r_hyper(int kin_node, int kout_g, int tm, int degree_node) {
	

	//cout<<kin_node<<" "<<kout_g<<" "<<tm<<" "<<degree_node<<endl;
	return LOG_TABLE.right_cumulative_function(degree_node, kout_g, tm, kin_node);

}*/



inline double hyper_table(int kin_node, int kout_g, int tm, int degree_node) {
	
	return LOG_TABLE.hyper(kin_node, kout_g, tm, degree_node);
	

}




/*double compute_global_fitness(int kin_node, int kout_g, int tm, int degree_node, double minus_log_total, int number_of_neighs, int Nstar) {
	
	// so, it must be completely true that kin_node is referred to the node and not to the module
	
	
	double topologic=  LOG_TABLE.right_cumulative_function(degree_node, kout_g, tm, kin_node);
	
	double mu=  (number_of_neighs +1.)  / (Nstar+1.);
	double weight_part= log_together(minus_log_total, kin_node)  *  mu;

	
	
	return log_together(-log(topologic) - log(weight_part), 2);


}
*/

inline double r_help(int degree_node_out, int kout_g_out, int tmout, int kin_node_out, double minus_log_total_out, double mu) {

	//cout<<"pok"<<endl;
	

	double b2out=LOG_TABLE.right_cumulative_function(degree_node_out, kout_g_out, tmout, kin_node_out+1);
	
	
	//cout<<"pok2 "<<b2out<<" "<<(hyper_table(kin_node_out, kout_g_out, tmout,  degree_node_out))<<endl;
	
	double topologic_out=  b2out + ran4() * (hyper_table(kin_node_out, kout_g_out, tmout,  degree_node_out));

	topologic_out /= mu;
	
	//cout<<mu<<endl;
	if(topologic_out>1)
		topologic_out=1;


	//cout<<topologic_out<<"<-< indeg: "<<kin_node_out<<" / "<<degree_node_out<<" tm: "<<tmout<<endl;
	
	
	
	double weight_part_out= log_together(minus_log_total_out, kin_node_out);
	//cout<<"weighted_part "<<weight_part_out<<" "<<log_together(-log(topologic_out) - log(weight_part_out), 2)<<" ...  "<<minus_log_total_out<<" "<<kin_node_out<<endl;
	//cherr();
	
	if(topologic_out<=0 || weight_part_out<=0 )
		return -1;

	//cout<<degree_node_out<<" "<<kin_node_out<<" "<<tmout<<" "<<kout_g_out<<" "<<minus_log_total_out<<" degree-kp-tm-kout_g-mtlw   "<<-log(topologic_out) - log(weight_part_out)<<endl;

	return -log(topologic_out) - log(weight_part_out);

}


double compute_global_fitness_randomized(int kin_node_in, int kout_g_in, int kin_node_out, int kout_g_out, int tmin, int tmout, int degree_node_in, 
										int degree_node_out, double minus_log_total_in, double minus_log_total_out, int number_of_neighs, int Nstar) {

	// so, it must be completely true that kin_node is referred to the node and not to the module
	
	

	/*cout<<"compute_global_fitness_randomized"<<endl;
	cout<<"degree: "<<degree_node_in<<" "<<degree_node_out<<" tmin "<<tmin<<" "<<tmout;
	cout<<" kout_g "<<kout_g_in<<" "<< kout_g_out <<" kin_node "<<kin_node_in<<" "<< kin_node_out <<" Nstar "<<Nstar<<" nn "<<number_of_neighs<<" minus_log_total "<<minus_log_total_in<<" "<< minus_log_total_out <<endl;
	*/
	
	
	// degree_node_in should be here *************
	
	if(paras.weighted==false) {
		
		if(degree_node_in==0) {
			
			double b2out=LOG_TABLE.right_cumulative_function(degree_node_out, kout_g_in, tmout, kin_node_out+1);
			double topologic_out=  b2out + ran4() * (hyper_table(kin_node_out, kout_g_in, tmout,  degree_node_out));
			if(topologic_out<=0)
				return 1e-200;
			
			return topologic_out;
			
		}
		
		
		if(degree_node_out==0) {
			
			double b2out=LOG_TABLE.right_cumulative_function(degree_node_in, kout_g_out, tmin, kin_node_in+1);
			double topologic_out=  b2out + ran4() * (hyper_table(kin_node_in, kout_g_out, tmin,  degree_node_in));
			if(topologic_out<=0)
				return 1e-200;
			
			return topologic_out;
			
		}

		
		
		double b2out=LOG_TABLE.right_cumulative_function(degree_node_out, kout_g_in, tmout, kin_node_out+1);
		double topologic_out=  b2out + ran4() * (hyper_table(kin_node_out, kout_g_in, tmout,  degree_node_out));
		if(topologic_out<=0)
			return 1e-200;
		
		
		double b2in=LOG_TABLE.right_cumulative_function(degree_node_in, kout_g_out, tmin, kin_node_in+1);
		double topologic_in=  b2out + ran4() * (hyper_table(kin_node_in, kout_g_out, tmin,  degree_node_in));
		if(topologic_in<=0)
			return 1e-200;
		
		
		return log_together( -log(topologic_in) - log(topologic_out), 2);

	
		
	
	
	}
	
	
	
	
	
	
	
	double mu=  (number_of_neighs +1.)  / (Nstar+1.);

	/*cout<<"*********************************"<<endl;
	cout<<degree_node_out<<" "<<kout_g_in<<" "<<tmout<<" "<<kin_node_out<<" "<<endl;
	cout<<degree_node_in<<" "<<kout_g_out<<" "<<tmin<<" "<<kin_node_in<<" "<<endl;*/

	
	if(degree_node_in==0) {
		
		
		
		double x= r_help(degree_node_out, kout_g_in, tmout, kin_node_out, minus_log_total_out, mu);
		if(x<0)
			return 1e-200;
		
		return log_together(x, 2);
		
	}
	
	if(degree_node_out==0) {
		
		double x= r_help(degree_node_in, kout_g_out, tmin, kin_node_in, minus_log_total_in, mu);
		if(x<0)
			return 1e-200;
		
		return log_together(x, 2);
		
	}

	double xin= r_help(degree_node_in, kout_g_out, tmin, kin_node_in, minus_log_total_in, mu);
	double xout= r_help(degree_node_out, kout_g_in, tmout, kin_node_out, minus_log_total_out, mu);

	
	if(xin<0 || xout<0)
		return 1e-200;
	
	
	return log_together(xin + xout, 4);

	
	
	//double tog= log_together(-log(topologic_in) - log(weight_part_in) -log(topologic_out) - log(weight_part_out), 4);
	//return min(weight_part , tog);
	//return min(min(topologic, weight_part) , tog);
	//return tog;




}


/*

double compute_global_fitness_randomized(int kin_node_in, int kout_g_in, int kin_node_out, int kout_g_out, int tmin, int tmout, int degree_node_in, 
										int degree_node_out, double minus_log_total_in, double minus_log_total_out, int number_of_neighs, int Nstar, bool dummy) {

	// so, it must be completely true that kin_node is referred to the node and not to the module
	
	
	cout<<"compute_global_fitness_randomized"<<endl;
	cout<<"degree: "<<degree_node_in<<" "<<degree_node_out<<" tmin "<<tmin<<" "<<tmout;
	cout<<" kout_g "<<kout_g_in<<" "<< kout_g_out <<" kin_node "<<kin_node_in<<" "<< kin_node_out <<" Nstar "<<Nstar<<" nn "<<number_of_neighs<<" minus_log_total "<<minus_log_total_in<<" "<< minus_log_total_out <<endl;
	
	double mu=  (number_of_neighs +1.)  / (Nstar+1.);

	
	
	if(kin_node_in==0) {
		
		double x= r_help(degree_node_out, kout_g_out, tmout, kin_node_out, minus_log_total_out, mu);
		if(x<0)
			return 1e-200;
		
		return log_together(x, 2);
		
	}
	
	if(kin_node_out==0) {
		
		double x= r_help(degree_node_in, kout_g_in, tmin, kin_node_in, minus_log_total_in, mu);
		if(x<0)
			return 1e-200;
		
		return log_together(x, 2);
		
	}

	double xout= r_help(degree_node_out, kout_g_out, tmout, kin_node_out, minus_log_total_out, mu);
	double xin= r_help(degree_node_in, kout_g_in, tmin, kin_node_in, minus_log_total_in, mu);
	
	if(xin<0 || xout<0)
		return 1e-200;
	
	
	return log_together(xin + xout, 4);

	
	
	//double tog= log_together(-log(topologic_in) - log(weight_part_in) -log(topologic_out) - log(weight_part_out), 4);
	//return min(weight_part , tog);
	//return min(min(topologic, weight_part) , tog);
	//return tog;




}

//*/
/*

double compute_global_fitness_randomized(int kin_node, int kout_g, int tm, int degree_node, double minus_log_total) {
	
	// so, it must be completely true that kin_node is referred to the node and not to the module
	
	
	//cout<<"compute_global_fitness_randomized"<<endl;
	//cout<<"degree: "<<degree_node<<" tm "<<tm<<" kout_g "<<kout_g<<" kin_node "<<kin_node<<" Nstar "<<0<<" nn "<<0<<" minus_log_total "<<minus_log_total<<endl;
	
	
	double b2=LOG_TABLE.right_cumulative_function(degree_node, kout_g, tm, kin_node+1);
	//cout<<"b2: "<<b2<<endl;
	
	double topologic=  b2 + ran4() * (hyper_table(kin_node, kout_g, tm,  degree_node));

	
	
	//cout<<"topologic "<<topologic<<" mu "<<1.<<" t/mu "<<topologic/1<<" --- --- ---- "<<b2<<" "<<(hyper_table(kin_node, kout_g, tm,  degree_node));
	
	
	
	
	double weight_part= log_together(minus_log_total, kin_node);
	
	//cout<<" weight_part "<<log_together(minus_log_total, kin_node);
	//cout<<" total f:::  "<<log_together(-log(topologic) - log(weight_part), 2)<<" kin: "<<kin_node<<" / "<<degree_node<<" mtlw "<<minus_log_total<<endl;
	//return weight_part;
	if(topologic<=0 || weight_part<=0)
		return 1e-200;
	
	return min(topologic, weight_part);




}
*/




double compute_global_fitness_randomized_short(int kin_node, int kout_g, int tm, int degree_node, double minus_log_total) {
	
	// so, it must be completely true that kin_node is referred to the node and not to the module
	
	
	//cout<<"compute_global_fitness_randomized"<<endl;
	//cout<<"degree: "<<degree_node<<" tm "<<tm<<" kout_g "<<kout_g<<" kin_node "<<kin_node<<" Nstar "<<0<<" nn "<<0<<" minus_log_total "<<minus_log_total<<endl;
	
	
	double b2=LOG_TABLE.right_cumulative_function(degree_node, kout_g, tm, kin_node+1);
	//cout<<"b2: "<<b2<<endl;
	
	double topologic=  b2 + ran4() * (hyper_table(kin_node, kout_g, tm,  degree_node));

	
	
	//cout<<"topologic "<<topologic<<" mu "<<1.<<" t/mu "<<topologic/1<<" --- --- ---- "<<b2<<" "<<(hyper_table(kin_node, kout_g, tm,  degree_node));
	
	
	
	
	double weight_part=2;
	
	if(paras.weighted)
		weight_part=log_together(minus_log_total, kin_node);
	
	//cout<<" weight_part "<<log_together(minus_log_total, kin_node);
	//cout<<" total f:::  "<<log_together(-log(topologic) - log(weight_part), 2)<<" kin: "<<kin_node<<" / "<<degree_node<<" mtlw "<<minus_log_total<<endl;
	//return weight_part;
	if(topologic<=0 || weight_part<=0)
		return 1e-200;
	
	return min(topologic, weight_part);




}







class facts {




	public:
	
		
		facts(int in_id, int in_od, double in_iw, double in_ow, multimap<double, int>::iterator c, int ide, int ode)  {
			
			
			fitness_iterator=c; 
			
			internal_indegree=in_id; 
			minus_log_total_wrin=in_iw; 
			
			internal_outdegree=in_od;
			minus_log_total_wrout=in_ow; 

			indegree=ide;
			outdegree=ode;
		
		
		
		};
		
		
		~facts(){};
		
		
		int indegree;
		int outdegree;
		
		int internal_indegree;						
		int internal_outdegree;
		
		
		double minus_log_total_wrin;								// wr is the right part of the exponential for the weights, this is the sum over the internal stubs of that
		double minus_log_total_wrout;								// wr is the right part of the exponential for the weights, this is the sum over the internal stubs of that
		multimap<double, int>::iterator fitness_iterator;
		
		
};




class weighted_tabdeg {


	public:
				
		weighted_tabdeg(){};
		~weighted_tabdeg(){};
		
		void _set_(weighted_tabdeg &);
		
		
		void clear();
		void edinsert(int a, int kpin, int kpout, int ktin, int ktout, double mtlwin, double mtlwout, double fit);		// this function inserts element a (or edit it if it was already inserted)
		bool erase(int a);
		void set_deque(deque<int> & );


		int size() {return lab_facts.size();};
		
		
		void print_nodes(ostream &);
		bool update_group(int a, int delta_deg_in, int delta_deg_out, double delta_win, double delta_wout, int nstar, int nn, int kout_g_in, int kout_g_out, int tmin, int tmout, int ktin, int ktout, deque<int> & to_be_erased);
		bool update_neighs(int a, int delta_deg_in, int delta_deg_out, double delta_win, double delta_wout, int nstar, int kout_g_in, int kout_g_out, int tmin, int tmout, int ktin, int ktout);
		//int internal_degree_of(int a);
		int best_node(int & lab, double & best_fitness, int kout_g_in, int kout_g_out, int Nstar, int nneighs, int tmin, int tmout);
		int worst_node(int & lab, double & worst_fitness, int kout_g_in, int kout_g_out, int Nstar, int nneighs, int tmin, int tmout);
		int external_xs(int a, int b, int c, int d, int Nstar, int nneighs, deque<pair<double, int> > & logrs_lab, const double & max_r_one, const double & maxr_two);
		int external_xs(int a, int b, int c, int d, int Nstar, int nneighs, deque<pair<double, int> > & logrs_lab, const double & max_r_one, const double & maxr_two, set<int> & ONL);
		bool is_internal(int a);
		
		
		map<int, facts> lab_facts;													// maps the label into the facts
		multimap<double, int> fitness_lab;											// maps the fitness into the label  (this can be optimized)
	
	
	

};


void weighted_tabdeg::clear() {

	lab_facts.clear();
	fitness_lab.clear();


}


/*int weighted_tabdeg::internal_degree_of(int a) {
	
	
	map<int, facts>::iterator itm= lab_facts.find(a);
	if(itm!=lab_facts.end()) {
		
		
		return itm->second.internal_degree;
		
	}
	
	return -1;


}*/

void weighted_tabdeg::edinsert(int a, int kpin, int kpout, int ktin, int ktout, double mtlwin, double mtlwout, double fit) {		// this function inserts element a (or edit it if it was already inserted)
	
	
	erase(a);
	
	multimap<double, int>::iterator	fiit = fitness_lab.insert(make_pair(fit, a));
	facts F(kpin, kpout, mtlwin, mtlwout, fiit, ktin, ktout);
	
	
	
	lab_facts.insert(make_pair(a, F));
	

}


bool weighted_tabdeg::erase(int a) {		// this function erases element a if exists (and returns true)
	
	
	
	map<int, facts>::iterator itm= lab_facts.find(a);
	if(itm!=lab_facts.end()) {
		
		
		fitness_lab.erase(itm->second.fitness_iterator);
		lab_facts.erase(itm);
		
		return true;
		
	}
	
	return false;
	
}



bool weighted_tabdeg::is_internal(int a) {
	
	
	map<int, facts>::iterator itm= lab_facts.find(a);
	if(itm==lab_facts.end())
		return false;
	
	
	return true;
	
}




void weighted_tabdeg::set_deque(deque<int> & vv) {
	
	vv.clear();
	
	for(map<int, facts>::iterator itm= lab_facts.begin(); itm!=lab_facts.end(); itm++)
		vv.push_back(itm->first);

}


void weighted_tabdeg::print_nodes(ostream & outb) {
	
	
	cout<<"printing nodes:.. (lab intk mtlw fitness degree) "<<size()<<endl;
	
	for(map<int, facts>::iterator itm= lab_facts.begin(); itm!=lab_facts.end(); itm++)
		cout<<itm->first<<" "<<"indegrees(in - out):: "<<itm->second.internal_indegree<<" "<<itm->second.internal_outdegree<<" weights "<<
		itm->second.minus_log_total_wrin<<" "<<itm->second.minus_log_total_wrin<<" "<<(itm->second.fitness_iterator)->first<<" "<<itm->second.indegree<<endl;
	
	
	
}


int weighted_tabdeg::worst_node(int & lab, double & worst_fitness, int kout_g_in, int kout_g_out, int Nstar, int nneighs, int tmin, int tmout) {
	
	
	
	//cout<<"worst_node fitness - lab - (cgroup)"<<endl;
	//prints(fitness_lab);
	
	// I can try to compute the fitness here
	
	lab=-1;
	worst_fitness=-1;
	
	multimap<double, int>:: iterator bit= fitness_lab.end();
	if(bit==fitness_lab.begin())
		return -1;
	
	
	int stopper=0;
	while(bit!=fitness_lab.begin()) {
		
		bit--;
		map<int, facts> :: iterator itm= lab_facts.find(bit->second);
		
		
		
		int kout_g_out_prime= kout_g_out + itm->second.internal_indegree + itm->second.internal_outdegree - itm->second.outdegree;
		int kout_g_in_prime= kout_g_in + itm->second.internal_indegree + itm->second.internal_outdegree - itm->second.indegree;
		
		
		//cout<<"uno"<<endl;
		
		double F= compute_global_fitness_randomized(itm->second.internal_indegree, kout_g_in_prime,	itm->second.internal_outdegree, kout_g_out_prime,
													tmin + itm->second.indegree, tmout + itm->second.outdegree, itm->second.indegree, itm->second.outdegree,
													itm->second.minus_log_total_wrin, itm->second.minus_log_total_wrout, nneighs+1, Nstar+1);
		
		
		
		//cout<<itm->first<<" "<<F<<" worstnode... lab-fit "<<endl;
		
		
		if(F>worst_fitness) {
			
			worst_fitness= F;
			lab=itm->first;
		}
		
		stopper++;
		if(stopper==num_up_to)
			break;
	
	}
	
	
	
	
	
	
	return 0;

}








int weighted_tabdeg::best_node(int & lab, double & best_fitness, int kout_g_in, int kout_g_out, int Nstar, int nneighs, int tmin, int tmout) {
	
	
	// I can try to compute the fitness here
	
	/*cout<<"NE BEST NODE "<<endl;
	
	cout<<"fitness_lab  "<<endl;
	prints(fitness_lab);*/

	
	
	lab=-1;
	best_fitness=1;
	
	multimap<double, int>:: iterator bit= fitness_lab.begin();
	if(bit==fitness_lab.end()) {
		return -1;
	}
	
	int stopper=0;
	while(bit!=fitness_lab.end()) {
		
		map<int, facts> :: iterator itm= lab_facts.find(bit->second);
		
		//cout<<"due"<<endl;
		double F= compute_global_fitness_randomized(itm->second.internal_indegree, kout_g_in, itm->second.internal_outdegree, kout_g_out, tmin, tmout, itm->second.indegree, itm->second.outdegree,
													itm->second.minus_log_total_wrin, itm->second.minus_log_total_wrout, nneighs, Nstar);
		
		
		
		
		//cout<<itm->first<<" "<<F<<" ... node-fit"<<endl;
		
		
		if(F<best_fitness) {
			
			best_fitness= F;
			lab=itm->first;
		}
		
		stopper++;
		if(stopper==num_up_to)
			break;
		
		bit++;


	}
	
	
	
	
	
	
	return 0;

}






int weighted_tabdeg::external_xs(int kout_g_in, int kout_g_out, int tmin, int tmout, int Nstar, int nneighs, deque<pair<double, int> > & logrs_lab, const double & max_r_one, const double & maxr_two, set<int> & ONL) {
	
		
	double max_r=min(max_r_one, maxr_two);
	
	//cout<<"---> "<<max_r<<endl;

	logrs_lab.clear();

	multimap<double, int>:: iterator bit= fitness_lab.begin();
	if(bit==fitness_lab.end()) {
		return -1;
	}
	

	int counter=0;
	while(bit!=fitness_lab.end()) {
		
		
		map<int, facts> :: iterator itm= lab_facts.find(bit->second);
		
		
		if(ONL.find(itm->first)!= ONL.end()) {
			
			double F;
			/*if(new_print)
				F= compute_global_fitness_randomized_print(itm->second.internal_degree, kout_g, tm, itm->second.degree, itm->second.minus_log_total_wr, nneighs, Nstar);
			else*/
			
			//cout<<"tre"<<endl;
			F= compute_global_fitness_randomized(itm->second.internal_indegree, kout_g_in, itm->second.internal_outdegree, kout_g_out, tmin, tmout, itm->second.indegree, itm->second.outdegree,
													itm->second.minus_log_total_wrin, itm->second.minus_log_total_wrout, nneighs, Nstar);
			
			
			//cout<<"ext_nodes:::: "<<itm->first<<" "<<F<<" <-new FFIITT old-> "<<bit->first<<" "<<max_r<<endl;
		
			if(F>max_r && counter>num_up_to)
				break;
		
			logrs_lab.push_back(make_pair(-F, itm->first));
		}
		
		bit++;
		counter++;

	}
	
	
	/*
	if(new_print) {
		cout<<"EXTERNAL stuff"<<endl;
		
			//*
			for(int i=0; i<logrs_lab.size(); i++) {
				
				cout<<"lab-fit-ext: "<<logrs_lab[i].first<<" "<<logrs_lab[i].second<<endl;
			
			
			}
			
			cout<<"end .... EXTERNAL stuff"<<endl;
	}//*/


	if(logrs_lab.empty())
		return -1;
	
	sort(logrs_lab.begin(), logrs_lab.end());
	
		
	
	return 0;

}




int weighted_tabdeg::external_xs(int kout_g_in, int kout_g_out, int tmin, int tmout, int Nstar, int nneighs, deque<pair<double, int> > & logrs_lab, const double & max_r_one, const double & maxr_two) {
	
		
	double max_r=min(max_r_one, maxr_two);
	
	//cout<<"---> "<<max_r<<endl;

	logrs_lab.clear();

	multimap<double, int>:: iterator bit= fitness_lab.begin();
	if(bit==fitness_lab.end()) {
		return -1;
	}
	

	int counter=0;
	while(bit!=fitness_lab.end()) {
		
		
		map<int, facts> :: iterator itm= lab_facts.find(bit->second);
				
		//double F= compute_global_fitness_randomized(itm->second.internal_degree, kout_g, tm, itm->second.degree, itm->second.minus_log_total_wr, nneighs, Nstar);
		//cout<<"++++++++++++++++++++++++++++++"<<endl;
		
		//cout<<"quatt"<<endl;
		double F= compute_global_fitness_randomized(itm->second.internal_indegree, kout_g_in, itm->second.internal_outdegree, kout_g_out, tmin, tmout, itm->second.indegree, itm->second.outdegree,
											itm->second.minus_log_total_wrin, itm->second.minus_log_total_wrout, nneighs, Nstar);
		
		
		//cout<<"ext_nodes:::: "<<itm->first<<" "<<F<<" <-new FFIITT old-> "<<bit->first<<" "<<max_r<<endl;
		//cherr();
		
		if(F>max_r && counter>num_up_to)
			break;
		
		logrs_lab.push_back(make_pair(-F, itm->first));
		
		bit++;
		counter++;
	}
	
	
	/*
	cout<<"EXTERNAL stuff"<<endl;
	
	
	for(int i=0; i<logrs_lab.size(); i++) {
		
		cout<<"lab-fit-ext: "<<logrs_lab[i].first<<" "<<logrs_lab[i].second<<endl;
	
	
	}
	
	cout<<"end .... EXTERNAL stuff"<<endl;
	//*/
	


	if(logrs_lab.empty())
		return -1;
	
	sort(logrs_lab.begin(), logrs_lab.end());
	
		
	
	return 0;

}



void weighted_tabdeg::_set_(weighted_tabdeg & one)  {

	
	clear();
	for(map<int, facts>::iterator itm= one.lab_facts.begin(); itm!=one.lab_facts.end(); itm++)
		edinsert(itm->first, itm->second.internal_indegree, itm->second.internal_outdegree, itm->second.indegree, itm->second.outdegree, 
		itm->second.minus_log_total_wrin, itm->second.minus_log_total_wrout, (itm->second.fitness_iterator)->first);
		
		//(itm->first, itm->second.internal_degree, itm->second.degree, itm->second.minus_log_total_wr, (itm->second.fitness_iterator)->first);
	
	
	

}



bool weighted_tabdeg::update_group(int a, int delta_degree_out, int delta_degree_in, double delta_mtlwout, double delta_mtlwin, 
											int nstar, int nn, int kout_g_in, int kout_g_out, int tmin, int tmout, int ktin, int ktout, deque<int> & to_be_erased) {

	
	
	// this function is to change the internal degree and mtlw of a certain node (to insert it or erase if necessary)
	
	
	map<int, facts>::iterator itm= lab_facts.find(a);
	if(itm==lab_facts.end())
		return false;
	
	
	itm->second.minus_log_total_wrin+=delta_mtlwin;
	itm->second.minus_log_total_wrout+=delta_mtlwout;
	
	itm->second.internal_indegree+=delta_degree_in;
	itm->second.internal_outdegree+=delta_degree_out;
	
	
	if( (itm->second.internal_indegree + itm->second.internal_outdegree)==0 && size()>1) {
		to_be_erased.push_back(a);
		return true;
	}
	
	
	
	
	//cout<<"UPdating... group "<<a<<" delta_deg_in "<<delta_degree_in<<" "<<endl;
	
	
	

	int kout_g_out_prime= kout_g_out + itm->second.internal_indegree + itm->second.internal_outdegree - itm->second.outdegree;
	int kout_g_in_prime= kout_g_in + itm->second.internal_indegree + itm->second.internal_outdegree - itm->second.indegree;
		
	//cout<<"cin"<<endl;
	double fit= compute_global_fitness_randomized(itm->second.internal_indegree, kout_g_in_prime,	itm->second.internal_outdegree, kout_g_out_prime,
													tmin + itm->second.indegree, tmout + itm->second.outdegree, itm->second.indegree, itm->second.outdegree,
													itm->second.minus_log_total_wrin, itm->second.minus_log_total_wrout, nn+1, nstar+1);
	

	fitness_lab.erase(itm->second.fitness_iterator);
	multimap<double, int>::iterator	fiit = fitness_lab.insert(make_pair(fit, a));
	itm->second.fitness_iterator=fiit;
	

	return true;


}







bool weighted_tabdeg::update_neighs(int a, int delta_deg_out, int delta_deg_in, double delta_wout, double delta_win, int nstar, int kout_g_in, int kout_g_out, int tmin, int tmout, int ktin, int ktout) {
	
	
	
	// this function is to change the internal degree and mtlw of a certain node (to insert it or erase if necessary)
	
	//cout<<"UPdating... neighs "<<a<<" "<<kt<<endl;


	map<int, facts>::iterator itm= lab_facts.find(a);
	if(itm==lab_facts.end()) {
		
		edinsert(a, 0, 0, ktin, ktout, 0, 0, 1);
		itm= lab_facts.find(a);
	}
	
	
	itm->second.internal_indegree+=delta_deg_in;
	itm->second.internal_outdegree+=delta_deg_out;


	if((itm->second.internal_indegree + itm->second.internal_outdegree)==0) {
		//cout<<"erased from neigh update "<<a<<endl;
		erase(a);
		return true;
	}
	
	itm->second.minus_log_total_wrin+=delta_win;
	itm->second.minus_log_total_wrout+=delta_wout;
	
	
	
	//cout<<"sei"<<endl;
	double fit= compute_global_fitness_randomized(itm->second.internal_indegree, kout_g_in,	itm->second.internal_outdegree, kout_g_out,
													tmin, tmout, itm->second.indegree, itm->second.outdegree,
													itm->second.minus_log_total_wrin, itm->second.minus_log_total_wrout, size(), nstar);

	
	
	fitness_lab.erase(itm->second.fitness_iterator);
	multimap<double, int>::iterator	fiit = fitness_lab.insert(make_pair(fit, a));
	itm->second.fitness_iterator=fiit;
	
	
	return true;
	
	
}



