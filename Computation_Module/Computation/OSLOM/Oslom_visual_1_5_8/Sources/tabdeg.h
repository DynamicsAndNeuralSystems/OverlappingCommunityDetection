
// external_xs can be implemented better
// best_node and worst_node should be stochastic (?) I don't know. it's cumbersome

#include "./standard_package/standard_include.cpp"
//#include "border_score.h"



//#include <gsl/gsl_cdf.h>
//#include <gsl/gsl_randist.h>
//#include <gsl/gsl_sf_gamma.h>
//#define bin_vs_hyper 0.2

/*
double compute_r_hyper_andrea(int kin_node, int kout_g, int tm, int degree_node) {


	double sum=0;
	for(int i=kin_node; i<=degree_node; i++)
		sum+=compute_hypergeometric(i, degree_node, kout_g, tm);
	
	if(sum<0)
		sum=0;
	
	return sum;



}
//*/


inline double right_tail_order_exponential(int N, int pos, double x) {


	// N is the total number of variables, pos is from 1 to N. 1 is the smallest.
	// the routine computes the probality c_pos=  p(X_pos >= x)
	
	
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


inline double compute_r_hyper(int kin_node, int kout_g, int tm, int degree_node) {
	

	//cout<<kin_node<<" "<<kout_g<<" "<<tm<<" "<<degree_node<<endl;
	return LOG_TABLE.right_cumulative_function(degree_node, kout_g, tm, kin_node);

}



inline double hyper_table(int kin_node, int kout_g, int tm, int degree_node) {
	
	return LOG_TABLE.hyper(kin_node, kout_g, tm, degree_node);
	

}



typedef multiset<pair<int, int> >  muspi;		// muspi is a multiset of internal degree - label


class tabdeg {


	public:
				
		tabdeg(){};
		~tabdeg(){};
		
		void _set_(tabdeg &);
		bool is_internal(int );
		void edinsert(int , int , int );
		bool erase(int, int);
		int indegof(int );
		int size() {return nodes_deginl.size();};
		void print_indeg(ostream &);
		void print_nodes(ostream &);
		int totin();
		int best_node(int, int, double &, int &, int &);
		int worst_node(int, int, double &, int &, int &);
		void set_deque(deque<int> & );
		void clear();
		int external_xs(int kout_g, int tm, deque<pair<double, int> > & logrs_lab, const double & max_r, const double & maxr_two);
		//int external_xs(int kout_g, int tm, deque<pair<double, int> > & logrs_lab, double max_r, set<int> & erased);
	
	
	
	private:

		map<int, muspi> deg_inlab;								// maps the degree into the pair internal_degree - label
		map <int, muspi::iterator> nodes_deginl;				// maps the label into an iterator to the pair internal_degree - label
	
	
	
	

};


void tabdeg::clear() {

	deg_inlab.clear();
	nodes_deginl.clear();


}

bool tabdeg::is_internal(int a) {

	
	map<int, muspi::iterator >::iterator itm= nodes_deginl.find(a);
	if(itm==nodes_deginl.end())
		return false;
	else
		return true;




}

void tabdeg::edinsert(int a, int kp, int kt) {		// this function inserts element a (or edit it if it was already inserted)
	
	erase(a, kt);
	
	map<int, muspi >::iterator itf=deg_inlab.find(kt);
	
	if (itf==deg_inlab.end()) {
		
		muspi first;
		itf = deg_inlab.insert(make_pair(kt, first)).first;
		
	}
	
		
	muspi::iterator itms= itf->second.insert(make_pair(kp, a));
	nodes_deginl.insert(make_pair(a, itms));
	
	


}


bool tabdeg::erase(int a, int kt) {		// this function erases element a if exists (and returns true)
	
	
	
	map<int, muspi::iterator >::iterator itm= nodes_deginl.find(a);
	if(itm!=nodes_deginl.end()) {
		
		deg_inlab[kt].erase(itm->second);
		nodes_deginl.erase(itm);
		return true;
		
	}
		
	
	return false;
	

}


int tabdeg::indegof(int a) {		// return the internal degree of a, 0 if it's not internal
	
	map<int, muspi::iterator >::iterator itm= nodes_deginl.find(a);
	if(itm!=nodes_deginl.end())
		return itm->second->first;
	else
		return 0;
		


}




void tabdeg::set_deque(deque<int> & vv) {
	
	
	vv.clear();
	
	for(map<int, muspi::iterator >::iterator itm= nodes_deginl.begin(); itm!=nodes_deginl.end(); itm++)
		vv.push_back(itm->first);



}


void tabdeg::print_nodes(ostream & outb) {
	
	for(map<int, muspi::iterator >::iterator itm= nodes_deginl.begin(); itm!=nodes_deginl.end(); itm++)
		outb<<itm->first<<"\t"<<itm->second->first<<endl;
	
	
	
}

void tabdeg::print_indeg(ostream & outb) {
	
	for(map<int, muspi>::iterator itm= deg_inlab.begin(); itm!=deg_inlab.end(); itm++) {
		outb<<itm->first<<":\t";
		for(muspi::iterator itmu=itm->second.begin(); itmu!=itm->second.end(); itmu++)
			outb<<itmu->first<<" ";
		outb<<endl;
		
	}
	
	
}


int tabdeg::totin() {
	
	int kkin=0;
	
	for(map<int, muspi >::iterator itm= deg_inlab.begin(); itm!=deg_inlab.end(); itm++) {
		
		for(muspi::iterator itmm=itm->second.begin(); itmm!=itm->second.end(); itmm++)
			kkin+=itmm->first;
			
	}
	
	
	return kkin;
	
}

int tabdeg::best_node(int kout_g, int tm, double & lowest_r, int & benode_kp, int & degree_ref) {
	
	// this works only for external nodes
	//cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;

	benode_kp=0;
	lowest_r=1.1;
	deque<int> best_nodes;


	for(map<int, muspi>::iterator itm= deg_inlab.begin(); itm!=deg_inlab.end(); itm++) if(!itm->second.empty()) {
		
		
		
		muspi::iterator itb= itm->second.end();
		itb--;
		
		double r=1;
		
		//cout<<"itb->first "<<itb->first<<" "<<benode_kp<<endl;
		if(itb->first>benode_kp)
			r=compute_r_hyper(itb->first, kout_g, tm, itm->first);
		
		if(r<lowest_r) {
			
			
			degree_ref=itm->first;
			benode_kp=itb->first;
			lowest_r=r;
			best_nodes.clear();
			
			while(itb->first==benode_kp) {
			
				best_nodes.push_back(itb->second);
				if(itb==itm->second.begin())
					break;
				
				itb--;
			
			}
		
		
		
		}
		
		
		
	}
	
	if(best_nodes.empty()) {
		//cerr<<"tabdeg is empty!"<<endl;
		benode_kp=0;
		degree_ref=0;
		lowest_r=1;
		return -1;
	
	}
	
	
	//prints(best_nodes);
	return best_nodes[irand(best_nodes.size()-1)];



}








int tabdeg::worst_node(int kout_g, int tm, double & biggest_r, int & wnode_kp, int & wnode_kt) {
	
	// this works only for internal nodes
	
	if(nodes_deginl.empty()) {

		wnode_kp=0;
		wnode_kt=0;
		biggest_r=1;
		return -1;
	}

	
	
	biggest_r=-1;
	deque<int> worst_nodes;

	map<int, muspi>::iterator itmuspi= deg_inlab.end();
	itmuspi--;		// it starts from the end
	wnode_kp=itmuspi->first +1;
	itmuspi++;

	
	while(itmuspi!=deg_inlab.begin()) {
	
		itmuspi--;
		
		if(!itmuspi->second.empty()) {
			
			
			
			muspi::iterator itb= itmuspi->second.begin();
			
			double r=-2;
			
			//cout<<"itb->first "<<itb->first<<" "<<wnode_kp<<" on "<<itmuspi->first<<endl;
			if(itb->first<wnode_kp) {
				r=compute_r_hyper(itb->first, kout_g + 2 * itb->first - itmuspi->first, tm + itmuspi->first, itmuspi->first);
				//cout<<"yes : "<<r<<" "<<deg_inlab.size()<<endl;
			}
			
			/*
			else {
			
				r=compute_r_hyper(itb->first, kout_g + 2 * itb->first - itmuspi->first, tm + itmuspi->first, itmuspi->first);
				cout<<"no : "<<r<<endl;
			
			
			} 
			//*/
			
			if(r>biggest_r) {
				
				wnode_kp=itb->first;
				wnode_kt=itmuspi->first;
				biggest_r=r;

				worst_nodes.clear();
				
				while(itb->first==wnode_kp) {
				
					worst_nodes.push_back(itb->second);
					itb++;
					
					if(itb==itmuspi->second.end())
						break;
					
					
				
				}
			
			
			
			}
			
			
			
		}
	}
	
	
	
	//cout<<"worst nodes "<<endl;
	//prints(worst_nodes);
	return worst_nodes[irand(worst_nodes.size()-1)];



}



/*
double tabdeg::r_of_this_internal(int kout_g, int tm, int kin_node, int degree) {

	// it works only for internal nodes
	
	
	double r=compute_r_hyper(kin_node, kout_g + 2 * kin_node - degree, tm + degree, degree);
	cout<<"r_of_this_internal: "<<kin_node<<" "<<degree<<" "<<r<<endl;

	return r;




}
*/








int tabdeg::external_xs(int kout_g, int tm, deque<pair<double, int> > & logrs_lab, const double & max_r_one, const double & maxr_two) {
	
	// this works only for external nodes
	// it computes all the -log(r) for the neighbors such that r<max_r 
	
	double max_r=min(max_r_one, maxr_two);
	//cout<<"max_r: "<<max_r<<endl;
	
	
	//cout<<"external size: "<<size()<<endl;

	logrs_lab.clear();


	for(map<int, muspi>::iterator itm= deg_inlab.begin(); itm!=deg_inlab.end(); itm++) if(!itm->second.empty()) {
		
		
		
		muspi::iterator itb= itm->second.end();
		itb--;
		
		while(true) {
			
			
			double b2=compute_r_hyper(itb->first +1, kout_g, tm,  itm->first);

			if(b2>max_r)
				break;
			
			
			double r= b2 + ran4() * (hyper_table(itb->first, kout_g, tm,  itm->first));
			
			//cout<<itb->first<<" "<<kout_g<<" "<<tm<<" "<<itm->first<<" "<<r<<endl;
			
			
			/*if(itb->first==1 && itm->first==1)
				cout<<r<<" "<<kout_g<<" "<<tm<<endl;*/
			
			
			//if(itb->first>0 && (r<max_r || itb->first==itm->first)) 
			if(itb->first>0 && r<max_r) 
				logrs_lab.push_back(make_pair(-r, itb->second));
			
			
			if(itb==itm->second.begin())
				break;
			itb--;

		
		}
		
	}
	
	//cout<<"..."<<endl;
	
	if(logrs_lab.empty())
		return -1;
	
	
	sort(logrs_lab.begin(), logrs_lab.end());
		
	
	return 0;



}



void tabdeg::_set_(tabdeg & one)  {

	
	clear();
	for(map<int, muspi>::iterator itm= one.deg_inlab.begin(); itm!=one.deg_inlab.end(); itm++) {
		
		for(muspi::iterator itmu=itm->second.begin(); itmu!=itm->second.end(); itmu++)
			edinsert(itmu->second, itmu->first, itm->first);
		
	}
	


}



