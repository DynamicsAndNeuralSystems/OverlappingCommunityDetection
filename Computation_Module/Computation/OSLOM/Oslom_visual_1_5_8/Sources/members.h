#include "./standard_package/standard_include.cpp"
# define equal_up_to 0.995



// this just tells you if you already found a module like that or not

class ne_modules {		// 'ne' stays for 'not equal'


	public:
				
		ne_modules(int a){ clear(a); };
		~ne_modules(){};
		
		
		void erase(const deque<int> &);
		bool insert(const deque<int> &);
		bool insert(const deque<int> &, int & );
		void clear(int);
		deque<deque<int> > memberships;		// on each row the label of a module

		int size(){return sizes.size(); };
	private:
		
		deque<int> sizes;						// the sizes of the modules (label-size)
		
};






void ne_modules::clear(int dim) {
	
	
	sizes.clear();
	memberships.clear();
	deque<int> first;
	for(int i=0; i<dim; i++)
		memberships.push_back(first);


}


void ne_modules::erase(const deque<int> & c) {



	map<int, int> com_ol;		// it maps the index of the comms into the overlap
	
	for(int i=0; i<c.size(); i++) {
		
		for(int j=0; j<memberships[c[i]].size(); j++)
			int_histogram(memberships[c[i]][j], com_ol);
	
	}
	
	
	for(map<int, int>::iterator itm=com_ol.begin(); itm!=com_ol.end(); itm++)
		if(itm->second==c.size() && itm->second==sizes[itm->first]) {		// if they are equal
			
			for(int i=0; i<c.size(); i++)
				memberships[c[i]].erase(remove(memberships[c[i]].begin(), memberships[c[i]].end(), itm->first), memberships[c[i]].end());
			
		
		}

}

bool ne_modules::insert(const deque<int> & c) {

	int gm;
	
	return insert(c, gm);

}

bool ne_modules::insert(const deque<int> & c, int & gm) {

	
	gm=0;
	
	map<int, int> com_ol;		// it maps the index of the comms into the overlap
	
	for(int i=0; i<c.size(); i++) {
		
		for(int j=0; j<memberships[c[i]].size(); j++)
			int_histogram(memberships[c[i]][j], com_ol);
	
	}
	
	
	for(map<int, int>::iterator itm=com_ol.begin(); itm!=com_ol.end(); itm++)
		if(double(itm->second) >= equal_up_to * c.size() && double(itm->second) >= equal_up_to * sizes[itm->first]) {
			
			gm=itm->first;
			return false;
		
		}
	
	
	int new_name=sizes.size();
	sizes.push_back(c.size());
	
	for(int i=0; i<c.size(); i++)
		memberships[c[i]].push_back(new_name);
	
	gm=new_name;
	
	return true;

}









