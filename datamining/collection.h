#ifndef COLLECTION_H_
#define COLLECTION_H_

#include <set>
#include <vector>
#include <climits>
#include "../objects/object.h"

using namespace  std;

class Collection
{

public:
	vector<Object*> objects;
	vector<int> labels;

	set<int> classes;

	// number of variables
	int D;

	// the shortest object  and the longest object
	int min_length, max_length;


	Collection();
	Collection(vector<Object*> _objects, vector<int> _labels);
	virtual ~Collection();
	void display();
	void display_all();


	// method to read time series from the UCR time series archive.
	// each row is a time series
	/*
	 * class  observations...
	 * ----------------------
	 * 3  v1 v2 v3 v4 v5 ....
	 * 1  v1 v2 v3 v4 v5 ....
	 */
	static Collection* read_UTS(string name_file, bool normalize=false, int limit = INT_MAX);


	
	// method to read time series from the UCI machine learning repository and UCR
	// 
	/*
	 * id_serie time_index class  attrib_1  attrib_2 ... attrib_D
	 * ----------------------------------------------------------
	 * 	 
	 */
	static Collection* read_MTS(string name_file, bool normalize=false);
	static Collection* read_MTS_inv(string name_file, bool normalize=false);


};

#endif /* COLLECTION_H_ */
