#ifndef KNNCLASSIFIER_H_
#define KNNCLASSIFIER_H_

#include "../indexing/index.h"
#include "../datamining/collection.h"

class KnnClassifier
{
public:
	static double* knn_test(Collection* train, Collection* test, Distance* dist, int k=1);

	static double* nn_test(Collection* train, Collection* test, Distance* dist);

	static double* knn_test(Index* ind_train, vector<int> train_labels, Collection* test, int k=1);
	
	static double* nn_test(Index* ind_train, vector<int> train_labels, Collection* test);
	
	static double* cv_test(Collection* all_data, Distance* dist, int n=10, int k=1);

	static double* leaveone_cv_test(Collection* all_data, Distance* dist, int k=1);


	static int predict(Collection* train, Object* q, Distance*  dist, int k);	
	static int predict(Index* ind_train, vector<int> labels_train, Object* q, int k);

	static int nn_predict(Collection* train, Object* Q, Distance* dist);
	static int nn_predict(Index* ind_train, vector<int> labels_train, Object* q);

};

class count_class{
public:
	// counter & position of the first object 
	int count, pos_first; 
	count_class()
	{
		count =0;  pos_first = -1;
	};
	
	count_class(int _count, int _pos_first=-1)
	{
		count = _count; pos_first =_pos_first;
	};
};

#endif /* KNNCLASSIFIER_H_ */
