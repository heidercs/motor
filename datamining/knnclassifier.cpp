#include "knnclassifier.h"
#include "../distances/distance.h"
#include "../util/functions.h"
#include "../indexing/indexlineal.h"


double* KnnClassifier::nn_test(Collection* train, Collection* test, Distance* dist)
{
	int class_pred, correct=0;

    double sum_time=0;   
	
	//testing
	int i=0;

    for(Object* Q : test->objects)
    {
        auto start = getTime();
        class_pred = nn_predict(train, Q, dist);
		sum_time += diffTime(start, getTime());
				
        //cout<<"["<<Q->get_id_file()<<"] \tclass:"<<labels_test[i]<<", class pred:"<<class_pred<<endl;

        if(class_pred == test->labels[i])     correct++;

        if((i+1) % 10 == 0)
        	cout<<"\r["<<int(i*100.0/test->objects.size()) <<"%]  "<<flush;

		i++;        
    }    
    cout<<"\r[100%]   \n";

    double error=(test->objects.size() - correct)/(test->objects.size()*1.0);

    double *metrics = new double[2];
    metrics[0] = error;
    metrics[1] = sum_time;
	
	return metrics;
}

double* KnnClassifier::knn_test(Collection* train, Collection* test, Distance* dist, int k)
{
	if(k == 1)
		return nn_test(train, test, dist);
	
	Index* ind_train = new IndexLineal(dist);
	ind_train->insert_all(train->objects);	

	int class_pred, correct=0;        
    double  sum_time=0;

	 //testing
	int i=0;
    for(Object* Q : test->objects)
    {
        auto start = getTime();
        class_pred = predict(ind_train, train->labels, Q, k);
		sum_time += diffTime(start, getTime());
				
        //cout<<"["<<Q->get_id_file()<<"] \tclass:"<<labels_test[i]<<", class pred:"<<class_pred<<endl;

        if(class_pred == test->labels[i])     correct++;

        if((i+1) % 10 == 0)
        	cout<<"\r["<<int(i*100.0/test->objects.size()) <<"%]   "<<flush;

		i++;        
    }    
    cout<<"\r[100%]   \n";

    double error=(test->objects.size() - correct)/(test->objects.size()*1.0);

    double *metrics = new double[2];
    metrics[0] = error;
    metrics[1] = sum_time;

	delete ind_train;

	return metrics;
}

double* KnnClassifier::nn_test(Index* ind_train, vector<int> train_labels, Collection* test)
{
	int class_pred, correct=0;

    double sum_time=0;   
	
	//testing
	int i=0;

    for(Object* Q : test->objects)
    {
        auto start = getTime();
        class_pred = nn_predict(ind_train, train_labels, Q);
		sum_time += diffTime(start, getTime());
				
        //cout<<"["<<Q->get_id_file()<<"] \tclass:"<<labels_test[i]<<", class pred:"<<class_pred<<endl;

        if(class_pred == test->labels[i])     correct++;

        if((i+1) % 10 == 0)
        	cout<<"\r["<<int(i*100.0/test->objects.size()) <<"%]  "<<flush;

		i++;        
    }    
    cout<<"\r[100%]   \n";

    double error=(test->objects.size() - correct)/(test->objects.size()*1.0);

    double *metrics = new double[2];
    metrics[0] = error;
    metrics[1] = sum_time;
	
	return metrics;
}

double* KnnClassifier::knn_test(Index* ind_train, vector<int> train_labels, Collection* test, int k)
{
	int class_pred, correct=0;

    double sum_time=0;   
	
	//testing
	int i=0;

    for(Object* Q : test->objects)
    {
        auto start = getTime();
        class_pred = predict(ind_train, train_labels, Q, k);
		sum_time += diffTime(start, getTime());
				
        //cout<<"["<<Q->get_id_file()<<"] \tclass:"<<labels_test[i]<<", class pred:"<<class_pred<<endl;

        if(class_pred == test->labels[i])     correct++;

        if((i+1) % 10 == 0)
        	cout<<"\r["<<int(i*100.0/test->objects.size()) <<"%]  "<<flush;

		i++;        
    }    
    cout<<"\r[100%]   \n";

    double error=(test->objects.size() - correct)/(test->objects.size()*1.0);

    double *metrics = new double[2];
    metrics[0] = error;
    metrics[1] = sum_time;
	
	return metrics;
}

double* KnnClassifier::cv_test(Collection* all_data, Distance* dist, int n, int k)
{
	double *metrics = new double[2];
	double sum_error=0, sum_time=0;	
	int init, end, N;

	Collection *train = new Collection(), *test = new Collection();

    N = all_data->objects.size();
	init = 0;

	for(int i = 0 ; i <  n; i++)
	{		
		end = std::min(N, (int)round( (i+1.0)*N/(1.0*n) ) ) ;

        for (int j = 0; j < init; ++j)        
        {
            test->objects.push_back(all_data->objects[j]);
            test->labels.push_back(all_data->labels[j]);
        }
                    
        for (int j = init; j < end; ++j)
        {
            train->objects.push_back(all_data->objects[j]);
            train->labels.push_back(all_data->labels[j]);
        }
        
        for (int j = end; j < N; ++j)
        {
            test->objects.push_back(all_data->objects[j]);
            test->labels.push_back(all_data->labels[j]);
        }           

        if(k == 1)
			metrics = nn_test(train, test, dist);
		else
			metrics = knn_test(train, test, dist, k);

        cout<<"i="<<i+1<<", train="<<train->objects.size()<<", test="<<test->objects.size()<<endl;
        cout<<"\t--> error="<<metrics[0]<<endl;

		sum_error += metrics[0];
		sum_time += metrics[1];

        init = end + 1;

		train->objects.clear();
		test->objects.clear();
        train->labels.clear();
        test->labels.clear();		
	}

    metrics[0] = sum_error / n; // avg. error
    metrics[1] = sum_time; // total time

    return metrics;
}

double* KnnClassifier::leaveone_cv_test(Collection* all_data, Distance* dist, int k)
{
	double sum_time=0;	
	int label, q_label, N = all_data->objects.size(), correct = 0;
	Object* q;

	for (int i = 0; i < N; ++i)
	{
		q = all_data->objects[0];
		q_label = all_data->labels[0];

		all_data->objects.erase(all_data->objects.begin());
		all_data->labels.erase(all_data->labels.begin());

        auto start = getTime();
		
		if(k == 1)
			label = nn_predict(all_data, q, dist);
		else
			label = predict(all_data, q, dist, k);

		sum_time += diffTime(start, getTime());

		if(label == q_label)	correct++;

		all_data->objects.push_back(q);
		all_data->labels.push_back(q_label);
	}


	double *metrics = new double[2];
    metrics[0] = (N - correct)/(N*1.0); // error
    metrics[1] = sum_time;

    return metrics;
}


int KnnClassifier::predict(Collection* train, Object* q, Distance*  dist, int k)
{
	if(k == 1)
		return nn_predict(train, q, dist);

	Index* ind_train = new IndexLineal(dist);
	ind_train->insert_all(train->objects);

	int best_class = predict(ind_train, train->labels, q, k);

	delete ind_train;

	return best_class;
}	

int KnnClassifier::nn_predict(Collection* train, Object* Q, Distance* dist)
{
    double d, best_so_far=DBL_MAX;
    int pos_nn = -1;

    int i=0;
	for(Object* P : train->objects)
    {
		d = dist->d(P, Q) ;
		
        if(d < best_so_far)
        {
            pos_nn = i;
            best_so_far = d;
        }

		i++;
    }

    return train->labels[pos_nn];
}

int KnnClassifier::predict(Index* ind_train, vector<int> labels_train, Object* q, int k)
{
	int label;
	count_class cc, cc_max(0);

    unordered_map<int /* label */, count_class /* <count, position> */ > hash;

	vector<PairObjectDistance> res = ind_train->knn(q, k);

	for(int i = 0; i < (int) res.size(); i++)
	{
		label = labels_train[ res[i].object->get_id_file() ];
		if(hash.find(label) == hash.end())
			hash[label] = count_class(0, i);
		else
			hash[label].count = hash[label].count + 1;
	}

	int best_class = -1;
	for(auto it = hash.begin(); it != hash.end(); it++)
	{
		cc = it->second;
		if(cc.count > cc_max.count)
		{
			cc_max = cc;
			best_class = it->first;
		}
		else if(cc.count == cc_max.count)
		{
			if(cc.pos_first < cc_max.pos_first)
			{
				cc_max = cc;
				best_class = it->first;
			}
		}
	}

	return best_class;
}


int KnnClassifier::nn_predict(Index* ind_train, vector<int> labels_train, Object* q)
{
	vector<PairObjectDistance> res = ind_train->knn(q, 1);
	return labels_train[res[0].object->get_id_file()];
}