#include <fstream>
#include <cstring>

#include "collection.h"
#include "../objects/vector.h"
#include "../objects/matrix.h"


Collection::Collection()
{
	this->min_length = 1000000000;
	this->max_length = 0;
	this->D = 0;
}

Collection::Collection(vector<Object*> _objects, vector<int> _labels)
{
	this->objects = _objects;
	this->labels = _labels;
	this->min_length = 1000000000;
	this->max_length = 0;
	this->D = 0;

	for(auto lb : labels)
	{
		classes.insert(lb);
	}

	for(Object* P : objects)
	{
		min_length = std::min(min_length, P->get_n());
		max_length = std::max(max_length, P->get_n());
	}

	if(dynamic_cast<Vector*>(objects[0]))
	{// univariate time series
		this->D = 1;
	}
	else if(dynamic_cast<Matrix*>(objects[0]))
	{// multivariate time series
		this->D = ((Matrix*)objects[0])->get_D();
	}
}

Collection::~Collection()
{
	// freee memory
	vector<Object*>().swap(objects);
	vector<int>().swap(labels);
	set<int>().swap(classes);
}



void Collection::display()
{
	cout<<"# objects: "<<objects.size()<<endl;	
	cout<<"# variables: "<<D<<endl;
	cout<<"# classes: "<<classes.size()<<endl;	
	cout<<"Min & Max length: ["<<min_length<<" - "<<max_length<<"]"<<endl;
}

void Collection::display_all()
{
	int i=0;
	for(Object* P : objects)
	{
		cout<<"----------------------------------\n";
		P->display(); cout<<endl;
		cout<<"label:"<<labels[i++]<<"\n";				
	}
}

Collection* Collection::read_UTS(string name_file, bool normalize, int limit)
{
	ifstream input(name_file, ios::in);

	if(!input.good())
	{
		cout<<"The file "<<name_file<<" was not found \n";
		return NULL;
	}

	char buffer[32194];
	char *tokens;
	vector<double> vec;
	double *arr;

	vector<Object*> objects;
	vector<int> labels;

	int index=0;
	while(!input.eof() && index < limit)
	{
		input.getline(buffer,32194);
		if(strlen(buffer)==0) continue;

		tokens=strtok(buffer," \t"); // class
		labels.push_back(atof(tokens));

		tokens=strtok(NULL," \t");
		while(tokens != NULL)
		{
			vec.push_back(atof(tokens));
			tokens = strtok(NULL, " \t");
		}

		arr = new double[vec.size()];
		std::copy(vec.begin(), vec.end(), arr);		

		objects.push_back(new Vector(arr, vec.size(), index, normalize));

		vec.clear();
		index++;		
	}

	input.close();

	return new Collection(objects, labels);
}


Collection* Collection::read_MTS(string name_file, bool normalize)
{
    ifstream input(name_file,ios::in);

    if(!input.good())
	{
		cout<<"The file "<<name_file<<" was not found \n";
		return NULL;
	}

    char buffer[32194];
    char *token;
    int id_ant=-1,lb_ant=-1;
	int id,lb;
	int n, D, i, j;
	double** ma;

	vector<double> obs;
	vector< vector<double> > mts;

	vector<Object*> objects;
	vector<int> labels;

    while(!input.eof())
    {
        input.getline(buffer,32194);
        if(strlen(buffer)==0) continue;

        token=strtok(buffer," ,\t");//serie
		id=atoi(token);
		token=strtok(NULL," \t"); // index time
		token=strtok(NULL," \t"); // class label
        lb=atoi(token);
        token=strtok(NULL," \t"); // observacion
		obs.clear();
        while(token != NULL) // elemnts
        {
            obs.push_back(atof(token));
            token = strtok(NULL," \t");
        }

        if(id != id_ant && mts.size() > 0)
		{
			n = mts.size();
			D = mts[0].size();
			
			ma = new double*[n];
			for( i = 0; i < n; ++i)
			{
				ma[i] = new double[D];
				for ( j = 0; j < D; ++j)				
					ma[i][j] = mts[i][j];				
			}

			objects.push_back(new Matrix(ma, n, D, id_ant, normalize));

			labels.push_back(lb_ant);			

			mts.clear();
		}

		mts.push_back(obs);

		id_ant = id;	lb_ant = lb;
    }

    if(mts.size() > 0)
	{
		n = mts.size();
		D = mts[0].size();
		
		ma = new double*[n];
		for( i = 0; i < n; ++i)
		{
			ma[i] = new double[D];
			for ( j = 0; j < D; ++j)				
				ma[i][j] = mts[i][j];				
		}

		objects.push_back(new Matrix(ma, n, D, id_ant));

		labels.push_back(lb_ant);	
	}

    input.close();

    return new Collection(objects, labels);
}


Collection* Collection::read_MTS_inv(string name_file, bool normalize)
{
    ifstream input(name_file,ios::in);

    if(!input.good())
	{
		cout<<"The file "<<name_file<<" was not found \n";
		return NULL;
	}

    char buffer[32194];
    char *token;
    int id_ant=-1,lb_ant=-1;
	int id,lb;
	int n, D, i, j;
	double** ma;

	vector<double> obs;
	vector< vector<double> > mts;

	vector<Object*> objects;
	vector<int> labels;


    while(!input.eof())
    {
        input.getline(buffer,32194);
        if(strlen(buffer)==0) continue;


        token=strtok(buffer," ,\t");//serie
		id=atoi(token);
		token=strtok(NULL," \t"); // index time
		token=strtok(NULL," \t"); // class label
        lb=atoi(token);
        token=strtok(NULL," \t"); // observacion
		obs.clear();
        while(token != NULL) // elemnts
        {
            obs.push_back(atof(token));
            token = strtok(NULL," \t");
        }

        if(id != id_ant && mts.size() > 0)
		{
			n = mts.size();
			D = mts[0].size();
			// ma: matriz inversa de mts
			ma = new double*[D];
			for( j = 0; j < D; ++j)
			{
				ma[j] = new double[n];
				for ( i = 0; i < n; ++i)
				{
					ma[j][i] = mts[i][j];
				}
			}

			objects.push_back(new Matrix(ma, D, n, id_ant, normalize, true));
			labels.push_back(lb_ant);			

			mts.clear();
		}

		mts.push_back(obs);

		id_ant = id;	lb_ant = lb;
    }

    if(mts.size() > 0)
	{
		n = mts.size();
		D = mts[0].size();
		// ma: matriz inversa de mts
		ma = new double*[D];
		for( j = 0; j < D; ++j)
		{
			ma[j] = new double[n];
			for ( i = 0; i < n; ++i)
			{
				ma[j][i] = mts[i][j];
			}
		}

		objects.push_back(new Matrix(ma, D, n, id_ant) );

		labels.push_back(lb_ant);			

		mts.clear();
	}

    input.close();

    return new Collection(objects, labels);
}


