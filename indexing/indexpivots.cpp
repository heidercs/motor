#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <math.h>
#include <float.h>
#include <set>

#include "indexpivots.h"
#include "../util/functions.h"


IndexPivots::IndexPivots(Distance *_dist, string _class_object) 
{
	path_ind = "";
    type_store = IN_RAM;
    class_object = _class_object;

    this->dist=_dist;
    this->n_pivots  = -1;        
}

IndexPivots::IndexPivots(Distance *_dist, string _path_ind, string _class_object) 
{
	path_ind = _path_ind;
    type_store = IN_RAM;
    class_object = _class_object;

    this->dist=_dist;
    this->n_pivots = -1;        
}

IndexPivots::~IndexPivots()
{    
	if(!matrix.empty())
	{
		for(int i=0;i<(int)matrix.size();i++)
			delete[] matrix[i];
	}
	
	vector<Object*>().swap(objects);
	vector<Object*>().swap(pivots);
	vector<double*>().swap(matrix);
}

void IndexPivots::set_params(int _n_pivots)
{
    this->n_pivots = _n_pivots;
    clean_index();
}

void IndexPivots::clean_index()
{	
	if(type_store == IN_DISK)
	{
		string file = path_ind + "/pivots";
        remove(file.c_str());
    }
    vector<Object*>().swap(objects); // limpiar data
	vector<Object*>().swap(pivots);
   	vector<double*>().swap(matrix);
}

bool IndexPivots::insert(Object *P)
{
	if(P->name() != class_object)
    {
        std::cerr<<"This data type is not support\n";
        return false;
    }
    
    double* dist_pivots = new double[n_pivots];
	for(int j=0;j<n_pivots;j++)
	{
		dist_pivots[j] = dist->d(P, pivots[j]);
	}

	this->matrix.push_back(dist_pivots);
    write_object(P);

    return true; 
}


void IndexPivots::insert_all(vector<Object*> _objects)
{	
	//seleccionar n_pivots  
	select_pivots_semirand(_objects);
	//select_pivots_rand(_objects);

	//construir matriz de distancias	
	matrix = vector<double*>((int)objects.size());
	double* dist_pivots;
	int i=0, j=0;
	for(Object* P : objects) 
    {		
		dist_pivots = new double[n_pivots];
		for(j=0;j<n_pivots;j++)
		{
			dist_pivots[j] = dist->d(P, pivots[j]);
		}
		matrix[i] = dist_pivots;
		i++;
	}
}

	
vector<PairObjectDistance> IndexPivots::knn(Object* Q, int k)
{
    double dist_max_k = DBL_MAX, LB=0, d=0;
    List<PairObjectDistance> *heap_max=new List<PairObjectDistance>();
    PairObjectDistance par;
    Object* P;    

	
   	// search on pivots	
	double* matrix_q=new double[n_pivots];	
	for(int j=0; j<n_pivots; j++)
	{
		d = dist->d(pivots[j], Q); CONT_DIST++;
		matrix_q[j] = d;		

		if(heap_max->get_num_nodes() < k)
		{
			heap_max->add_heap_max(PairObjectDistance(pivots[j], d));			
		}
		else
		{
			if(d <= dist_max_k)
			{
                heap_max->add_heap_max(PairObjectDistance(pivots[j], d));
            	heap_max->del_first();
        		dist_max_k = heap_max->get_first().dist; 
			}
		}
	}
	
	
	//preparacion de la cola de prioridad (en este caso, set ordenado)
	vector<PairObjectDistance> heap_min;
	for(int i=0; i<(int)this->objects.size(); i++) 
	{
		LB = DBL_MAX;
		for(int j=0; j<n_pivots; j++)
		{
			d = fabs(matrix_q[j] - this->matrix[i][j]);			
			if(d < LB)
			{
				LB = d;
			}			
		}
		//d(q, x) >= LB
		//si LB > r => d(q, x) > r => descartado (y parar)
		if(LB <= dist_max_k)
			heap_min.push_back(PairObjectDistance(this->objects[i], LB));
	}
	std::sort(heap_min.begin(), heap_min.end());
		

	//revision de la cola de prioridad
	for(auto &par : heap_min)
	{
		if(par.dist > dist_max_k && heap_max->get_num_nodes() == k)    break;
		
		P =  read_object(par.object);
		d = dist->d(P, Q); CONT_DIST++;

		if(heap_max->get_num_nodes()<k)
		{
			heap_max->add_heap_max(PairObjectDistance(P, d));
		}//if... set no lleno (no hay radio para descartar)
		else
		{
			if(d < dist_max_k)
			{
				heap_max->add_heap_max(PairObjectDistance(P, d));
                heap_max->del_first();
                dist_max_k=heap_max->get_first().get_distance();
			}
			else if(type_store == IN_DISK)
                delete P;
			
		}//else... set lleno			
	}//while... cada candidato
	
	 // colocar de menor a mayor 
    vector<PairObjectDistance> results = heap_max->toVectorInv();
        
    delete heap_max;
    heap_min.clear();

	return results;
}

void IndexPivots::select_pivots_semirand(vector<Object*> _objects, int seed)
{
	int desp=_objects.size()/n_pivots;//Tamanio del desplazamiento
    int candidato = 0, i = 0, j = 0;

	this->pivots = vector<Object*>(this->n_pivots);

	//seleccionar n_pivots  
	vector<int> agregados(this->n_pivots);
	srand (time(NULL) * seed);
    for(j=0; j < this->n_pivots;  j++)
    {
        candidato = (int) (rand() % desp) ;        
        agregados[j] = j*desp + candidato;
    }

	//separar pivotes de datos
	i=0; j=0;
	for(Object* P : _objects)
    {
		if(i == agregados[j])
			pivots[j++] = P;
		else
			this->objects.push_back(P);
		i++;
    }
}

void IndexPivots::select_pivots_rand(vector<Object*> _objects, int seed)
{
	int N = (int) _objects.size();
    int candidato=0, i = 0, j = 0;

	this->pivots = vector<Object*>(this->n_pivots);

	//seleccionar n_pivots  
	std::set<int>  agregados;	
	srand (time(NULL) * seed);
	for(j=0; j<n_pivots; j++)
	{
		do
		{
			candidato = (int) (rand() % N);
		}
		while(agregados.find(candidato)!=agregados.end());
		
		agregados.insert(candidato);		
	}

	//separar pivotes de datos
	i=0; j=0;
	for(Object* P : _objects)
    {
		if(agregados.find(i) != agregados.end())
			pivots[j++] = P;
		else
			this->objects.push_back(P);
		i++;
    }
}

Object* IndexPivots::read_object(Object* p)
{	
	if(type_store == IN_RAM)
		return p;
	else
	{		
		return read(path_ind + "/pivots", class_object, ((ItemDisk*)p)->get_position_file());
	}
}

void IndexPivots::write_object(Object* p)
{
	if(type_store == IN_RAM)
		objects.push_back(p);	
	else // IN_DISK
	{
		long pos_in_file=write(path_ind+"/ind_lineal", p);		
		objects.push_back(new ItemDisk(pos_in_file));			
	}	
}

void IndexPivots::display()
{
	cout<<"- number of pivots:"<<this->n_pivots<<endl;
	for(Object* P : pivots)
	{
		cout<<"\t "<<P->get_id_file()<<", ";
	}
	cout<<endl;

	cout<<"- size data:"<<objects.size()<<endl;
	Object* P;
	for(Object* p : objects)
	{
		P = read_object(p); 		        
        std::cout<<"\t "<<P->get_id_file()<<", ";
		if(type_store == IN_DISK)
			delete P;
    }
	cout<<endl;
}

int IndexPivots::size()
{
    return (int)objects.size();
}

void IndexPivots::reset_counters()
{
    CONT_DIST=0;
    CONT_MINDISTS=0;
}

int IndexPivots::num_access()
{
    return CONT_DIST;
}

int IndexPivots::num_dists()
{
    return CONT_DIST;
}

int IndexPivots::num_mindists()
{
    return CONT_MINDISTS;
}
