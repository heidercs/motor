#ifndef INDEXPIVOTS_H
#define INDEXPIVOTS_H

#include <iostream>
#include <vector>

#include "index.h"
#include "../distances/distance.h"
#include "../util/list.h"

class IndexPivots: public Index
{
protected:
    //store all the pivots, size =  K
    vector<Object*> pivots;

    //store the rest of objects, size = N - K
    vector<Object*> objects;	

	//number of pivots
	int n_pivots;

    //Distance 
    Distance* dist;

    //Matriz de distancias con los Pivots
    vector<double*> matrix;	


    //guarda el numero de buckets accesados y distancias computadas por la ultima consulta efectuada
    long int CONT_DIST, CONT_MINDISTS;

public:
    IndexPivots(Distance *_dist, string _class_object); // IN_RAM
    IndexPivots(Distance *_dist, string _path_ind, string _class_object); //IN_DISK
    virtual ~IndexPivots();   
    
     void set_params(int _n_pivots);

    void display();
    bool insert(Object *q);    
    int size();    
    vector<PairObjectDistance> knn(Object* q, int k);     
    void clean_index();    
    void insert_all(vector<Object*> _objects);
    // retorna el numero de buckets creados
    int num_buckets(){return 0;};
    // retorna el numero de accesos y distancias computadas zde la ultima consulta
    int num_access(); int num_dists();  int num_mindists();
    void reset_counters();

protected:      

    // escribe el object en disco/RAM
    virtual void write_object(Object* p);   
    // lee los objetos desde disco/RAM
    Object* read_object(Object* p);    

    void select_pivots_semirand(vector<Object*> _objects, int seed = 1);
    void select_pivots_rand(vector<Object*> _objects, int seed = 1);
};

#endif // IndexLINEAL_H
