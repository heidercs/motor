#ifndef INDEXSNAKEHCSTREAM_H
#define INDEXSNAKEHCSTREAM_H

#include "../indexing/pairs.h"
#include "../indexing/index.h"
#include "stream.h"   
#include "detectevent.h"  
#include <list> 

using namespace std::tr1;

// snake table using FIFO Compact 
class IndexSnakeHCStream: public Index, public DetectEvent
{
private:
    //object collection
    vector<Object*> objects;
    //Metric Measure
    Distance* dist;

    //pointer to data stream
    Stream* stream;
    // size of the window
    int window;
    // type of repleacement a cell. 
    int type_rep;    
	//snake table
	vector<PairObjectDistance> *table;	
	// number of queries
	int n_query;	

	// the most discord subsequence
    PairPosDistance far_nn;

    long int CONT_BUCKETS, CONT_DIST, CONT_MINDISTS;

public:
    IndexSnakeHCStream(Stream* _stream, Distance* _dist);
    IndexSnakeHCStream(Stream* _stream, Distance* _dist,  string _path_ind);

    virtual ~IndexSnakeHCStream();

    void set_params(int _n_query);

public: //overwrites
	// actualiza el vecino mas cercano de todas las subsecuencias: llamar la primera vez
	void update_all(int _window);
		
	// obtiene la subsecuencia mas discordiante	
    PairPosDistance run_the_most_discord();
	    
    //retornar los K subsecuencias mas discordiantes y no solapadas    
    std::vector<PairPosDistance> get_topk_discords(int k);
    
    // verificar si es la mas discordiante
    bool is_the_most_discord(int s_init);

    // actualiza la estructura interna para los ultimos datos insertados en el stream
	bool fix_structure();
	
	// computa el vecino mas cercano de s_init. Actualiza el best_nn para todas las subsecuencias
    // s_init: posicion inicial de la subsecuencia de consulta. Asume que el vector best_nn esta completo 
    PairPosDistance get_nn(int s_init, double far_dist=0);

    int get_window() { return window; }
	
	// llama a la funcion de la clase padre para limpiar el indice
    void clean_index();

    void clean_table();

    void reset_counters();

    int num_buckets();
       // retorna el numero de accesos y distancias computadas de la ultima consulta
    int num_access(); int num_dists();  int num_mindists();

    bool insert(Object* p);

    void display();

    int size();

    vector<PairObjectDistance> knn(Object* q, int k) { return vector<PairObjectDistance>(); }
 
private:		
   
     
    // compute the nearest neigbor by snake table technique using the repleacement heuristics FIFO Sparce
    PairPosDistance get_nn(Object *Q, double far_dist);  
		
	
	// verifica si existe overlap entre la subsecuencia "par" y alguna otra de la "cola"
	bool overlap(int pos, std::vector<PairPosDistance> res);   
		
    //inicialize the table
    PairPosDistance inicialize_table(int index_q);

    Object* read_object(Object* p);
};

#endif // INDICELINEAL_H
