#ifndef INDEXPIVOTSSTREAM_H
#define INDEXPIVOTSSTREAM_H

#include "../distances/lineal.h"   
#include "../indexing/pairs.h"
#include "../indexing/indexpivots.h"
#include "../util/list.h"   
#include "streamuni.h"  
#include "streammulti.h"  
#include "detectevent.h"   


class IndexPivotsStream: public IndexPivots, public DetectEvent
{

private:   	
    //pointer to data stream
    Stream* stream;    
    // size of the window
    int window;

    // the most discord subsequence
    PairPosDistance far_nn;
    
public:
    IndexPivotsStream(Stream* _stream, Distance* _dist);
    IndexPivotsStream(Stream* _stream, Distance* _dist, string _path_ind);
	virtual ~IndexPivotsStream();
   
	// actualiza el vecino mas cercano de todas las subsecuencias: llamar la primera vez
	void update_all(int _window);
		
	// obtiene la subsecuencia mas discordiante	
    PairPosDistance run_the_most_discord();
	    
    //retornar los K subsecuencias mas discordiantes y no solapadas        
    vector<PairPosDistance>  get_topk_discords(int k);
    
    // verificar si es la mas discordiante
    bool is_the_most_discord(int s_init);
    
    // actualiza la estructura interna para los ultimos datos insertados en el stream
	bool fix_structure(); 
	
	// computa el vecino mas cercano de s_init. aplicando  non-self match
    PairPosDistance get_nn(int s_init, double far_dist=0);
    PairPosDistance get_nn(Object* Q, double far_dist=0);


    int get_window() { return window; }
    	  	
 
private:
	
	// genera un vector de todos los buckets, 
    // tal que el bucket con la minima cantidad de elementos esta al inicio del vector
    vector<Object*> generate_outerloop();

    //search the discord into the set of pivots
    PairPosDistance get_discord_pivot();
    
};

#endif // INDICELINEAL_H
