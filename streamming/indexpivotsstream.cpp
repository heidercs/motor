#include <float.h>
#include <math.h>
#include <stdlib.h> 
#include <limits.h> 

#include "indexpivotsstream.h"
#include "../util/functions.h"


IndexPivotsStream::IndexPivotsStream(Stream* _stream, Distance* _dist):
    IndexPivots(_dist, _stream->get_class_object())
{    
	this->stream=_stream;	
	this->window=0;	
    this->far_nn.dist = 0;
    this->far_nn.pos = -1;
}

IndexPivotsStream::IndexPivotsStream(Stream* _stream, Distance* _dist, string _path_ind):
    IndexPivots(_dist, _path_ind, _stream->get_class_object())
{    
	this->stream=_stream;	
	this->window=0;	
    this->far_nn.dist = 0;
    this->far_nn.pos = -1;
}

IndexPivotsStream::~IndexPivotsStream() { }

 
bool IndexPivotsStream::is_the_most_discord(int s_init)
{	
    if(far_nn.dist == 0)
    {
        run_the_most_discord();
    }
    else
    {
        int n_ant=size();
        fix_structure();
        int n_now=size();

        PairPosDistance nn;

        for(int i=n_ant; i < n_now; i++)
        {
            nn =get_nn(i, far_nn.dist);

            if(nn.pos != -1 && nn.dist >= far_nn.dist)
            {
                far_nn.dist = nn.dist;
                far_nn.pos = i;
            }            
        }
    }

    return far_nn.pos == s_init;
}


void IndexPivotsStream::update_all(int _window)
{	
    if(_window == this->window)
    {
        fix_structure();
        return;
    }

    clean_index();

	this->window=_window;	
    int n_len=stream->get_n();

    
    std::vector<Object*> objects;
    for(int i=0;i< n_len - window + 1;i++)
    {
        objects.push_back(stream->get_sub(i, window));
    }

    insert_all(objects);//bulk    
}


bool IndexPivotsStream::fix_structure()
{
	bool fixed=false;
	int n_len=stream->get_n();
	int n=size();
		
	// verificar si el hash contiene a todas las subsecuencias del stream
	if(n < n_len - window + 1)
	{
		int from=n;
		
		for(int i=from; i <n_len - window + 1; i++) 
		{
			insert(stream->get_sub(i, window));		
		}
		
		fixed=true;							
	}		
				
	return fixed;	
}   

PairPosDistance IndexPivotsStream::run_the_most_discord()
{	
	vector<Object*> out_elems;
	PairPosDistance nn;
    Object* Q;

	// comprobar el indice 
	fix_structure();
	
    // busca la mas discordiante primero en los pivots
    PairPosDistance far = get_discord_pivot();
    
    // iterador externo 
	// ordena los objects de mayor a menor respecto la distancia a todos los pivots
    vector<Object*> ind_outerloop = generate_outerloop();
	//vector<Object*> ind_outerloop = this->objects;
   
    for(Object* q : ind_outerloop)
    {
		// ordena los buckets respecto a MINDIST entre la palabra  w_out y los keys de cada bucket		
        Q = read_object(q); 

        nn = get_nn(Q, far.dist); //Exact search

        if(nn.pos != -1 && nn.dist >= far.dist)         // update far
        {
            far.dist = nn.dist;
            far.pos =  Q->get_id_file();
        }         

		if(type_store == IN_DISK)
            delete Q; //free memory
    }

    this->far_nn = far;
	
	return far;
}

/** Version ineficiente, hay que mejorar **/
vector<PairPosDistance>  IndexPivotsStream::get_topk_discords(int k)
{
	return vector<PairPosDistance>();
}


PairPosDistance IndexPivotsStream::get_discord_pivot()
{
    PairPosDistance far(-1, 0);
    double d, nn_dist;
    int i=0,j=0;

    for(Object* Q : pivots)
    { 
        nn_dist = DBL_MAX;

        //buscar en los demas pivots
        for(Object* P : pivots)
        {
            if(abs(P->get_id_file() - Q->get_id_file()) >= window) // non-self match
            {
                d = dist->d(P, Q); CONT_DIST++;

                if(d <= nn_dist)
                {
                    nn_dist = d;
                }
            }
        }

        //buscar en la matriz de similitudes
        i = 0;
        for(Object *P : objects)
        {
            if(abs(P->get_id_file() - Q->get_id_file()) >= window) // non-self match
            {
                d = matrix[i][j];

                if(d <= nn_dist)
                {
                    nn_dist = d;
                }
            }
            i++;
        }

        if(nn_dist >= far.dist)
        {
            far.dist = nn_dist;
            far.pos =  Q->get_id_file();
        }

        j++;
    }

    return far;
}

PairPosDistance IndexPivotsStream::get_nn(int s_init, double far_dist)
{
    Object* Q = stream->get_sub(s_init, window);

    PairPosDistance nn = get_nn(Q,  far_dist);

    delete Q;

    return nn;
}

PairPosDistance IndexPivotsStream::get_nn(Object* Q, double far_dist)
{
    PairPosDistance  nn;
    Object* P;
    double d, LB;
    bool break_soon = false;

    nn.dist=DBL_MAX;
    nn.pos=-1;

    // search on pivots 
    double* matrix_q=new double[n_pivots];  
    for(int j=0; j<n_pivots; j++)
    {
        P = pivots[j];      

        if(abs(P->get_id_file() - Q->get_id_file()) >= window) // non-self match
        {
            d = dist->d(P, Q); CONT_DIST++;
            matrix_q[j] = d;     

            if(d <= nn.dist)
            {
                nn.dist = d;
                nn.pos = P->get_id_file();
            }
        }
        else
            matrix_q[j] = -1; 
    }

    //preparacion de la cola de prioridad (en este caso, set ordenado)
    vector<PairObjectDistance> heap_min;
    for(int i=0; i<(int)this->objects.size(); i++) 
    {
        LB = DBL_MAX;
        for(int j=0; j<n_pivots; j++)
        {
            if(matrix_q[j] == -1) continue;
            d = fabs(matrix_q[j] - this->matrix[i][j]);         
            if(d < LB)
            {
                LB = d;
            }           
        }
        //d(q, x) >= LB
        //si LB > r => d(q, x) > r => descartado (y parar)
        
        if(LB == DBL_MAX)  LB = 0;
        
        heap_min.push_back(PairObjectDistance(this->objects[i], LB));
    }
    std::sort(heap_min.begin(), heap_min.end());


    //revision de la cola de prioridad
    for(auto &par : heap_min)
    {
        if(par.dist > nn.dist)    break;
        
        P =  read_object(par.object);
        
        if(abs(P->get_id_file() - Q->get_id_file()) >= window) // non-self match
        {
            d = dist->d(P, Q); CONT_DIST++;
   
            if(d <= nn.dist)
            {
                nn.dist = d;
                nn.pos = P->get_id_file();
            }
            else if(type_store == IN_DISK)
                delete P;
                
            break_soon = (nn.pos != -1 && nn.dist < far_dist); // break soon
        }

        if(break_soon)     break;
    }

    vector<PairObjectDistance>().swap(heap_min);
    delete[] matrix_q;
    
    return nn;
}


vector<Object*> IndexPivotsStream::generate_outerloop()
{		
    int N = (int)objects.size();
    vector<PairPosDistance> score_object(N);
    vector<Object*> outer(N);
    int i=0, j=0;
    double sum = 0;

	for(i=0; i<N; i++)  
	{        
        sum = 0;
        for (j = 0; j < n_pivots; ++j)
        {
            sum += matrix[i][j];
        }
        score_object[i] =  PairPosDistance(i, sum);       
	}

    //order descendingly
    std::sort(score_object.begin(), score_object.end(), 
              [](PairPosDistance left, PairPosDistance right)
                 {return left.dist > right.dist;}
              );

    //copy to outer vector
    i=0;
    for(PairPosDistance &par : score_object)
    {
        outer[i++] = objects[par.pos];
    }

    vector<PairPosDistance>().swap(score_object);
        
	return outer;    
}


