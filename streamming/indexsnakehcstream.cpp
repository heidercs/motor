#include "indexsnakehcstream.h"
#include <float.h>
#include <math.h>
#include <stdlib.h> 


IndexSnakeHCStream::IndexSnakeHCStream(Stream* _stream, Distance* _dist)
{
    path_ind = "";
    type_store = IN_RAM;
    class_object = _stream->get_class_object();
    dist=_dist;

	this->stream=_stream;		
    this->n_query = 5;
    this->table = NULL;    
}

IndexSnakeHCStream::IndexSnakeHCStream(Stream* _stream, Distance* _dist,  string _path_ind)
{
    path_ind = _path_ind;
    type_store = IN_DISK;
    class_object = _stream->get_class_object();
    dist=_dist;

	this->stream=_stream;		
    this->n_query = 5;
    this->table = NULL;    
}

IndexSnakeHCStream::~IndexSnakeHCStream()
{	    
}

void IndexSnakeHCStream::set_params(int _n_query)
{
    this->n_query = _n_query;
}

void IndexSnakeHCStream::update_all(int _window)
{
    if(_window == this->window)
    {
        fix_structure();
        return;
    }

    clean_index();

    this->window=_window;
    int n_len=stream->get_n();
    for(int i=0;i<n_len - window + 1; i++)
    {
        insert(stream->get_sub(i, window));
    }
}
 
bool IndexSnakeHCStream::is_the_most_discord(int s_init)
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
            nn = get_nn(i, far_nn.dist);

            if(nn.pos != -1 && nn.dist >= far_nn.dist)
            {
                far_nn.dist = nn.dist;
                far_nn.pos = i;
            }
        }
    }

    return far_nn.pos == s_init;
}

bool IndexSnakeHCStream::fix_structure()
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

PairPosDistance IndexSnakeHCStream::inicialize_table(int index_q)
{
    clean_table();    

    Object *P;
    double d = 0;
    PairPosDistance nn(-1, DBL_MAX);
    Object* Q = read_object(objects[index_q]);    

    int n = size();    
    this->table = new vector<PairObjectDistance>[n];

    for(int j = 0; j< n ;j++)
    {        
        if(abs(j - index_q) >= window) // non-self match
        {
            P = read_object(objects[j]);
            d = dist->d(Q, P); CONT_DIST++;

            if(d < nn.dist)
            {
                nn.dist = d;
                nn.pos = j;
            }

            table[j].push_back(PairObjectDistance(Q, d));
        }
    }    

    return nn;
}


PairPosDistance IndexSnakeHCStream::run_the_most_discord()
{	
	fix_structure();//verificar que todas las subsecuencias esten completas
	
    PairPosDistance nn, far(-1, 0);

    Object* Q = NULL;    

    far = inicialize_table(0);  //inicialize with the first element 

    // compute with the rest of elements
    for(int j = 1; j< size() ;j++)
	{
        Q = read_object(objects[j]);
        
        nn = get_nn(Q, far.dist);             
        							
        if(nn.pos != -1 && nn.dist >= far.dist)
		{
            far.dist = nn.dist;
			far.pos = j;
		}

        if(type_store == IN_DISK)
            delete Q; //free memory         
    }

    clean_table();
		
	return far;
}


std::vector<PairPosDistance> IndexSnakeHCStream::get_topk_discords(int k)
{
		
	// guarda las subsecuencia en orden de la mas discordiante: vecino mas cercano (nn) con mayor distancia    
    vector<PairPosDistance> all_best_nn;
    PairPosDistance nn,inn_max;
    
    Object* Q = NULL;
    nn = inicialize_table(0);
    nn.pos = 0;
    all_best_nn.push_back(nn);  
    
    for(int j = 1; j< size() ;j++)
    {
        Q = read_object(objects[j]);
       
        nn = get_nn(Q, 0);

        nn.pos = Q->get_id_file();                              
        all_best_nn.push_back(nn);                                         
        
        if(type_store == IN_DISK)
            delete Q; //free memory    
    }
    
            
    // buscar los topk sin solapamiento
    vector<PairPosDistance>  res;   
    for(int j=0;j <k; j++)
    {
        inn_max=PairPosDistance(-1,-1);
        for(auto& nn : all_best_nn)
        {
            if(nn.dist > inn_max.dist)
            {
                if(!overlap(nn.pos, res))
                {
                    inn_max.dist = nn.dist;
                    inn_max.pos = nn.pos;
                }
            }   
        }       
        res.push_back(inn_max);     
    }
    vector<PairPosDistance>().swap(all_best_nn);

    clean_table();
    
    return res;
}


PairPosDistance IndexSnakeHCStream::get_nn(int s_init, double far_dist)
{
    Object* Q = stream->get_sub(s_init, window);

    PairPosDistance nn;

    nn =  get_nn(Q, far_dist);

    return nn;
}


PairPosDistance IndexSnakeHCStream::get_nn(Object* Q, double far_dist)
{               
    int n=size(), pos;
    double d=0, lb;
    Object *P;
    PairPosDistance nn(-1, DBL_MAX);
    bool  break_soon = false;

    int index_q = Q->get_id_file(),  i, j;    
    
    double dist_min;  int pos_min;
    
    // cache of queries 
    double* cache = new double[n];
    for(j =0; j<n; j++)    cache[j] = -1;

    for(j = 0; j<n && !break_soon; j++)
    {       
        if(abs(j - index_q) >= window) // non-self match
        {               

            pos_min=0; dist_min = DBL_MAX; lb = 0;
            for(i=0; i<(int)table[j].size())
            {
                pos = table[j][i].object->get_id_file();

                if(cache[pos] == -1)
                {
                    d = dist->d(Q,  table[j][i].object); CONT_DIST++;
                    cache[pos] = d;
                }
                else
                    d = cache[pos];

                if(d < dist_min) {dist_min=d; pos_min = i;}                
                
                lb = fabs(d - table[j][i].dist);

                if(lb > nn.dist) break;
            }

            if(lb <= nn.dist)
            {                       
                P = read_object(objects[j]);                
                d = dist->d(Q, P); CONT_DIST++;
                cache[j] = d;

                if(d < nn.dist)
                {
                    nn.dist =d;
                    nn.pos  =j;
                }

                break_soon = (nn.dist!=-1 && nn.dist < far_dist);
                
                if(type_store == IN_DISK)
                    delete P;      

                //Highest/Compact
                if((int)table[j].size() == this->n_query)                
                    table[j][pos_min] = PairObjectDistance(Q, d);
                else
                    table[j].push_back(PairObjectDistance(Q, d));
            }            
        }               
    }

    delete[] cache;
 
    return nn;
}


bool IndexSnakeHCStream::overlap(int pos, std::vector<PairPosDistance> res)
{	
	bool has_overlap=false;	
	PairPosDistance par;
		
	for(int i=0;i<(int)res.size() && !has_overlap;i++)
	{
		par=res[i];
		has_overlap=abs(pos - par.pos) < window;
	}
	
	return has_overlap;
}



//////////////////////////////////////////////////////////
bool IndexSnakeHCStream::insert(Object* p)
{
    if(class_object.empty())
        class_object = p->name();
    else if(class_object.compare(p->name()) != 0)
    {
        std::cout<<"** Incompatible Value **\n";
        return false;
    }   
    
    if(type_store == IN_RAM)
        objects.push_back(p);   
    else // IN_DISK
    {
        long pos_in_file=write(path_ind+"/ind_snake", p);      
        objects.push_back(new ItemDisk(pos_in_file));           
    }

    return true;
}   

Object* IndexSnakeHCStream::read_object(Object* p)
{   
    if(type_store == IN_RAM)
        return p;
    else
    {       
        return read(path_ind + "/ind_snake", class_object, ((ItemDisk*)p)->get_position_file());
    }
}

void IndexSnakeHCStream::clean_index()
{
    vector<Object*>().swap(this->objects);
    if(type_store == IN_DISK)
    {
        string file=path_ind + "/ind_snake";
        remove(file.c_str());
    }
    
}

void IndexSnakeHCStream::clean_table()
{
    if(table != NULL)
    {       
        delete[] table;        
    }    
    
}

void IndexSnakeHCStream::display()
{
    Object* object;
    for(Object* p : objects)
    {
        object =  read_object(p);
        object->display();
        std::cout<<std::endl;
        if(type_store == IN_DISK)
            delete object;
    }
}

int IndexSnakeHCStream::size()
{
    return (int)objects.size();
}

void IndexSnakeHCStream::reset_counters()
{
    CONT_BUCKETS=CONT_DIST=0;
    CONT_MINDISTS=0;
}

int IndexSnakeHCStream::num_buckets()
{
    return 0;
}


int IndexSnakeHCStream::num_access()
{
    return CONT_BUCKETS;
}

int IndexSnakeHCStream::num_dists()
{
    return CONT_DIST;
}

int IndexSnakeHCStream::num_mindists()
{
    return CONT_MINDISTS;
}
