#include <math.h>
#include "pip.h"


/**
 * Return the position of the farthest point regarding the line 
 * between the position left and right of the time series
 * @param  p          [data vector]
 * @param  left       [position left]
 * @param  right      [position right]
 * @return            [posicion PIP]
 */
int Pip::found_pip(double* p, int left, int right)
{
    int pc;
    double dist,distmay;
    double fact = (p[right] - p[left])/(right-left);

    distmay=0;
    pc=(left + right)/2.0;

    for(int i=left+1; i<=right-1; i++)
    {
        dist = p[left] + (i-left) * fact  - p[i];                
        dist = fabs(dist);

        if(dist > distmay)
        {
            distmay=dist;
            pc=i;
        }

    }

    return pc;
}


/**
 * Compute recursively new_n perceptual important points from p
 */
std::vector<Pip> Pip::get_rec_pip(double* p, int N, int new_n)
{
	std::vector<Pip> coefs(new_n);
    int i;
	
	if(new_n >= N)
	{
		for (i = 0; i < N; i++)		
			coefs[i] = Pip(i, p[i]);

        return coefs;
	}
	
    List<Pip> *l=new List<Pip>();
    l->add_heap_min(Pip(0,p[0]));
    get_rec_pip(p,l,0,N-1,new_n-2);
    l->add_heap_min(Pip(N-1,p[N-1]));
    
    
    l->reset();
    i = 0;
    while(l->has_next())
		coefs[i++] = (l->next());	
	
	delete l;	
   
    return coefs;
}

void Pip::get_rec_pip(double* p, List<Pip>* l, int izq, int der, int n)
{
    if(n > 0 && der > izq)
    {
        int n_izq, pc;
        pc = found_pip(p,izq,der);
        l->add_heap_min(Pip(pc, p[pc]));
        n--;
        
        // porpotionally distributing on both sides
        n_izq=round((pc-izq-1)*n/((der-izq-1)*1.0));

        // find the rest 
        get_rec_pip(p,l,izq,pc,n_izq);
        get_rec_pip(p,l,pc,der,n-n_izq);
    }
}


/**
 * Compute new_n perceptual important points by segmentation strategy 
 */
double* Pip::get_seg_pip(double* p, int N, int new_n)
{
    double *coefs = new double[new_n];

    if(new_n >= N)
    {
        for (int i = 0; i < N; i++)     
            coefs[i] = p[i];

        return coefs;
    }

    double w = 1.0*N/new_n;
    int pos, izq = 0, der = 0;

    for (int i = 0; i < new_n; ++i)
    {
        der = std::min((int) round((i+1) * w), N) - 1;

        pos = found_pip(p, izq, der);

        coefs[i] = p[pos];    

        izq = der;    
    }

    return coefs;
}

double* Pip::get_rec_pip_double(double* p, int N, int new_n)
{
    std::vector<Pip> vec = get_rec_pip(p, N, new_n);

    double *coefs   = new double[(int)vec.size()];
    
    for(int i=0; i<(int)vec.size(); i++)
        coefs[i] = vec[i].value;

    return coefs;
}

/**
 * return all points ordering by their pip priority. 
 */
std::vector<Pip> Pip::get_all_pip(double* p, int N)
{    
    std::vector<Pip> coefs(N);

    List<Pip> *l=new List<Pip>();
    l->add_heap_min(Pip(0,p[0]));    
    l->add_heap_min(Pip(N-1,p[N-1]));  

    int pos, der, izq, ant, j, i = 0;
    coefs[i++] = Pip(0,p[0]);
    coefs[i++] = Pip(N-1,p[N-1]);


    while(i < N) 
    {    
        ant = i - 1;
        l->reset();   
        izq = l->next().index;     

        while(l->has_next())
        {        
            der = l->next().index;

            if(der > izq + 1)
            {    
                pos = found_pip(p, izq, der);
                coefs[i++] = Pip(pos, p[pos]);
            }

            izq = der;
        }

        for(j = ant; j < i; j++)
            l->add_heap_min(coefs[j]);
    }

    delete l;

    return coefs;
}   

std::ostream& operator<<(std::ostream& output, Pip p )
{
    output <<p.index + 1<<":"<<p.value<<", ";
    return output;
}

bool operator >(Pip p1, Pip p2)
{
    return p1.index > p2.index;
}

bool operator <(Pip p1, Pip p2)
{
    return p1.index < p2.index;
}
