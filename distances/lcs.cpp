#include "lcs.h"
#include "../objects/vector.h"
#include "../objects/matrix.h"
#include <iostream>
#include <math.h>

LCS::LCS(int _type_lcs, double _fraccion_b, double _epsilon)
{
    this->type_lcs=_type_lcs;
    this->fraccion_b=0.25;
    this->epsilon=0.05;
    
    if(type_lcs==1)
	{
		this->actualizar_b=true;
		if(_fraccion_b>0 && _fraccion_b<=1.0)
			this->fraccion_b=_fraccion_b;
		if(_epsilon>0 && _epsilon<=1.0)
			this->epsilon=_epsilon;			
	}
	else
		this->actualizar_b=false;
	
}

LCS::~LCS()
{
}

bool LCS::cumple_condicion(Object* p, Object* q, int i, int j)
{
    bool correct=true;
    double d;
    
    //verificar que los indices esten dentro de la banda    
    correct=fabs(i-j)<=band;

    if(correct && dynamic_cast<Vector*>(p) && dynamic_cast<Vector*>(q))
    {                
        if(type_lcs==0)//si es LCS sobre secuencias discretas, verificar igualdad 
			correct=((Vector*)p)->get(i)  ==  ((Vector*)q)->get(j);		
		else
		{//diferencia absoluta entre elemeto i,j de p,q
			d=fabs(((Vector*)p)->get(i)  -  ((Vector*)q)->get(j));			
			correct=d<epsilon;
		}	        
		
    }
    else if(correct && dynamic_cast<Matrix*>(p) && dynamic_cast<Matrix*>(q))
    {        
        double* pi=((Matrix*)p)->get(i);
        double* qj=((Matrix*)q)->get(j);
        //a partir de la dimension mas pequenia
        int n=std::min(((Matrix*)p)->get_D(),((Matrix*)q)->get_D());
        d=0;
		//Aplicar distance euclidiana, mas estable que aplicar umbral a cada elemento del vector
        for(int i=0; i<n; i++)
            d+=pow(pi[i]-qj[i],2);
        d=sqrt(d);
        correct=d<epsilon;        
    }    
        
    return correct;
}

double LCS::d(Object *p,Object *q)
{
    if(!(dynamic_cast<Vector*>(p) && dynamic_cast<Vector*>(q)) &&
            !(dynamic_cast<Matrix*>(p) && dynamic_cast<Matrix*>(q)) )
    {        
        std::cerr<<"*** Esta distance no soporta este tipo de object *** \n";
        return -1;
    }

    ////////////////////////////////////////////////////////////////////////////////////
    
    double d;

    int nq=q->get_n();
    int np=p->get_n();

    int i,j;
    
    int** m=new int*[nq+1];
    for(i=0; i<nq+1; i++)
        m[i]=new int[np+1];
	
	if(type_lcs==1 && actualizar_b)
    {
        band=round(fraccion_b*np);
    }
    else 
		band=std::max(nq,np);
		

    //Se inicializa la primer fila y la primera columna con 0    
    for(i=0; i<=nq; i++)
    {
		m[i][0]=0;        
    }
    
    for(j=0; j<=np; j++)
	{
		m[0][j]=0;
	}
  

    for(i=1; i<=nq; i++)
    {  
        for(j=1; j<=np; j++)
        {			
            if(cumple_condicion(q,p,i-1,j-1) ) 
            {
				m[i][j]=m[i-1][j-1]+1;
			}   
			else
			{
				m[i][j]=std::max(m[i-1][j], m[i][j-1]);
			}	
			std::cout<<m[i][j]<<"\t";		
        }        
        std::cout<<"\n";
    }
    d=m[nq][np];
           

    for(i=0; i<=nq; i++)
        delete[] m[i];
    delete[] m;
        

    return 1.0 - d/(1.0*std::min(nq,np)); //transformar LCSS (similitud) a distance (dis-similitud)
}


