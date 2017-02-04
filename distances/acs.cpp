#include "acs.h"
#include "../objects/vector.h"
#include "../objects/matrix.h"
#include <iostream>
#include <math.h>

ACS::ACS(int _type_acs, double _epsilon)
{
    this->type_acs=_type_acs;    
    this->epsilon=0.05;
    
    if(type_acs==1)
	{	
		if(_epsilon>0 && _epsilon<1.0)
			this->epsilon=_epsilon;			
	}
}

ACS::~ACS()
{
}

bool ACS::cumple_condicion(Object* p, Object* q, int i, int j)
{
    double d;
    bool correct=true;

    if(dynamic_cast<Vector*>(p) && dynamic_cast<Vector*>(q))
    {                
        if(type_acs==0)//si es LCS sobre secuencias discretas, verificar igualdad 		
			correct=((Vector*)p)->get(i)  ==  ((Vector*)q)->get(j);
				
		else//sobre numeros reales, aplicar umbral epsilon
		{//diferencia absoluta entre elemeto i,j de p,q
			d=fabs(((Vector*)p)->get(i)  -  ((Vector*)q)->get(j));
			correct=d<epsilon;				
		}	        
    }
    else if(dynamic_cast<Matrix*>(p) && dynamic_cast<Matrix*>(q))
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

double ACS::d(Object *p,Object *q)
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
    
    double** m=new double*[nq+1];
    for(i=0; i<nq+1; i++)
        m[i]=new double[np+1];
		
    //Se inicializa la primer fila y la primera columna con 0    
    for(i=0; i<=nq; i++)
    {
		m[i][0]=1;        
    }
    
    for(j=0; j<=np; j++)
	{
		m[0][j]=1;
	}  

    for(i=1; i<=nq; i++)
    {  
        for(j=1; j<=np; j++)
        {	
            
            if(cumple_condicion(q,p,i-1,j-1) ) 
            {
				m[i][j]=2*m[i-1][j-1];
			}   
			else
			{
				m[i][j]=m[i-1][j] + m[i][j-1] - m[i-1][j-1];
			}			
        }        
    }
    d=m[nq][np];
        

    for(i=0; i<=nq; i++)
        delete[] m[i];
    delete[] m;
    

	//transformar ACS (similitud) a distance (dis-similitud),
	//normalizando con el total de subsecuencias pow(2,n), incluyendo la vacia	 
    return 1 - d/pow(2.0,std::min(nq,np)); 
}


