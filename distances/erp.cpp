#include "erp.h"
#include "../objects/vector.h"
#include "../objects/matrix.h"
#include "../util/functions.h"
#include <iostream>
#include <math.h>

ERP::ERP(double _gap)
{    
	if(_gap >= 0)
		this->gap=_gap;	
	else
		this->gap=0;				
}

ERP::~ERP()
{
	
}

double ERP::dist(Object* p, Object* q, int i, int j)
{
    double d=-1;

    if(dynamic_cast<Vector*>(p) && dynamic_cast<Vector*>(q))
    {        
		//Manhattan
		d=fabs(((Vector*)p)->get(i)  -  ((Vector*)q)->get(j));		
    }
    else if(dynamic_cast<Matrix*>(p) && dynamic_cast<Matrix*>(q))
    {
        //distance Manhattan entre elemeto i,j de p,q
        double* pi=((Matrix*)p)->get(i);
        double* qj=((Matrix*)q)->get(j);
        //a partir de la dimension mas pequenia
        int n=std::min(((Matrix*)p)->get_D(),((Matrix*)q)->get_D());
        d=0;
        for(int i=0; i<n; i++)
        {
			d+=fabs(pi[i]-qj[i]);
		}		
    }

    return d;
}

double ERP::dist(Object* p, int i, int gap)
{
    double d=-1;

    if(dynamic_cast<Vector*>(p))
    {        
		//Manhattan
		d=fabs(((Vector*)p)->get(i)  -  gap);		
    }
    else if(dynamic_cast<Matrix*>(p))
    {
        //distance Manhattan entre elemeto i,j de p,q
        double* pi=((Matrix*)p)->get(i);        
        //a partir de la dimension mas pequenia
        int n=((Matrix*)p)->get_D();
        d=0;
        for(int i=0; i<n; i++)
        {
			d+=fabs(pi[i]-gap);
		}		
    }

    return d;
}


double ERP::d(Object *p,Object *q)
{
    if(!(dynamic_cast<Vector*>(p) && dynamic_cast<Vector*>(q)) &&
            !(dynamic_cast<Matrix*>(p) && dynamic_cast<Matrix*>(q)) )
    {        
        std::cerr<<"*** Esta distance no soporta este tipo de object *** \n";
        return -1;
    }

    ////////////////////////////////////////////////////////////////////////////////////
    
    double d,costo;

    int nq=q->get_n();
    int np=p->get_n();

    int i,j;
    
    double** m=new double*[nq+1];
    for(i=0; i<nq+1; i++)
        m[i]=new double[np+1];
        
	m[0][0]=0;	
    //Se inicializa la primer fila y la primera columna  
    costo = 0;
    for(i=1; i<=nq; i++)
    {
		costo+=dist(q,i-1,this->gap);		
		m[i][0]=costo;        
    }
    
    costo=0; 
    for(j=1; j<=np; j++)
	{
		costo+=dist(p,j-1,this->gap);		
		m[0][j]=costo;        
	}  

    for(i=1; i<=nq; i++)
    {  
        for(j=1; j<=np; j++)
        {	
            m[i][j]=min3(	m[i-1][j-1] + dist(p,q,j-1, i-1),
							m[ i ][j-1] + dist(p,j-1, this->gap),
							m[i-1][ j ] + dist(q,i-1, this->gap) );
						
        }        
    }
    d=m[nq][np];
        

    for(i=0; i<=nq; i++)
        delete[] m[i];
    delete[] m;
    
    return d; 
}

double ERP::lower_bound(Object *p,Object *q)
{
	if(!(dynamic_cast<Vector*>(p) && dynamic_cast<Vector*>(q)) &&
            !(dynamic_cast<Matrix*>(p) && dynamic_cast<Matrix*>(q)) )
    {        
        std::cerr<<"*** Esta distance no soporta este tipo de object *** \n";
        return -1;
    }
    
    double sumq=0;
    for(int i=0;i<q->get_n();i++)
		sumq+=dist(q,i,this->gap); 
		
	double sump=0;
    for(int i=0;i<p->get_n();i++)
		sump+=dist(p,i,this->gap);
		
	return fabs(sumq - sump); //diferencia de areas
}


