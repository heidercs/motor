#include "banddtw.h"
#include "../objects/vector.h"
#include "../objects/matrix.h"
#include "../objects/vectorpair.h"
#include "../objects/matrixpair.h"
#include "../objects/signature.h"
#include "../descriptors/paa.h"
#include "../util/functions.h"
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <float.h>

BandDTW::BandDTW(int _type_band, double _fraction_r, int _type_lineal)
{
    type_band=_type_band;    
    dist_lineal = new Lineal(_type_lineal);
    
    //solo aplica a DTW con banda: type_band=1,2
        r_max=0; 
        fraction_r=_fraction_r;  //el radio de la banda se actualiza a una fracion de la dimension
}

double BandDTW::d(Object* p,Object* q)
{    
    double d = -1;
    int np = p->get_n();
    int nq = q->get_n();

    if(dynamic_cast<Vector*>(p) && dynamic_cast<Vector*>(q))
    {
        double* _p=((Vector*)p)->get_elems();
        double* _q=((Vector*)q)->get_elems();
        d = dist(_p, _q, np, nq);
    }
    else if(dynamic_cast<Matrix*>(p) && dynamic_cast<Matrix*>(q))
    {        
        int Dp = ((Matrix*)p)->get_D();
        int Dq = ((Matrix*)q)->get_D();
        double** _p=((Matrix*)p)->get_elems();
        double** _q=((Matrix*)q)->get_elems();
        bool inverse = ((Matrix*)p)->is_inverse() && ((Matrix*)q)->is_inverse();

        if(inverse)
        // for us, this option indicate that
        // the distribution of each variable is different
        // therefore, we combine the matches between variables
        {
            int n=std::min(np, nq);

            d = 1;
            for(int i=0;i<n;i++)
                d *= 1 + dist(_p[i], _q[i], Dp, Dq);
            d = d - 1;
        }
        // in the normal case, we apply the multidimensional distance
        else
        {
            int D = std::min(Dp, Dq);
            d = dist(_p, _q, np, nq, D);
        }

    }
    else if(dynamic_cast<Signature*>(p) && dynamic_cast<Signature*>(q))
    {
        int D = std::min(((Signature*)p)->get_D(), ((Signature*)q)->get_D()); 
        double **_p =   ((Signature*)p)->get_centroids();            
        double **_q =   ((Signature*)q)->get_centroids();        
        bool by_levels = ((Signature*)q)->is_by_levels() && ((Signature*)p)->is_by_levels();

        if(by_levels) 
            d = dist_by_levels(_p, _q, np, nq, D);   
        else        
        {
            double *_wp =  ((Signature*)p)->get_weights();
            double *_wq =  ((Signature*)q)->get_weights();

            if(_wp != NULL && _wq != NULL)
                d = dist(_p, _q, _wp, _wq, np, nq, D);           
            else
                d = dist(_p, _q, np, nq, D);
        }
       
    }
    else if(dynamic_cast<VectorPair*>(p) && dynamic_cast<VectorPair*>(q))
    {
        double* _vp=((VectorPair*)p)->get_value();
        double* _vq=((VectorPair*)q)->get_value();
        double* _sp=((VectorPair*)p)->get_slope();
        double* _sq=((VectorPair*)q)->get_slope();
        bool by_levels = ((VectorPair*)q)->is_by_levels() && ((VectorPair*)p)->is_by_levels();

        if(by_levels) 
            d = dist_by_levels(_vp, _vq, _sp, _sq, np, nq);   
        else
            d = dist(_vp, _vq, _sp, _sq, np, nq);  

    }
    else if(dynamic_cast<MatrixPair*>(p) && dynamic_cast<MatrixPair*>(q))
    {
        int D = std::min(((MatrixPair*)p)->get_D(), ((MatrixPair*)q)->get_D());
        double** _vp=((MatrixPair*)p)->get_value();
        double** _vq=((MatrixPair*)q)->get_value();
        double** _sp=((MatrixPair*)p)->get_slope();
        double** _sq=((MatrixPair*)q)->get_slope();
        bool by_levels = ((MatrixPair*)q)->is_by_levels() && ((MatrixPair*)p)->is_by_levels();
        
        if(by_levels) 
            d = dist_by_levels(_vp, _vq, _sp, _sq, np, nq, D);   
        else               
            d = dist(_vp, _vq, _sp, _sq, np, nq, D);        
    }
    else
    {        
        std::cerr<<"*** this distance does not support this type of object *** \n";
    }
   

    if(dist_lineal->get_type_lineal() == 2)
        d = sqrt(d);
    
    return d;    
}

double BandDTW::d(Object *p, Object *q, double th_dist)
{
    if(dist_lineal->get_type_lineal() == 2)
        th_dist *= th_dist;

    double d = -1;

    int np = p->get_n();
    int nq = q->get_n();

    if(dynamic_cast<VectorPair*>(p) && dynamic_cast<VectorPair*>(q))
    {
        bool by_levels = ((VectorPair*)q)->is_by_levels() && ((VectorPair*)p)->is_by_levels();
        if(by_levels) 
        {    
            double* _vp=((VectorPair*)p)->get_value();
            double* _vq=((VectorPair*)q)->get_value();
            double* _sp=((VectorPair*)p)->get_slope();
            double* _sq=((VectorPair*)q)->get_slope();            
            
            d = dist_by_levels(_vp, _vq, _sp, _sq, np, nq, th_dist);              
        }        
            
    }
    else if(dynamic_cast<MatrixPair*>(p) && dynamic_cast<MatrixPair*>(q))
    {
        bool by_levels = ((MatrixPair*)q)->is_by_levels() && ((MatrixPair*)p)->is_by_levels();
        if(by_levels) 
        {
            int D = std::min(((MatrixPair*)p)->get_D(), ((MatrixPair*)q)->get_D());
            double** _vp=((MatrixPair*)p)->get_value();
            double** _vq=((MatrixPair*)q)->get_value();
            double** _sp=((MatrixPair*)p)->get_slope();
            double** _sq=((MatrixPair*)q)->get_slope();

            d = dist_by_levels(_vp, _vq, _sp, _sq, np, nq, D, th_dist); 
        }
    }

    if(d == -1)        
        std::cerr<<" *** The threslhod-distance does not support this type of object *** \n";    


    if(dist_lineal->get_type_lineal() == 2)
        d = sqrt(d);

    return d;
}

/////////////////
void BandDTW::set_band(int _type_band, double _fraction_r)
{
    type_band=_type_band;
    //r_max actualizable
    r_max=0;
    fraction_r=_fraction_r;        
}

void BandDTW::set_band(int _type_band, int _r_max)
{
    type_band=_type_band;
    //r_max fijo
    r_max = _r_max;
    fraction_r = -1;        
}


double BandDTW::dist(double *p, double *q,  int np, int nq)
{
    double d;

    int i,j,b_min,b_max,pos;

    double** m=new double*[nq+1];
    for(i=0; i<nq+1; i++)
        m[i]=new double[np+1];

    //fijar el radio de la band r_max
    fix_r_max(np);

    //Por ahora, uso una inicializacion del interior de la matrix
    //Esto es porque si las bands cambian, hay que limpiar partes de la matrix
    //Esto puede hacerse mejor, sobretodo si las bands no cambian
    for(i=0; i<=nq; i++)
    {
        for(j=0; j<=np; j++)
        {
            m[i][j]=DBL_MAX;
        }
    }
    // Punto de inicio en 0.0
    m[0][0]=0;

    for(i=1; i<=nq; i++)
    {
        pos=((i-1.0)/(nq-1.0))*(np-1.0);
        b_min=bandMin(pos, np);
        b_max=bandMax(pos, np);

        //notar que la band es inclusiva, por lo que sumo 1 en ambos lados
        for(j=1+b_min; j<=1+b_max; j++)
        {
            //insertar, borrar, reemplazar
            m[i][j]=dist_lineal->dist(p[j-1],q[i-1]) 
                    + min3(m[i-1][j], m[i][j-1], m[i-1][j-1]); // simetry property
        }
    }

    d=m[nq][np];//DTW, resultado de la alineacion

    // Time-Normalized Distance Measure (solve the scaling time)    
    // symmetric distance
    d = d/(nq + np); 


    for(i=0; i<=nq; i++)
        delete[] m[i];
    delete[] m;
    
    return d;
}

double BandDTW::dist(double **p, double **q,  int np, int nq, int D)
{
    double d;

    int i,j,b_min,b_max,pos;

    double** m=new double*[nq+1];
    for(i=0; i<nq+1; i++)
        m[i]=new double[np+1];

    //fijar el radio de la band r_max
    fix_r_max(np);

    for(i=0; i<=nq; i++)
    {
        for(j=0; j<=np; j++)
        {
            m[i][j]=DBL_MAX;
        }
    }
    // Punto de inicio en 0.0
    m[0][0]=0;

    for(i=1; i<=nq; i++)
    {
        pos=((i-1.0)/(nq-1.0))*(np-1.0);
        b_min=bandMin(pos, np);
        b_max=bandMax(pos, np);

        //notar que la band es inclusiva, por lo que sumo 1 en ambos lados
        for(j=1+b_min; j<=1+b_max; j++)
        {
            //insertar, borrar, reemplazar
            m[i][j]=dist_lineal->dist(p[j-1], q[i-1], D) 
                    + min3(m[i-1][j], m[i][j-1], m[i-1][j-1]);
        }
    }

    d=m[nq][np];//DTW, resultado de la alineacion

    // Time-Normalized Distance Measure (solve the scaling time)    
    // symmetric distance
    d = d/(nq + np); 

    for(i=0; i<=nq; i++)
        delete[] m[i];
    delete[] m;

    
    return d;
}

double BandDTW::dist(double **p, double **q, double *wp, double *wq, int np, int nq, int D)
{
    double d;

    int i,j,b_min,b_max,pos;

    double** m=new double*[nq+1];
    for(i=0; i<nq+1; i++)
        m[i]=new double[np+1];

    //fijar el radio de la band r_max
    fix_r_max(np);

    for(i=0; i<=nq; i++)
    {
        for(j=0; j<=np; j++)
        {
            m[i][j]=DBL_MAX;
        }
    }
    // Punto de inicio en 0.0
    m[0][0]=0;

    for(i=1; i<=nq; i++)
    {
        pos=((i-1.0)/(nq-1.0))*(np-1.0); // su equivalente en P
        b_min=bandMin(pos, np);
        b_max=bandMax(pos, np);

        //notar que la band es inclusiva, por lo que sumo 1 en ambos lados
        for(j=1+b_min; j<=1+b_max; j++)
        {
            //insertar, borrar, reemplazar
            m[i][j]=dist_lineal->dist(p[j-1], q[i-1], wp[j-1], wq[i-1], D) 
                    + min3(m[i-1][j], m[i][j-1], m[i-1][j-1]);
        }
    }

    d=m[nq][np];//DTW, resultado de la alineacion

    // Time-Normalized Distance Measure (solve the scaling time)    
    // symmetric distance
    d = d/(nq + np); 

    for(i=0; i<=nq; i++)
        delete[] m[i];
    delete[] m;

    return d;
}


double BandDTW::dist_by_levels(double **p, double **q,  int np, int nq, int D, double th_dist)
{
    double d = 0;

    int n = std::min(np, nq);

    int izq, der, max_level = log2(n + 1) - 1;

    for (int l = 0; l <= max_level && d <= th_dist; ++l)        
    {
        izq = pow(2, l);
        der = pow(2, l + 1) - 1;

        d += partial_dist(p, q, izq - 1, der - 1, D);
    }
    
    return d/n;
    //return d;
}

double BandDTW::partial_dist(double **p, double **q,  int izq, int der, int D)
{
    double d;

    int i,j,b_min,b_max,pos;

    int n = der - izq + 1;

    double** m=new double*[n+1];
    for(i=0; i<n + 1; i++)
        m[i] = new double[n+1];

    //fijar el radio de la band r_max
    fix_r_max(n);

    for(i=0; i<=n; i++)
    {
        for(j=0; j<=n; j++)
        {
            m[i][j]=DBL_MAX;
        }
    }
    // Punto de inicio en 0.0
    m[0][0]=0;

    for(i=1; i<=n; i++)
    {
        pos = (i - 1.0); 
        b_min = bandMin(pos, n);
        b_max = bandMax(pos, n);

        //notar que la band es inclusiva, por lo que sumo 1 en ambos lados
        for(j=1+b_min; j<=1+b_max; j++)
        {
            //insertar, borrar, reemplazar
            m[i][j] = dist_lineal->dist(p[izq + j-1], q[izq + i-1], D) 
                    + min3(m[i-1][j], m[i][j-1], m[i-1][j-1]);
        }
    }

    d=m[n][n];//DTW, resultado de la alineacion

    for(i=0; i<=n; i++)
        delete[] m[i];
    delete[] m;

    return d;
}



double BandDTW::dist(double *vp, double *vq, double *sp, double *sq, int np, int nq)
{
    double d;

    int i,j,b_min,b_max,pos;

    double** m=new double*[nq+1];
    for(i=0; i<nq+1; i++)
        m[i]=new double[np+1];

    //fijar el radio de la band r_max
    fix_r_max(np);

    for(i=0; i<=nq; i++)
    {
        for(j=0; j<=np; j++)
        {
            m[i][j]=DBL_MAX;
        }
    }
    // Punto de inicio en 0.0
    m[0][0]=0;

    for(i=1; i<=nq; i++)
    {
        pos=((i-1.0)/(nq-1.0))*(np-1.0); // su equivalente en P
        b_min=bandMin(pos, np);
        b_max=bandMax(pos, np);

        //notar que la band es inclusiva, por lo que sumo 1 en ambos lados
        for(j=1+b_min; j<=1+b_max; j++)
        {
            //insertar, borrar, reemplazar
            m[i][j]=dist_lineal->dist(vp[j-1], vq[i-1], sp[j-1], sq[i-1]) 
                    + min3(m[i-1][j], m[i][j-1], m[i-1][j-1]);
        }
    }

    d=m[nq][np];//DTW, resultado de la alineacion

    // Time-Normalized Distance Measure (solve the scaling time)    
    // symmetric distance
    d = d/(nq + np); 

    for(i=0; i<=nq; i++)
        delete[] m[i];
    delete[] m;

    return d;
}


double BandDTW::dist(double **vp, double **vq, double **sp, double **sq, int np, int nq, int D)
{
    double d;

    int i,j,b_min,b_max,pos;

    double** m=new double*[nq+1];
    for(i=0; i<nq+1; i++)
        m[i]=new double[np+1];

    //fijar el radio de la band r_max
    fix_r_max(np);

    for(i=0; i<=nq; i++)
    {
        for(j=0; j<=np; j++)
        {
            m[i][j]=DBL_MAX;
        }
    }
    // Punto de inicio en 0.0
    m[0][0]=0;

    for(i=1; i<=nq; i++)
    {
        pos=((i-1.0)/(nq-1.0))*(np-1.0); // su equivalente en P
        b_min=bandMin(pos, np);
        b_max=bandMax(pos, np);

        //notar que la band es inclusiva, por lo que sumo 1 en ambos lados
        for(j=1+b_min; j<=1+b_max; j++)
        {
            //insertar, borrar, reemplazar
            m[i][j]=dist_lineal->dist(vp[j-1], vq[i-1], sp[j-1], sq[i-1], D) 
                    + min3(m[i-1][j], m[i][j-1], m[i-1][j-1]);
        }
    }

    d=m[nq][np];//DTW, resultado de la alineacion

    // Time-Normalized Distance Measure (solve the scaling time)    
    // symmetric distance
    d = d/(nq + np); 

    for(i=0; i<=nq; i++)
        delete[] m[i];
    delete[] m;

    return d;
}


double BandDTW::dist_by_levels(double* vp, double* vq, double* sp, double* sq,  int np, int nq, double th_dist)
{
    double d = 0;

    int n = std::min(np, nq);

    int izq, der, max_level = log2(n + 1) - 1;

    for (int l = 0; l <= max_level && d <= th_dist; ++l)
    {
        izq = pow(2, l);
        der = pow(2, l + 1) - 1;

        d += partial_dist(vp, vq, sp, sq, izq - 1, der - 1);        
    }

    return d/n;
    //return d;
}


double BandDTW::partial_dist(double *vp, double *vq, double *sp, double *sq, int izq, int der)
{
    double d;

    int i,j,b_min,b_max,pos;

    int n = der - izq + 1;

    double** m=new double*[n+1];
    for(i=0; i<n + 1; i++)
        m[i]=new double[n+1];

    //fijar el radio de la band r_max
    fix_r_max(n);

    for(i=0; i<=n; i++)
    {
        for(j=0; j<=n; j++)
        {
            m[i][j]=DBL_MAX;
        }
    }
    // Punto de inicio en 0.0
    m[0][0]=0;

    for(i=1; i<=n; i++)
    {
        pos = (i - 1.0); 
        b_min=bandMin(pos, n);
        b_max=bandMax(pos, n);

        //notar que la band es inclusiva, por lo que sumo 1 en ambos lados
        for(j=1+b_min; j<=1+b_max; j++)
        {
            //insertar, borrar, reemplazar
            m[i][j]=dist_lineal->dist(vp[izq + j-1], vq[izq + i-1], sp[izq + j-1], sq[izq + i-1]) 
                    + min3(m[i-1][j], m[i][j-1], m[i-1][j-1]);
        }
    }

    d=m[n][n];//DTW, resultado de la alineacion

    for(i=0; i<=n; i++)
        delete[] m[i];
    delete[] m;

    return d;
}


double BandDTW::dist_by_levels(double** vp, double** vq, double** sp, double** sq,  int np, int nq, int D, double th_dist)
{
    double  d = 0;

    int n = std::min(np, nq);

    int izq, der, max_level = log2(n + 1) - 1;

    for (int l = 0; l <=max_level && d <= th_dist; ++l)
    {
        izq = pow(2, l);
        der = pow(2, l + 1) - 1;

        d += partial_dist(vp, vq, sp, sq, izq - 1, der - 1, D) ;
        //d += partial_dist(vp, vq, sp, sq, izq - 1, der - 1, D) / (1 + l);
    }

    return d/n;
    //return d;
}


double BandDTW::partial_dist(double** vp, double** vq, double** sp, double** sq, int izq, int der, int D)
{
    double d;

    int i,j,b_min,b_max,pos;

    int n = der - izq + 1;

    double** m=new double*[n+1];
    for(i=0; i<n + 1; i++)
        m[i]=new double[n+1];

    //fijar el radio de la band r_max
    fix_r_max(n);

    for(i=0; i<=n; i++)
    {
        for(j=0; j<=n; j++)
        {
            m[i][j]=DBL_MAX;
        }
    }
    // Punto de inicio en 0.0
    m[0][0]=0;

    for(i=1; i<=n; i++)
    {
        pos = (i - 1.0); 
        b_min=bandMin(pos, n);
        b_max=bandMax(pos, n);

        //notar que la band es inclusiva, por lo que sumo 1 en ambos lados
        for(j=1+b_min; j<=1+b_max; j++)
        {
            //insertar, borrar, reemplazar
            m[i][j]=dist_lineal->dist(vp[izq + j-1], vq[izq + i-1], sp[izq + j-1], sq[izq + i-1], D) 
                    + min3(m[i-1][j], m[i][j-1], m[i-1][j-1]);
        }
    }

    d=m[n][n];//DTW, resultado de la alineacion

    for(i=0; i<=n; i++)
        delete[] m[i];
    delete[] m;

    return d;
}



int BandDTW::bandMin(int pos, int dimension)
{
	if(type_band == NON_BAND)
	{
		return 0;
	}
    else if(type_band == SAKOE_CHIBA)
    {
        return std::max(0, pos - r_max);
    }
    else if(type_band == ITAKURA)
    {
        int r_local=0;
        if(pos < dimension/2)
        {
            r_local=1+round((((float)r_max*2)/dimension)*pos);
        }
        else
        {
            r_local=1+round((((float)r_max*2)/dimension)*(dimension-pos-1));
        }
        return std::max(0,pos-r_local);
    }
    
    return -1;
    
}

int BandDTW::bandMax(int pos, int dimension)
{
	if(type_band == NON_BAND)
	{
        return dimension-1;
    }
	else if(type_band == SAKOE_CHIBA)
    {
        return std::min(pos+r_max,dimension-1);
    }
    else if(type_band == ITAKURA)
    {
        int r_local=0;
        if(pos < dimension/2)
        {
            //return (int)(c1+c2*pos);
            r_local=1+round((r_max*2.0/dimension)*pos);            
        }
        else
        {
            //return (int)(c1+c2*(q.dimension()-pos));
            r_local=1+round((r_max*2.0/dimension)*(dimension-pos-1));
        }        
        return std::min(pos+r_local,dimension-1);
    }
    
    return dimension;    
}

double BandDTW::lower_bound(Object *p, Object *q)
{
    double* _p=NULL,*_q=NULL;    
    int np=p->get_n();
    int nq=q->get_n();
    int n;
    double lb=0.0,val=0;
    
    _p=((Vector*)p)->get_elems();
	_q=((Vector*)q)->get_elems();
	n=std::min(np,nq);
            
	// igualar longitud de ambas series con Piecewise Aggregate Approximation (PAA)	
    if(np > n)
    {
		_p=Paa::get_paa(_p,np,n);
	}
	
	if(nq > n)
	{
		_q=Paa::get_paa(_q,nq,n);
	}
			
    double* U=get_U(_q,n);
	double* L=get_L(_q,n);
    
    for(int i=0; i<n; i++)
    {		
        val = 0.0;     
        if(_p[i] > U[i])
        {
            val=fabs(_p[i] - U[i]);                
        }
        else if(_p[i] < L[i])
        {
            val=fabs(_p[i] - L[i]);                        
        }
        
        if(dist_lineal->get_type_lineal() == 1)
            lb += val;
        else
			lb += pow(val,2);	
    }

    lb=lb/(nq+np); // Time-Normalized Distance Measure (similarly to DTW dist)
    
    //liberar memoria si se creo una representacion PAA
    if(np > n)
	{		
		delete[] _p;
	}
	if(nq > n)
	{		
		delete[] _q;
	}
	
	delete[] U;
	delete[] L;
	
    if(dist_lineal->get_type_lineal() == 1)
		return lb;
	else    
		return sqrt(lb);
}

double BandDTW::lower_bound(Object *p, double *L, double *U, int n)
{
    double* _p=NULL;    
    int np=p->get_n();    
    
    double lb=0.0,val=0;
    
    _p=((Vector*)p)->get_elems();
	        
	// igualar longitud de ambas series con Piecewise Aggregate Approximation (PAA)	
    if(np > n)
    {
		_p=Paa::get_paa(_p, np, n);
	}	
    
    for(int i=0; i<n; i++)
    {		 
		val = 0.0;             
        if(_p[i] > U[i])
        {
            val=fabs(_p[i] - U[i]);                
        }
        else if(_p[i] < L[i])
        {
            val=fabs(_p[i] - L[i]);                        
        }
        
        if(dist_lineal->get_type_lineal() == 1)
            lb+=val;
        else
			lb+=pow(val,2);	
    }

    lb=lb/(n+np); // Time-Normalized Distance Measure (similarly to DTW dist)
    
    //liberar memoria si se creo una representacion PAA
    if(np > n)
	{		
		delete[] _p;
	}
		
    if(dist_lineal->get_type_lineal() == 1)
		return lb;
	else    
		return sqrt(lb);
}

double BandDTW::upper_bound(Object *p,Object *q)
{
    double* _p=NULL, *_q=NULL;
    int np=p->get_n();
    int nq=q->get_n();
	int n;
	double ub=0.0,val=0;
        
    _p=((Vector*)p)->get_elems();
	_q=((Vector*)q)->get_elems();
	n=std::min(np,nq);
            
	// igualar longitud de ambas series con Piecewise Aggregate Approximation (PAA)	
    if(np>n)
    {
		_p=Paa::get_paa(_p,np,n);
	}
	
	if(nq>n)
	{
		_q=Paa::get_paa(_q,nq,n);
	}
		   
	double* U=get_U(_q,n);
	double* L=get_L(_q,n);	   
	
    for(int i=0; i<n; i++)
    {       
        val=std::max(fabs(U[i] -_p[i]), fabs(_p[i] - L[i]));
        if(dist_lineal->get_type_lineal() == 1)
			ub+=val;
		else
			ub+=pow(val,2);	
    }
    
    ub=(2*r_max+1.0)*ub;

    ub=ub/(nq+np); // Time-Normalized Distance Measure (similarly to DTW dist)
        
    
    //liberar memoria si se cre\'o una representacion PAA
	if(np>n)
	{		
		delete[] _p;
	}
	if(nq>n)
	{		
		delete[] _q;
	}
	
	delete[] U;
	delete[] L;
		
    if(dist_lineal->get_type_lineal() == 1)
		return ub;
	else    
		return sqrt(ub);
}

double BandDTW::upper_bound(Object *p, double *L, double *U, int n)
{
    double* _p=NULL;
    int np=p->get_n();    
	double ub=0.0,val=0;
        
    _p=((Vector*)p)->get_elems();
	        
	// igualar longitud de ambas series con Piecewise Aggregate Approximation (PAA)	
    if(np>n)
    {
		_p=Paa::get_paa(_p,np,n);
	}
	
	
    for(int i=0; i<n; i++)
    {       
        val=std::max(fabs(U[i] -_p[i]), fabs(_p[i] - L[i]));
        if(dist_lineal->get_type_lineal() == 1)
			ub+=val;
		else
			ub+=pow(val,2);	
    }
    ub=(2*r_max+1.0)*ub;

    ub=ub/(n+np); // Time-Normalized Distance Measure (similarly to DTW dist)
        
    
    //liberar memoria si se cre\'o una representacion PAA
	if(np>n)
	{		
		delete[] _p;
	}
	
		
    if(dist_lineal->get_type_lineal() == 1)
		return ub;
	else    
		return sqrt(ub);
}

double* BandDTW::get_U(double* q,int nq)
{    
	//fijar el radio de la band r_max 
	fix_r_max(nq);
	
	double* U=new double[nq];
	if(type_band == NON_BAND)
	{//si no se utiliza band, usar el maximo
		double d_max=0.0;
		for(int pos=0;pos<nq;pos++)
			d_max=std::max(q[pos],d_max);
			
		for(int pos=0;pos<nq;pos++)
			U[pos]=d_max;	
	}
	else
	{    
		for(int pos=0;pos<nq;pos++)
		{
			U[pos]=0.0;
			for(int i=bandMin(pos, nq); i<=bandMax(pos, nq); i++)
			{
				if(q[i] > U[pos])
				{
					U[pos]=q[i];
				}
			}
		}
	}
    return U;
}

double* BandDTW::get_L(double* q,int nq)
{    
	//fijar el radio de la band r_max 
	fix_r_max(nq);
	 
    double* L=new double[nq];
    if(type_band == NON_BAND)
     {//si no se utiliza band, usar el minimo 
		double d_min=DBL_MAX;
		for(int pos=0;pos<nq;pos++)
			d_min=std::min(q[pos],d_min);
			
		for(int pos=0;pos<nq;pos++)
			L[pos]=d_min;	
	 }
	 else
	 {    
		for(int pos=0;pos<nq;pos++)
		{
			L[pos]=DBL_MAX;
			for(int i=bandMin(pos, nq); i<=bandMax(pos, nq); i++)
			{
				if(q[i] < L[pos])
				{
					L[pos]=q[i];
				}
			}
		}
	}
    return L;
}

void BandDTW::fix_r_max(int n)
{
	if(type_band == NON_BAND)
		r_max = round(n/2.0);
	else if(fraction_r > -1)    
        r_max = round(fraction_r * n);        
}

int BandDTW::get_type_lineal()
{
    return dist_lineal->get_type_lineal();
}

