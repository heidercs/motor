#include <iostream>
#include <math.h>
#include "controlband.h"
#include "../util/functions.h"


double* ControlBand::suavizar(double* p, int n, int k)
{
    double* suave=new double[n];
    int init, end;
    for(int i=0; i<n; i++)
    {
        init=std::max(0,i-(k/2));
        end=std::min(init + k, n);
        suave[i]=media(p,init,end);
    }
    return suave;
}



//retorna las posiciones de los puntos accidentales
// k reduce, x aumenta
bool* ControlBand::get_control_band(double* p,int n, double fact_k, double fact_x)
{
    bool* cb=new bool[n];
    int k=std::max(3.0,fact_k*n);
    double* suave=suavizar(p,n,k);
    double rad=fact_x*desviacion(p,n);

    for(int i=0; i<n; i++)
    {
        cb[i]=(fabs(p[i] - suave[i])>rad);
    }
    delete[] suave;
    return cb;
}
