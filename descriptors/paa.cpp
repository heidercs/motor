#include "paa.h"
#include <math.h>
#include <iostream>

int* Paa::get_paa(int* p, int n, int new_N)
{
	if(new_N >= n) return p;
	
    int* paa=new int[new_N];
    int curr_dim, init=0;
    double frac=n/(new_N*1.0);
    double res=0;

    for(int i=0; i<new_N; i++)
    {
        curr_dim=(int)(frac + res);
        res=(frac + res) - curr_dim;
        paa[i]=0;        
        for(int j=init; j<init+curr_dim; j++)
        {
            paa[i]+=p[j];
        }
        paa[i]=paa[i]/(curr_dim*1.0);
        init+=curr_dim;
    }

    return paa;
}

double* Paa::get_paa(double* p, int n, int new_N)
{
	double* paa=new double[new_N];
	
	if(new_N >= n) 
	{
		for (int i = 0; i < n; i++)		
			paa[i] = p[i];
		return paa;
	}
	
    int curr_dim, init=0;
    double frac=n/(new_N*1.0);
    double res=0;

    for(int i=0; i<new_N; i++)
    {
        curr_dim=(int)(frac + res);
        res=(frac + res) - curr_dim;
        paa[i]=0;        
        for(int j=init; j<init+curr_dim; j++)
        {
            paa[i]+=p[j];
        }
        paa[i]=paa[i]/(curr_dim*1.0);
        init+=curr_dim;
    }

    return paa;
}
