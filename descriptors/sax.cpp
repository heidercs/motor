#include <math.h>
#include <float.h>
#include <iostream>
#include <stdlib.h>

#include "sax.h"

#include "paa.h"
#include "../util/functions.h"
#include "../util/normaldistribution.h"

////////////////  static constants /////////////////////////
constexpr double Sax::quantil_2[];
constexpr double Sax::quantil_3[];
constexpr double Sax::quantil_4[];
constexpr double Sax::quantil_5[];
constexpr double Sax::quantil_6[];
constexpr double Sax::quantil_7[];
constexpr double Sax::quantil_8[];
constexpr double Sax::quantil_9[];
constexpr double Sax::quantil_10[];
constexpr double Sax::quantil_16[];
////////////////////////////////////////////////////////////


Sax::Sax(int _card, double _mean, double _sd)
{
    card=_card;
    tam_bits = ceil(log2(card));
    mean=_mean;
    sd=_sd;
    
    if(card<2 || (card>10 && card!=16))
    {
        std::cout<<"This number of regions is not considered\n";
        exit(0);
    }
    lookup_table=NULL;

    compute_breakpoints();
}



double Sax::get_quantil(int p)
{
    if(p <= 0)      return -DBL_MAX; // first breakpoint = - INF
    if(p >= card)   return +DBL_MAX; // last break point = + INF

    return break_points[p-1];
}

void Sax::compute_breakpoints()
{   
    //numero de breakpoints: card - 1, 
    if(mean == 0 && sd == 1) 
    {    
        const double *sel;
        switch(card)
        {
            case 2: sel = quantil_2; break;
            case 3: sel = quantil_3; break;
            case 4: sel = quantil_4; break;
            case 5: sel = quantil_5; break;
            case 6: sel = quantil_6; break;
            case 7: sel = quantil_7; break;
            case 8: sel = quantil_8; break;
            case 9: sel = quantil_9; break;
            case 10: sel = quantil_10;   break;
            case 16: sel = quantil_16;   break;        
        }

        break_points = new double[card - 1];   
        std::copy(sel, sel + card -1, break_points);
    }
    else
    {
        break_points = new double[card - 1];
        for(int i = 1;i <= card-1; i++)
        {
            break_points[i-1] = qnorm(1.0*i/card, mean, sd);
        }
    }    
}

int* Sax::run_sax(double*p, int n,  int new_n)
{
    double* paa=Paa::get_paa(p, n, new_n);
    int* sax=new int[new_n];

    int i,j;

    for(i=0; i<new_n; i++)
    {
        for(j=0; j<card-1; j++)
        {
            if(paa[i]<break_points[j])
            {
                sax[i]=j;
                break;
            }
        }
        if(j==card-1)
            sax[i]=j;
    }

    delete[] paa;

    return sax;
}

int* Sax::run_sax_from_paa(double* paa, int new_n)
{
    int* sax=new int[new_n];

    int i,j;

    for(i=0; i<new_n; i++)
    {
        for(j=0; j<card-1; j++)
        {
            if(paa[i]<break_points[j])
            {
                sax[i]=j;
                break;
            }
        }
        if(j==card-1)
            sax[i]=j;
    }

    return sax;
}


bool* Sax::get_sax_bits(double*p, int n,  int new_n)
{    
    int* sax_int = run_sax(p, n, new_n);    
    bool *bits=convert_bits(sax_int, new_n);    
    
    delete[] sax_int;

    return bits;
}


bool* Sax::convert_bits(int* sax, int n)
{
    int new_n = tam_bits * n;

    bool* bits=new bool[new_n];

    for(int i=0; i<new_n; i++)
        bits[i]=0;

    int k=0,num;

    for(int i=0; i<new_n; i+=tam_bits)
    {
        num=sax[i/tam_bits];
        k=tam_bits-1;
        while(num>0)
        {
            bits[i + k--] = num % 2;
            num=num/2;
        }
    }
    return bits;
}

double Sax::cost_sax(int s1, int s2)
{
    return lookup_table[s1][s2];
}


void Sax::construir_lookup()
{
    //construir la tabla de similitud entre simbolos SAX
    lookup_table=new double*[card];
    
    for(int i=0;i<card;i++)
    {
        lookup_table[i]=new double[card];
        for(int j=0;j<card;j++)
        {
            lookup_table[i][j]=0;
            if(abs(i-j) > 1)
                lookup_table[i][j]=get_quantil(std::max(i,j)) - get_quantil(std::min(i,j)+1) ;            
        }        
    }
}


int Sax::get_tambits()
{
    return this->tam_bits;
}


int Sax::get_card()
{
    return this->card;
}

Sax::~Sax()
{
    if(lookup_table!=NULL)
    {
        for(int i=0;i<card;i++)
            delete[] lookup_table[i];

        delete[] lookup_table;
    }

    delete[] break_points;
}

///////////////////////////////////////////////////////////


int* Sax::get_sax(double*p, int n, int new_n, int card)
{   
    Sax sax(card);
    return sax.run_sax(p, n, new_n);
}

double* Sax::get_sax_double(double*p, int n, int new_n, int card)
{    
    Sax sax(card);
    int* s = sax.run_sax(p, n, new_n);   
    double *r = new double[new_n];
    std::copy(s, s + new_n, r);

    delete[] s;

    return r;
}


///////////////////////////////////////////////////////////////

SaxFile::SaxFile()
{
    this->serie=NULL;
    this->repr_sax=NULL;
    this->n_sax=0;
}

SaxFile::SaxFile(Vector *_serie, int* _repr_sax, int _n_sax)
{
    this->serie=_serie;
    this->repr_sax=_repr_sax;
    this->n_sax=_n_sax;
}

SaxFile::SaxFile(int* _repr_sax,int _n_sax)
{
    this->repr_sax=_repr_sax;
    this->n_sax=_n_sax;
}


SaxFile::SaxFile(Vector *_serie)
{
    this->serie=_serie;
}

