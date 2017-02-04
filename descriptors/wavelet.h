#ifndef WAVELET_H
#define WAVELET_H

#include <iostream>
#include <cmath>

class Wavelet
{
public:
    /** Haar Transform **/
    static double* get_haar(double *vec, int n,int w=0);
    static double** get_haar2D(double **matrix, int n, int m);
};


#endif
