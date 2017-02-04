#ifndef FOURIER_H
#define FOURIER_H

#include "../util/complex.h"

class Fourier
{

public:
    static double* get_coef_fourier(double* v, int n,int num_coef);
    static double** get_coef_fourier2D(double** ma, int n,int m,int num_coef);

};

#endif // FURIER_H
